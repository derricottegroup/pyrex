"""
An Implementation of Euler Integration for reaction path following.

References:
[1] C. Gonzalez and H. B. Schlegel, J. Chem. Phys. 90(4):2154 (An Improved Algorithm for Reaction-Path Following)
"""

__authors__  = "Wallace D. Derricotte"
__credits__  = ["Wallace D. Derricotte"]

__copyright__  = "(c) 2018-2019, Derricotte Research Group"
__license__    = "MIT License"
__date__       = "2019-01-02"

import numpy as np
import psi4
import os
import sys
import json
from pyscf import gto, scf, dft, grad, solvent, mp 
from pyscf.solvent import ddcosmo, ddcosmo_grad

#####################
# Read in Parameters#
#####################

class Params():
    def __init__(self):
        """
            Initialize .json file provided by user in command line, read input and store variables.
        """
        self.grace_period = 50
        self.e_conv = 1e-5
        self.qm_program = 'psi4'
        self.nthreads = 1
        self.do_solvent = False
        #self.damp = -0.1
        json_input = sys.argv[1]
        self.read_input(json_input)            
        
    def json2xyz(self, input_params):
        """
            Geometries in the .json format (MOLSSI standard) are given as an array. This functions
            takes the given geometry and converts it to a string with a format similar to xyz files.
        """
        charge = int(input_params['molecule']['molecular_charge'])
        mult = int(input_params['molecule']['molecular_multiplicity'])
        geom = ''
        geom += '\n%d %d\n' %(charge, mult)
        symbols = input_params['molecule']['symbols']
        geometry = input_params['molecule']['geometry']
        for i in range(len(symbols)):
            geom += "%s  " %symbols[i]
            for j in range(3):
                if(j==2):
                    geom += "   %f   \n" %geometry[(i + 2*i) + j]
                else:
                    geom += "   %f   " %geometry[(i + 2*i) + j]
        return geom

    def read_input(self,json_input):
        """ 
            Reads Input from .json input file. (date: 01/02/19 - WDD)

            General Parameters Read in from Input:
            -------------------------------------
                geometry(array) -- Geometry array is read in and converted from json->xyz format
                basis(string) -- Specifies orbital basis set for QM calculations
                method(string) -- Specifies type of QM calculations
                symbols(array) -- Array of symbols specifying each atom, should correspond with rows
                in the geometry specification.
                natoms(int) -- Takes the length of the symbols array to get number of atoms.
                keywords(map) -- Any valid psi4 keywords can be placed here for QM calculations.

            IRC Specific Parameters Read in from Input:
            ------------------------------------------
                direction(string) -- specify IRC to perform a "forward" or "backward" descent.
                step_size(float) -- step size to be taken by IRC algorithm in units amu^(1/2)*bohr.
                mode(int) -- integer corresponding to normal mode you want to follow.
                normal_mode_file(string) -- file that contains normal mode output (currently only
                compatible with psi4 normal mode writer files)
        """
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        if 'molecule' in input_params:
            if 'geometry' in input_params['molecule']:
                self.geometry = self.json2xyz(input_params)
            if 'symbols' in input_params['molecule']:
                self.symbols = input_params['molecule']['symbols']
                self.natoms = len(input_params['molecule']['symbols'])
            if 'molecular_charge' in input_params['molecule']:
                self.charge = int(input_params['molecule']['molecular_charge'])
            if 'molecular_multiplicity' in input_params['molecule']:
                self.mult = int(input_params['molecule']['molecular_multiplicity'])
        if 'model' in input_params:
            if 'basis' in input_params['model']:
                self.basis = input_params['model']['basis']
            if 'method' in input_params['model']:
                self.method = input_params['model']['method']
        if 'keywords' in input_params:
            self.keywords = input_params['keywords']
        if 'irc' in input_params:    
            if 'direction' in input_params['irc']:
                self.direction = input_params['irc']['direction']
            if 'step_size' in input_params['irc']:
                self.step_size = input_params['irc']['step_size']
            if 'mode' in input_params['irc']:
                self.mode = input_params['irc']['mode']
            if 'normal_mode_file' in input_params['irc']:
                self.normal_mode_file = input_params['irc']['normal_mode_file']
                self.ts_vec = self.normal_mode_reader()
            if 'grace_period' in input_params['irc']:
                self.grade_period = input_params['irc']['grace_period']
            if 'e_conv' in input_params['irc']:
                self.e_conv = input_params['irc']['e_conv']
        if 'pyrex' in input_params:
            if 'qm_program' in input_params['pyrex']:
                self.qm_program = input_params['pyrex']['qm_program']
            if 'xc_functional' in input_params['pyrex']:
                self.xc_functional = input_params['pyrex']['xc_functional']
            if 'eps' in input_params['pyrex']:
                self.eps = input_params['pyrex']['eps']
            if 'do_solvent' in input_params['pyrex']:
                self.do_solvent = input_params['pyrex']['do_solvent']
            if 'nthreads' in input_params['pyrex']:
                self.nthreads = input_params['pyrex']['nthreads'] 
    def normal_mode_reader(self):
        """
            Reads the file containing the normal modes of vibration. This function currently
            only works with the Psi4 normal mode output. In order to produce a compatible file
            run a frequency calculation in Psi4 with the option "normal_mode_writer" set True.
            
            Parameters:
            ----------
                self(self) -- contains all shared parameters.
            Returns:
                ts_vec(array) -- Array containing the normal mode vector.
        """
        natoms = self.natoms
        direction = self.direction
        normal_mode_file = open(self.normal_mode_file, 'r')
        ts_vec = []
        for line in normal_mode_file:
            if("vibration %d" %self.mode in line):
                line = next(normal_mode_file)
                for i in range(natoms):
                    trj = line.split()
                    trj = list(map(float, trj))
                    np_trj = np.array(trj)
                    if(direction=='backward'):
                        ts_vec.append(-1.0*np_trj)
                    else:
                        ts_vec.append(np_trj)
                    line = next(normal_mode_file)
        return ts_vec

########################
## Gradient Functions ##
########################

def energy_calc(params, current_geom, mol):
    energy = 0.0
    if(params.qm_program=='pyscf'):
        pymol = gto.Mole()
        pymol.verbose = 0
        geom_vec = []
        for i in range(params.natoms):
            atom = [params.symbols[i],]
            atom_coords = []
            for j in range(3):
                atom_coords.append(current_geom[i][j])
            atom_coords = tuple(atom_coords)
            atom.append(atom_coords)
            geom_vec.append(atom)
        #print(geom_vec)
        pymol.atom = geom_vec
        pymol.unit = 'Bohr'
        pymol.basis = params.basis
        pymol.charge = params.charge
        pymol.spin = params.mult - 1
        pymol.build()
        if(params.method == "scf"):
            scf_obj = scf.RHF(pymol)
        if(params.method == "mp2"): #TODO Doesn't work yet. Few things left to figure out. 
            mf = scf.RHF(pymol).kernel()
            scf_obj = mp.MP2(mf)
        #if(params.method == "dft"):
        #    scf_obj = dft.RKS(mol)
        #    scf_obj.xc = params.xc_functional
        if(params.do_solvent):
            solv_obj = ddcosmo.ddcosmo_for_scf(scf_obj)
            solv_obj.with_solvent.eps = params.eps
            solv_cosmo = solvent.ddCOSMO(solv_obj)
            energy = solv_cosmo.kernel()
        else:
            energy = scf_obj.kernel()

    if(params.qm_program=='psi4'):
        mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        grad_method = "%s/%s" %(params.method,params.basis)
        psi4.core.set_output_file("psi4_out.dat", False)
        psi4.set_options(params.keywords)
        psi4.set_num_threads(params.nthreads)
        energy = psi4.energy(grad_method) 
    return energy


def grad_calc(params,current_geom, mol):
    """
        Uses Psi4 to calculate the energy gradient and returns the mass-weighted
        gradient and energy. Here any of the keywords the user provides in the 
        .json input are used to set the options for the energy calculation.

        Parameters:
        ----------
            params(self) -- contains initialized shared parameters.
            current_geom(np array) -- Matrix of size natoms x 3 containing the geometry.
            mol(psi4.Molecule) -- Psi4 molecule object containing the current molecule.
        Returns:
        -------
            grad_mw(np array) -- Mass weighted gradient matrix of size natoms x 3.
            E(float) -- single-point energy from Psi4 calculation. 
    """
    if(params.qm_program=='pyscf'):
        pymol = gto.Mole()
        pymol.verbose = 0
        geom_vec = []
        for i in range(params.natoms):
            atom = [params.symbols[i],]
            atom_coords = []
            for j in range(3):
                atom_coords.append(current_geom[i][j])
            atom_coords = tuple(atom_coords)
            atom.append(atom_coords)
            geom_vec.append(atom)
        #print(geom_vec)
        pymol.atom = geom_vec
        pymol.unit = 'Bohr'
        pymol.basis = params.basis
        pymol.charge = params.charge
        pymol.spin = params.mult - 1
        pymol.build()
        scf_obj = scf.RHF(pymol)
        if(params.do_solvent):
            solv_obj = ddcosmo.ddcosmo_for_scf(scf_obj)
            solv_obj.with_solvent.eps = params.eps
            solv_cosmo = solvent.ddCOSMO(solv_obj).run()
            E = solv_cosmo.scf()
            grad = solv_obj.nuc_grad_method().kernel()
        else:
            E = scf_obj.kernel()
            grad = scf_obj.nuc_grad_method().kernel()
        #print(grad)
        grad_mw = mass_weight(params.natoms, grad, mol)        
    if(params.qm_program=='psi4'):
        mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        grad_method = "%s/%s" %(params.method,params.basis)
        psi4.core.set_output_file("psi4_out.dat", False)
        psi4.set_options(params.keywords)
        psi4.set_num_threads(params.nthreads)
        E, wfn = psi4.energy(grad_method,return_wfn=True)
        psi4.set_num_threads(params.nthreads)
        grad = np.asarray(psi4.gradient(grad_method,ref_wfn=wfn))
        #print(grad)
        grad_mw = mass_weight(params.natoms, grad, mol)
    return grad_mw, E


def parabolic_fit(xs, ys):
    fit = np.polyfit(xs, ys, deg=2)
    fit = np.poly1d(fit)
    minima = fit.deriv().r
    real_minima = minima[minima.imag==0].real
    return real_minima


def ishida_morokuma(output_file):
    """
        This function runs the Ishida-Morokuma irc procedure
    """
    max_steps = 1000
    params = Params()
    line_step_size = 0.3333*params.step_size
    #line_step_size = 0.025*params.step_size
    current_geom = params.geometry
    mol = psi4.geometry(params.geometry)
    starting_vec = np.asarray(params.ts_vec)
    grad_method = "%s/%s" %(params.method,params.basis)
    steps = 0
    E = 0.0
    previous_E = 0.0
    del_E = 0.0
    energies = []
    current_geom = np.asarray(mol.geometry())
    #print(current_geom)
    output = open(output_file, "a")
    output.write('\n\n--Intrinsic Reaction Coordinate (%s)--\n' %(params.direction))
    output.write('\n--------------------------------------------------------------------------------------\n')
    output.write('\n{:>20} {:>20} {:>20} {:>20}\n'.format('Coordinate', 'E', 'Delta E', 'Gradient Norm'))
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    last_energy = None
    while (steps <= max_steps):
        if(steps==0):
            grad_0 = mass_weight(params.natoms, starting_vec, mol)
            E_0 = energy_calc(params, current_geom, mol)
        else:
            grad_0, E_0  = grad_calc(params, current_geom, mol)
        
        if(last_energy):
            del_E = E_0 - last_energy
        if((last_energy and (E_0 > last_energy) and steps > params.grace_period)):
            output = open(output_file, "a") 
            print("\nIRC Energy Has Increased! You're Likely Near a Minimum!\n")
            output.close()
            break
        if(last_energy and np.abs(del_E)<params.e_conv and steps > params.grace_period):
            output = open(output_file, "a")
            print("\nIRC Has Converged!\n")
            output.close()
            break
        mol.save_xyz_file('imk_step_'+str(steps)+'.xyz',False)
        coords_1 = euler_step(params.natoms, current_geom, grad_0,params.step_size,mol)
        current_geom = coords_1
        #mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        grad_1, E_1 = grad_calc(params, current_geom, mol) 
        grad_0_norm = np.linalg.norm(grad_0)
        grad_1_norm = np.linalg.norm(grad_1)

        # Calculate Bisector (Eq. 6)
        D = grad_0/grad_0_norm - grad_1/grad_1_norm
        D_normed = D/np.linalg.norm(D)

        line_xs = [0,]
        line_energies = [E_1,]

        line_step_size_thresh = 1.5*line_step_size
        #line_step_size_thresh = 2.0*line_step_size
        
        # Find useful point by projecting grad_1 on D
        grad_1_normed = grad_1/grad_1_norm
        step_D1 = grad_1*D_normed*D_normed*line_step_size
        step_D1_norm = np.linalg.norm(step_D1)
       # if step_D1_norm < line_step_size_thresh:
       #     coords_1 = mass_weight_geom(params.natoms, coords_1, mol)
       #     current_geom = coords_1 + step_D1
       #     current_geom = un_mass_weight_geom(params.natoms, current_geom, mol) 
       #     mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
       #     step_D1_E = psi4.energy(grad_method)
       #     line_xs.append(step_D1_norm)
       #     line_energies.append(step_D1_E)
        # Otherwise take a step along D
       # else:
        step_D2 = line_step_size*D_normed
        step_D2_norm = np.linalg.norm(step_D2)
        coords_1 = mass_weight_geom(params.natoms, coords_1, mol)
        current_geom = coords_1 + step_D2
        current_geom = un_mass_weight_geom(params.natoms, coords_1, mol)
        #mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        #step_D2_E = psi4.energy(grad_method)
        step_D2_E = energy_calc(params, current_geom, mol)
        #line_xs.append(step_D2_norm)
        line_xs.append(step_D2_norm)
        line_energies.append(step_D2_E)

        # Calculate 3rd point by taking a half step size
        if(line_energies[1] >= line_energies[0]):
            step_D3 = 0.5*line_step_size*D_normed # Half Step Size
            #new_del = 0.5*line_step_size
        else:
            step_D3 = 2.0*line_step_size*D_normed #Double Step Size
            #new_del = 2.0*line_step_size
        
        step_D3_norm = np.linalg.norm(step_D3)
        #coords_1 = mass_weight_geom(params.natoms, coords_1, mol)
        current_geom = coords_1 + step_D3
        current_geom = un_mass_weight_geom(params.natoms, current_geom, mol)
        #mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        #step_D3_E = psi4.energy(grad_method)
        step_D3_E = energy_calc(params, current_geom, mol)
        line_xs.append(step_D3_norm)
        line_energies.append(step_D3_E) 

        real_minimum = parabolic_fit(line_xs, line_energies)
        
        current_geom = coords_1 + (real_minimum*D_normed)
        current_geom = un_mass_weight_geom(params.natoms, current_geom, mol)
        
        mol.set_geometry(psi4.core.Matrix.from_array(current_geom))

        last_energy = E_0
        
        #params.damp -= 0.5

        if(params.direction=="backward"):
            coord = -1*steps*params.step_size
        else:
            coord = steps*params.step_size
        print_step(output_file,coord, E_0, del_E, grad_0_norm)
        steps = steps+1
    output = open(output_file, "a")
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    with open('irc_%s.xyz' %params.direction,'w') as outfile:
        for i in range(steps):
            with open('imk_step_'+str(i)+'.xyz')as infile:
                outfile.write(infile.read())
        os.system('rm imk_step*')

            


def mass_weight(natoms,grad, mol):
    """
        Mass weights the given gradient

        Parameters:
        ----------
            natoms(int) -- number of atoms in the molecule.
            grad(np array) -- matrix of size natoms x 3 containing the gradients.
        Returns:
        -------
            grad_mw(np array) -- Mass weighted gradient matrix of size natoms x 3.
    """
    grad_mw = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            grad_mw[i][j] = grad[i][j]/np.sqrt(mol.mass(i))
    return grad_mw

def mass_weight_geom(natoms,geom,mol):
    """
        Mass weights the given coordinates

        Parameters:
        ----------
            natoms(int) -- number of atoms in the molecule.
            geom(np array) -- matrix of size natoms x 3 containing the cartesian geometry (bohr).
        Returns:
        -------
            coord_mw(np array) -- Mass weighted coordinate matrix of size natoms x 3 in units
            amu^(1/2)*bohr.
    """
    coord_mw = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            coord_mw[i][j] = geom[i][j]*np.sqrt(mol.mass(i))
    return coord_mw

def un_mass_weight_geom(natoms,geom_mw,mol):
    """
        Un-Mass weights the given coordinates. Displacements are done along the mass weighted
        coordinate and then un-mass-weighted prior to gradient calculations.

        Parameters:
        ----------
            natoms(int) -- number of atoms in the molecule.
            geom_mw(np array) -- matrix of size natoms x 3 containing the cartesian geometry
            in units amu^(1/2)*bohr
        Returns:
        -------
            coord(np array) -- Coordinate matrix of size natoms x 3 in units bohr.
    """
    coord = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            coord[i][j] = geom_mw[i][j]/np.sqrt(mol.mass(i))
    return coord

def euler_step(natoms,current_geom,grad,step_size,mol):
    """
        Take a single Euler step along the gradient. The coordinates are first mass weighted prior
        to the Euler step and then un-mass-weighted for printing/gradient calculation. This 
        procedure is adapted from Equation 2 in Ref [1].

        Parameters:
        ----------
            natoms(int) -- number of atoms in the molecule.
            current_geom(np array) -- matrix of size natoms x 3 containing the cartesian 
            geometry in units bohr.
            grad(np array) -- mass-weighted gradient matrix of size natoms x 3 containing
            the gradients in units amu^(1/2)*bohr.
            step_size(float) -- user provided IRC step size in units amu^(1/2)*bohr. 
    """
    current_geom = mass_weight_geom(natoms, current_geom,mol)
    grad_norm = np.asarray(np.linalg.norm(grad))
    current_geom -= step_size*(grad/grad_norm)
    #real_step = step_size/grad_norm
    #current_geom -= step_size*(grad)
    #current_geom -= real_step*(grad)
    current_geom = un_mass_weight_geom(natoms, current_geom,mol)
    return current_geom

def irc(output_file):
    """
        This function runs the irc procedure
    """
    max_steps = 1000
    params = Params()
    
    mol = psi4.geometry(params.geometry)
    print(mol)
    starting_vec = np.asarray(params.ts_vec)

    steps = 0
    E = 0.0
    previous_E = 0.0
    energies = []
    current_geom = np.asarray(mol.geometry())
    #print(current_geom)
    output = open(output_file, "a")
    output.write('\n\n--Intrinsic Reaction Coordinate (%s)--\n' %(params.direction))
    output.write('\n-------------------------------------------------------------------------------------')
    output.write('\n{:>20} {:>20} {:>20} {:>20}\n'.format('Coordinate', 'E', 'Delta E', 'Gradient Norm'))
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    while (steps <= max_steps):
        mol.save_xyz_file('euler_step_'+str(steps)+'.xyz',False)
        if(steps==0):
            grad = mass_weight(params.natoms, starting_vec, mol)
        else:
            grad, E  = grad_calc(params, current_geom, mol)
    
        current_geom = euler_step(params.natoms, current_geom, grad,params.step_size,mol)
        if(steps > 20):
            if(E>previous_E):
                #print("pyREX: New energy is greater! Likely near a minimum!")
                break
        steps = steps + 1
        del_E = E - previous_E
        previous_E = E
        if(params.direction=="backward"):
            coord = -1*steps*params.step_size
        else:
            coord = steps*params.step_size
        print_step(output_file,coord, E, del_E, grad)
    output = open(output_file, "a")
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    with open('irc_%s.xyz' %params.direction,'w') as outfile:
        for i in range(steps):
            with open('euler_step_'+str(i)+'.xyz')as infile:
                outfile.write(infile.read())
        os.system('rm euler_step*')

def print_step(output_file, coord, energy, del_E, grad):
    output = open(output_file, "a")
    grad_norm = np.linalg.norm(grad)
    output.write('\n{:>20.4f} {:>20.7f} {:>20.10f} {:>20.10f}\n'.format(coord, energy, del_E, grad_norm)) 
    output.close()  
