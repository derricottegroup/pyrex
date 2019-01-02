#!/usr/bin/env python

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

#####################
# Read in Parameters#
#####################

class Params():
    def __init__(self):
        """
            Initialize .json file provided by user in command line, read input and store variables.
        """
        json_input = sys.argv[1]
        self.read_input(json_input)            
        
    def json2xyz(self, input_params):
        """
            Geometries in the .json format (MOLSSI standard) are given as an array. This functions
            takes the given geometry and converts it to a string with a format similar to xyz files.
        """
        charge = input_params['molecule']['molecular_charge']
        mult = input_params['molecule']['molecular_multiplicity']
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
        if input_params['molecule']['geometry']:
            self.geometry = self.json2xyz(input_params)
        if input_params['model']['basis']:
            self.basis = input_params['model']['basis']
        if input_params['model']['method']:
            self.method = input_params['model']['method']
        if input_params['molecule']['symbols']:
            self.symbols = input_params['molecule']['symbols']
            self.natoms = len(input_params['molecule']['symbols'])
        if input_params['keywords']:
            self.keywords = input_params['keywords']    
        if input_params['irc']['direction']:
            self.direction = input_params['irc']['direction']
        if input_params['irc']['step_size']:
            self.step_size = input_params['irc']['step_size']
        if input_params['irc']['mode']:
            self.mode = input_params['irc']['mode']
        if input_params['irc']['normal_mode_file']:
            self.normal_mode_file = input_params['irc']['normal_mode_file']
            self.ts_vec = self.normal_mode_reader()
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
    """
    mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
    grad_method = "%s/%s" %(params.method,params.basis)
    psi4.core.set_output_file("psi4_out.dat", False)
    psi4.set_options(params.keywords)
    E, wfn = psi4.energy(grad_method,return_wfn=True)
    grad = np.asarray(psi4.gradient(grad_method,ref_wfn=wfn))
    grad_mw = mass_weight(params.natoms, grad)
    return grad_mw, E

def mass_weight(natoms,grad):
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

def mass_weight_geom(natoms,geom):
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

def un_mass_weight_geom(natoms,geom_mw):
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

def euler_step(natoms,current_geom,grad,step_size):
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
    current_geom = mass_weight_geom(natoms, current_geom)
    grad_norm = np.asarray(np.linalg.norm(grad))
    current_geom -= step_size*(grad/grad_norm)
    current_geom = un_mass_weight_geom(natoms, current_geom)
    return current_geom

max_steps = 1000
params = Params()

mol = psi4.geometry(params.geometry)
#grad, E = grad_calc(params)
starting_vec = np.asarray(params.ts_vec)
#print(starting_vec)
steps = 0
E = 0.0
previous_E = 0.0
energies = []
current_geom = np.asarray(mol.geometry())
#print(current_geom)
while (steps <= max_steps):
    mol.save_xyz_file('euler_step_'+str(steps)+'.xyz',False)
    if(steps==0):
        grad = mass_weight(params.natoms, starting_vec)
        #grad = starting_vec
    else:
        grad, E  = grad_calc(params, current_geom, mol)

    current_geom = euler_step(params.natoms, current_geom, grad,params.step_size)
    if(steps > 100):
        if(E>previous_E):
            print("pyREX: New energy is greater! Likely near a minimum!")
            break
    steps = steps + 1
    previous_E = E
with open('irc_%s.xyz' %params.direction,'w') as outfile:
    for i in range(steps):
        with open('euler_step_'+str(i)+'.xyz')as infile:
            outfile.write(infile.read())
    os.system('rm euler_step*')
