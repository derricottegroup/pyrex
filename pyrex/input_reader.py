import sys
import os
import json
import numpy as np
from .header import *
from pydantic import BaseModel, Field
from typing import Dict, List, Optional 

class Params(object):
    def __init__(self,json_data):
        """
            Initialize .json file provided by user in command line, read input and store variables.
        """
        self.do_irc = False
        self.do_solvent = False
        self.pcm_solvent = "Water"
        self.do_energy = False
        self.active_site = False
        self.do_frag = False
        self.do_sapt = False
        self.do_fsapt = False
        self.fsapt_analyze = False
        self.do_atomic = False
        self.eps = 80.4 # Water by default
        self.do_polarization = False
        self.do_eda = False
        self.set_memory = False
        self.memory_allocation =  ""
        self.do_conceptualdft = False
        self.do_fragility_spec = False
        self.single_connectivity_matrix = False
        self.do_supermolecular = False
        self.do_surf_scan = False
        self.scan_type = 'relaxed'
        self.surf_scan_mode = False
        self.constraint_type = None
        self.constrained_atoms = None
        self.constrained_values = None
        self.do_rexplot = False
        self.energy_file = None
        self.sapt_file = None
        self.orca_header = ""
        self.orca_block = ""
        self.qm_program = "psi4"
        self.grace_period = 50
        self.e_conv = 1e-5
        self.nthreads = 1
        #if(sys.argv[0]!="pytest"):
        #    json_input = sys.argv[1]
        #else:
        #json_input = "input.json" 
        self.read_input(json_data)
        # Load Output file
        output_filename = "pyrex_output.dat"
        #json_data=open(json_input).read()
        header(output_filename, json_data)
        self.outfile = output_filename
    def read_input(self, json_data):
        #json_data=open(json_input).read()
        input_params = json.loads(json_data)
        """
            Molecule Block: Defines the properties of the molecular system.
           
            Paramters:
            ---------
                symbols(list) -- list containing the symbol for each atom in order.
                molecular_charge(str) -- charge of the total molecular system
                molecular_multiplicity(int) -- multiplicity of the total molecular system
                fragments(list of lists) -- list of lists containing the integers defining fragments
                                            of the system. Currently only accepts two fragments.
                fragment_charges(list) -- defining the charge on each fragment. 
                fragment_multiplicities -- defining the multiplicity of each fragment.
                geometry(list) -- defines the xyz coordinates for the molecule (numbers only!) 
        """
        if 'molecule' in input_params:
            if 'symbols' in input_params['molecule']:
                self.symbols = input_params['molecule']['symbols']
                self.natoms = len(input_params['molecule']['symbols'])
            if 'molecular_charge' in input_params["molecule"]:
                self.molecular_charge = int(input_params["molecule"]["molecular_charge"])
            if 'molecular_multiplicity' in input_params["molecule"]:
                self.molecular_multiplicity = input_params["molecule"]["molecular_multiplicity"]
            if 'fragments' in input_params["molecule"]:
                self.frag_A = input_params["molecule"]["fragments"][0]
                self.natoms_A = len(self.frag_A)
                self.frag_B = input_params["molecule"]["fragments"][1]
                self.natoms_B = len(self.frag_B)
                self.fraglist = input_params["molecule"]["fragments"]
                self.do_frag = True
            if 'fragment_charges' in input_params["molecule"]:
                self.charge_A = input_params["molecule"]["fragment_charges"][0]
                self.charge_B = input_params["molecule"]["fragment_charges"][1]
                self.frag_charge = input_params["molecule"]["fragment_charges"]
            if 'fragment_multiplicities' in input_params["molecule"]:
                self.mult_A = input_params["molecule"]["fragment_multiplicities"][0]
                self.mult_B = input_params["molecule"]["fragment_multiplicities"][1]
                self.frag_mult = input_params["molecule"]["fragment_multiplicities"]
            if 'geometry' in input_params['molecule']:
                self.geometry = self.json2xyz(input_params)
        """
            Model Block: Defines the level of theory for the calculations. 

            Paramters:
            ---------
                basis(str) -- sets the basis set for the calculation
                method(str) -- specify the method you want to use (scf, mp2, dft, etc.). 
                               CAUTION: If using DFT with PySCF, you must say "dft" here
                               and specify the xc functional in the Pyrex block. 
        """
        if 'model' in input_params:
            if 'basis' in input_params['model']:
                self.basis = input_params['model']['basis']
            if 'method' in input_params['model']:
                self.method = input_params['model']['method']
        """
            Keywords Block: Only used for PSI4, any valid PSI4 keyword for global options can 
                            be placed here and used for all modules.  
        """
        if 'keywords' in input_params:
            self.keywords = input_params['keywords']
        else:
            self.keywords = {} 
        """
            IRC Block: Options for intrinsic reaction coordinate calculations. 

            Paramters:
            ---------
                direction(str) -- specify if the irc will go in the forward(positive gradient) or 
                                  backward(negative gradient) direction. 
                step_size(float) -- specify the irc step size in units au amu^(1/2). CAUTION: the 
                                    IRC algorithm currently implemented is very unstable at any 
                                    step size higher than 0.01 au amu^(1/2). 
                mode(int) -- specify which vibrational mode in your normal mode file you will be following. 
                normal_mode_file(str) -- the name of the file containing the normal mode you want to follow. 
                grace_period(int) -- specify how many iterations may pass before the energy check 
                                     ends the IRC due to small energy change or energy increase. This
                                     can be helpful for really flat PESes. 
                e_conv(float) -- specify the energy convergence threshold. 
        """
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
            self.do_irc = True
        """
            Pyrex Block: Options relevant for Pyrex related tasks.

            Paramters:
            ---------
                qm_program(str) -- name of the quantum chemistry program you want to use (PSI4 or PYSCF)
                do_energy(bool) -- do you want to calculate the energy?
                xc_functional(str) -- specify an exchange-correlation functional for DFT calculations. 
                nthreads(int) -- specify the number of threads to use for your calculation. 
                do_solvent(bool) -- do you want to use an implicit solvent model?
                eps(float) -- sets the vacuum permitivitty constant for implicit solvation. 
                do_conceptualdft(bool) -- do you want to calculate properties from conceptual DFT?
                                          this includes chemical potential, reaction electronic flux,
                                          chemical hardness/softness, etc. 
                do_fragility_spec(bool) -- do you want to calculate atomic fragility spectrum. 
                                           CAUTION: this requires calculation of Hessian at 
                                           every point and can be expensive even for relatively 
                                           small systems.  
                do_supermolecular(bool) -- do you want to calculate the interaction energy using 
                                           a supermolecular approach?
                energy_read(str) -- name of file that contains pre-computed reaction energies 
                force_max(float) -- value of the coordinate (in units au amu^(1/2)) that 
                                    corresponds to the force maximum structure. 
                force_min(float) -- value of the coordinate (in units au amu^(1/2)) that
                                    corresponds to the force minimum structure. 
                do_sapt(bool) -- do you want to calculate interaction energies using
                                 symmetry-adapted perturbation theory?
                sapt_read(str) -- name of file that contains pre-computed SAPT interaction energies
                sapt_method(str) -- specify the level of SAPT you would like to use. 
                do_atomic(bool) -- do you want to calculate atomic force contributions?
                do_polarization(bool) -- do you want to decompose reaction electronic flux
                                         into its polarization and transfer components? 
                irc_stepsize(float) -- specify the stepsize that the IRC has. 
                irc_filename(str) -- name of the file that contains the IRC geometries. 
        """
        if 'pyrex' in input_params:
            if 'qm_program' in input_params['pyrex']:
                self.qm_program = input_params['pyrex']['qm_program']
            if 'do_energy' in input_params['pyrex']:
                self.do_energy = bool(input_params['pyrex']['do_energy'])
            if 'xc_functional' in input_params['pyrex']:
                self.xc_functional = input_params['pyrex']['xc_functional']
            if 'nthreads' in input_params['pyrex']:
                self.nthreads = input_params['pyrex']['nthreads']
            if 'set_memory' in input_params['pyrex']:
                self.set_memory = True
                self.memory_allocation = str(input_params['pyrex']['set_memory'])
            if 'do_solvent' in input_params['pyrex']:
                self.do_solvent = bool(input_params['pyrex']['do_solvent'])
            if 'pcm_solvent' in input_params['pyrex']:
                self.pcm_solvent = str(input_params['pyrex']['pcm_solvent'])
            if 'eps' in input_params['pyrex']:
                self.eps = input_params['pyrex']['eps']
            if 'do_conceptualdft' in input_params['pyrex']:
                self.do_conceptualdft = bool(input_params['pyrex']['do_conceptualdft'])
            if 'do_fragility_spec' in input_params['pyrex']:
                self.do_fragility_spec = bool(input_params['pyrex']['do_fragility_spec'])
            if 'single_connectivity_matrix' in input_params['pyrex']:
                self.single_connectivity_matrix = bool(input_params['pyrex']['single_connectivity_matrix'])
            if 'do_supermolecular' in input_params['pyrex']:
                self.do_supermolecular = bool(input_params['pyrex']['do_supermolecular'])
            if 'energy_read' in input_params['pyrex']:
                self.energy_file = input_params['pyrex']['energy_read']
            if 'force_max' in input_params['pyrex']:
                self.force_max = input_params['pyrex']['force_max']
            if 'force_min' in input_params['pyrex']:
                self.force_min = input_params['pyrex']['force_min']
            if 'restart' in input_params['pyrex']:
                self.restart = bool(input_params['pyrex']['restart']) #TODO Implement this functionality
            if 'do_sapt' in input_params['pyrex']:
                self.do_sapt = bool(input_params['pyrex']['do_sapt'])
            if 'sapt_read' in input_params['pyrex']:
                self.sapt_file = input_params['pyrex']['sapt_read']
            if 'sapt_method' in input_params["pyrex"]:
                self.sapt_method = input_params["pyrex"]["sapt_method"]
            if 'do_atomic' in input_params['pyrex']:
                self.do_atomic = bool(input_params['pyrex']['do_atomic'])
            if 'do_polarization' in input_params['pyrex']:
                self.do_polarization = bool(input_params['pyrex']['do_polarization'])
            if 'irc_stepsize' in input_params['pyrex']:
                self.irc_stepsize = input_params['pyrex']['irc_stepsize']
            if 'irc_filename' in input_params['pyrex']:
                self.irc_filename = input_params['pyrex']['irc_filename']
                self.irc_grab()
            if 'orca_header' in input_params['pyrex']:
                self.orca_header = input_params['pyrex']['orca_header']
            if 'orca_block' in input_params['pyrex']:
                self.orca_block = input_params['pyrex']['orca_block']
            if 'active_site' in input_params['pyrex']:
                self.active_site_indices = input_params['pyrex']['active_site']
                self.active_site = True
        """
            Surf Scan Block: Options relevant for surface scans in Pyrex.
                             (Currently only compatible with PSI4)

            Paramters:
            ---------
                scan_type(str) -- relaxed or unrelaxed surface scan?
                constraint_type(str) -- bond, angle, or dihedral constraints? 
                constrained_atoms(list) -- atoms involved in defining constrained bond, angle, etc. 
                constrained_values(list) -- range of values you wish to scan from. 
        """
        if 'surf_scan' in input_params:
            self.do_surf_scan = True
            if 'scan_type' in input_params['surf_scan']:
                self.scan_type = input_params['surf_scan']['scan_type']
            if 'constraint_type' in input_params['surf_scan']:
                self.constraint_type = input_params['surf_scan']['constraint_type']
            if 'constrained_atoms' in input_params['surf_scan']:
                atoms_array = input_params['surf_scan']['constrained_atoms']
                atoms_string = ""
                for atom in atoms_array:
                    atoms_string = atoms_string + "%s " %atom
                self.constrained_atoms = atoms_string
            if 'constrained_values' in input_params['surf_scan']:
                const_inp = input_params['surf_scan']['constrained_values']
                if(const_inp[0]>const_inp[1]):
                    tmp_array = np.arange(const_inp[1],const_inp[0],const_inp[2])
                    tmp_array = np.flip(tmp_array,0)
                else:
                    tmp_array = np.arange(const_inp[0],const_inp[1],const_inp[2])
                self.constrained_values = tmp_array
        #TODO Finish adding documentation. Start adding tests/samples/examples
        if 'fsapt' in input_params:
            self.do_fsapt = True
            self.monomer_A_frags = input_params['fsapt']['monomer_A_frags']
            self.monomer_B_frags = input_params['fsapt']['monomer_B_frags']
            self.monomer_A_labels = input_params['fsapt']['monomer_A_labels']
            self.monomer_B_labels = input_params['fsapt']['monomer_B_labels']
        if 'rexplot' in input_params:
            self.do_rexplot = True
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

    def irc_grab(self):
        irc = []
        geometries = []
        coordinates = []
        full_irc = open(self.irc_filename, "r")
        self.surf_scan_mode = False
        # Grab and store geometries from the IRC
        for line in full_irc:
            if("Full IRC Point" in line): #TODO Make this compatible with surface scan xyz files
                geom = []
                irc_num_line = line.split()
                irc_num = int(irc_num_line[3])
                for i in range(self.natoms):
                    line = next(full_irc)
                    geom.append(line.lstrip())
                irc.append((irc_num, geom))
                geometries.append(geom)
                coordinate = irc_num*self.irc_stepsize
                coordinates.append(round(coordinate,3))
            if("surface scan with" in line):
                self.surf_scan_mode = True
                geom = []
                irc_num_line = line.split()
                irc_num = float(irc_num_line[7])
                for i in range(self.natoms):
                    line = next(full_irc)
                    geom.append(line.lstrip())
                irc.append((irc_num, geom))
                geometries.append(geom)
                coordinate = irc_num
                coordinates.append(round(coordinate,3))
        full_irc.close()
        if(irc==[]):
            with open(self.irc_filename) as f:
                irc_num = 0
                for line in f:
                    geom = []
                    line = next(f)
                    for i in range(self.natoms):
                        line = next(f)
                        geom.append(line.lstrip())
                    irc_num = irc_num + 1
                    irc.append((irc_num, geom))
                    geometries.append(geom)
                    coordinate = irc_num
                    coordinates.append(round(coordinate,3))
        self.irc = irc
        self.geometries = geometries
        self.coordinates = coordinates
    def normal_mode_reader(self):
        """
            Reads the file containing the normal modes of vibration. This function currently
            only works with the Psi4 normal mode output. In order to produce a compatible file
            run a frequency calculation in Psi4 with the option "normal_mode_writer" set True.

            Parameters:
            ----------
                self(self) -- contains all shared parameters.
            Returns:
                ts_vec(array) -- Array containing the normal mode vector
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


