import sys
import os
import json




class Params(object):
    def __init__(self):
        """
            Initialize .json file provided by user in command line, read input and store variables.
        """
        self.do_irc = False
        self.do_solvent = False
        self.do_energy = False
        self.do_frag = False
        self.do_sapt = False
        self.do_atomic = False
        self.eps = 80.4 # Water by default
        self.do_polarization = False
        self.do_eda = False
        self.do_conceptualdft = False
        self.do_fragility_spec = False
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
        self.qm_program = "psi4"
        json_input = sys.argv[1]
        self.read_input(json_input)
        # Load Output file
        output_filename = "pyrex_output.dat"
        json_data=open(json_input).read()
        header(output_filename, json_data)
    def read_input(self, json_input):
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        if 'molecule' in input_params:
            if 'molecular_charge' in input_params['molecule']:
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
                self.do_frag = True
            if 'fragment_charges' in input_params["molecule"]:
                self.charge_A = input_params["molecule"]["fragment_charges"][0]
                self.charge_B = input_params["molecule"]["fragment_charges"][1]
            if 'fragment_multiplicities' in input_params["molecule"]:
                self.mult_A = input_params["molecule"]["fragment_multiplicities"][0]
                self.mult_B = input_params["molecule"]["fragment_multiplicities"][1]
            if 'geometry' in input_params['molecule']:
                self.geometry = self.json2xyz(input_params)
        if 'model' in input_params:
            if 'basis' in input_params['model']:
                self.basis = input_params['model']['basis']
            if 'method' in input_params['model']:
                self.method = input_params['model']['method']
        if 'keywords' in input_params:
            self.keywords = input_params['keywords']
        if 'irc' in input_params:
            self.do_irc = True
        if 'pyrex' in input_params:
            if 'qm_program' in input_params['pyrex']:
                self.qm_program = input_params['pyrex']['qm_program']
            if 'do_energy' in input_params['pyrex']:
                self.do_energy = bool(input_params['pyrex']['do_energy'])
            if 'xc_functional' in input_params['pyrex']:
                self.xc_functional = input_params['pyrex']['xc_functional']
            if 'nthreads' in input_params['pyrex']:
                self.nthreads = input_params['pyrex']['nthreads']
            if 'do_solvent' in input_params['pyrex']:
                self.do_solvent = input_params['pyrex']['do_solvent']
            if 'eps' in input_params['pyrex']:
                self.eps = input_params['pyrex']['eps']
            if 'do_conceptualdft' in input_params['pyrex']:
                self.do_conceptualdft = bool(input_params['pyrex']['do_conceptualdft'])
            if 'do_fragility_spec' in input_params['pyrex']:
                self.do_fragility_spec = bool(input_params['pyrex']['do_fragility_spec'])
            if 'do_supermolecular' in input_params['pyrex']:
                self.do_supermolecular = bool(input_params['pyrex']['do_supermolecular'])
            if 'energy_read' in input_params['pyrex']:
                self.energy_file = input_params['pyrex']['energy_read']
            if 'sapt_read' in input_params['pyrex']:
                self.sapt_file = input_params['pyrex']['sapt_read']
            if 'force_max' in input_params['pyrex']:
                self.force_max = input_params['pyrex']['force_max']
            if 'force_min' in input_params['pyrex']:
                self.force_min = input_params['pyrex']['force_min']
            if 'restart' in input_params['pyrex']:
                self.restart = bool(input_params['pyrex']['restart']) #TODO Implement this functionality
            if 'do_sapt' in input_params['pyrex']:
                self.do_sapt = bool(input_params['pyrex']['do_sapt'])
            if 'do_atomic' in input_params['pyrex']:
                self.do_atomic = bool(input_params['pyrex']['do_atomic'])
            if 'do_polarization' in input_params['pyrex']:
                self.do_polarization = bool(input_params['pyrex']['do_polarization'])
            if 'sapt_method' in input_params["pyrex"]:
                self.sapt_method = input_params["pyrex"]["sapt_method"]
            if 'irc_stepsize' in input_params['pyrex']:
                self.irc_stepsize = input_params['pyrex']['irc_stepsize']
            if 'irc_filename' in input_params['pyrex']:
                self.irc_filename = input_params['pyrex']['irc_filename']
                self.irc_grab()
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
        if 'fsapt' in input_params:
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
        self.irc = irc
        self.geometries = geometries
        self.coordinates = coordinates
