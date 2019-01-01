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
        json_input = sys.argv[1]
        self.read_input(json_input)            
        
    def json2xyz(self, input_params):
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
        # Read in input from JSON input file
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        #print("DATA READ")
        if input_params['molecule']['geometry']:
            self.geometry = self.json2xyz(input_params)
        if input_params['model']['basis']:
            self.basis = input_params['model']['basis']
            #self.level_of_theory = "%s/%s" %(self.method,self.basis)
        if input_params['model']['method']:
            self.method = input_params['model']['method']
            #self.level_of_theory = "%s/%s" %(self.method,self.basis)
        if input_params['molecule']['symbols']:
            self.symbols = input_params['molecule']['symbols']
            self.natoms = len(input_params['molecule']['symbols'])    
        if input_params['dvv']['cons_vel']:
            self.cons_vel = input_params['dvv']['cons_vel']
        if input_params['dvv']['err_tol']:
            self.err_tol = input_params['dvv']['err_tol']
        if input_params['dvv']['e_conv']:
            self.e_conv = input_params['dvv']['e_conv']
        if input_params['dvv']['direction']:
            self.direction = input_params['dvv']['direction']
        if input_params['dvv']['mode']:
            self.mode = input_params['dvv']['mode']
        if input_params['dvv']['mode_freq']:
            self.mode_freq = input_params['dvv']['mode_freq']
        if input_params['dvv']['normal_mode_file']:
            self.normal_mode_file = input_params['dvv']['normal_mode_file']
            self.ts_vec = self.normal_mode_reader()
    def normal_mode_reader(self):
        natoms = self.natoms
        direction = self.direction
        normal_mode_file = open(self.normal_mode_file, 'r')
        ts_vec = []
        for line in normal_mode_file:
            if("vibration %d" %self.mode in line):
                #print(line)
                line = next(normal_mode_file)
                #print(line)
                for i in range(natoms):
                    trj = line.split()
                    #print(trj)
                    trj = list(map(float, trj))
                    np_trj = np.array(trj)
                    #print(np_trj)
                    if(direction=='backward'):
                        ts_vec.append(-1.0*np_trj)
                    else:
                        ts_vec.append(np_trj)
                    line = next(normal_mode_file)
        #print(ts_vec)
        return ts_vec

########################
## Gradient Functions ##
########################

def grad_calc(params,current_geom, mol):
    mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
    #print(np.asarray(mol.geometry))
    grad_method = "%s/%s" %(params.method,params.basis)
    psi4.core.set_output_file("psi4_out.dat", False)
    E, wfn = psi4.energy(grad_method,return_wfn=True)
    grad = np.asarray(psi4.gradient(grad_method,ref_wfn=wfn))
    grad_mw = mass_weight(params, grad)
    return grad_mw, E

def mass_weight(params,grad):
    #mol = psi4.geometry(params.geometry)
    grad_mw = np.zeros((params.natoms, 3))
    for i in range(params.natoms):
        for j in range(3):
            grad_mw[i][j] = grad[i][j]*np.sqrt(mol.mass(i))
    return grad_mw

def euler_step(current_geom,grad):
    grad_norm = np.asarray(np.linalg.norm(grad))
    #print(grad_norm)
    #print(grad)
    #ratio = grad/grad_norm
    current_geom -= 0.01*(grad/grad_norm)
    return current_geom

max_steps = 480
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
        grad = mass_weight(params, starting_vec)
    else:
        grad, E  = grad_calc(params, current_geom, mol)

    current_geom = euler_step(current_geom, grad)
    if(steps > 25):
        if(E>previous_E):
            print(E)
            print("New energy is greater!")
            break
    print(current_geom)
    steps = steps + 1
    print(E)
    previous_E = E
with open('irc_%s.xyz' %params.direction,'w') as outfile:
    for i in range(steps):
        with open('euler_step_'+str(i)+'.xyz')as infile:
            outfile.write(infile.read())
    os.system('rm euler_step*')
