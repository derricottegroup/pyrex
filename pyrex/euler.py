#!/usr/bin/env python

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
        if input_params['molecule']['geometry']:
            self.geometry = self.json2xyz(input_params)
        if input_params['model']['basis']:
            self.basis = input_params['model']['basis']
        if input_params['model']['method']:
            self.method = input_params['model']['method']
        if input_params['molecule']['symbols']:
            self.symbols = input_params['molecule']['symbols']
            self.natoms = len(input_params['molecule']['symbols'])    
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
    mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
    grad_method = "%s/%s" %(params.method,params.basis)
    psi4.core.set_output_file("psi4_out.dat", False)
    E, wfn = psi4.energy(grad_method,return_wfn=True)
    grad = np.asarray(psi4.gradient(grad_method,ref_wfn=wfn))
    grad_mw = mass_weight(params.natoms, grad)
    return grad_mw, E

def mass_weight(natoms,grad):
    grad_mw = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            grad_mw[i][j] = grad[i][j]/np.sqrt(mol.mass(i))
    return grad_mw

def mass_weight_geom(natoms,geom):
    coord_mw = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            coord_mw[i][j] = geom[i][j]*np.sqrt(mol.mass(i))
    return coord_mw

def un_mass_weight_geom(natoms,geom_mw):
    coord = np.zeros((natoms, 3))
    for i in range(natoms):
        for j in range(3):
            coord[i][j] = geom_mw[i][j]/np.sqrt(mol.mass(i))
    return coord

def euler_step(natoms,current_geom,grad,step_size):
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
    if(steps > 25):
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
