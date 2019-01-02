"""
    Implementation of Reaction Fragility Spectrum
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
from geomparser import *


class Params():
    def __init__(self):
        """
            Initialize .json file provided by user in command line, read input and store variables.
        """
        json_input = sys.argv[1]
        self.read_input(json_input)
        self.geometry_builder()
    def read_input(self, json_input):
        json_data=open(json_input).read()
        input_params = json.loads(json_data)
        #print(input_params['molecule']['molecular_charge'])
        #TODO Add explanation of read in variables similar to euler.py
        if input_params['molecule']['symbols']:
            self.symbols = input_params['molecule']['symbols']
        if input_params['molecule']['molecular_charge']:
            self.molecular_charge = int(input_params['molecule']['molecular_charge'])
        if input_params['molecule']['molecular_multiplicity']:
            self.molecular_multiplicity = input_params['molecule']['molecular_multiplicity']
        if input_params['model']['basis']:
            self.basis = input_params['model']['basis']
        if input_params['model']['method']:
            self.method = input_params['model']['method']
            self.natoms = len(input_params['molecule']['symbols'])
        if input_params['keywords']:
            self.keywords = input_params['keywords']
        if input_params['pyrex']['irc_filename']:
            self.irc_filename = input_params['pyrex']['irc_filename']
        if input_params['pyrex']['irc_stepsize']:
            self.irc_stepsize = input_params['pyrex']['irc_stepsize']

    def geometry_builder(self):
        full_irc = open(self.irc_filename, "r")
        natoms = self.natoms
        irc = []
        geometries = []
        coordinates = []
        for line in full_irc:
            if "Full IRC Point" in line:
                geom = []
                irc_num_line = line.split()
                irc_num = int(irc_num_line[3])
                for i in range(natoms):
                    line = next(full_irc)
                    geom.append(line.lstrip())
                irc.append((irc_num, geom))
                geometries.append(geom)
                coordinates.append(irc_num*self.irc_stepsize)
        geomparser = Geomparser(natoms, self.molecular_charge, self.molecular_multiplicity, geometries, coordinates)
        self.geoms = geomparser.geombuilder()   
        self.irc = irc
        self.coordinates = coordinates


params = Params()
geoms = params.geoms
def compute_hessian(params, geom):
    # Use PSI4 to calculate Hessian matrix
    psi4.core.set_output_file("hessian.out", False)
    mol = psi4.geometry(geom)
    #print(self.keywords)
    psi4.set_options(params.keywords)
    H = psi4.hessian(params.method)
    Hess = np.array(H)
    return Hess

csv_output = open("fragility_spectrum.csv", "w")
count = 0
for geom in geoms:
    print("Calculating fragility for IRC Point %.3f" %params.coordinates[count])
    H = compute_hessian(params,geom)
    #Grab Square blocks of Hessian
    increment = 0
    hess_traces = []
    for z in range(params.natoms):
        H_atom = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
              H_atom[i][j] = H[i+increment][j+increment]
        print("Reaction Fragility for %s atom" %params.symbols[z])
        hess_trace = np.trace(H_atom)
        print(hess_trace)
        increment += 3 
    count += 1 

#TODO This is almost done. Just have to take the numerical derivatives of the traces in order to
#     get the fragility spectrum. This process is very slow, wonder if Hessian updating can 
#     increase efficiency? (Potential Chem. Phys. Lett./IJQC paper idea!)
