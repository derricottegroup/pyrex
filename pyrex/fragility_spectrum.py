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


#class Params():
#    def __init__(self):
#        """
#            Initialize .json file provided by user in command line, read input and store variables.
#        """
#        json_input = sys.argv[1]
#        self.read_input(json_input)
#        self.geometry_builder()
#        self.oldgrad = 0.0
#        self.oldcoord = 0.0
#    def read_input(self, json_input):
#        json_data=open(json_input).read()
#        input_params = json.loads(json_data)
#        #print(input_params['molecule']['molecular_charge'])
#        #TODO Add explanation of read in variables similar to euler.py
#        if 'symbols' in input_params['molecule']:
#            self.symbols = input_params['molecule']['symbols']
#            self.natoms = len(input_params['molecule']['symbols'])
#        if 'molecular_charge' in input_params['molecule']:
#            self.molecular_charge = int(input_params['molecule']['molecular_charge'])
#        if 'molecular_multiplicity' in input_params['molecule']:
#            self.molecular_multiplicity = input_params['molecule']['molecular_multiplicity']
#        if 'basis' in input_params['model']:
#            self.basis = input_params['model']['basis']
#        if 'method' in input_params['model']:
#            self.method = input_params['model']['method']
#        if 'keywords' in input_params:
#            self.keywords = input_params['keywords']
#        if 'irc_filename' in input_params['pyrex']:
#            self.irc_filename = input_params['pyrex']['irc_filename']
#        if 'irc_stepsize' in input_params['pyrex']:
#            self.irc_stepsize = input_params['pyrex']['irc_stepsize']
#
#    def geometry_builder(self):
#        full_irc = open(self.irc_filename, "r")
#        natoms = self.natoms
#        irc = []
#        geometries = []
#        coordinates = []
#        for line in full_irc:
#            if "Full IRC Point" in line:
#                geom = []
#                irc_num_line = line.split()
#                irc_num = int(irc_num_line[3])
#                for i in range(natoms):
#                    line = next(full_irc)
#                    geom.append(line.lstrip())
#                irc.append((irc_num, geom))
#                geometries.append(geom)
#                coordinates.append(irc_num*self.irc_stepsize)
#        geomparser = Geomparser(natoms, self.molecular_charge, self.molecular_multiplicity, geometries, coordinates)
#        self.geoms = geomparser.geombuilder()   
#        #self.irc = irc_filename

def compute_hessian(params, geom):
    # Use PSI4 to calculate Hessian matrix
    psi4.core.set_output_file("hessian.out", False)
    mol = psi4.geometry(geom)
    psi4.set_options(params.keywords)
    H = psi4.hessian("%s/%s" %(params.method,params.basis))
    Hess = np.array(H)
    return Hess

def vectorize(coords, natoms):
    # Function to transform gradient and coordinate matrices into vectors
    increment = 0
    vectorized_coords = np.zeros((natoms*3, 1))
    for i in range(natoms):
        for j in range(3):
            vectorized_coords[j+increment][0] = coords[i][j]
        increment += 3
    return vectorized_coords

def murtagh_sargent_hess_update(hess, del_g, del_x):
    hess_disp = del_g - np.dot(H,del_x) #Eq 15 J. Chem. Theory Comput. 2005, 1, 61-69
    del_hess = (np.dot(hess_disp,hess_disp.T))/(np.dot(hess_disp.T,del_x))
    return del_hess

def powell_symmetric_broyden_hess_update(hess, del_g, del_x):
    hess_disp = del_g - np.dot(H,del_x) #Eq 16 J. Chem. Theory Comput. 2005, 1, 61-69
    trans = np.dot(del_x, del_x.T)
    first_term = (np.dot(hess_disp,del_x.T) + np.dot(del_x,hess_disp.T))/(np.dot(del_x.T,del_x))
    second_term = (np.dot(del_x.T,hess_disp)*trans)/(np.dot(del_x.T,del_x)*np.dot(del_x.T,del_x))
    del_hess = np.subtract(first_term, second_term)
    return del_hess

def bofill_hess_update(hess, del_g, del_x):
    hess_disp = del_g - np.dot(H,del_x) #Eqs 17 and 18 J. Chem. Theory Comput. 2005, 1, 61-69
    print(hess_disp)
    coeff = np.square(np.dot(del_x.T,hess_disp))/(np.square(del_x)*np.square(hess_disp))
    h_ms = murtagh_sargent_hess_update(hess, del_g, del_x)
    h_psb = powell_symmetric_broyden_hess_update(hess, del_g, del_x)
    hess_disp = coeff*h_ms + (1.0 - coeff)*h_psb
    return hess_disp

def fragility_spec(params,geoms,output_file):
    #geoms = params.geometries
    #print(geoms) 
    output = open(output_file, "a")
    output.write("\nReaction Fragility Spectrum")
    output.write("\n-----------------------------------")
    output.close()
    count = 0
    hess_update_track = 0
    hess_traces = []
    grad_method = "%s/%s" %(params.method,params.basis)
    c_matrices = []
    for geom in geoms:
        output = open(output_file, "a")
        output.write("\nCalculating Hessian Trace for IRC Point %.3f\n" %params.coordinates[count])
        output.close()
       # if(hess_update_track==0 or params.coordinates[count]==0.0):
        #print("prior to Hessian calc")
        H = compute_hessian(params,geom)
        #print(H)
        #mol = psi4.geometry(geom)
        #params.oldgrad = np.asarray(psi4.gradient(grad_method))
        #params.oldcoords = np.asarray(mol.geometry())
       # else:
       #     mol = psi4.geometry(geom)
       #     grad = np.asarray(psi4.gradient(grad_method))
       #     coords = np.asarray(mol.geometry())
       #     print(params.oldcoords)
       #     print(coords)
       #     print(coords-params.oldcoords)
       #     #eigvals, del_x = np.linalg.eig(H)
       #     del_x = np.resize((coords - params.oldcoords),3*params.natoms)
       #     print(del_x)
       #     del_g = np.reshape((grad - params.oldgrad),3*params.natoms)
       #     print(del_g)
       #     hess_diff = powell_symmetric_broyden_hess_update(H, del_g, del_x)
       #     params.oldcoords = coords
       #     params.oldgrad = grad
       #     H += hess_diff
        #Grab Square blocks of Hessian
        connectivity_matrix = np.zeros((params.natoms, params.natoms)) # Building Blank Connectivity Matrix
        hess_trace_bond = []
        increment_i = 0
        for i in range(params.natoms):
            increment_j = 0
            for j in range(params.natoms):
                output = open(output_file, "a")
                H_bond = np.zeros((3,3))
                #print(increment_i)
                #print(increment_j)
                for k in range(3):
                    for l in range(3):
                        H_bond[k][l] = H[k+increment_i][l+increment_j]
                output.write("\nHessian Trace for %s%d-%s%d Bond\n" %(params.symbols[i],i, params.symbols[j],j))
                hess_trace = np.trace(H_bond)
                hess_trace_bond.append(hess_trace)
                connectivity_matrix[i][j] = hess_trace 
                output.write("%.5f\n" %hess_trace)
                increment_j += 3 
                hess_traces.append(hess_trace_bond)
                #count += 1
                output.close()
            increment_i += 3
        c_matrices.append(connectivity_matrix) #Store connectivity matrix for current geometry
        count +=1    
        #TODO: Continue Working on Connectivity Matrix Code, implement all of the techniques from the paper
        #hess_update_track +=1
        #if(hess_update_track==2):
        #    hess_update_track = 0
    
    atom_ordered_traces = [] # Order such that traces for each atom are in same array
    for i in range(params.natoms):
        atom_traces = []
        for j in range(len(geoms)):
            atom_traces.append(hess_traces[j][i])
        atom_ordered_traces.append(atom_traces)

    # Grab Divergences from Connectivity Matrices
    divergences = [] # Empty list for divergences
    for i in range(params.natoms):
        for j in range(params.natoms):
            divergence_elements = []
            for k in range(len(geoms)):
                divergence_elements.append(c_matrices[k][i][j])
            divergences.append(("%s%d-%s%d" %(params.symbols[i],i, params.symbols[j],j),divergence_elements))


    # Calculate Fragility Spectrum for Each Divergence
    fragilities = [] # Empty list for fragilities
    for i in range(len(divergences)):
        bond_title = divergences[i][0]
        divergence = divergences[i][1]
        fragile = -1.0*np.gradient(divergence, params.irc_stepsize)
        fragilities.append((bond_title,fragile))
            
            
    # Take derivative of trace to get fragility spectrum for each atom
    # fragilities = []
    # for i in range(params.natoms):
    #     fragility = np.gradient(atom_ordered_traces[i],params.irc_stepsize)
    #     fragilities.append(fragility)
    #     #print("Fragility for atom %s%d" %(params.symbols[i], i))
    #     #print(fragility)
    
    csv_output = open("frag_spec.csv", "w")
    csv_output.write("Coordinate,")
    for i in range(len(fragilities)):
        csv_output.write("%s," %(fragilities[i][0]))
    csv_output.write("\n")
    for i in range(len(params.geoms)):
        csv_output.write("%.4f," %params.coordinates[i])
        for j in range(len(fragilities)):
            csv_output.write("%.8f," %fragilities[j][1][i])  
        csv_output.write("\n")     
    # for i in range(params.natoms):
    #     csv_output.write("%s%d," %(params.symbols[i], i))
    # csv_output.write("\n")
    # for i in range(len(params.geoms)):
    #     csv_output.write("%.4f," %params.coordinates[i])
    #     for j in range(params.natoms):
    #         csv_output.write("%.8f," %fragilities[j][i])
    #     csv_output.write("\n")
