"""
    Implementation of Reaction Fragility Spectrum and Related Methods

    References:
        (1) Komorowski et al. Phys. Chem. Chem. Phys, 2016, 18, 32658
        (2) Zaklika et al. J. Phys. Chem. A 2019, 123, 4274-4283
        (3) Ordon et al. J. Phys. Chem. A 2020, 124, 1076-1086

    Equations from these papers are referenced in comments throughout the code
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
import sparrow_interface 

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
                hess_trace = np.trace(H_bond) # Eqn.9 and Eqn.10 Ref.(2)
                hess_trace_bond.append(hess_trace)
                connectivity_matrix[i][j] = hess_trace # Eqn.8 Ref.(2)
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
        fragile = -1.0*np.gradient(divergence, params.irc_stepsize) # Eqn.20 Ref.(2)
        fragilities.append((bond_title,fragile))
    
    # Calculate Atomization energies described in Ref 3
    relative_energies = []
    for i in range(params.natoms):
        for j in range(params.natoms):
            bond_title = "%s%d-%s%d" %(params.symbols[i],i, params.symbols[j],j)
            atomization_energies = []
            for k in range(len(geoms)):
                atomization_energy = 0
                c_matrix_trace = np.trace(c_matrices[k])
                if(i==j):
                    atomization_energy = c_matrices[k][i][j]/c_matrix_trace # Eqn.17 Ref.(3)
                else:
                    atomization_energy = -2.0*c_matrices[k][i][j]/c_matrix_trace # Eqn.18 Ref.(3)
                atomization_energies.append(atomization_energy)
            relative_energies.append((bond_title,atomization_energies))

    relative_energy_variations = []
    for i in range(len(relative_energies)):
        bond_title = relative_energies[i][0]
        relative_energy = relative_energies[i][1]
        variation = np.gradient(relative_energy, params.irc_stepsize)
        relative_energy_variations.append((bond_title,variation))       
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

    csv_output = open("rel_energy_variations.csv", "w")
    csv_output.write("Coordinate,")
    for i in range(len(relative_energy_variations)):
        csv_output.write("%s," %(relative_energy_variations[i][0]))
    csv_output.write("\n")
    for i in range(len(params.geoms)):
        csv_output.write("%.4f," %params.coordinates[i])
        for j in range(len(relative_energy_variations)):
            csv_output.write("%.8f," %relative_energy_variations[j][1][i])  
        csv_output.write("\n")    
    # for i in range(params.natoms):
    #     csv_output.write("%s%d," %(params.symbols[i], i))
    # csv_output.write("\n")
    # for i in range(len(params.geoms)):
    #     csv_output.write("%.4f," %params.coordinates[i])
    #     for j in range(params.natoms):
    #         csv_output.write("%.8f," %fragilities[j][i])
    #     csv_output.write("\n")
