#!/usr/bin/env python
"""
Main pyREX interface
"""
__authors__ = "Wallace D. Derricotte"
__credits__ = ["Wallace D. Derricotte"]

__copyright__ = "(c) 2018, Derricotte Research Group"
__license__ = "BSD-3-Clause"
__date__ = "2018-5-29"
import psi4
import numpy as np
import os
from header import *
from input_parser import *
import calctools
import random
import re
import datetime
import json
import euler
import input_reader
import fragility_spectrum
import quotes
import csvread
import pandas as pd
from decimal import Decimal
from concept_dft import *
from scf_class import *
import surface_scan
from geomparser import *
from atomic import atomic_decomp
from sapt_class import *
from prettytable import PrettyTable

import sys

if(os.path.isdir('psi4_output')):
    pass
else:
    os.makedirs("psi4_output")

# Load Output file
output_filename = "pyrex_output.dat"
#header(output_filename, json_data)

# Load User Parameters
params = Params()
#level_of_theory = "%s/%s" %(params.method,params.basis) # Level of Theory for Total Energies

#########
## IRC ##
#########

if(params.do_irc):
    euler.ishida_morokuma(output_filename)

##################
## Surface Scan ##
##################

if(params.do_surf_scan):
    surface_scan.surf_psi4(params, output_filename)

#############################
# Initialize Common Classes #
#############################
if(params.do_energy or params.do_sapt):
    geomparser = Geomparser(params.natoms, params.molecular_charge, params.molecular_multiplicity, params.geometries, params.coordinates)

    scf_instance = scf_class(params, output_filename)

####################################################
# Build Geometries and Print Geometric Information #
####################################################
if(params.do_energy or params.do_sapt):
    geoms = geomparser.geombuilder()
    geomparser.atomic_distances()

if(params.do_energy and params.qm_program == 'pyscf'):
    pyscf_geoms = geomparser.pyscf_geombuilder()
    #print(pyscf_geoms)

##########################
# Fragment Optimizations #
##########################
if(params.do_frag==True):
    iso_frag_A = geomparser.iso_frag(params.charge_A, params.mult_A, params.frag_A)
    iso_frag_B = geomparser.iso_frag(params.charge_B, params.mult_B, params.frag_B)
    e_A = scf_instance.opt("A", params.natoms_A, iso_frag_A[0])
    e_B = scf_instance.opt("B", params.natoms_B, iso_frag_B[0])

###########################
# Run Energy Calculations #
###########################

if(params.do_energy and params.energy_file==None):
    if(params.qm_program == 'psi4'):
        energies, wavefunctions = scf_instance.psi4_scf(geoms)
        nelec = wavefunctions[0].nalpha()

    if(params.qm_program == 'pyscf' and params.method == 'scf'):
        energies,frontier_orb_energies = scf_instance.pyscf_scf(pyscf_geoms)
     
    if(params.qm_program == 'pyscf' and params.method == 'dft'):
        energies,frontier_orb_energies = scf_instance.pyscf_dft(pyscf_geoms, params.xc_functional)

# Store energies in .csv file
    energy_csv = open("energy.csv", "w+")
    energy_csv.write("Coordinate,Energy\n")
    for i in range(len(params.coordinates)):
        energy_csv.write("%f, %f\n" %(params.coordinates[i], energies[i]))
    energy_csv.close()

# Read in energy values if given

if(params.energy_file!=None):
    data = pd.read_csv(params.energy_file)
    energies = data['Energy'].values

###########################
# Reaction Force Analysis #
###########################
if(params.do_energy or params.energy_file!=None):
    coordinates = params.coordinates
    reaction_force_values = -1.0*np.gradient(energies,params.irc_stepsize)
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Force Analysis--\n')
    output.write('\n-------------------------------------------------------------------------------------')
    output.write('\n{:>20} {:>20} {:>20}\n'.format('IRC Point', 'F (Hartree)', 'F (kcal)'))
    output.write('-------------------------------------------------------------------------------------\n')
    
    force_count = 0
    for reaction_force in reaction_force_values:
        output.write('\n{:>20} {:>20.7f} {:>20.7f}\n'.format(params.coordinates[force_count], reaction_force, reaction_force*627.51))
        force_count = force_count+1
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    
    force_csv = open("force.csv", "w+")
    force_csv.write("Coordinate,Force\n")
    for i in range(len(coordinates)):
        force_csv.write("%f, %f\n" %(coordinates[i], reaction_force_values[i]))
    force_csv.close()
    
    index_min = np.argmin(np.asarray(reaction_force_values))
    #print(reaction_force_values[index_min])
    if(params.surf_scan_mode):
        index_ts = np.argmax(np.asarray(energies))
    else:
        index_ts  = coordinates.index(0.0000)
    index_max = np.argmax(np.asarray(reaction_force_values))
    
    r_struct = geoms[0]
    force_min_struct = geoms[index_min]
    ts_struct = geoms[index_ts]
    force_max_struct = geoms[index_max]
    p_struct = geoms[-1]
    
    #key_structs = [r_struct, force_min_struct, ts_struct, force_max_struct, p_struct]
    
    #scf_instance.psi4_molden(key_structs) 
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Partitioning--\n')
    output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], coordinates[index_min]))
    output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(coordinates[index_min], coordinates[index_max]))
    output.write("\nProduct Region:            %.3f ------> %.3f\n" %(coordinates[index_max], coordinates[-1]))

############################
# Conceptual DFT Analysis  #
############################

# Grabbing chemical potential and calculating REF using re_flux class
if(params.do_conceptualdft):
    c_dft = concept_dft()
    if(params.molecular_multiplicity==1):
        if(params.qm_program == 'psi4'):
            potentials = np.array(c_dft.potential(wavefunctions))
        if(params.qm_program == 'pyscf'):
            potentials = np.array(c_dft.potential_pyscf(frontier_orb_energies))
    else:
        potentials = np.array(c_dft.potential_open_shell(wavefunctions))
    if(params.qm_program == 'psi4'):
        hardness = np.array(c_dft.hardness(wavefunctions))
    if(params.qm_program == 'pyscf'):
        hardness = np.array(c_dft.hardness_pyscf(frontier_orb_energies))
    reaction_electronic_flux = np.array(c_dft.electronic_flux(params.irc_stepsize))
    electrophilicity = []
    for i in range(len(hardness)):
        electrophil = (potentials[i]*potentials[i])/(2.0*hardness[i])
        electrophilicity.append(electrophil)

    flux_csv = open("conceptual_dft.csv", "w+")
    flux_csv.write("Coordinate,Chemical Potential,Hardness,Electrophilicity,Reaction Electronic Flux\n")
    for i in range(len(coordinates)):
        flux_csv.write("%f,%f,%f,%f,%f\n" %(coordinates[i], potentials[i], hardness[i], electrophilicity[i], reaction_electronic_flux[i]))
    flux_csv.close()


###########################
# Reaction Work Integrals #
###########################
if(params.do_energy or params.energy_file!=None):
    W_1 = -1.0*np.trapz(reaction_force_values[:index_min+1], dx=params.irc_stepsize)
    W_2 = -1.0*np.trapz(reaction_force_values[index_min:index_ts+1], dx=params.irc_stepsize)
    W_3 = -1.0*np.trapz(reaction_force_values[index_ts:index_max+1], dx=params.irc_stepsize)
    W_4 = -1.0*np.trapz(reaction_force_values[index_max:], dx=params.irc_stepsize)

    output.write('\n\n--Reaction Works--\n')
    output.write('\n-------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Units', 'W_1', 'W_2', 'W_3', 'W_4'))
    output.write('\n-------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}\n'.format('Hartree', W_1, W_2, W_3, W_4))
    output.write('\n{:>15} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}\n'.format('Kcal/mol', W_1*627.51, W_2*627.51, W_3*627.51, W_4*627.51))
    output.write('\n-------------------------------------------------------------------------------------')
    output.close()

#strain_energies = []
#chemical_potentials = []
#chemical_potentials_A = []
#chemical_potentials_B = []
#reaction_electronic_flux = []
#reaction_electronic_flux_A = []
#reaction_electronic_flux_B = []
#polarization_flux = []
#transfer_flux = []
#int_energies = []
#electrostatics = []
#exchange = []
#induction = []
#dispersion = []

if(params.do_sapt==True or params.do_supermolecular==True):
    sapt_geometries = geomparser.sapt_geombuilder(params.charge_A, params.mult_A, params.charge_B, params.mult_B, params.frag_A, params.frag_B)
    #p_sapt_geometries = geomparser.sapt_geombuilder(p_charge_A, p_mult_A, p_charge_B, p_mult_B, product_frag_A, product_frag_B)

#if(params.do_polarization==True):
    #frag_A_geoms = geomparser.frag_ghost(charge_A, mult_A, frag_A_atom_list)
    #frag_B_geoms = geomparser.frag_ghost(charge_B, mult_B, frag_B_atom_list)
    #energies_A, wavefunctions_A = scf_instance.psi4_scf(frag_A_geoms)
    #energies_B, wavefunctions_B = scf_instance.psi4_scf(frag_B_geoms)
    #nelec_A = wavefunctions_A[0].nalpha()
    #nelec_B = wavefunctions_B[0].nalpha()

#######################################
# Atomic Force And Work Decomposition #
#######################################

if(params.do_atomic==True):
    atomic_forces = atomic_decomp(params, output_filename, geoms)

    W_1_conts = [] # Array to store the contributions to Work in Region 1
    for i in range(params.natoms):
        W_1_conts.append(np.sum(atomic_forces[i][:index_min]))
    W_2_conts = [] # Array to store the contributions to Work in Region 2
    for i in range(params.natoms):
        W_2_conts.append(np.sum(atomic_forces[i][index_min:index_ts]))
    W_3_conts = [] # Array to store the contributions to Work in Region 3
    for i in range(params.natoms):
        W_3_conts.append(np.sum(atomic_forces[i][index_ts:index_max]))
    W_4_conts = [] # Array to store the contributions to Work in Region 4
    for i in range(params.natoms):
        W_4_conts.append(np.sum(atomic_forces[i][index_max:]))
    output = open(output_filename, "a")
    output.write('\n\n--Atomic Decomposition of Reaction Works (kcal/mol)--\n')
    output.write('\n-------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Atom', 'W_1', 'W_2', 'W_3', 'W_4'))
    output.write('\n-------------------------------------------------------------------------------------')
    for i in range(params.natoms): 
        output.write('\n{:>15} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}\n'.format('%s%d' %(params.symbols[i],i+1), W_1_conts[i]*627.51, W_2_conts[i]*627.51, W_3_conts[i]*627.51, W_4_conts[i]*627.51))
    output.write('\n-------------------------------------------------------------------------------------')
    output.close()


######################################################
# Symmetry-Adapted Perturbation Theory Decomposition #
######################################################

if(params.energy_file!=None and (params.do_sapt==True or params.do_supermolecular==True)):
    #print(params.coordinates)
    index_min = params.coordinates.index(params.force_min)
    index_max = params.coordinates.index(params.force_max)
    index_ts = params.coordinates.index(0.0000)

if(params.do_supermolecular==True):
    irc_step_size = params.irc_stepsize
    coordinates = params.coordinates
    #sapt_method = params.sapt_method
    basis = params.basis
    del_E = [(energy - (e_A + e_B)) for energy in energies]
    sapt_ = sapt(sapt_geometries[:index_ts+1], params, params.method, basis, output_filename)
    int_  = sapt_.psi4_super()
    strain_e = []
    for i in range(len(energies[:index_ts+1])):
        strain_e.append((del_E[i] - int_[i]))
    reaction_force_strain_ = -1.0*np.gradient(strain_e,irc_step_size)
    reaction_force_int_ = -1.0*np.gradient(int_,irc_step_size)
    sapt_f_csv = open("super_force.csv", "w+")
    sapt_f_csv.write("Coordinate,F_strain,F_int\n")
    for i in range(len(coordinates[:index_ts+1])):
        sapt_f_csv.write("%f,%f,%f\n" %(coordinates[i], reaction_force_strain_[i], reaction_force_int_[i]))
    sapt_f_csv.close()
    sapt_e_csv = open("super_energy.csv", "w+")
    sapt_e_csv.write("Coordinate,E_strain,E_int\n")
    for i in range(len(coordinates[:index_ts+1])):
        sapt_e_csv.write("%f,%f,%f\n" %(coordinates[i], strain_e[i], int_[i]))
    sapt_e_csv.close()

    # Reaction Works for Region 1
    W_1_strain = -1.0*np.trapz(reaction_force_strain_[:index_min+1], dx=irc_step_size)
    W_1_int = -1.0*np.trapz(reaction_force_int_[:index_min+1], dx=irc_step_size)

    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 1)--\n')
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)'))
    output.write('-------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f}\n'.format(W_1_strain*627.509, W_1_int*627.509))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()

    # Reaction Works for Region 2
    W_2_strain = -1.0*np.trapz(reaction_force_strain_[index_min:index_ts+1], dx=irc_step_size)
    W_2_int = -1.0*np.trapz(reaction_force_int_[index_min:index_ts+1], dx=irc_step_size)

    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 2)--\n')
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)'))
    output.write('-------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f}\n'.format(W_2_strain*627.509, W_2_int*627.509))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()




if(params.do_sapt==True):
    irc_step_size = params.irc_stepsize
    coordinates = params.coordinates
    sapt_method = params.sapt_method
    basis = params.basis
    del_E = [(energy - (e_A + e_B)) for energy in energies]
    print(del_E)
    print(energies)
    output = open(output_filename, "a")
    output.write('\n\n--Symmetry-Adapted Perturbation Theory Energy Decomposition--\n')
    output.close()
    if(params.sapt_file!=None):
        sapt_data = pd.read_csv(params.sapt_file)
        int_ = sapt_data['E_int'].values
        elst_ = sapt_data['E_elst'].values
        exch_ = sapt_data['E_exch'].values
        ind_ = sapt_data['E_ind'].values
        disp_ = sapt_data['E_disp'].values
    else:
        sapt_ = sapt(sapt_geometries[:index_ts+1], params, sapt_method, basis, output_filename)
        int_, elst_, exch_, ind_, disp_  = sapt_.psi4_sapt()
    strain_e = []
    for i in range(len(energies[:index_ts+1])):
        strain_e.append((del_E[i] - int_[i]))
    reaction_force_strain_ = -1.0*np.gradient(strain_e,irc_step_size)
    reaction_force_int_ = -1.0*np.gradient(int_,irc_step_size)
    reaction_force_elst_ = -1.0*np.gradient(elst_,irc_step_size)
    reaction_force_exch_ = -1.0*np.gradient(exch_,irc_step_size)
    reaction_force_ind_ = -1.0*np.gradient(ind_,irc_step_size)
    reaction_force_disp_ = -1.0*np.gradient(disp_,irc_step_size)

   # reaction_force_int_limit = -1.0*np.gradient(int_[:index_ts+1],irc_step_size)
   # output = open(output_filename, "a")
   # output.write('\n\n--SAPT Energy Decomposition (Region 1)--\n')
   # output.close()
   # del_E = [(energy - (e_A + e_B)) for energy in energies]
   # sapt_1 = sapt(sapt_geometries[0:index_min+1], sapt_method, basis, output_filename)
   # int_1, elst_1, exch_1, ind_1, disp_1  = sapt_1.psi4_sapt()
   # strain_e_1 = []
   # for i in range(len(energies[0:index_min+1])):
   #     strain_e_1.append(del_E[i] - int_1[i])
   # reaction_force_strain_1 = -1.0*np.gradient(strain_e_1,irc_step_size)
   # reaction_force_int_1 = -1.0*np.gradient(int_1,irc_step_size)
   # reaction_force_elst_1 = -1.0*np.gradient(elst_1,irc_step_size)
   # reaction_force_exch_1 = -1.0*np.gradient(exch_1,irc_step_size)
   # reaction_force_ind_1 = -1.0*np.gradient(ind_1,irc_step_size)
   # reaction_force_disp_1 = -1.0*np.gradient(disp_1,irc_step_size)
    sapt_f_csv = open("sapt_force.csv", "w+")
    sapt_f_csv.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    for i in range(len(coordinates[:index_ts+1])):
        sapt_f_csv.write("%f,%f,%f,%f,%f,%f,%f\n" %(coordinates[i], reaction_force_strain_[i], reaction_force_int_[i], reaction_force_elst_[i], reaction_force_exch_[i], reaction_force_ind_[i], reaction_force_disp_[i]))
    sapt_f_csv.close()
    sapt_e_csv = open("sapt_energy.csv", "w+")
    sapt_e_csv.write("Coordinate,E_strain,E_int,E_elst,E_exch,E_ind,E_disp\n")
    for i in range(len(coordinates[:index_ts+1])):
        sapt_e_csv.write("%f,%f,%f,%f,%f,%f,%f\n" %(coordinates[i], strain_e[i], int_[i], elst_[i], exch_[i], ind_[i], disp_[i]))
    sapt_e_csv.close()

    # Reaction Works for Region 1
    W_1_strain = -1.0*np.trapz(reaction_force_strain_[:index_min+1], dx=irc_step_size)
    W_1_int = -1.0*np.trapz(reaction_force_int_[:index_min+1], dx=irc_step_size)
    W_1_elst = -1.0*np.trapz(reaction_force_elst_[:index_min+1], dx=irc_step_size)
    W_1_exch = -1.0*np.trapz(reaction_force_exch_[:index_min+1], dx=irc_step_size)
    W_1_ind = -1.0*np.trapz(reaction_force_ind_[:index_min+1], dx=irc_step_size)
    W_1_disp = -1.0*np.trapz(reaction_force_disp_[:index_min+1], dx=irc_step_size)
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 1)--\n')
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)', 'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    output.write('-------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_1_strain*627.509, W_1_int*627.509, W_1_elst*627.509, W_1_exch*627.509, W_1_ind*627.509, W_1_disp*627.509))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()

   # output = open(output_filename, "a")
   # output.write('\n\n--SAPT Energy Decomposition (Region 2)--\n')
   # output.close()
   # sapt_2 = sapt(sapt_geometries[index_min:index_ts+1], sapt_method, basis, output_filename)
   # int_2, elst_2, exch_2, ind_2, disp_2  = sapt_2.psi4_sapt()
   # strain_e_2 = []
   # for i in range(len(energies[index_min:index_ts+1])):
   #     strain_e_2.append(del_E[index_min+i] - int_2[i])
   # reaction_force_strain_2 = -1.0*np.gradient(strain_e_2,irc_step_size)
   # reaction_force_int_2 = -1.0*np.gradient(int_2,irc_step_size)
   # reaction_force_elst_2 = -1.0*np.gradient(elst_2,irc_step_size)
   # reaction_force_exch_2 = -1.0*np.gradient(exch_2,irc_step_size)
   # reaction_force_ind_2 = -1.0*np.gradient(ind_2,irc_step_size)
   # reaction_force_disp_2 = -1.0*np.gradient(disp_2,irc_step_size)
   # sapt_f_csv_2 = open("sapt_force_region_2.csv", "w+")
   # sapt_f_csv_2.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
   # for i in range(len(coordinates[index_min:index_ts+1])):
   #     sapt_f_csv_2.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[index_min+i], reaction_force_strain_2[i], reaction_force_int_2[i], reaction_force_elst_2[i], reaction_force_exch_2[i], reaction_force_ind_2[i], reaction_force_disp_2[i]))
   # sapt_f_csv_2.close()

    # Reaction Works for Region 2
    W_2_strain = -1.0*np.trapz(reaction_force_strain_[index_min:index_ts+1], dx=irc_step_size)
    W_2_int = -1.0*np.trapz(reaction_force_int_[index_min:index_ts+1], dx=irc_step_size)
    W_2_elst = -1.0*np.trapz(reaction_force_elst_[index_min:index_ts+1], dx=irc_step_size)
    W_2_exch = -1.0*np.trapz(reaction_force_exch_[index_min:index_ts+1], dx=irc_step_size)
    W_2_ind = -1.0*np.trapz(reaction_force_ind_[index_min:index_ts+1], dx=irc_step_size)
    W_2_disp = -1.0*np.trapz(reaction_force_disp_[index_min:index_ts+1], dx=irc_step_size)
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 2)--\n')
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)',
'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_2_strain*627.51, W_2_int*627.51, W_2_elst*627.51, W_2_exch*627.51, W_2_ind*627.51, W_2_disp*627.51))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()


   # # Reaction Works for Region 3
   # W_3_strain = -1.0*np.trapz(reaction_force_strain_[index_ts:index_max+1], dx=irc_step_size)
   # W_3_int = -1.0*np.trapz(reaction_force_int_[index_ts:index_max+1], dx=irc_step_size)
   # W_3_elst = -1.0*np.trapz(reaction_force_elst_[index_ts:index_max+1], dx=irc_step_size)
   # W_3_exch = -1.0*np.trapz(reaction_force_exch_[index_ts:index_max+1], dx=irc_step_size)
   # W_3_ind = -1.0*np.trapz(reaction_force_ind_[index_ts:index_max+1], dx=irc_step_size)
   # W_3_disp = -1.0*np.trapz(reaction_force_disp_[index_ts:index_max+1], dx=irc_step_size)
   # output = open(output_filename, "a")
   # output.write('\n\n--Reaction Work Decomposition (Region 3)--\n')
   # output.write('\n---------------------------------------------------------------------------------------------------')
   # output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)',
#'W_#elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
   # output.write('---------------------------------------------------------------------------------------------------\n')
   # output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_3_strain*627.51, W_3_int*627.51, W_3_elst*627.51, W_3_exch*627.51, W_3_ind*627.51, W_3_disp*627.51))
   # output.write('---------------------------------------------------------------------------------------------------\n')
   # output.close()

   # # Reaction Works for Region 4
   # W_4_strain = -1.0*np.trapz(reaction_force_strain_[index_max:], dx=irc_step_size)
   # W_4_int = -1.0*np.trapz(reaction_force_int_[index_max:], dx=irc_step_size)
   # W_4_elst = -1.0*np.trapz(reaction_force_elst_[index_max:], dx=irc_step_size)
   # W_4_exch = -1.0*np.trapz(reaction_force_exch_[index_max:], dx=irc_step_size)
   # W_4_ind = -1.0*np.trapz(reaction_force_ind_[index_max:], dx=irc_step_size)
   # W_4_disp = -1.0*np.trapz(reaction_force_disp_[index_max:], dx=irc_step_size)
   # output = open(output_filename, "a")
   # output.write('\n\n--Reaction Work Decomposition (Region 4)--\n')
   # output.write('\n-----------------------------------------------------------------------------------------------------')
   # output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)', 'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
   # output.write('---------------------------------------------------------------------------------------------------\n')
   # output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_4_strain*627.51, W_4_int*627.51, W_4_elst*627.51, W_4_exch*627.51, W_4_ind*627.51, W_4_disp*627.51))
   # output.write('---------------------------------------------------------------------------------------------------\n')
   # output.close()

   # output = open(output_filename, "a")
   # output.write('\n\n--SAPT Energy Decomposition (Region 3)--\n')
   # output.close()
   # sapt_3 = sapt(sapt_geometries[index_ts:index_max+1], sapt_method, basis, output_filename)
   # int_3, elst_3, exch_3, ind_3, disp_3  = sapt_3.psi4_sapt()
   # strain_e_3 = []
   # for i in range(len(energies[index_ts:index_max+1])):
   #     strain_e_3.append(del_E[index_ts+i] - int_3[i])
   # reaction_force_strain_3 = -1.0*np.gradient(strain_e_3,irc_step_size)
   # reaction_force_int_3 = -1.0*np.gradient(int_3,irc_step_size)
   # reaction_force_elst_3 = -1.0*np.gradient(elst_3,irc_step_size)
   # reaction_force_exch_3 = -1.0*np.gradient(exch_3,irc_step_size)
   # reaction_force_ind_3 = -1.0*np.gradient(ind_3,irc_step_size)
   # reaction_force_disp_3 = -1.0*np.gradient(disp_3,irc_step_size)
   # output = open(output_filename, "a")
   # output.write('\n\n--SAPT Energy Decomposition (Region 4)--\n')
   # output.close()
   # sapt_4 = sapt(sapt_geometries[index_max:], sapt_method, basis, output_filename)
   # int_4, elst_4, exch_4, ind_4, disp_4  = sapt_4.psi4_sapt()
   # strain_e_4 = []
   # for i in range(len(energies[index_ts:index_max+1])):
   #     strain_e_4.append(del_E[index_max+i] - int_4[i])
   # reaction_force_strain_4 = -1.0*np.gradient(strain_e_4,irc_step_size)
   # reaction_force_int_4 = -1.0*np.gradient(int_4,irc_step_size)
   # reaction_force_elst_4 = -1.0*np.gradient(elst_4,irc_step_size)
   # reaction_force_exch_4 = -1.0*np.gradient(exch_4,irc_step_size)
   # reaction_force_ind_4 = -1.0*np.gradient(ind_4,irc_step_size)
   # reaction_force_disp_4 = -1.0*np.gradient(disp_4,irc_step_size)

   # int_2.pop(0)
   # elst_2.pop(0) 
   # exch_2.pop(0) 
   # ind_2.pop(0) 
   # disp_2.pop(0)

   # int_3.pop(0)
   # elst_3.pop(0)
   # exch_3.pop(0)
   # ind_3.pop(0)
   # disp_3.pop(0)

   # int_4.pop(0)
   # elst_4.pop(0)
   # exch_4.pop(0)
   # ind_4.pop(0)
   # disp_4.pop(0)

   # int_all = int_1 + int_2 + int_3 + int_4
   # elst_all = elst_1 + elst_2 + elst_3 + elst_4
   # exch_all = exch_1 + exch_2 + exch_3 + exch_4
   # ind_all = ind_1 + ind_2 + ind_3 + ind_4
   # disp_all = disp_1 + disp_2 + disp_3 + disp_4
   # strain_e_all = []
   # for i in range(len(energies)):
   #     strain_e_all.append(del_E[i] - int_all[i])
   # reaction_force_strain_all = -1.0*np.gradient(strain_e_all,irc_step_size)
   # reaction_force_int_all = -1.0*np.gradient(int_all,irc_step_size)
   # reaction_force_elst_all = -1.0*np.gradient(elst_all,irc_step_size)
   # reaction_force_exch_all = -1.0*np.gradient(exch_all,irc_step_size)
   # reaction_force_ind_all = -1.0*np.gradient(ind_all,irc_step_size)
   # reaction_force_disp_all = -1.0*np.gradient(disp_all,irc_step_size)
   # sapt_f_csv_all = open("sapt_force.csv", "w+")
   # sapt_f_csv_all.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
   # for i in range(len(coordinates)):
   #     sapt_f_csv_all.write("%f,%f,%f,%f,%f,%f,%f\n" %(coordinates[i], reaction_force_strain_all[i], reaction_force_int_all[i], reaction_force_elst_all[i], reaction_force_exch_all[i], reaction_force_ind_all[i], reaction_force_disp_all[i]))
   # sapt_f_csv_all.close()
   # sapt_e_csv_all = open("sapt_energy.csv", "w+")
   # sapt_e_csv_all.write("Coordinate,E_strain,E_int,E_elst,E_exch,E_ind,E_disp\n")
   # for i in range(len(coordinates)):
   #     sapt_e_csv_all.write("%f,%f,%f,%f,%f,%f,%f\n" %(coordinates[i], strain_e_all[i], int_all[i], elst_all[i], exch_all[i], ind_all[i], disp_all[i]))
   # sapt_e_csv_all.close()

###############################
# Reaction Fragility Spectrum #
###############################

if(params.do_fragility_spec):
    fragility_spectrum.fragility_spec(output_filename)

############################
# Rexplot Plotting Utility #
############################

if(params.do_rexplot):
    csvread.plot()

if(params.do_polarization==True):
    potentials_A = np.array(ref.potential(wavefunctions_A))
    reaction_electronic_flux_A = (nelec_A/nelec)*np.array(ref.electronic_flux(irc_step_size))
    potential_diff_A = (nelec_A/nelec)*(potentials_A - potentials)
    transfer_flux_A = np.gradient(potential_diff_A)
    potentials_B = np.array(ref.potential(wavefunctions_B))
    reaction_electronic_flux_B = (nelec_B/nelec)*np.array(ref.electronic_flux(irc_step_size))
    potential_diff_B = (nelec_B/nelec)*(potentials_B - potentials)
    transfer_flux_B = np.gradient(potential_diff_B)
    transfer_flux = transfer_flux_A + transfer_flux_B 
    polarization_flux = reaction_electronic_flux - transfer_flux

output = open(output_filename, "a")
output.write("\n***pyREX Exiting Successfully***\n")
rand_int = random.randint(0, len(quotes.quotes)-1)
#print(rand_int)
output.write(quotes.quotes[rand_int])
