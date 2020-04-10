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
import random
import re
import datetime
import json
import euler
import input_reader
import input_gen
import fragility_spectrum
import quotes
import csvread
import pandas as pd
from decimal import Decimal
from concept_dft import *
from scf_class import *
from input_pydantic import *
import surface_scan
from geomparser import *
from atomic import atomic_decomp
from sapt_class import *
import argparse
import sys


# Load Output file
output_filename = "pyrex_output.dat"
#header(output_filename, json_data)
pyrex_help = """PYREX is an open-source python toolkit for intrinsic reactivity analysis"""
# Read in User Parameters
def read_params(args):
    if(os.path.isdir('psi4_output')):
        pass
    else:
        os.makedirs("psi4_output")
    if(args["input"]):
        json_input = args["input"]
    else:
        json_input = "input.json"
    json_data=open(json_input).read()
    params = input_reader.Params(json_data)
    return params

def params_read(PyrexInput):
    inp_data = build_inp(PyrexInput)
    json_input = json.dumps(inp_data)
    params = input_reader.Params(json_input)
    return params 



#TODO: Finish adding this type of structure. Use VMD_Cube from Evangelista Lab as a guide.
# parser = argparse.ArgumentParser(description=pyrex_help)
# def read_options():
#     parser.add_argument('--input', type=str, help="Specify input filename",default=None)
#     #parser.add_argument('--output', type=str, help="Specify output filename",default='pyrex_output.dat')
#     parser.add_argument('--rexplot', action="store_true", help="Use PYREX plotting utility")
#     parser.add_argument('--join_irc', action="store_true", help='Combine forward and backward IRCs. (string, default = false)')
#     parser.add_argument('--input_gen', help="Use Built in PYREX input generator", default=None)
#     parser.add_argument('--connectivity_matrix', action="store_true", help="Produce Connectivity Matrix for a single structure")
#     parser.add_argument('--fsapt_analyze', action="store_true", help="Utility to perform reaction F-SAPT decomposition of reaction force")

    #output_filename = args.output 
#def pyrex(output_filename):
# Load User Parameters
#params = input_reader.Params()
#level_of_theory = "%s/%s" %(params.method,params.basis) # Level of Theory for Total Energies

#########
## IRC ##
#########
def calculate_irc(params):
    if(params.do_irc):
        euler.ishida_morokuma(params,output_filename)

##################
## Surface Scan ##
##################
def calculate_surface_scan(params):
    if(params.do_surf_scan):
        surface_scan.surf_psi4(params, output_filename)

#############################
# Initialize Common Classes #
#############################

def geometry_builder(params):
    if(params.do_energy or params.do_sapt):
        geomparser = Geomparser(params.natoms, params.molecular_charge, params.molecular_multiplicity, params.geometries, params.coordinates)

####################################################
# Build Geometries and Print Geometric Information #
####################################################
    if(params.do_energy or params.do_sapt):
        geoms = geomparser.geombuilder()
        params.geoms = geoms 
        geomparser.atomic_distances()

    if(params.do_polarization):
        frag_geoms_irc = []
        for i in range(len(params.fraglist)):
            frag_geoms = geomparser.frag_ghost(params.frag_charge[i], params.frag_mult[i], params.fraglist[i])
            frag_geoms_irc.append(frag_geoms) 
        params.frag_geoms_irc = frag_geoms_irc 

    if(params.do_energy and params.qm_program == 'pyscf'):
        pyscf_geoms = geomparser.pyscf_geombuilder()
        params.geoms = pyscf_geoms  

    if(params.do_energy and params.qm_program == 'orca'):
        orca_geoms = geomparser.orca_geombuilder()
        params.geoms = orca_geoms

##########################
# Fragment Optimizations #
##########################
    if(params.do_frag==True and params.fsapt_analyze==False):
        scf_instance = scf_class(params, output_filename)
        iso_frag_A = geomparser.iso_frag(params.charge_A, params.mult_A, params.frag_A)
        iso_frag_B = geomparser.iso_frag(params.charge_B, params.mult_B, params.frag_B)
        params.e_A = scf_instance.opt("A", params.natoms_A, iso_frag_A[0])
        params.e_B = scf_instance.opt("B", params.natoms_B, iso_frag_B[0])

###########################
# Run Energy Calculations #
###########################

def energy_calculations(params):
    if(params.do_energy and params.energy_file==None):
        scf_instance = scf_class(params, output_filename)
        if(params.qm_program == 'psi4'):
            params.energies, params.wavefunctions = scf_instance.psi4_scf(params.geoms)
            nelec = params.wavefunctions[0].nalpha()
        if(params.do_polarization):
            frag_energies_irc = []
            frag_wavefunctions_irc = []
            for i in range(len(params.fraglist)):
                output = open(output_filename, "a")
                output.write('\n\n-- Calculating the Energy of Fragment %d--\n' %(i+1))
                output.close()
                frag_energies, frag_wavefunctions = scf_instance.psi4_scf(params.frag_geoms_irc[i])
                frag_energies_irc.append(frag_energies)
                frag_wavefunctions_irc.append(frag_wavefunctions)
            params.frag_energies_irc = frag_energies_irc
            params.frag_wavefunctions_irc = frag_wavefunctions_irc 

        if(params.qm_program == 'pyscf' and params.method == 'scf'):
            params.energies,params.frontier_orb_energies = scf_instance.pyscf_scf(params.geoms)
   
        if(params.qm_program == 'pyscf' and params.method == 'dft'):
            params.energies,params.frontier_orb_energies = scf_instance.pyscf_dft(params.geoms, params.xc_functional)

        if(params.qm_program == 'orca'):
            params.energies,params.frontier_orb_energies = scf_instance.orca_scf(params.geoms)

        if(params.qm_program == 'sparrow'):
            params.energies = scf_instance.sparrow_scf(params.geoms)

# Store energies in .csv file
        energy_csv = open("energy.csv", "w+")
        energy_csv.write("Coordinate,Energy\n")
        for i in range(len(params.coordinates)):
            energy_csv.write("%f, %f\n" %(params.coordinates[i], params.energies[i]))
        energy_csv.close()

# Read in energy values if given

    if(params.energy_file!=None):
        data = pd.read_csv(params.energy_file)
        params.energies = data['Energy'].values

###########################
# Reaction Force Analysis #
###########################

def reaction_force_analysis(params):
    if(params.do_energy or params.energy_file!=None):
        coordinates = params.coordinates
        reaction_force_values = -1.0*np.gradient(params.energies,params.irc_stepsize)
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
        
        params.index_min = np.argmin(np.asarray(reaction_force_values))
        ##print(reaction_force_values[index_min])
        #if(params.surf_scan_mode):
        #    params.index_ts = np.argmax(np.asarray(params.energies))
        #else:
        #    params.index_ts  = coordinates.index(0.0000)
        params.index_max = np.argmax(np.asarray(reaction_force_values))
        
        #r_struct = params.geoms[0]
        #force_min_struct = params.geoms[params.index_min]
        #ts_struct = params.geoms[params.index_ts]
        #force_max_struct = params.geoms[params.index_max]
        #p_struct = params.geoms[-1]
        
        #key_structs = [r_struct, force_min_struct, ts_struct, force_max_struct, p_struct]
        
        #scf_instance.psi4_molden(key_structs) 
        output = open(output_filename, "a")
        output.write('\n\n--Reaction Partitioning--\n')
        output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], coordinates[params.index_min]))
        output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(coordinates[params.index_min], coordinates[params.index_max]))
        output.write("\nProduct Region:            %.3f ------> %.3f\n" %(coordinates[params.index_max], coordinates[-1]))
        params.reaction_force_values = reaction_force_values 

###########################
# Reaction Force Constant #
###########################

def reaction_force_constant(params):
    if(params.do_energy or params.energy_file!=None):
        coordinates = params.coordinates
        reaction_force_constant_values = -1.0*np.gradient(params.reaction_force_values,params.irc_stepsize)
        output = open(output_filename, "a")
        output.write('\n\n--Reaction Force Constant--\n')
        output.write('\n-------------------------------------------------------------------------------------')
        output.write('\n{:>20} {:>20} {:>20}\n'.format('IRC Point', 'k (Hartree)', 'k (kcal)'))
        output.write('-------------------------------------------------------------------------------------\n')

        force_count = 0
        for reaction_force_constant in reaction_force_constant_values:
            output.write('\n{:>20} {:>20.7f} {:>20.7f}\n'.format(params.coordinates[force_count], reaction_force_constant, reaction_force_constant*627.51))
            force_count = force_count+1
        output.write('-------------------------------------------------------------------------------------\n')
        output.close()
        
        force_constant_csv = open("force_constant.csv", "w+")
        force_constant_csv.write("Coordinate,Force Constant\n")
        for i in range(len(coordinates)):
            force_constant_csv.write("%f, %f\n" %(coordinates[i], reaction_force_constant_values[i]))
        force_constant_csv.close()

##############################
# Activation Strain Analysis #
##############################

# Calculating interaction energies for activation strain analysis
def activation_strain_analysis(params):
    if(params.do_polarization):
        e_int = []
        del_e = []
        strain_e = []
        for i in range(len(params.coordinates)):
            frag_energy_sum = 0.0
            for j in range(len(params.fraglist)):
                frag_energy_sum += params.frag_energies_irc[j][i]
            current_e_int = params.energies[i] - frag_energy_sum
            current_del_e = params.energies[i] - (params.e_A + params.e_B)
            current_strain_e = current_del_e - current_e_int 
            del_e.append(current_del_e)
            e_int.append(current_e_int)
            strain_e.append(current_strain_e)
        act_strain = open("activation_strain.csv", "w+")
        act_strain.write("Coordinate,Energy,Interaction Energy,Strain Energy\n")
        for i in range(len(params.coordinates)):
            act_strain.write("%f,%f,%f,%f\n" %(params.coordinates[i], del_e[i], e_int[i], strain_e[i]))
        act_strain.close()
        


############################
# Conceptual DFT Analysis  #
############################

# Grabbing chemical potential and calculating REF using re_flux class
def conceptual_dft_analysis(params):
    if(params.do_conceptualdft):
        c_dft = concept_dft()
        if(params.molecular_multiplicity==1):
            if(params.qm_program == 'psi4'):
                potentials = np.array(c_dft.potential(params.wavefunctions))
            if(params.do_polarization):
                frag_potentials_irc = []
                for i in range(len(params.fraglist)):
                    frag_potentials = np.array(c_dft.potential(params.frag_wavefunctions_irc[i]))
                    frag_potentials_irc.append(frag_potentials)
            if(params.qm_program == 'pyscf' or params.qm_program == 'orca'):
                potentials = np.array(c_dft.potential_pyscf(params.frontier_orb_energies))
        else:
            potentials = np.array(c_dft.potential_open_shell(params.wavefunctions))
        if(params.qm_program == 'psi4'):
            hardness = np.array(c_dft.hardness(params.wavefunctions))
        if(params.qm_program == 'pyscf' or params.qm_program == 'orca'):
            hardness = np.array(c_dft.hardness_pyscf(params.frontier_orb_energies))

        reaction_electronic_flux = np.array(c_dft.electronic_flux(potentials, params.irc_stepsize))
        
        if(params.do_polarization):
            potential_diffs = c_dft.potential_diff_frag(len(params.fraglist),potentials,frag_potentials_irc)
            transfer_flux = np.array(c_dft.electronic_flux(potential_diffs, params.irc_stepsize))
            polarization_flux = np.array(reaction_electronic_flux - transfer_flux)
        electrophilicity = []
        for i in range(len(hardness)):
            electrophil = (potentials[i]*potentials[i])/(2.0*hardness[i])
            electrophilicity.append(electrophil)
    
        flux_csv = open("conceptual_dft.csv", "w+")
        flux_csv.write("Coordinate,Chemical Potential,Hardness,Electrophilicity,Reaction Electronic Flux\n")
        for i in range(len(params.coordinates)):
            flux_csv.write("%f,%f,%f,%f,%f\n" %(params.coordinates[i], potentials[i], hardness[i], electrophilicity[i], reaction_electronic_flux[i]))
        flux_csv.close()

        if(params.do_polarization):
            ref_decomp = open("flux_decomposition.csv", "w+")
            ref_decomp.write("Coordinate,Reaction Electronic Flux,Polarization Flux,Transfer Flux\n")
            for i in range(len(params.coordinates)):
                ref_decomp.write("%f,%f,%f,%f\n" %(params.coordinates[i], reaction_electronic_flux[i], polarization_flux[i], transfer_flux[i]))
            ref_decomp.close()


###########################
# Reaction Work Integrals #
###########################
def total_reaction_work_integrals(params):
    if(params.do_energy or params.energy_file!=None):
        no_ts_skip = False
        reaction_force_values = params.reaction_force_values
        #params.index_min = np.argmin(np.asarray(reaction_force_values))
        #params.index_max = np.argmax(np.asarray(reaction_force_values))
        if(params.surf_scan_mode):
            params.index_ts = np.argmax(np.asarray(params.energies))
        else:
            try:
                params.index_ts  = params.coordinates.index(0.0000)
            except ValueError:
                no_ts_skip = True
        if(no_ts_skip==False):
            W_1 = -1.0*np.trapz(reaction_force_values[:params.index_min+1], dx=params.irc_stepsize)
            W_2 = -1.0*np.trapz(reaction_force_values[params.index_min:params.index_ts+1], dx=params.irc_stepsize)
            W_3 = -1.0*np.trapz(reaction_force_values[params.index_ts:params.index_max+1], dx=params.irc_stepsize)
            W_4 = -1.0*np.trapz(reaction_force_values[params.index_max:], dx=params.irc_stepsize)
            output = open(output_filename, "a") 
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

def sapt_geometry_build(params):
    if(params.do_sapt==True or params.do_supermolecular==True):
        geomparser = Geomparser(params.natoms, params.molecular_charge, params.molecular_multiplicity, params.geometries, params.coordinates)
        params.sapt_geometries = geomparser.sapt_geombuilder(params.charge_A, params.mult_A, params.charge_B, params.mult_B, params.frag_A, params.frag_B)
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

def atomic_force_decomposition(params):
    if(params.do_atomic==True):
        atomic_forces = atomic_decomp(params, output_filename, params.geoms)
    
        W_1_conts = [] # Array to store the contributions to Work in Region 1
        for i in range(params.natoms):
            W_1_conts.append(np.sum(atomic_forces[i][:params.index_min]))
        W_2_conts = [] # Array to store the contributions to Work in Region 2
        for i in range(params.natoms):
            W_2_conts.append(np.sum(atomic_forces[i][params.index_min:params.index_ts]))
        W_3_conts = [] # Array to store the contributions to Work in Region 3
        for i in range(params.natoms):
            W_3_conts.append(np.sum(atomic_forces[i][params.index_ts:params.index_max]))
        W_4_conts = [] # Array to store the contributions to Work in Region 4
        for i in range(params.natoms):
            W_4_conts.append(np.sum(atomic_forces[i][params.index_max:]))
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

def sapt_energy_calc(params):
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
        del_E = [(energy - (params.e_A + params.e_B)) for energy in params.energies]
        sapt_ = sapt(params.sapt_geometries[:index_ts+1], params, params.method, basis, output_filename)
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
        #index_min = params.index_min
        #index_max = params.index_max
        #index_ts = params.index_ts 
        del_E = [(energy - (params.e_A + params.e_B)) for energy in params.energies]
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
            sapt_ = sapt(params.sapt_geometries[:index_ts+1], params, sapt_method, basis, output_filename)
            int_, elst_, exch_, ind_, disp_  = sapt_.psi4_sapt()
        strain_e = []
        for i in range(len(params.energies[:index_ts+1])):
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

def fsapt_analysis(params):
    sapt_ = sapt(params.sapt_geometries, params, params.method, params.basis, output_filename)
    sapt_.fsapt_post_analysis()
    
###############################
# Reaction Fragility Spectrum #
###############################
def reaction_fragility_spectrum(params):
    if(params.do_fragility_spec):
        fragility_spectrum.fragility_spec(params,params.geoms,output_filename)

############################
# Rexplot Plotting Utility #
############################

def rexplot_utility(json_input):
    csvread.plot(json_input)

#def flux_polarization(params):
#    if(params.do_polarization==True):
#        potentials_A = np.array(ref.potential(wavefunctions_A))
#        reaction_electronic_flux_A = (nelec_A/nelec)*np.array(ref.electronic_flux(irc_step_size))
#        potential_diff_A = (nelec_A/nelec)*(potentials_A - potentials)
#        transfer_flux_A = np.gradient(potential_diff_A)
#        potentials_B = np.array(ref.potential(wavefunctions_B))
#        reaction_electronic_flux_B = (nelec_B/nelec)*np.array(ref.electronic_flux(irc_step_size))
#        potential_diff_B = (nelec_B/nelec)*(potentials_B - potentials)
#        transfer_flux_B = np.gradient(potential_diff_B)
#        transfer_flux = transfer_flux_A + transfer_flux_B 
#        polarization_flux = reaction_electronic_flux - transfer_flux

def success(params):
    output = open(output_filename, "a")
    output.write("\n***pyREX Exiting Successfully***\n")
    rand_int = random.randint(0, len(quotes.quotes)-1)
    output.write(quotes.quotes[rand_int])

# def main():
#     read_options()
#     args = vars(parser.parse_args())
#     # These are all methods that need to skip the normal pyrex stuff
#     if(args["input_gen"]=="irc"):
#         input_gen.build_irc()
#     elif(args["input_gen"]=="energy"):
#         input_gen.build_standard()
#     elif(args["input_gen"]=="frag"):
#         input_gen.build_frag()
#     elif(args["join_irc"]==True):
#         input_gen.join_irc()
#     elif(args["rexplot"]==True):
#         json_input = args["input"]
#         rexplot_utility(json_input)
#     elif(args["fsapt_analyze"]==True):
#         params = read_params(args)
#         params.fsapt_analyze = True
#         geometry_builder(params)
#         sapt_geometry_build(params)
#         fsapt_analysis(params)
#     elif(args["connectivity_matrix"]==True):
#         params = read_params(args)
#         fragility_spectrum.single_connectivity_matrix(params,output_filename)

#     # Normal pyrex stuff
#     else:
#         params = read_params(args)
#         calculate_irc(params)
#         calculate_surface_scan(params)
#         geometry_builder(params)
#         energy_calculations(params)
#         reaction_force_analysis(params)
#         reaction_force_constant(params)
#         activation_strain_analysis(params)
#         conceptual_dft_analysis(params)
#         total_reaction_work_integrals(params)
#         sapt_geometry_build(params)
#         atomic_force_decomposition(params)
#         sapt_energy_calc(params)
#         reaction_fragility_spectrum(params)
#         #rexplot_utility(params)
#         #flux_polarization(params)
#         success(params)

# if __name__ == '__main__':
#     main()
