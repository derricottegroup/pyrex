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
import re
import datetime
import json
import euler
from re_flux import *
from scf_class import *
from geomparser import *
from sapt_class import *
from prettytable import PrettyTable

import sys

input_file = ""

#if(data["molecule"]["fragments"]):
#	print("THIS LOGIC WORKS!!!")

if len(sys.argv) == 1:
    input_file = "pyrex_input.dat"
if len(sys.argv) > 1:
    input_file = sys.argv[1]

json_data=open(input_file).read()

data = json.loads(json_data)

if(data["irc"]):
    euler.irc()

print(input_file)

user_values = input_parser(input_file)

if(os.path.isdir('psi4_output')):
    pass
else:
    os.makedirs("psi4_output")

output_filename = "pyrex_output.dat"
header(output_filename)
irc_filename = user_values['irc_filename']
full_irc = open(irc_filename, "r")
#output = open(output_filename, "w+")
irc = []

atom_symbols = []

do_frag = data["pyrex"]["do_frag"]
do_polarization = data["pyrex"]["do_polarization"]
do_sapt = data["pyrex"]["do_sapt"]
do_eda = data["pyrex"]["do_eda"]

print(do_polarization)
#charge_A = 0 #Specify total charge on Monomer B
if(do_frag==True):
    r_charge_A = data["molecule"]["fragment_charges"][0]
    r_mult_A = data["molecule"]["fragment_multiplicities"][0] #Multiplicity on Monomer A
    p_charge_A = data["molecule"]["fragment_charges"][0]
    p_mult_A = data["molecule"]["fragment_multiplicities"][0] #Multiplicity on Monomer A
    reactant_frag_A = data["molecule"]["fragments"][0]
    product_frag_A = data["molecule"]["fragments"][0]
    r_charge_B = data["molecule"]["fragment_charges"][1] #Specify total charge on Monomer B
    r_mult_B = data["molecule"]["fragment_multiplicities"][1] #Multiplicity on Monomer B
    p_charge_B = data["molecule"]["fragment_charges"][1] #Specify total charge on Monomer B
    p_mult_B = data["molecule"]["fragment_multiplicities"][1]#Multiplicity on Monomer B
    reactant_frag_B = data["molecule"]["fragments"][1]
    product_frag_B = data["molecule"]["fragments"][1]
    r_natoms_A = len(reactant_frag_A)
    r_natoms_B = len(reactant_frag_B)
    p_natoms_A = len(product_frag_A)
    p_natoms_B = len(product_frag_B)
else:
    charge_A = None
    mult_A = None
    frag_A_atom_list = None
    charge_B = None
    mult_B = None
    frag_B_atom_list = None
    natoms_A = None
    natoms_B = None

charge_dimer = data["molecule"]["molecular_charge"] #Specify total charge on the supermolecular complex
mult_dimer = data["molecule"]["molecular_multiplicity"] #Specify multiplicity of the supermolecular complex

#charge_mult = [charge_dimer,mult_dimer,charge_A,mult_A,charge_B,mult_B]
irc_step_size = data["pyrex"]["irc_stepsize"] #in units au*amu^(1/2), Psi4 default is 0.2
method = data["model"]["method"]
sapt_method = data["pyrex"]["sapt_method"]
basis = data["model"]["basis"]

level_of_theory = "%s/%s" %(method,basis) # Level of Theory for Total Energies


# Grab number of atoms (natoms) from the top of the XYZ file.
natoms = int(full_irc.readline())
coordinates = []
geometries = []
# Grab and store geometries from the IRC
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
        coordinates.append(irc_num*irc_step_size)
#TODO Store all coordinates in tuple array --> (,(dimer_geom,frag_a,frag_b),) prior to scf calls 
#irc_energies = []
reaction_force = []
reaction_force_strain = []
reaction_force_int = []
reaction_force_elst = []
reaction_force_exch = []
reaction_force_ind = []
reaction_force_disp = []

strain_energies = []
chemical_potentials = []
chemical_potentials_A = []
chemical_potentials_B = []
reaction_electronic_flux = []
reaction_electronic_flux_A = []
reaction_electronic_flux_B = []
polarization_flux = []
transfer_flux = []
#int_energies = []
#electrostatics = []
#exchange = []
#induction = []
#dispersion = []


e_A = 0.0
e_B = 0.0

table_header = ['IRC Point', 'E', 'Delta E', 'Potential']

t_pol_header = ['IRC Point', 'Del E', 'E_int', 'E_strain' ,'E_elst', 'E_pauli', 'E_orb']

t_sapt_header = ['IRC Point', 'E_int', 'E_elst', 'E_exch', 'E_ind', 'E_disp']

#TODO Incorporate the functions below into the main pyREX interface!
#table_header = ['IRC Point', 'E', 'Delta E', 'Potential', 'Potential A', 'Potential B']
t = PrettyTable(table_header)
t_pol = PrettyTable(t_pol_header)
t_sapt = PrettyTable(t_sapt_header)


geomparser = Geomparser(natoms, charge_dimer, mult_dimer, geometries, coordinates)

scf_instance = scf_class(data, output_filename)

geoms = geomparser.geombuilder()
geomparser.atomic_distances()

#frag_geoms_no_ghost = geomparser.frag_no_ghost(charge_A, mult_A, frag_A_atom_list)

if(do_sapt==True):
    r_sapt_geometries = geomparser.sapt_geombuilder(r_charge_A, r_mult_A, r_charge_B, r_mult_B, reactant_frag_A, reactant_frag_B)
    p_sapt_geometries = geomparser.sapt_geombuilder(p_charge_A, p_mult_A, p_charge_B, p_mult_B, product_frag_A, product_frag_B)

if(do_frag==True):
    r_iso_frag_A = geomparser.iso_frag(r_charge_A, r_mult_A, reactant_frag_A)
    r_iso_frag_B = geomparser.iso_frag(r_charge_B, r_mult_B, reactant_frag_B)
    p_iso_frag_A = geomparser.iso_frag(p_charge_A, p_mult_A, product_frag_A)
    p_iso_frag_B = geomparser.iso_frag(p_charge_B, p_mult_B, product_frag_B)
    r_e_A = scf_instance.opt("A", r_natoms_A, r_iso_frag_A[0])
    r_e_B = scf_instance.opt("B", r_natoms_B, r_iso_frag_B[0])
    #p_e_A = scf_instance.opt("A", p_natoms_A, p_iso_frag_A[-1])
    #p_e_B = scf_instance.opt("B", p_natoms_B, p_iso_frag_B[-1]) 

energies, wavefunctions = scf_instance.psi4_scf(geoms)
nelec = wavefunctions[0].nalpha()

energy_csv = open("energy.csv", "w+")
energy_csv.write("Coordinate,Energy\n")
for i in range(len(coordinates)):
    energy_csv.write("%f, %f\n" %(coordinates[i], energies[i]))
energy_csv.close()

# Calculate Reaction Forces directly after scf
reaction_force_values = -1.0*np.gradient(energies,irc_step_size)
output = open(output_filename, "a")
output.write('\n\n--Reaction Force Analysis--\n')
output.write('\n-------------------------------------------------------------------------------------')
output.write('\n{:>20} {:>20} {:>20}\n'.format('IRC Point', 'F (Hartree)', 'F (kcal)'))
output.write('-------------------------------------------------------------------------------------\n')

force_count = 0
for reaction_force in reaction_force_values:
    output.write('\n{:>20} {:>20.7f} {:>20.7f}\n'.format(force_count, reaction_force, reaction_force*627.51))
    force_count = force_count+1
output.write('-------------------------------------------------------------------------------------\n')
output.close()

force_csv = open("force.csv", "w+")
force_csv.write("Coordinate,Force\n")
for i in range(len(coordinates)):
    force_csv.write("%f, %f\n" %(coordinates[i], reaction_force_values[i]))
force_csv.close()

index_min = np.argmin(np.asarray(reaction_force_values))
print(reaction_force_values[index_min])
index_ts  = coordinates.index(0.0000)
index_max = np.argmax(np.asarray(reaction_force_values))

output = open(output_filename, "a")
output.write('\n\n--Reaction Partitioning--\n')
output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], coordinates[index_min]))
output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(coordinates[index_min], coordinates[index_max]))
output.write("\nProduct Region:            %.3f ------> %.3f\n" %(coordinates[index_max], coordinates[-1
]))

#Calculate Reaction Work For Each Region
W_1 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, 0, index_min)
W_2 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_min, index_ts)
W_3 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_ts, index_max)
W_4 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_max, len(reaction_force_values)-1)

output = open(output_filename, "a")
output.write('\n\n--Reaction Works--\n')
output.write('\n-------------------------------------------------------------------------------------')
output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('Units', 'W_1', 'W_2', 'W_3', 'W_4'))
output.write('\n-------------------------------------------------------------------------------------')
output.write('\n{:>15} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}\n'.format('Hartree', W_1, W_2, W_3, W_4))
output.write('\n{:>15} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}\n'.format('Kcal/mol', W_1*627.51, W_2*627.51, W_3*627.51, W_4*627.51))
output.write('\n-------------------------------------------------------------------------------------')
output.close()

if(do_polarization==True):
    frag_A_geoms = geomparser.frag_ghost(charge_A, mult_A, frag_A_atom_list)
    frag_B_geoms = geomparser.frag_ghost(charge_B, mult_B, frag_B_atom_list)
    energies_A, wavefunctions_A = scf_instance.psi4_scf(frag_A_geoms)
    energies_B, wavefunctions_B = scf_instance.psi4_scf(frag_B_geoms)
    nelec_A = wavefunctions_A[0].nalpha()
    nelec_B = wavefunctions_B[0].nalpha()

r_del_E = [(energy - (r_e_A + r_e_B)) for energy in energies]
#p_del_E = [(energy - (p_e_A + p_e_B)) for energy in energies]




if(do_sapt==True):
    output = open(output_filename, "a")
    output.write('\n\n--SAPT Energy Decomposition (Region 1)--\n')
    output.close()
    sapt_1 = sapt(r_sapt_geometries[0:index_min+1], sapt_method, basis, output_filename)
    int_1, elst_1, exch_1, ind_1, disp_1  = sapt_1.psi4_sapt()
    r_strain_e_1 = []
    for i in range(len(energies[0:index_min+1])):
        r_strain_e_1.append(r_del_E[i] - int_1[i])
    reaction_force_strain_1 = -1.0*np.gradient(r_strain_e_1,irc_step_size)
    reaction_force_int_1 = -1.0*np.gradient(int_1,irc_step_size)
    reaction_force_elst_1 = -1.0*np.gradient(elst_1,irc_step_size)
    reaction_force_exch_1 = -1.0*np.gradient(exch_1,irc_step_size)
    reaction_force_ind_1 = -1.0*np.gradient(ind_1,irc_step_size)
    reaction_force_disp_1 = -1.0*np.gradient(disp_1,irc_step_size)
    sapt_f_csv_1 = open("sapt_force_region_1.csv", "w+")
    sapt_f_csv_1.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    for i in range(len(coordinates[0:index_min+1])):
        sapt_f_csv_1.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[i], reaction_force_strain_1[i], reaction_force_int_1[i], reaction_force_elst_1[i], reaction_force_exch_1[i], reaction_force_ind_1[i], reaction_force_disp_1[i]))
    sapt_f_csv_1.close()
    W_1_strain = -1.0*np.trapz(reaction_force_strain_1, dx=irc_step_size)
    W_1_int = -1.0*np.trapz(reaction_force_int_1, dx=irc_step_size)
    W_1_elst = -1.0*np.trapz(reaction_force_elst_1, dx=irc_step_size)
    W_1_exch = -1.0*np.trapz(reaction_force_exch_1, dx=irc_step_size)
    W_1_ind = -1.0*np.trapz(reaction_force_ind_1, dx=irc_step_size)
    W_1_disp = -1.0*np.trapz(reaction_force_disp_1, dx=irc_step_size)
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 1 : %.4fkcal/mol)--\n' %(W_1*627.51))
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)', 'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    output.write('-------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_1_strain*627.51, W_1_int*627.51, W_1_elst*627.51, W_1_exch*627.51, W_1_ind*627.51, W_1_disp*627.51))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()

    output = open(output_filename, "a")
    output.write('\n\n--SAPT Energy Decomposition (Region 2)--\n')
    output.close()
    sapt_2 = sapt(r_sapt_geometries[index_min:index_ts+1], sapt_method, basis, output_filename)
    int_2, elst_2, exch_2, ind_2, disp_2  = sapt_2.psi4_sapt()
    r_strain_e_2 = []
    for i in range(len(energies[index_min:index_ts+1])):
        r_strain_e_2.append(r_del_E[index_min+i] - int_2[i])
    reaction_force_strain_2 = -1.0*np.gradient(r_strain_e_2,irc_step_size)
    reaction_force_int_2 = -1.0*np.gradient(int_2,irc_step_size)
    reaction_force_elst_2 = -1.0*np.gradient(elst_2,irc_step_size)
    reaction_force_exch_2 = -1.0*np.gradient(exch_2,irc_step_size)
    reaction_force_ind_2 = -1.0*np.gradient(ind_2,irc_step_size)
    reaction_force_disp_2 = -1.0*np.gradient(disp_2,irc_step_size)
    sapt_f_csv_2 = open("sapt_force_region_2.csv", "w+")
    sapt_f_csv_2.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    for i in range(len(coordinates[index_min:index_ts+1])):
        sapt_f_csv_2.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[index_min+i], reaction_force_strain_2[i], reaction_force_int_2[i], reaction_force_elst_2[i], reaction_force_exch_2[i], reaction_force_ind_2[i], reaction_force_disp_2[i]))
    sapt_f_csv_2.close()
    W_2_strain = -1.0*np.trapz(reaction_force_strain_2, dx=irc_step_size)
    W_2_int = -1.0*np.trapz(reaction_force_int_2, dx=irc_step_size)
    W_2_elst = -1.0*np.trapz(reaction_force_elst_2, dx=irc_step_size)
    W_2_exch = -1.0*np.trapz(reaction_force_exch_2, dx=irc_step_size)
    W_2_ind = -1.0*np.trapz(reaction_force_ind_2, dx=irc_step_size)
    W_2_disp = -1.0*np.trapz(reaction_force_disp_2, dx=irc_step_size)
    output = open(output_filename, "a")
    output.write('\n\n--Reaction Work Decomposition (Region 2 : %.4fkcal/mol)--\n' %(W_2*627.51))
    output.write('\n-------------------------------------------------------------------------------------------------')
    output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)',
'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_2_strain*627.51, W_2_int*627.51, W_2_elst*627.51, W_2_exch*627.51, W_2_ind*627.51, W_2_disp*627.51))
    output.write('--------------------------------------------------------------------------------------------------\n')
    output.close()
    output = open(output_filename, "a")
    output.write('\n\n--SAPT Energy Decomposition (Region 3)--\n')
    output.close()
    sapt_3 = sapt(r_sapt_geometries[index_ts:index_max+1], sapt_method, basis, output_filename)
    int_3, elst_3, exch_3, ind_3, disp_3  = sapt_3.psi4_sapt()
    r_strain_e_3 = []
    for i in range(len(energies[index_ts:index_max+1])):
        r_strain_e_3.append(r_del_E[index_ts+i] - int_3[i])
    reaction_force_strain_3 = -1.0*np.gradient(r_strain_e_3,irc_step_size)
    reaction_force_int_3 = -1.0*np.gradient(int_3,irc_step_size)
    reaction_force_elst_3 = -1.0*np.gradient(elst_3,irc_step_size)
    reaction_force_exch_3 = -1.0*np.gradient(exch_3,irc_step_size)
    reaction_force_ind_3 = -1.0*np.gradient(ind_3,irc_step_size)
    reaction_force_disp_3 = -1.0*np.gradient(disp_3,irc_step_size)
    output = open(output_filename, "a")
    output.write('\n\n--SAPT Energy Decomposition (Region 4)--\n')
    output.close()
    sapt_4 = sapt(r_sapt_geometries[index_max:], sapt_method, basis, output_filename)
    int_4, elst_4, exch_4, ind_4, disp_4  = sapt_4.psi4_sapt()
    r_strain_e_4 = []
    for i in range(len(energies[index_ts:index_max+1])):
        r_strain_e_4.append(r_del_E[index_max+i] - int_4[i])
    reaction_force_strain_4 = -1.0*np.gradient(r_strain_e_4,irc_step_size)
    reaction_force_int_4 = -1.0*np.gradient(int_4,irc_step_size)
    reaction_force_elst_4 = -1.0*np.gradient(elst_4,irc_step_size)
    reaction_force_exch_4 = -1.0*np.gradient(exch_4,irc_step_size)
    reaction_force_ind_4 = -1.0*np.gradient(ind_4,irc_step_size)
    reaction_force_disp_4 = -1.0*np.gradient(disp_4,irc_step_size)

    int_2.pop(0)
    elst_2.pop(0) 
    exch_2.pop(0) 
    ind_2.pop(0) 
    disp_2.pop(0)

    int_3.pop(0)
    elst_3.pop(0)
    exch_3.pop(0)
    ind_3.pop(0)
    disp_3.pop(0)

    int_4.pop(0)
    elst_4.pop(0)
    exch_4.pop(0)
    ind_4.pop(0)
    disp_4.pop(0)

    int_all = int_1 + int_2 + int_3 + int_4
    elst_all = elst_1 + elst_2 + elst_3 + elst_4
    exch_all = exch_1 + exch_2 + exch_3 + exch_4
    ind_all = ind_1 + ind_2 + ind_3 + ind_4
    disp_all = disp_1 + disp_2 + disp_3 + disp_4
    r_strain_e_all = []
    for i in range(len(energies)):
        r_strain_e_all.append(r_del_E[i] - int_all[i])
    reaction_force_strain_all = -1.0*np.gradient(r_strain_e_all,irc_step_size)
    reaction_force_int_all = -1.0*np.gradient(int_all,irc_step_size)
    reaction_force_elst_all = -1.0*np.gradient(elst_all,irc_step_size)
    reaction_force_exch_all = -1.0*np.gradient(exch_all,irc_step_size)
    reaction_force_ind_all = -1.0*np.gradient(ind_all,irc_step_size)
    reaction_force_disp_all = -1.0*np.gradient(disp_all,irc_step_size)
    sapt_f_csv_all = open("sapt_force.csv", "w+")
    sapt_f_csv_all.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    for i in range(len(coordinates)):
        sapt_f_csv_all.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[i], reaction_force_strain_all[i], reaction_force_int_all[i], reaction_force_elst_all[i], reaction_force_exch_all[i], reaction_force_ind_all[i], reaction_force_disp_all[i]))
    sapt_f_csv_all.close()
    sapt_e_csv_all = open("sapt_energy.csv", "w+")
    sapt_e_csv_all.write("Coordinate,E_strain,E_int,E_elst,E_exch,E_ind,E_disp\n")
    for i in range(len(coordinates)):
        sapt_e_csv_all.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[i], r_strain_e_all[i], int_all[i], elst_all[i], exch_all[i], ind_all[i], disp_all[i]))
    sapt_e_csv_all.close()
    #output = open(output_filename, "a")
    #output.write('\n\n--SAPT Energy Decomposition (Region 3)--\n')
    #output.close()
    #sapt_3 = sapt(p_sapt_geometries[index_ts:index_max+1], sapt_method, basis, output_filename)
    #int_3, elst_3, exch_3, ind_3, disp_3  = sapt_3.psi4_sapt()
    #p_strain_e_3 = []
    #for i in range(len(energies[index_ts:index_max+1])):
    #    p_strain_e_3.append(p_del_E[index_ts+i] - int_3[i])
    #reaction_force_strain_3 = -1.0*np.gradient(p_strain_e_3,irc_step_size)
    #reaction_force_int_3 = -1.0*np.gradient(int_3,irc_step_size)
    #reaction_force_elst_3 = -1.0*np.gradient(elst_3,irc_step_size)
    #reaction_force_exch_3 = -1.0*np.gradient(exch_3,irc_step_size)
    #reaction_force_ind_3 = -1.0*np.gradient(ind_3,irc_step_size)
    #reaction_force_disp_3 = -1.0*np.gradient(disp_3,irc_step_size)
    #sapt_f_csv_3 = open("sapt_force_region_3.csv", "w+")
    #sapt_f_csv_3.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    #for i in range(len(coordinates[index_ts:index_max+1])):
    #    sapt_f_csv_3.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[index_ts+i], reaction_force_strain_3[i], reaction_force_int_3[i], reaction_force_elst_3[i], reaction_force_exch_3[i], reaction_force_ind_3[i], reaction_force_disp_3[i]))
    #sapt_f_csv_3.close()
    #W_3_strain = -1.0*np.trapz(reaction_force_strain_3, dx=irc_step_size)
    #W_3_int = -1.0*np.trapz(reaction_force_int_3, dx=irc_step_size)
    #W_3_elst = -1.0*np.trapz(reaction_force_elst_3, dx=irc_step_size)
    #W_3_exch = -1.0*np.trapz(reaction_force_exch_3, dx=irc_step_size)
    #W_3_ind = -1.0*np.trapz(reaction_force_ind_3, dx=irc_step_size)
    #W_3_disp = -1.0*np.trapz(reaction_force_disp_3, dx=irc_step_size)
    #output = open(output_filename, "a")
    #output.write('\n\n--Reaction Work Decomposition (Region 3 : %.4fkcal/mol)--\n' %(W_3*627.51))
    #output.write('\n--------------------------------------------------------------------------------------------------')
    #output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)','W_e#lst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    #output.write('---------------------------------------------------------------------------------------------------\n')
    #output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_3_strain*627.51, W_3_int*627.51, W_3_elst*627.51, W_3_exch*627.51, W_3_ind*627.51, W_3_disp*627.51))
    #output.write('---------------------------------------------------------------------------------------------------\n')
    #output.close()
    #output = open(output_filename, "a")
    #output.write('\n\n--SAPT Energy Decomposition (Region 4)--\n')
    #output.close()
    #sapt_4 = sapt(p_sapt_geometries[index_max:], sapt_method, basis, output_filename)
    #int_4, elst_4, exch_4, ind_4, disp_4  = sapt_4.psi4_sapt()
    #p_strain_e_4 = []
    #for i in range(len(energies[index_max:])):
    #    p_strain_e_4.append(p_del_E[index_max+i] - int_4[i])
    #reaction_force_strain_4 = -1.0*np.gradient(p_strain_e_4,irc_step_size)
    #reaction_force_int_4 = -1.0*np.gradient(int_4,irc_step_size)
    #reaction_force_elst_4 = -1.0*np.gradient(elst_4,irc_step_size)
    #reaction_force_exch_4 = -1.0*np.gradient(exch_4,irc_step_size)
    #reaction_force_ind_4 = -1.0*np.gradient(ind_4,irc_step_size)
    #reaction_force_disp_4 = -1.0*np.gradient(disp_4,irc_step_size)
    #sapt_f_csv_4 = open("sapt_force_region_4.csv", "w+")
    #sapt_f_csv_4.write("Coordinate,F_strain,F_int,F_elst,F_exch,F_ind,F_disp\n")
    #for i in range(len(coordinates[index_max:])):
    #    sapt_f_csv_4.write("%f, %f, %f, %f, %f, %f, %f\n" %(coordinates[index_max+i], reaction_force_strain_4[i], reaction_force_int_4[i], reaction_force_elst_4[i], reaction_force_exch_4[i], reaction_force_ind_4[i], reaction_force_disp_4[i]))
    #sapt_f_csv_4.close()
    #W_4_strain = -1.0*np.trapz(reaction_force_strain_4, dx=irc_step_size)
    #W_4_int = -1.0*np.trapz(reaction_force_int_4, dx=irc_step_size)
    #W_4_elst = -1.0*np.trapz(reaction_force_elst_4, dx=irc_step_size)
    #W_4_exch = -1.0*np.trapz(reaction_force_exch_4, dx=irc_step_size)
    #W_4_ind = -1.0*np.trapz(reaction_force_ind_4, dx=irc_step_size)
    #W_4_disp = -1.0*np.trapz(reaction_force_disp_4, dx=irc_step_size)
    #output = open(output_filename, "a")
    #output.write('\n\n--Reaction Work Decomposition (Region 4 : %.4fkcal/mol)--\n' %(W_4*627.51))
    #output.write('\n----------------------------------------------------------------------------------------------------')
    #output.write('\n{:>15} {:>15} {:>15} {:>15} {:>15} {:>15}\n'.format('W_strain(kcal)', 'W_int(kcal)', 'W_elst(kcal)','W_exch(kcal)', 'W_ind(kcal)', 'W_disp(kcal)'))
    #output.write('----------------------------------------------------------------------------------------------------\n')
    #output.write('\n{:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f} {:>15.5f}\n'.format(W_4_strain*627.51, W_4_int*627.51, W_4_elst*627.51, W_4_exch*627.51, W_4_ind*627.51, W_4_disp*627.51))
    #output.write('---------------------------------------------------------------------------------------------------\n')
    #output.close()
# Grabbing chemical potential and calculating REF using re_flux class
#ref = re_flux()
#potentials = np.array(ref.potential(wavefunctions))
#reaction_electronic_flux = np.array(ref.electronic_flux(irc_step_size))

if(do_polarization==True):
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

# List Comprehensions!!! WOOT! 
#irc_energies = [energy[0] for energy in energies] # SCF Total Energies
#chemical_potentials = [potential[0] for potential in potentials] # Chemical Potentials
#del_E = [(energy - (e_A + e_B)) for energy in energies]  #Energy Change relative to optimized fragment
#for i in range(len(energies)):
#    t.add_row([i+1,"%.7f" %energies[i],"%.7f" %del_E[i], "%.7f" %potentials[i]])
#output = open(output_filename, "a")
#output.write("\n\n--Reaction Energy Analysis--\n\n")
#output.write("%s\n" %t.get_string())
#output.close()
#
#for i in range(len(energies)):
#    if(do_frag==True):
#        strain_energies.append(del_E[i] - int_[i])
#    
#    t.add_row([i+1,"%.7f" %energies[i],"%.7f" %del_E[i], "%.7f" %potentials[i]])
#    #if(do_eda==True):
#        #t_pol.add_row([i+1,"%.7f" %del_E[i], "%.7f" %interaction_energies[i][0], "%.7f" %(del_E[i] - interaction_energies[i][0]),"%.7f" %interaction_energies[i][1], "%.7f" %interaction_energies[i][2],  "%.7f" %interaction_energies[i][3]])
#    if(do_sapt==True):
#        t_sapt.add_row([i+1, "%.7f" %int_[i], "%.7f" %elst_[i], "%.7f" %exch_[i], "%.7f" %ind_[i], "%.7f" %disp_[i]])
#
#
#output = open(output_filename, "a")
#output.write("\n\n--Reaction Energy Analysis--\n\n")
#output.write("%s\n" %t.get_string())
#output.close()

if(do_eda==True):
    output = open(output_filename, "a")
    output.write("\n\n--Energy Decomposition Analysis--\n\n")
    output.write("%s\n" %t_pol.get_string())
    output.close()

#if(do_sapt==True):
#    output = open(output_filename, "a")
#    output.write("\n\n--Symmetry Adapted Perturbation Theory--\n\n")
#    output.write("%s\n" %t_sapt.get_string())
#    output.close()


#force_coordinates = coordinates[2:len(coordinates)-2]
#reaction_force_values = -1.0*np.gradient(energies,irc_step_size)
#reaction_electronic_flux = -1.0*np.gradient(potentials,irc_step_size)
#reaction_electronic_flux_A = -1.0*np.gradient(chemical_potentials_A,irc_step_size)
#reaction_electronic_flux_B = -1.0*np.gradient(chemical_potentials_B,irc_step_size)

#print(reaction_force_values)
#if(do_sapt==True):
#    reaction_force_strain = -1.0*np.gradient(strain_energies,irc_step_size)
#    reaction_force_int = -1.0*np.gradient(int_,irc_step_size)
#    reaction_force_elst = -1.0*np.gradient(elst_,irc_step_size)
#    reaction_force_exch = -1.0*np.gradient(exch_,irc_step_size)
#    reaction_force_ind = -1.0*np.gradient(ind_,irc_step_size)
#    reaction_force_disp = -1.0*np.gradient(disp_,irc_step_size)


#output = open(output_filename, "a")
#output.write("\n\n--IRC Force Partitioning--\n\n")
#output.write("Minimum Force =  %.10f\n" %(min(reaction_force_values)))
#output.write("Maximum Force =   %.10f\n" %(max(reaction_force_values)))


#for i in range(len(reaction_force_values)):
#    reaction_force.append((coordinates[i],reaction_force_values[i]))
#print(reaction_force)

#index_min = np.argmin(np.asarray(reaction_force_values))
#index_ts  = coordinates.index(0.0000)
#index_max = np.argmax(np.asarray(reaction_force_values))

#output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], coordinates[index_min]))
#output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(coordinates[index_min], coordinates[index_max]))
#output.write("\nProduct Region:            %.3f ------> %.3f\n" %(coordinates[index_max], coordinates[-1]))

# Calculate Work in Reactant Region
#W_1 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, 0, index_min)
#if(do_sapt==True):
#    W_1_strain = -1.0*calctools.num_integrate(coordinates, reaction_force_strain, 0, index_min)
#    W_1_int = -1.0*calctools.num_integrate(coordinates, reaction_force_int, 0, index_min)
#    W_1_elst = -1.0*calctools.num_integrate(coordinates, reaction_force_elst, 0, index_min)
#    W_1_exch = -1.0*calctools.num_integrate(coordinates, reaction_force_exch, 0, index_min)
#    W_1_ind = -1.0*calctools.num_integrate(coordinates, reaction_force_ind, 0, index_min)
#    W_1_disp = -1.0*calctools.num_integrate(coordinates, reaction_force_disp, 0, index_min)

#Calculate Work in Transition State Region 1
#W_2 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_min, index_ts)
#if(do_sapt==True):
#    W_2_strain = -1.0*calctools.num_integrate(coordinates,reaction_force_strain,index_min,index_ts)
#    W_2_int = -1.0*calctools.num_integrate(coordinates, reaction_force_int,index_min, index_ts)
#    W_2_elst = -1.0*calctools.num_integrate(coordinates, reaction_force_elst, index_min, index_ts)
#    W_2_exch = -1.0*calctools.num_integrate(coordinates, reaction_force_exch, index_min, index_ts)
#    W_2_ind = -1.0*calctools.num_integrate(coordinates, reaction_force_ind, index_min, index_ts)
#    W_2_disp = -1.0*calctools.num_integrate(coordinates, reaction_force_disp, index_min, index_ts)


#Calculate Work in Transition State Region 2
#W_3 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_ts, index_max)
#if(do_sapt==True):
#    W_3_strain = -1.0*calctools.num_integrate(force_coordinates, reaction_force_strain, index_ts, index_max)
#    W_3_int = -1.0*calctools.num_integrate(force_coordinates, reaction_force_int,index_ts, index_max)
#    W_3_elst = -1.0*calctools.num_integrate(force_coordinates, reaction_force_elst, index_ts, index_max)
#    W_3_exch = -1.0*calctools.num_integrate(force_coordinates, reaction_force_exch, index_ts, index_max)
#    W_3_ind = -1.0*calctools.num_integrate(force_coordinates, reaction_force_ind, index_ts, index_max)
#    W_3_disp = -1.0*calctools.num_integrate(force_coordinates, reaction_force_disp, index_ts, index_max)

#Calculate Work in Product Region
#W_4 = -1.0*calctools.num_integrate(coordinates, reaction_force_values, index_max, len(reaction_force)-1)
#if(do_sapt==True):
#    W_4_strain = -1.0*calctools.num_integrate(force_coordinates, reaction_force_strain, index_max, len(reaction_force)-1)
#    W_4_int = -1.0*calctools.num_integrate(force_coordinates, reaction_force_int, index_max, len(reaction_force)-1)
#    W_4_elst = -1.0*calctools.num_integrate(force_coordinates, reaction_force_elst, index_max, len(reaction_force)-1)
#    W_4_exch = -1.0*calctools.num_integrate(force_coordinates, reaction_force_exch, index_max, len(reaction_force)-1)
#    W_4_ind = -1.0*calctools.num_integrate(force_coordinates, reaction_force_ind, index_max, len(reaction_force)-1)
#    W_4_disp = -1.0*calctools.num_integrate(force_coordinates, reaction_force_disp, index_max, len(reaction_force)-1)

#output.write("\n\n--Work Integrals--\n\n")
#t_work = PrettyTable(["Unit","W_1", "W_2", "W_3", "W_4", "E_act", "E_react"])
#t_work.add_row(["Hartree", "%.7f" %W_1, "%.7f" %W_2, "%.7f" %W_3, "%.7f" %W_4,"%.7f" %(W_1+W_2), "%.7f" %(W_1+W_2+W_3+W_4)])
#t_work.add_row(["kcal/mol","%.7f" %(W_1*627.51), "%.7f" %(W_2*627.51), "%.7f" %W_3,"%.7f" %W_4,"%.7f" %((W_1+W_2)*627.51), "%.7f" %((W_1+W_2+W_3+W_4)*627.51)])
#output.write(t_work.get_string())
#
#if(do_sapt==True):
#    output.write("\n\n--Decomposition of Work Integrals Using SAPT--\n\n")
#    output.write("***DISCLAIMER: The sum of the components below will only be equal to the total work integral if you used a method that includes dispersion (i.e. MPn,CCSD, etc.) If the work was calculated with SCF or dispersionless DFT then W = W_strain + W_elst + W_exch + W_ind.")
#    output.write("\n\n--W_1 Decomposition (%.2f kcal/mol)--\n\n" %(W_1*627.51))
#    t_w1sapt = PrettyTable(["W_strain", "W_elst", "W_exch" , "W_ind", "W_disp"])
#    t_w1sapt.add_row(["%.7f" %(W_1_strain*627.51), "%.7f" %(W_1_elst*627.51), "%.7f" %(W_1_exch*627.51), "%.7f" %(W_1_ind*627.51),"%.7f" %(W_1_disp*627.51)])
#    output.write(t_w1sapt.get_string())
#
#if(do_sapt==True):
#    output.write("\n\n--W_2 Decomposition (%.2f kcal/mol)--\n\n" %(W_2*627.51))
#    t_w2sapt = PrettyTable(["W_strain", "W_elst", "W_exch" , "W_ind", "W_disp"])
#    t_w2sapt.add_row(["%.7f" %(W_2_strain*627.51), "%.7f" %(W_2_elst*627.51), "%.7f" %(W_2_exch*627.51), "%.7f" %(W_2_ind*627.51),"%.7f" %(W_2_disp*627.51)])
#    output.write(t_w2sapt.get_string())
#
##if(do_sapt==True):
##    output.write("\n\n--W_3 Decomposition (%.2f kcal/mol)--\n\n" %(W_3*627.51))
##    t_w3sapt = PrettyTable(["W_strain", "W_elst", "W_exch" , "W_ind", "W_disp"])
##    t_w3sapt.add_row(["%.7f" %(W_3_strain*627.51), "%.7f" %(W_3_elst*627.51), "%.7f" %(W_3_exch*627.51), "%.7f" %(W_3_ind*627.51),"%.7f" %(W_3_disp*627.51)])
##    output.write(t_w3sapt.get_string())
##if(do_sapt==True):
##    output.write("\n\n--W_4 Decomposition (%.2f kcal/mol)--\n\n" %(W_4*627.51))
##    t_w4sapt = PrettyTable(["W_strain", "W_elst", "W_exch" , "W_ind", "W_disp"])
##    t_w4sapt.add_row(["%.7f" %(W_4_strain*627.51), "%.7f" %(W_4_elst*627.51), "%.7f" %(W_4_exch*627.51), "%.7f" %(W_4_ind*627.51),"%.7f" %(W_4_disp*627.51)])
##    output.write(t_w4sapt.get_string())

# Calculate Energy Difference
#min_energy_index = np.argmin(np.asarray(energies))
#Del_E_raw = []
#for i in range(len(energies)):
#    Del_E_raw.append(energies[i] - energies[min_energy_index])
#
#Del_E = Del_E_raw

#Create CSV File

#output.write("----------------------------------------------------------------------------------------\n")
#output.write("   Coordinate(au amu^(1/2))            Energy_diff(au)                  Force(au)    \n")
#output.write("----------------------------------------------------------------------------------------\n")
#for i in range(len(reaction_force_values)):
#    output.write("       %.3f                           %.8f                             %.8f\n" %(force_coordinates[i], Del_E[i], reaction_force_values[i]))
#output.write("\n\n--Reaction Force and Electronic Flux Analysis--\n\n")
#t = PrettyTable(['Coordinate(au amu^(1/2))', 'Delta E', 'Force', 'Electronic Flux'])
##t.title = "pyREX Reaction Force Analysis Along Reaction Coordinate"
#for i in range(len(reaction_force_values)):
#    t.add_row(["%.2f" %coordinates[i], "%.7f" %Del_E[i], "%.7f" %reaction_force_values[i], "%.7f" %reaction_electronic_flux[i]])
#output.write(t.get_string())

#potentials_truncated = chemical_potentials[2:len(coordinates)-2]
#ref_csv = open("force_flux.csv","w+")
#ref_csv.write("Coordinate,DeltaE,Force,Chemical Potential,Reaction Electronic Flux\n")
#for i in range(len(reaction_force_values)):
#    ref_csv.write("%f,%f,%f,%f,%f\n" %(coordinates[i], Del_E[i], reaction_force_values[i], potentials[i],reaction_electronic_flux[i]))
#
#if(do_sapt==True):
#    sapt_e_csv = open("sapt_energy.csv" , "w+")
#    sapt_e_csv.write("Coordinate, E_elst, E_exch, E_ind, E_disp\n")
#    for i in range(len(coordinates)):
#        sapt_e_csv.write("%f, %f, %f, %f, %f\n" %(coordinates[i], elst_[i], exch_[i], ind_[i], disp_[i]))
#
#    sapt_f_csv = open("sapt_force.csv", "w+")
#    sapt_f_csv.write("Coordinate, F_elst, F_exch, F_ind, F_disp\n")
#    for i in range(len(coordinates)): 
#        sapt_f_csv.write("%f, %f, %f, %f, %f\n" %(coordinates[i], reaction_force_elst[i], reaction_force_exch[i], reaction_force_ind[i], reaction_force_disp[i]))   

output = open(output_filename, "a")
output.write("\n\n**pyREX Has Exited Successfully!**\n")
output.close()
