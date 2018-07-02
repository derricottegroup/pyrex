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
import scf
import sapt
import geomtools
import wfn
import re
import datetime
from prettytable import PrettyTable

import sys

input_file = ""

if len(sys.argv) == 1:
    input_file = "pyrex_input.dat"
if len(sys.argv) > 1:
    input_file = sys.arv[1]

print(input_file)

user_values = input_parser(input_file)

if(os.path.isdir('psi4_output')):
    pass
else:
    os.makedirs("psi4_output")

output_filename = "pyrex_output.dat"
header(output_filename)
irc_filename = user_values['irc_filename']
full_irc = open("full_irc.xyz", "r")
#output = open(output_filename, "w+")
csv_file = open("raw_data.csv","w+")
irc = []

atom_symbols = []

#TODO Add the ability to read in these options from an input file
# Input File Should Include: Charge (dimer and fragments), fragment atom list, step size of irc, irc_filename, level of theory, psi4 options. 
do_frag = user_values['do_frag']
do_polarization = user_values['do_polarization']
do_sapt = user_values['do_sapt']
print(do_polarization)
#charge_A = 0 #Specify total charge on Monomer B
charge_A = user_values['charge_A']
mult_A = user_values['mult_A'] #Specify multiplicity on Monomer B
frag_A_atom_list = user_values['atom_list_Frag_A']
charge_B = user_values['charge_B'] #Specify total charge on Monomer B
mult_B = user_values['mult_B'] #Specify multiplicity on Monomer B
frag_B_atom_list = user_values['atom_list_Frag_B']

natoms_A = len(frag_A_atom_list)
natoms_B = len(frag_B_atom_list)

charge_dimer = user_values['charge_dimer'] #Specify total charge on the supermolecular complex
mult_dimer = user_values['mult_dimer'] #Specify multiplicity of the supermolecular complex

charge_mult = [charge_dimer,mult_dimer,charge_A,mult_A,charge_B,mult_B]
irc_step_size = user_values['irc_step_size'] #in units au*amu^(1/2), Psi4 default is 0.2
method = user_values['method']
sapt_method = user_values['sapt_method']
basis = user_values['basis']

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
irc_energies = []
reaction_force = []
strain_energies = []
chemical_potentials = []
chemical_potentials_A = []
chemical_potentials_B = []
reaction_electronic_flux = []
reaction_electronic_flux_A = []
reaction_electronic_flux_B = []

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
 
psi_geometries = geomtools.geombuilder_array(natoms,charge_mult,geometries, frag_A_atom_list, frag_B_atom_list)

sapt_geometries = geomtools.saptbuilder(natoms,charge_mult,geometries, frag_A_atom_list, frag_B_atom_list)

e_A , e_B = scf.frag_opt_new(psi_geometries, level_of_theory, output_filename, natoms_A, natoms_B)

energies, wavefunctions, interaction_energies = scf.psi4_scf(psi_geometries, level_of_theory, pol=bool(do_polarization))

sapt_contributions = sapt.psi4_sapt(sapt_geometries, sapt_method, basis)

potentials = wfn.potential(wavefunctions, pol=do_polarization)

for i in range(len(energies)):
    del_E = 0.0
    irc_energies.append(energies[i][0])
    chemical_potentials.append(potentials[i][0])
    if(do_polarization==True):
        chemical_potentials_A.append(potentials[i][1])
        chemical_potentials_B.append(potentials[i][2])
    else:
        chemical_potentials_A.append(0.0)
        chemical_potentials_B.append(0.0)
    del_E = energies[i][0] - (e_A + e_B)
    t.add_row([i+1,"%.7f" %energies[i][0],"%.7f" %del_E, "%.7f" %potentials[i][0]])
    t_pol.add_row([i+1,"%.7f" %del_E, "%.7f" %interaction_energies[i][0], "%.7f" %(del_E - interaction_energies[i][0]),"%.7f" %interaction_energies[i][1], "%.7f" %interaction_energies[i][2],  "%.7f" %interaction_energies[i][3]])
    t_sapt.add_row([i+1, "%.7f" %sapt_contributions[i][0], "%.7f" %sapt_contributions[i][1], "%.7f" %sapt_contributions[i][2], "%.7f" %sapt_contributions[i][3], "%.7f" %sapt_contributions[i][4]])


output = open(output_filename, "a")
output.write("\n\n--Reaction Energy Analysis--\n\n")
output.write("%s\n" %t.get_string())
output.close()

if(do_polarization==True):
    output = open(output_filename, "a")
    output.write("\n\n--Energy Decomposition Analysis--\n\n")
    output.write("%s\n" %t_pol.get_string())
    output.close()

if(do_sapt==True):
    output = open(output_filename, "a")
    output.write("\n\n--Symmetry Adapted Perturbation Theory--\n\n")
    output.write("%s\n" %t_sapt.get_string())
    output.close()


force_coordinates = coordinates[2:len(coordinates)-2]
reaction_force_values = calctools.num_first_derivs(irc_energies,irc_step_size)
reaction_electronic_flux = calctools.num_first_derivs(chemical_potentials,irc_step_size)
reaction_electronic_flux_A = calctools.num_first_derivs(chemical_potentials_A,irc_step_size)
reaction_electronic_flux_B = calctools.num_first_derivs(chemical_potentials_B,irc_step_size)



output = open(output_filename, "a")
output.write("\n\n--IRC Force Partitioning--\n\n")
output.write("Minimum Force =  %.10f\n" %(min(reaction_force_values)))
output.write("Maximum Force =   %.10f\n" %(max(reaction_force_values)))


for i in range(len(reaction_force_values)):
    reaction_force.append((force_coordinates[i],reaction_force_values[i]))
#print(reaction_force)

index_min = np.argmin(np.asarray(reaction_force_values))
index_ts  = force_coordinates.index(0.0000)
index_max = np.argmax(np.asarray(reaction_force_values))

output.write("\nReactant Region:          %.3f ------> %.3f\n" %(coordinates[0], force_coordinates[index_min]))
output.write("\nTransition State Region:  %.3f ------> %.3f\n" %(force_coordinates[index_min], force_coordinates[index_max]))
output.write("\nProduct Region:            %.3f ------> %.3f\n" %(force_coordinates[index_max], coordinates[-1]))

# Calculate Work in Reactant Region
W_1 = -1.0*calctools.num_integrate(force_coordinates, reaction_force_values, 0, index_min)

#Calculate Work in Transition State Region 1
W_2 = -1.0*calctools.num_integrate(force_coordinates, reaction_force_values, index_min, index_ts)

#Calculate Work in Transition State Region 2
W_3 = -1.0*calctools.num_integrate(force_coordinates, reaction_force_values, index_ts, index_max)

#Calculate Work in Product Region
W_4 = -1.0*calctools.num_integrate(force_coordinates, reaction_force_values, index_max, len(reaction_force)-1)


output.write("\n\n--Work Integrals--\n\n")
t = PrettyTable(["Unit","W_1", "W_2", "W_3" , "W_4", "E_act", "E_react"])
t.add_row(["Hartree", "%.7f" %W_1, "%.7f" %W_2, "%.7f" %W_3, "%.7f" %W_4, "%.7f" %(W_1+W_2), "%.7f" %(W_1+W_2+ W_3 + W_4)])
t.add_row(["kcal/mol","%.7f" %(W_1*627.51), "%.7f" %(W_2*627.51), "%.7f" %(W_3*627.51), "%.7f" %(W_4*627.51),"%.7f" %((W_1+W_2)*627.51), "%.7f" %((W_1+W_2+ W_3 + W_4)*627.51)])
output.write(t.get_string())

# Calculate Energy Difference
min_energy_index = np.argmin(np.asarray(irc_energies))
Del_E_raw = []
for i in range(len(irc_energies)):
    Del_E_raw.append(irc_energies[i] - irc_energies[min_energy_index])

Del_E = Del_E_raw[2:len(coordinates)-2]

#Create CSV File

#output.write("----------------------------------------------------------------------------------------\n")
#output.write("   Coordinate(au amu^(1/2))            Energy_diff(au)                  Force(au)    \n")
#output.write("----------------------------------------------------------------------------------------\n")
#for i in range(len(reaction_force_values)):
#    output.write("       %.3f                           %.8f                             %.8f\n" %(force_coordinates[i], Del_E[i], reaction_force_values[i]))
output.write("\n\n--Reaction Force Analysis--\n\n")
t = PrettyTable(['Coordinate(au amu^(1/2))', 'Delta E', 'Force', 'Electronic Flux', 'Flux_A', 'Flux_B'])
#t.title = "pyREX Reaction Force Analysis Along Reaction Coordinate"
for i in range(len(reaction_force_values)):
    t.add_row(["%.2f" %force_coordinates[i], "%.7f" %Del_E[i], "%.7f" %reaction_force_values[i], "%.7f" %reaction_electronic_flux[i], "%.7f" %reaction_electronic_flux_A[i], "%.7f" %reaction_electronic_flux_B[i]])
output.write(t.get_string())

potentials_truncated = chemical_potentials[2:len(coordinates)-2]
csv_file.write("Coordinate,DeltaE,Force,Chemical Potential,Reaction Electronic Flux, Flux A, Flux B\n")
for i in range(len(reaction_force_values)):
    csv_file.write("%f,%f,%f,%f,%f,%f,%f\n" %(force_coordinates[i], Del_E[i], reaction_force_values[i], potentials_truncated[i],reaction_electronic_flux[i], reaction_electronic_flux_A[i], reaction_electronic_flux_B[i]))

output.write("\n\n**pyREX Has Exited Successfully!**\n")
output.close()
