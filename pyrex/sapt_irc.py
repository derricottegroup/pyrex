"""
Script for Calculating SAPT Energies Along IRC in Psi4
"""
__authors__ = "Wallace D. Derricotte"
__credits__ = ["Wallace D. Derricotte"]

__copyright__ = "(c) 2018, Derricotte Research Group"
__license__ = "BSD-3-Clause"
__date__ = "2018-5-29"

import psi4
import numpy as np
import os

def num_first_derivs(irc_energies,step_size):
    reaction_force = []
    for i in range(2,len(irc_energies)-2):
        current_force = -1.0*(-1.0*irc_energies[i+2] + 8.0*irc_energies[i+1] - 8.0*irc_energies[i-1] + 1.0*irc_energies[i-2])/(12.0*step_size)
        reaction_force.append(current_force)
    return reaction_force

def num_integrate(force_coordinates, reaction_force_values, lower_limit, upper_limit):
    integral_value = 0.0
    for i in range(lower_limit,upper_limit):
        integral_value += (force_coordinates[i+1] - force_coordinates[i])*((reaction_force_values[i] + reaction_force_values[i+1])/2.0)
    return integral_value

full_irc = open("full_irc.xyz", "r")
energy_values = open("energy_values.dat", "w+")
irc = []

charge_1 = 0 #Specify total charge on fragment 1
mult_1 = 1 #Specify multiplicity on fragment 1
frag_1_atom_list = [0,1,4]
charge_2 = 0 #Specify total charge on fragment 2
mult_2 = 1 #Specify multiplicity on fragment 2
frag_2_atom_list = [2,3,5]

fragment_line = 5 # Which line in your geometry divides the fragments

irc_step_size = 0.2 #in units au*amu^(1/2), Psi4 default is 0.2
level_of_theory = "sapt2" # Level of Theory for Total Energies

# Grab number of atoms (natoms) from the top of the XYZ file.
natoms = int(full_irc.readline())
coordinates = []

# Grab and store geometries from the IRC
for line in full_irc:
    if "Full IRC Point" in line:
        geom = []
        irc_num_line = line.split()
        irc_num = int(irc_num_line[3])
        for i in range(natoms):
            line = next(full_irc)
            geom.append(line)
        irc.append((irc_num, geom))
        coordinates.append(irc_num*irc_step_size)

irc_energies = []
irc_elst10_energies = []
irc_elst12_energies = []
irc_elst13_energies = []
irc_exch10_energies = []
irc_exch11_energies = []
irc_exch12_energies = []
irc_ind20_energies = []
irc_ind22_energies = []
irc_exchind20_energies = []
irc_exchind22_energies = []
irc_exchdisp20_energies = []
irc_disp20_energies = []
irc_disp21_energies = []
irc_disp22_energies = []
irc_disp30_energies = []
reaction_force = []
for i in range(len(irc)):
    geometry = ""
    geometry += "\n%d %d\n" %(charge_1, mult_1)
    for j in range(len(frag_1_atom_list)):
        geometry += irc[i][1][frag_1_atom_list[j]]
    geometry += "--"
    geometry += "\n%d %d\n" %(charge_2, mult_2)
    for j in range(len(frag_2_atom_list)):
        geometry += irc[i][1][frag_2_atom_list[j]]
    print(geometry)
    psi4.core.set_output_file("irc_%d.out" %irc[i][0], False)
    psi4.geometry(geometry)
    psi4.set_options({'reference': 'rhf', 'basis' : '6-311G**'})
    current_energy = psi4.energy(level_of_theory)
    elst10 = psi4.core.get_variable("SAPT ELST10,R ENERGY")
    elst12 = psi4.core.get_variable("SAPT ELST12,R ENERGY")
    elst13 = psi4.core.get_variable("SAPT ELST13,R ENERGY")
    exch10 = psi4.core.get_variable("SAPT EXCH10(S^2) ENERGY")
    exch11 = psi4.core.get_variable("SAPT EXCH11(S^2) ENERGY")
    exch12 = psi4.core.get_variable("SAPT EXCH12(S^2) ENERGY")
    ind20  = psi4.core.get_variable("SAPT IND20,R ENERGY")
    ind22  = psi4.core.get_variable("SAPT IND22 ENERGY")
    exchind20 = psi4.core.get_variable("SAPT EXCH-IND20,R ENERGY")
    exchind22 = psi4.core.get_variable("SAPT EXCH-IND22 ENERGY")
    exchdisp20 = psi4.core.get_variable("SAPT EXCH-DISP20,R ENERGY") 
    disp20 = psi4.core.get_variable("SAPT DISP20 ENERGY")
    disp21 = psi4.core.get_variable("SAPT DISP21 ENERGY")
    disp22 = psi4.core.get_variable("SAPT DISP22 ENERGY")
    disp30 = psi4.core.get_variable("SAPT DISP30 ENERGY")
    irc_energies.append(current_energy)
    irc_elst10_energies.append(elst10)
    irc_elst12_energies.append(elst12)
    irc_elst13_energies.append(elst13)
    irc_exch10_energies.append(exch10)
    irc_exch11_energies.append(exch11)
    irc_exch12_energies.append(exch12)
    irc_ind20_energies.append(ind20)
    irc_ind22_energies.append(ind22)
    irc_exchind20_energies.append(exchind20)
    irc_exchind22_energies.append(exchind22)
    irc_exchdisp20_energies.append(exchdisp20)
    irc_disp20_energies.append(disp20)
    irc_disp21_energies.append(disp21)
    irc_disp22_energies.append(disp22)
    irc_disp30_energies.append(disp30)
force_coordinates = coordinates[2:len(coordinates)-2]
reaction_force_values = num_first_derivs(irc_energies,irc_step_size)
reaction_force_values_elst10 = num_first_derivs(irc_elst10_energies,irc_step_size)
reaction_force_values_elst12 = num_first_derivs(irc_elst12_energies,irc_step_size)
reaction_force_values_elst13 = num_first_derivs(irc_elst13_energies,irc_step_size)
reaction_force_values_exch10 = num_first_derivs(irc_exch10_energies,irc_step_size)
reaction_force_values_exch11 = num_first_derivs(irc_exch11_energies,irc_step_size)
reaction_force_values_exch12 = num_first_derivs(irc_exch12_energies,irc_step_size)
reaction_force_values_ind20 = num_first_derivs(irc_ind20_energies,irc_step_size)
reaction_force_values_ind22 = num_first_derivs(irc_ind22_energies,irc_step_size)
reaction_force_values_exchind20 = num_first_derivs(irc_exchind20_energies,irc_step_size)
reaction_force_values_exchind22 = num_first_derivs(irc_exchind22_energies,irc_step_size)
reaction_force_values_exchdisp20 = num_first_derivs(irc_exchdisp20_energies,irc_step_size)
reaction_force_values_disp20 = num_first_derivs(irc_disp20_energies,irc_step_size)
reaction_force_values_disp21 = num_first_derivs(irc_disp21_energies,irc_step_size)
reaction_force_values_disp30 = num_first_derivs(irc_disp30_energies,irc_step_size)
for i in range(len(reaction_force_values)):
    reaction_force.append((force_coordinates[i],reaction_force_values[i]))
print(reaction_force)

index_min = force_coordinates.index(-1.4000000000000001)
index_ts  = force_coordinates.index(0.0000)
index_max = force_coordinates.index(1.0)

# Calculate Work in Reactant Region
W_1_int = -1.0*num_integrate(force_coordinates, reaction_force_values, 0, index_min)
W_1_elst10 = -1.0*num_integrate(force_coordinates, reaction_force_values_elst10, 0, index_min)
W_1_exch10 = -1.0*num_integrate(force_coordinates, reaction_force_values_exch10, 0, index_min)
W_1_exchind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchind20, 0, index_min)
W_1_ind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_ind20, 0, index_min)
W_1_exchdisp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchdisp20, 0, index_min)
W_1_disp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_disp20, 0, index_min)
#Calculate Work in Transition State Region 1
W_2_int = -1.0*num_integrate(force_coordinates, reaction_force_values, index_min, index_ts)
W_2_elst10 = -1.0*num_integrate(force_coordinates, reaction_force_values_elst10, index_min, index_ts)
W_2_exch10 = -1.0*num_integrate(force_coordinates, reaction_force_values_exch10, index_min, index_ts)
W_2_exchind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchind20, index_min, index_ts)
W_2_ind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_ind20, index_min, index_ts)
W_2_exchdisp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchdisp20, index_min, index_ts)
W_2_disp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_disp20, index_min, index_ts)
#Calculate Work in Transition State Region 2
W_3_int = -1.0*num_integrate(force_coordinates, reaction_force_values, index_ts, index_max)
W_3_elst10 = -1.0*num_integrate(force_coordinates, reaction_force_values_elst10, index_ts, index_max)
W_3_exch10 = -1.0*num_integrate(force_coordinates, reaction_force_values_exch10, index_ts, index_max)
W_3_exchind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchind20, index_ts, index_max)
W_3_ind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_ind20, index_ts, index_max)
W_3_exchdisp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchdisp20, index_ts, index_max)
W_3_disp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_disp20, index_ts, index_max)

#Calculate Work in Product Region
W_4_int = -1.0*num_integrate(force_coordinates, reaction_force_values, index_max, len(reaction_force)-1)
W_4_elst10 = -1.0*num_integrate(force_coordinates, reaction_force_values_elst10, index_max, len(reaction_force)-1)
W_4_exch10 = -1.0*num_integrate(force_coordinates, reaction_force_values_exch10, index_max, len(reaction_force)-1)
W_4_exchind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchind20, index_max, len(reaction_force)-1)
W_4_ind20 = -1.0*num_integrate(force_coordinates, reaction_force_values_ind20, index_max, len(reaction_force)-1)
W_4_exchdisp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_exchdisp20, index_max, len(reaction_force)-1)
W_4_disp20 = -1.0*num_integrate(force_coordinates, reaction_force_values_disp20, index_max, len(reaction_force)-1)

energy_values.write("W_1_int = %f kcal/mol\n" %(W_1_int*627.51))
energy_values.write("W_1_elst10 = %f kcal/mol\n" %(W_1_elst10*627.51))
energy_values.write("W_1_exch10 = %f kcal/mol\n" %(W_1_exch10*627.51))
energy_values.write("W_1_exchind20 = %f kcal/mol\n" %(W_1_exchind20*627.51))
energy_values.write("W_1_ind20 = %f kcal/mol\n" %(W_1_ind20*627.51))
energy_values.write("W_1_exchdisp20 = %f kcal/mol\n" %(W_1_exchdisp20*627.51))
energy_values.write("W_1_disp20 = %f kcal/mol\n" %(W_1_disp20*627.51))
energy_values.write("W_2_int = %f kcal/mol\n" %(W_2_int*627.51))
energy_values.write("W_2_elst10 = %f kcal/mol\n" %(W_2_elst10*627.51))
energy_values.write("W_2_exch10 = %f kcal/mol\n" %(W_2_exch10*627.51))
energy_values.write("W_2_exchind20 = %f kcal/mol\n" %(W_2_exchind20*627.51))
energy_values.write("W_2_ind20 = %f kcal/mol\n" %(W_2_ind20*627.51))
energy_values.write("W_2_exchdisp20 = %f kcal/mol\n" %(W_2_exchdisp20*627.51))
energy_values.write("W_2_disp20 = %f kcal/mol\n" %(W_2_disp20*627.51))
energy_values.write("W_3_int = %f kcal/mol\n" %(W_3_int*627.51))
energy_values.write("W_3_elst10 = %f kcal/mol\n" %(W_3_elst10*627.51))
energy_values.write("W_3_exch10 = %f kcal/mol\n" %(W_3_exch10*627.51))
energy_values.write("W_3_exchind20 = %f kcal/mol\n" %(W_3_exchind20*627.51))
energy_values.write("W_3_ind20 = %f kcal/mol\n" %(W_3_ind20*627.51))
energy_values.write("W_3_exchdisp20 = %f kcal/mol\n" %(W_3_exchdisp20*627.51))
energy_values.write("W_3_disp20 = %f kcal/mol\n" %(W_3_disp20*627.51))
energy_values.write("W_4_int = %f kcal/mol\n" %(W_4_int*627.51))
energy_values.write("W_4_elst10 = %f kcal/mol\n" %(W_4_elst10*627.51))
energy_values.write("W_4_exch10 = %f kcal/mol\n" %(W_4_exch10*627.51))
energy_values.write("W_4_exchind20 = %f kcal/mol\n" %(W_4_exchind20*627.51))
energy_values.write("W_4_ind20 = %f kcal/mol\n" %(W_4_ind20*627.51))
energy_values.write("W_4_exchdisp20 = %f kcal/mol\n" %(W_4_exchdisp20*627.51))
energy_values.write("W_4_disp20 = %f kcal/mol\n" %(W_4_disp20*627.51))
energy_values.write("\n\n")
energy_values.write("INTERACTION ENERGY\n")
energy_values.writelines(["%f," % energy  for energy in irc_energies])
energy_values.write("\n\n")
#irc_elst_energies
energy_values.write("INTERACTION ENERGY (Electrostatics)\n")
energy_values.writelines(["%f," % energy  for energy in irc_elst10_energies])
energy_values.write("\n\n")
energy_values.write("INTERACTION ENERGY (Exchange)\n")
energy_values.writelines(["%f," % energy  for energy in irc_exch10_energies])
energy_values.write("\n\n")
energy_values.write("INTERACTION ENERGY (Induction)\n")
energy_values.writelines(["%f," % energy  for energy in irc_ind20_energies])
energy_values.write("\n\n")
energy_values.write("INTERACTION ENERGY (Dispersion)\n")
energy_values.writelines(["%f," % energy  for energy in irc_disp20_energies])
energy_values.write("\n\n")
energy_values.write("COORDINATE\n")
energy_values.writelines(["%f," % coord  for coord in coordinates])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (ELECTROSTATICS 10)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_elst10])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (ELECTROSTATICS 12)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_elst12])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (ELECTROSTATICS 13)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_elst13])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (EXCHANGE-REPULSION 10)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_exch10])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (EXCHANGE-REPULSION 11)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_exch11])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (EXCHANGE-REPULSION 12)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_exch12])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (INDUCTION 20)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_ind20])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (Exchange-INDUCTION 20)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_exchind20])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (INDUCTION 22)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_ind22])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (Exchange-INDUCTION 22)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_exchind22])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (DISPERSION 20)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_disp20])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (DISPERSION 21)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_disp21])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (DISPERSION 22)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_disp22])
energy_values.write("\n\n")
energy_values.write("INTERACTION FORCE (DISPERSION 30)\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values_disp30])
energy_values.write("\n\n")
