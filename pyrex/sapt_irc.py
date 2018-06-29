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
charge_2 = -1 #Specify total charge on fragment 2
mult_2 = 1 #Specify multiplicity on fragment 2

fragment_line = 5 # Which line in your geometry divides the fragments

irc_step_size = 0.01 #in units au*amu^(1/2), Psi4 default is 0.2
level_of_theory = "sapt0/6-311G**" # Level of Theory for Total Energies

# Grab number of atoms (natoms) from the top of the XYZ file.
natoms = int(full_irc.readline())
coordinates = []

# Grab and store geometries from the IRC
for line in full_irc:
    if "IRC point" in line:
        geom = []
        irc_num_line = line.split()
        irc_num = int(irc_num_line[2])
        for i in range(natoms):
            line = next(full_irc)
            geom.append(line)
        irc.append((irc_num, geom))
        coordinates.append(irc_num*irc_step_size)

irc_energies = []
reaction_force = []
for i in range(len(irc)):
    geometry = ""
    geometry += "\n%d %d\n" %(charge_1, mult_1)
    for j in range(len(irc[i][1])):
        if j==(fragment_line):
            geometry += "--"
            geometry += "\n%d %d\n" %(charge_2, mult_2)
            geometry += irc[i][1][j]
        else:
            geometry += irc[i][1][j]
    #print(geometry)
    psi4.core.set_output_file("irc_%d.out" %irc[i][0], False)
    psi4.geometry(geometry)
    psi4.set_options({'reference': 'rhf'})
    current_energy = psi4.energy(level_of_theory)
    irc_energies.append(current_energy)
force_coordinates = coordinates[2:len(coordinates)-2]
reaction_force_values = num_first_derivs(irc_energies,irc_step_size)

for i in range(len(reaction_force_values)):
    reaction_force.append((force_coordinates[i],reaction_force_values[i]))
print(reaction_force)

index_min = np.argmin(np.asarray(reaction_force_values))
index_ts  = force_coordinates.index(1.37)
index_max = np.argmax(np.asarray(reaction_force_values))

# Calculate Work in Reactant Region
W_1 = -1.0*num_integrate(force_coordinates, reaction_force_values, 0, index_min)

#Calculate Work in Transition State Region 1
W_2 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_min, index_ts)

#Calculate Work in Transition State Region 2
W_3 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_ts, index_max)

#Calculate Work in Product Region
W_4 = -1.0*num_integrate(force_coordinates, reaction_force_values, index_max, len(reaction_force)-1)


energy_values.write("W_1 = %f kcal/mol\n" %(W_1*627.51))
energy_values.write("W_2 = %f kcal/mol\n" %(W_2*627.51))
energy_values.write("W_3 = %f kcal/mol\n" %(W_3*627.51))
energy_values.write("W_4 = %f kcal/mol\n" %(W_4*627.51))
energy_values.write("\n\n")
energy_values.writelines(["%f," % energy  for energy in irc_energies])
energy_values.write("\n\n")
energy_values.writelines(["%f," % coord  for coord in coordinates])
energy_values.write("\n\n")
energy_values.writelines(["%f," % force  for force in reaction_force_values])
