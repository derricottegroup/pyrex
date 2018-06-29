"""
Script for Calculating Energies Along IRC in Psi4
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

def ts_gradient(irc_energies, step_size):
    gradient_value = -1.0*(-1.0*irc_energies[4] + 8.0*irc_energies[3] - 8.0*irc_energies[1] + 1.0*irc_energies[0])/(12.0*step_size)
    return gradient_value

def num_integrate(force_coordinates, reaction_force_values, lower_limit, upper_limit):
    integral_value = 0.0
    for i in range(lower_limit,upper_limit):
        integral_value += (force_coordinates[i+1] - force_coordinates[i])*((reaction_force_values[i] + reaction_force_values[i+1])/2.0)
    return integral_value

full_irc = open("coords.xyz", "r")
#energy_values = open("energy_values.dat", "w+")
irc = []

charge_1 = 0 #Specify total charge on fragment 1
mult_1 = 1 #Specify multiplicity on fragment 1
charge_2 = -1 #Specify total charge on fragment 2
mult_2 = 1 #Specify multiplicity on fragment 2

fragment_line = 5 # Which line in your geometry divides the fragments

irc_step_size = 0.1 #in units au*amu^(1/2), Psi4 default is 0.2
level_of_theory = "ssapt0/6-311G**" # Level of Theory for Total Energies

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
    print(geometry)
    psi4.core.set_output_file("irc_%d.out" %irc[i][0], False)
    psi4.geometry(geometry)
    psi4.set_options({'reference': 'rhf'})
    current_energy = psi4.energy(level_of_theory)
    irc_energies.append(current_energy)
ts_gradient_value = ts_gradient(irc_energies,irc_step_size)
print(ts_gradient_value)
