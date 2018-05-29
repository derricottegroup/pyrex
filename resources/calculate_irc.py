"""
Script for Calculating Energies Along IRC in Psi4
"""
__authors__ = "Wallace D. Derricotte"
__credits__ = ["Wallace D. Derricotte"]

__copyright__ = "(c) 2018, Derricotte Research Group"
__license__ = "BSD-3-Clause"
__date__ = "2018-5-29"

import psi4
import numpy
import os

full_irc = open("full_irc.xyz", "r")
energy_values = open("energy_values.dat", "w+")
irc = []

charge = -1 #Specify total charge on your complex
mult = 1 #Specify multiplicity of your complex
irc_step_size = 0.2 #in units au*amu^(1/2), Psi4 default is 0.2
level_of_theory = "scf/sto-3g" # Level of Theory for Total Energies

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
for i in range(len(irc)):
    geometry = ""
    geometry += "\n%d %d\n" %(charge, mult)
    for j in range(len(irc[i][1])):
        geometry += irc[i][1][j]
    psi4.core.set_output_file("irc_%d.out" %irc[i][0], False)
    psi4.geometry(geometry)
    psi4.set_options({'reference': 'rhf'})
    current_energy = psi4.energy(level_of_theory)
    irc_energies.append(current_energy)
energy_values.writelines(["%f," % energy  for energy in irc_energies])
energy_values.write("\n\n")
energy_values.writelines(["%f," % coord  for coord in coordinates])
