"""
A Script for Building the Full IRC for a reaction from its forward and backward pieces in Psi4
"""
__authors__ = "Wallace D. Derricotte"
__credits__ = ["Wallace D. Derricotte"]

__copyright__ = "(c) 2018, Derricotte Research Group"
__license__ = "BSD-3-Clause"
__date__ = "2018-5-29"


irc_forward = open("irc_forward.xyz", "r") #opens forward direction of irc file from psi4
irc_backward = open("irc_backward.xyz", "r") #opens backward direction of irc file from psi4
full_irc_file = open("full_irc.xyz", "w+") #creates new file for full irc


# Make Empty list to store geometries along the IRC
full_irc = []

# Grab number of atoms (natoms) from the top of the XYZ file.
natoms = int(irc_forward.readline())

# Grab and store geometries from the forward IRC
for line in irc_forward:
    if "IRC point" in line:
        geom = []
        irc_num_line = line.split()
        irc_num = int(irc_num_line[2])
        for i in range(natoms):
            line = next(irc_forward)
            geom.append(line)
        full_irc.append((irc_num, geom))

# Grab and store geometries from the backward IRC
for line in irc_backward:
    if "IRC point" in line:
        geom = []
        irc_num_line = line.split()
        irc_num = int(irc_num_line[2])
        if (irc_num==0):
            pass
        else:
            for i in range(natoms):
                line = next(irc_backward)
                geom.append(line)
            full_irc.append((irc_num,geom))

# Sort the IRC
sorted_full_irc = sorted(full_irc)

# Write the Full IRC to an output file
for i in range(len(sorted_full_irc)):
    full_irc_file.write("%d\n" %natoms)
    full_irc_file.write("Full IRC Point %d\n" %sorted_full_irc[i][0])
    for j in range(len(sorted_full_irc[i][1])):
        full_irc_file.write("%s" %sorted_full_irc[i][1][j])
