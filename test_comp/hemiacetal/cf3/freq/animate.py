import os
import numpy as np



natoms = 13 
normal_mode_file = open("output.default.16063.molden_normal_modes", "r") 
atom_labels = []
all_vibs = []
coords = []
for line in normal_mode_file:
    if("[FR-COORD]" in line):
        line = next(normal_mode_file)
        for i in range(natoms):
            atom_position = line.split()
            atom_labels.append(atom_position[0])
            atom_position = list(map(float,atom_position[1:]))
            if(i!=(natoms-1)):
                line = next(normal_mode_file)
            else:
                pass
            coords.append(atom_position)
        coords = np.array(coords)*0.529177
    if("vibration" in line):
        #print(line)
        line = next(normal_mode_file)
        vib_gradient = []
        for i in range(natoms):
            #print(line)
            atom_gradient = line.split()
            atom_gradient = list(map(float,atom_gradient))
            #print(atom_gradient)
            if(i!=(natoms-1)):
                line = next(normal_mode_file)
            else:
                pass
            vib_gradient.append(atom_gradient)
        vib_gradient = np.array(vib_gradient)*0.529177
        all_vibs.append(vib_gradient)
#print(all_vibs)
#print(coords)


count = 1
frames = 25
for vib in all_vibs:
    backward_frames = []
    forward_frames = []
    for i in range(frames):
        backward_frames.append(coords + (float((frames-i)/frames))*vib)
    for i in range(frames):
        forward_frames.append(coords + (float(i/frames))*vib)
    outfile = open("frequency%d.xyz" %count, "w+")
    outfile.write("%d \n" %natoms)
    outfile.write("Vibrational Frequency frame \n")
    for i in range(natoms):
        outfile.write("%s" %atom_labels[i])
        for j in range(3):
            if(j==2):
                outfile.write("  %f \n" %coords[i][j])
            else:
                outfile.write("  %f " %coords[i][j])
    for forward in forward_frames:
        outfile.write("%d \n" %natoms)
        outfile.write("Vibrational Frequency frame \n")
        for i in range(natoms):
            outfile.write("%s" %atom_labels[i])
            for j in range(3):
                if(j==2):
                    outfile.write("  %f \n" %forward[i][j])
                else:
                    outfile.write("  %f " %forward[i][j])
    for backward in backward_frames:
        outfile.write("%d \n" %natoms)
        outfile.write("Vibrational Frequency frame \n")
        for i in range(natoms):
            outfile.write("%s" %atom_labels[i])
            for j in range(3):
                if(j==2):
                    outfile.write("  %f \n" %backward[i][j])
                else:
                    outfile.write("  %f " %backward[i][j])
    count = count+1
