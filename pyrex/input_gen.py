import os
import numpy as np
import json
from . import _make_iterencode


# Necessary Information for JSON File Encoding
json.encoder._make_iterencode = _make_iterencode._make_iterencode
indent = (2, None)

# Building a standard pyREX input file
def build_standard():
    irc_filename = input("What is the name of your coordinate file?[full_irc.xyz] ") or "full_irc.xyz"
    #irc_filename = "full_irc.xyz"
    input_file = {}
    input_file["molecule"] = {}
    input_file["model"] = {}
    input_file["pyrex"] = {}
    irc_file = open(irc_filename, "r")
    lines = irc_file.readlines()
    natoms = int(lines[0])
    symbols = []
    for i in range(natoms):
        line = lines[2+i].split()
        symbols.append(str(line[0]))
    # Build Molecule Block
    input_file["molecule"]["symbols"] = symbols
    input_file["molecule"]["molecular_charge"] = "0"
    input_file["molecule"]["molecular_multiplicity"] = 1
    # Build Model Block
    input_file["model"]["method"] = "scf"
    input_file["model"]["basis"] = "sto-3g"
    # Build Pyrex Block
    input_file["pyrex"]["nthreads"] = 1
    input_file["pyrex"]["irc_filename"] = irc_filename 
    input_file["pyrex"]["do_energy"] = True
    input_file["pyrex"]["do_conceptualdft"] = True
    input_file["pyrex"]["irc_stepsize"] = 0.2
    json_input = open("input.json", "w+")
    json.dump(input_file, json_input, indent=indent)
    json_input.close()
    return input_file

# Building a standard pyREX input file
def build_frag():
    irc_filename = input("What is the name of your coordinate file?[full_irc.xyz] ") or "full_irc.xyz"
    frag_A = [int(i) for i in input("Define the atoms that compose the first fragment:").split()]
    frag_B = [int(i) for i in input("Define the atoms that compose the second fragment:").split()]
    charge_A = int(input("What is the charge of the first fragment?")) 
    charge_B = int(input("What is the charge of the second fragment?"))
    mult_A = int(input("What is the multiplicity of the first fragment?"))
    mult_B = int(input("What is the multiplicity of the second fragment?"))

    input_file = {}
    input_file["molecule"] = {}
    input_file["model"] = {}
    input_file["pyrex"] = {}
    irc_file = open(irc_filename, "r")
    lines = irc_file.readlines()
    natoms = int(lines[0])
    symbols = []
    for i in range(natoms):
        line = lines[2+i].split()
        symbols.append(str(line[0]))
    # Build Molecule Block
    input_file["molecule"]["symbols"] = symbols
    input_file["molecule"]["molecular_charge"] = "0"
    input_file["molecule"]["molecular_multiplicity"] = 1
    input_file["molecule"]["fragments"] = [frag_A, frag_B]
    input_file["molecule"]["fragment_charges"] = [charge_A, charge_B]
    input_file["molecule"]["fragment_multiplicities"] = [mult_A, mult_B]
    # Build Model Block
    input_file["model"]["method"] = "scf"
    input_file["model"]["basis"] = "sto-3g"
    # Build Pyrex Block
    input_file["pyrex"]["nthreads"] = 1
    input_file["pyrex"]["irc_filename"] = irc_filename
    input_file["pyrex"]["do_energy"] = True
    input_file["pyrex"]["do_conceptualdft"] = True
    input_file["pyrex"]["do_polarization"] = True
    input_file["pyrex"]["irc_stepsize"] = 0.2
    json_input = open("input.json", "w+")
    json.dump(input_file, json_input, indent=indent)
    json_input.close()
    return input_file

def join_irc():
    irc_forward_file = input("What is the name of your forward irc file?[irc_forward.xyz]") or "irc_forward.xyz"
    irc_backward_file = input("What is the name of your backward irc file?[irc_backward.xyz]") or "irc_backward.xyz"
    full_irc_file = open("full_irc.xyz", "w+") #creates new file for full irc   
 
    # Make empty list to store geometries along the IRC
    full_irc = []

    natoms = 0

    # Grab and store geometries from the forward IRC
    with open(irc_forward_file) as f:
        irc_num = 0
        for line in f:
            natoms = int(line)
            geom = []
            line = next(f)
            for i in range(natoms):                                                        
                line = next(f)
                geom.append(line)
            full_irc.append((irc_num,geom))
            irc_num = irc_num + 1
    
    with open(irc_backward_file) as f:
        irc_num = 0
        for line in f: 
            natoms = int(line)
            geom = []
            line = next(f)
            for i in range(natoms):
                line = next(f)
                geom.append(line)
            if(irc_num!=0):
                full_irc.append((irc_num,geom))
            else:
                pass
            irc_num = irc_num - 1
    
    # Sort the IRC
    sorted_full_irc = sorted(full_irc)

    for i in range(len(sorted_full_irc)):
        full_irc_file.write("%d\n" %natoms)
        full_irc_file.write("Full IRC Point %d\n" %sorted_full_irc[i][0])
        for j in range(len(sorted_full_irc[i][1])):
            full_irc_file.write("%s" %sorted_full_irc[i][1][j])

def build_irc():
    struct_filename = input("What is the name of your xyz coordinate file?[struct.xyz] ") or "struct.xyz"
    input_file = {}
    input_file["molecule"] = {}
    input_file["model"] = {}
    input_file["pyrex"] = {}
    input_file["irc"] = {}
    struct_file = open(struct_filename, "r")
    lines = struct_file.readlines()
    natoms = int(lines[0])
    symbols = []
    structure = []
    for i in range(natoms):
        line = lines[2+i].split()
        symbols.append(str(line[0]))
        for i in range(3):
            structure.append(float(line[i+1]))
    # Build Molecule Block
    input_file["molecule"]["geometry"] = structure
    input_file["molecule"]["symbols"] = symbols
    input_file["molecule"]["molecular_charge"] = "0"
    input_file["molecule"]["molecular_multiplicity"] = 1
    # Build Model Block
    input_file["model"]["method"] = "scf"
    input_file["model"]["basis"] = "sto-3g"
    # Build Pyrex Block
    input_file["pyrex"]["nthreads"] = 1
    # Build IRC Block
    input_file["irc"]["direction"] = "forward"
    input_file["irc"]["step_size"] = 0.01
    input_file["irc"]["mode"] = 1
    input_file["irc"]["normal_mode_file"] = "normal_modes.dat"
    json_input = open("input.json", "w+")
    json.dump(input_file, json_input, indent=indent)
    json_input.close()
    format_geom("input.json")

def format_geom(json_input):
    json_file = open(json_input, "r")
    lines = json_file.readlines()
    for i in range(len(lines)):
        line = lines[i]
        if("geometry" in line):
            index = i
            coords = line.split(":")[1]
            newstring = "     \"geometry\" : "
            count = 0
            for character in coords:
                if(character==","):
                    count = count + 1
                    newstring += character 
                else:
                    newstring += character
                if(count==3):
                    count = 0
                    newstring += "\n                    "
            print(newstring)
    newstring = newstring.rstrip()
    newstring += "\n"
    lines[index] = newstring
    formatted_file = open("input.json", "w+")
    for line in lines:
        formatted_file.write(line)
