import os
import sys

def input_parser(input_file):
    user_values = {}
    infile = open(input_file, 'r').readlines()
    for line in infile:
        if line.startswith('Complex'):
            inp = line.split()
            user_values['charge_dimer'] = int(inp[1])
            user_values['mult_dimer']   = int(inp[2])
        if line.startswith('Frag_A'):
            inp = line.split()
            user_values['charge_A'] = int(inp[1])
            user_values['mult_A'] = int(inp[2])
        if line.startswith('Frag_B'):
            inp = line.split()
            user_values['charge_B'] = int(inp[1])
            user_values['mult_B'] = int(inp[2])
        if line.startswith('irc_step_size'):
            inp = line.split()
            user_values['irc_step_size'] = float(inp[1])
        if line.startswith('irc_filename'):
            inp = line.split()
            user_values['irc_filename'] = str(inp[1])
        if line.startswith('do_frag'):
            inp = line.split()
            user_values['do_frag'] = bool(inp[1])
        if line.startswith('method'):
            inp = line.split()
            user_values['method'] = str(inp[1])
        if line.startswith('basis'):
            inp = line.split()
            user_values['basis'] = str(inp[1])
        if line.startswith('do_polarization'):
            inp = line.split()
            user_values['do_polarization'] = str(inp[1])  
    return user_values