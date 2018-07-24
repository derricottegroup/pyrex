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
        if line.startswith('r_Frag_A'):
            inp = line.split()
            user_values['r_charge_A'] = int(inp[1])
            user_values['r_mult_A'] = int(inp[2])
        if line.startswith('r_Frag_B'):
            inp = line.split()
            user_values['r_charge_B'] = int(inp[1])
            user_values['r_mult_B'] = int(inp[2])
        if line.startswith('p_Frag_A'):
            inp = line.split()
            user_values['p_charge_A'] = int(inp[1])
            user_values['p_mult_A'] = int(inp[2])
        if line.startswith('p_Frag_B'):
            inp = line.split()
            user_values['p_charge_B'] = int(inp[1])
            user_values['p_mult_B'] = int(inp[2])
        if line.startswith('irc_step_size'):
            inp = line.split()
            user_values['irc_step_size'] = float(inp[1])
        if line.startswith('irc_filename'):
            inp = line.split()
            user_values['irc_filename'] = str(inp[1])
        if line.startswith('do_frag'):
            inp = line.split()
            user_values['do_frag'] = str_to_bool(inp[1])
        if line.startswith('do_eda'):
            inp = line.split()
            user_values['do_eda'] = str_to_bool(inp[1])
        if line.startswith('do_sapt'):
            inp = line.split()
            user_values['do_sapt'] = str_to_bool(inp[1])
        if line.startswith('method'):
            inp = line.split()
            user_values['method'] = str(inp[1])
        if line.startswith('sapt_method'):
            inp = line.split()
            user_values['sapt_method'] = str(inp[1])
        if line.startswith('basis'):
            inp = line.split()
            user_values['basis'] = str(inp[1])
        if line.startswith('do_polarization'):
            inp = line.split()
            user_values['do_polarization'] = str_to_bool(inp[1])
        if line.startswith('reactant_Frag_A'):
            inp = line.split('=')
            inp_array = inp[1].split()
            inp_array = [int(x) for x in inp_array]
            user_values['reactant_Frag_A'] = list(inp_array)
        if line.startswith('reactant_Frag_B'):
            inp = line.split('=')
            inp_array = inp[1].split()
            inp_array = [int(x) for x in inp_array]
            user_values['reactant_Frag_B'] = list(inp_array)
        if line.startswith('product_Frag_A'):
            inp = line.split('=')
            inp_array = inp[1].split()
            inp_array = [int(x) for x in inp_array]
            user_values['product_Frag_A'] = list(inp_array)
        if line.startswith('product_Frag_B'):
            inp = line.split('=')
            inp_array = inp[1].split()
            inp_array = [int(x) for x in inp_array]
            user_values['product_Frag_B'] = list(inp_array)
        else:
            pass 
    return user_values

def str_to_bool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Cannot convert %s to a boolean value" %string)
