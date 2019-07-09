"""
Class for Doing Potential Energy Scans along a Defined Geometric Coordinate
"""

__authors__  = "Wallace D. Derricotte"
__credits__  = ["Wallace D. Derricotte"]

__copyright__  = "(c) 2018-2019, Derricotte Research Group"
__license__    = "MIT License"
__date__       = "2019-01-02"

import numpy as np
import psi4
import os
import sys
import json


def surf_psi4(params,output_file):
    natoms = params.natoms
    geom = params.geometry
    geom += "symmetry c1\n"
    geom += "no_reorient\n"
    geom += "no_com"
    #print(geom)
    mol = psi4.geometry(geom)
    if(params.constraint_type=='angle'):
        coord_constraint = 'fixed_bend'
    if(params.constraint_type=='bond'):
        coord_constraint = 'fixed_distance'
    if(params.constraint_type=='dihedral'):
        coord_constraint = 'fixed_dihedral'
    output = open(output_file,"a")
    output.write('\n\n--%s surface scan for fixed %s of atoms %s--\n' %(params.scan_type, params.constraint_type, params.constrained_atoms))
    output.write('\n--------------------------------------------------------------------------------------\n')
    output.write('\n{:>20} {:>20}\n'.format('Coordinate', 'E'))
    output.write('-------------------------------------------------------------------------------------\n')
    output.close()
    surf_out = open("surface_scan.xyz", "w+")
    surf_out.close()
    for constrained_value in params.constrained_values:
        fixed_coord = params.constrained_atoms + str(constrained_value)
        #print(fixed_coord)
        psi4.set_options(params.keywords)
        psi4.set_module_options('Optking', {coord_constraint: fixed_coord})
        psi4.set_num_threads(params.nthreads)
        psi4.core.set_output_file("psi4_output/surf_%.2f.out" %constrained_value, False)
        surf_out = open("surface_scan.xyz", "a")
        surf_out.write("%s\n" %natoms)
        surf_out.write("%s surface scan with fixed %s of %f\n" %(params.scan_type, params.constraint_type, constrained_value))
        if(params.scan_type=='relaxed'):
            E = psi4.optimize(params.method)
            pre_string = mol.create_psi4_string_from_molecule()
            struct = pre_string.split("\n",5)[5]
            surf_out.write(struct)
        if(params.scan_type=='unrelaxed'):
            E = psi4.energy(params.method)
        output = open(output_file,"a")
        output.write('\n{:>20.4f} {:>20.7f}\n'.format(constrained_value, E))
        output.close()
