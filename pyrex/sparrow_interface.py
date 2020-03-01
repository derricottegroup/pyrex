import os
import time
import numpy as np

_log_output = False
def log(msg):
    global _log_output
    if _log_output:
        print(msg)

sparrow_exe = "sparrow"

def run_sparrow(xyzfile, charge, multiplicity, method, sparrow_output, additional_options):
    os.system("{} --structure {} --molecular_charge {} --spin_multiplicity {} --method {} {} > {}".format(
        sparrow_exe,xyzfile,charge,multiplicity, method, additional_options, sparrow_output))

def sparrow_energy(sparrow_output):
    output_file = open(sparrow_output, "r")
    e = 0.0
    for line in output_file:
        if "Energy [hartree]" in line:
            line = next(output_file)
            e = float(line)
    return e  

def sparrow_hessian(natoms, hessian_file):
    dimensions = 3*natoms
    hessian = np.zeros((dimensions,dimensions))
    hessian_output = open(hessian_file, "r")
    for line in hessian_output:
        if "Hessian" in line:
            line = next(hessian_output)
            for i in range(dimensions):
                row = line.split()
                for j in range(dimensions):
                    hessian[i][j] = float(row[j])
                if(i!=(dimensions-1)):
                    line = next(hessian_output)
                else:
                    pass
    hessian_output.close()
    return hessian

                        
        