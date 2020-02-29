import os
import time

_log_output = False
def log(msg):
    global _log_output
    if _log_output:
        print(msg)

sparrow_exe = "sparrow"

def run_sparrow(xyzfile, charge, multiplicity, method, sparrow_output):
    os.system("{} --structure {} --molecular_charge {} --spin_multiplicity {} --method {} > {}".format(sparrow_exe,xyzfile,charge,multiplicity, method, sparrow_output))

def sparrow_energy(sparrow_output):
    output_file = open(sparrow_output, "r")
    e = 0.0
    for line in output_file:
        if "Energy [hartree]" in line:
            line = next(output_file)
            e = float(line)
    return e  
    