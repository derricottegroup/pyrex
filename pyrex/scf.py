import numpy as np
import psi4

def frag_opt(geometry_frag, level_of_theory, outfile, label):
    output = open(outfile, "a")
    output.write("\n\n--Fragment %s Geometry Optimization--\n\n" %label)
    print("%s Geometry Optimization on Fragment %s" %(level_of_theory, label))
    output.write("\nInitial Geometry of Fragment:\n")
    output.write(geometry_frag[geometry_frag.find('\n')+4:])
    geometry_frag += "symmetry c1"
    psi4.set_options({'reference': 'rhf'})
    psi4.geometry(geometry_frag)
    psidump = "psi4_output/fragment_%s_opt.out" %label
    psi4.core.set_output_file(psidump, False)
    e_frag = psi4.optimize(level_of_theory)
    output.write("\nOutput Written To: %s\n" %psidump)
    output.write("Final Energy: %f\n" %e_frag)
    output.close()
    return e_frag
