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

def psi4_scf(geometries, level_of_theory, frag=False):
    energies = []
    wavefunctions = []
    for i in range(len(geometries)):
        psi4.core.set_output_file("psi4_output/irc_%d.out" %i, False)
        psi4.geometry(geometries[i][0])
        #TODO Give the user control over these SCF options
        psi4.set_options({'reference': 'rhf'})
        print("pyREX:Single Point Calculation on IRC Point %d" %(i))
        dimer_energy, dimer_wfn = psi4.energy(level_of_theory, return_wfn=True)
        if(frag):
            # Fragment A SCF
            psi4.core.set_output_file("psi4_output/irc_%d_A_scf.out" %i, False)
            psi4.geometry(geometries[i][1])
            psi4.set_options({'reference': 'rhf'})
            print("pyREX:Single Point Calculation on IRC Point %d (Fragment A)" %(i))
            frag_A_energy, frag_A_wfn = psi4.energy(level_of_theory, return_wfn=True)
            # Fragment B SCF
            psi4.core.set_output_file("psi4_output/irc_%d_B_scf.out" %i, False)
            psi4.geometry(geometries[i][1])
            psi4.set_options({'reference': 'rhf'})
            print("pyREX:Single Point Calculation on IRC Point %d (Fragment B)" %(i))
            frag_B_energy, frag_B_wfn = psi4.energy(level_of_theory, return_wfn=True)
        else:
            frag_A_energy = None
            frag_A_wfn = None
            frag_B_energy = None
            frag_B_wfn = None
        energies.append((dimer_energy,frag_A_energy,frag_B_energy))
        wavefunctions.append((dimer_wfn, frag_A_wfn, frag_B_wfn))
    return energies, wavefunctions
        
