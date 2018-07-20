import psi4
import numpy as np

def psi4_sapt(geometries, method, basis):
    sapt_contributions = []
    for i in range(len(geometries)):
        print("In SAPT Function")
        psi4.core.set_output_file("psi4_output/irc_%d_sapt.out" %i, False)
        geometry = geometries[i]
        geometry += "symmetry c1"
        psi4.geometry(geometry)
        #TODO Give the user control over these SCF options
        psi4.set_options({'reference': 'rhf', 'basis' : basis, 'exch_scale_alpha' : 'true', 'd_convergence' : 2})
        print("pyREX:SAPT Calculation on IRC Point %d" %(i))
        interaction_energy = psi4.energy(method)
        elst = psi4.core.get_variable("SAPT ELST ENERGY")
        exch = psi4.core.get_variable("SAPT EXCH ENERGY")
        ind  = psi4.core.get_variable("SAPT IND ENERGY")
        disp = psi4.core.get_variable("SAPT DISP ENERGY")
        sapt_contributions.append((interaction_energy, elst, exch, ind, disp))
    return sapt_contributions
