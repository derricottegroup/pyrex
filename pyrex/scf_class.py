import numpy as np
import psi4

class scf_class(object):

    def __init__(self,level_of_theory):
        self.level_of_theory = level_of_theory

    def psi4_scf(self, geometries):
        energies = []
        wavefunctions = []
        count = 0
        for geometry in geometries:
            psi4.core.set_output_file("psi4_output/irc_%d.out" %count, False)
            geom = geometry
            geom += "symmetry c1"
            psi4.geometry(geom)
            psi4.set_options({"reference" : "rhf"})
            print("pyREX:Single Point Calculation on IRC Point %d" %(count))
            energy, wfn = psi4.energy(self.level_of_theory, return_wfn=True)
            energies.append(energy)
            wavefunctions.append(wfn)
            count = count+1
        return energies, wavefunctions

    def opt(self, label, natoms, geom):
        frag = ""
        geom = geom.split('\n')[:(natoms+2)]
        for i in range(natoms+2):
            frag += "%s\n" %geom[i]
        frag += "symmetry c1"
        print("Geometry Optimization on Fragment %s" %label)
        psi4.set_options({"reference" : "rhf"})
        psi4.geometry(frag)
        psidump = "psi4_output/fragment_%s_opt.out" %label
        psi4.core.set_output_file(psidump, False)
        e = psi4.optimize(self.level_of_theory)
        return e
