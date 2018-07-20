import numpy as np
import psi4

class scf_class(object):

    def __init__(self,geometries,level_of_theory,scf_options)
        self.geometries = geometries
        self.level_of_theory = level_of_theory
        self.scf_options = scf_options

    def psi4_scf(self)
        energies = []
        wavefunctions = []
        count = 0
        for geometry in self.geometries:
            psi4.core.set_output_file("psi4_output/irc_%d.out" %count, False)
            geom = geometry[0]
            geom += "symmetry c1"
            psi4.geometry(geom)
            psi4.set_options(scf_options)
            print("pyREX:Single Point Calculation on IRC Point %d" %(count))
            energy, wfn = psi4.energy(level_of_theory, return_wfn=True)
            energies.append(energy)
            wavefunctions.append(wfn)
        return energy, wfn

    def frag_opt(self, natoms_A, natoms_B, opt_options):
        frag_A = ""
        geom = self.geometries[0][1]
        geom = geom.split('\n')[:(natoms_A+2)]
        for i in range(natoms_A+2):
            frag_A += "%s\n" %geom[i]
        frag_A += "symmetry c1"
        print("%s Geometry Optimization on Fragment A")
        psi4.set_options(opt_options)
        psi4.geometry(frag_A_geom)
        psidump = "psi4_output/fragment_A_opt.out"
        psi4.core.set_output_file(psidump, False)
        e_A = psi4.optimize(self.level_of_theory)
        frag_B = ""
        geom = self.geometries[0][2]
        geom = geom.split('\n')[:(natoms_B+2)]
        for i in range(natoms_B+2):
            frag_B += "%s\n" %geom[i]
        frag_B += "symmetry c1"
        print("%s Geometry Optimization on Fragment B")
        psi4.set_options(opt_options)
        psi4.geometry(frag_B_geom)
        psidump = "psi4_output/fragment_A_opt.out"
        psi4.core.set_output_file(psidump, False)
        e_B = psi4.optimize(self.level_of_theory)
        return e_A, e_B
