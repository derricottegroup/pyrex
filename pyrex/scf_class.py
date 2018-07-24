import numpy as np
import psi4

class scf_class(object):

    def __init__(self,level_of_theory,outfile):
        self.level_of_theory = level_of_theory
        self.outfile = outfile

    def psi4_scf(self, geometries):
        output = open(self.outfile, "a")
        output.write('\n\n--Reaction Energy--\n')
        output.write('\n-------------------------------------------------------------------------------------')
        output.write('\n{:>20} {:>20} {:>20} {:>20}\n'.format('IRC Point', 'E (Hartree)', 'HOMO (a.u.)','LUMO (a.u.)'))
        output.write('-------------------------------------------------------------------------------------\n')
        output.close()
        energies = []
        wavefunctions = []
        count = 0
        for geometry in geometries:
            output = open(self.outfile, "a")
            psi4.core.set_output_file("psi4_output/irc_%d.out" %count, False)
            geom = geometry
            geom += "symmetry c1"
            psi4.geometry(geom)
            psi4.set_options({"reference" : "rhf"})
            #print("pyREX:Single Point Calculation on IRC Point %d" %(count))
            energy, wfn = psi4.energy(self.level_of_theory, return_wfn=True)
            ndocc = wfn.doccpi()[0]
            eps = np.array(wfn.epsilon_a())
            homo_energy = eps[ndocc - 1]
            lumo_energy = eps[ndocc]
            energies.append(energy)
            wavefunctions.append(wfn)
            output.write('{:>20} {:>20.4f} {:>20.4f} {:>20.4f}\n'.format(count, energy, homo_energy, lumo_energy))
            count = count+1
            output.close()
        output = open(self.outfile, "a")
        output.write('-------------------------------------------------------------------------------------\n')
        output.close()
        return energies, wavefunctions

    def opt(self, label, natoms, geom):
        output = open(self.outfile, "a")
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
