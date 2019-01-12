"""
Base class to run SCF calculations
"""

import os
import numpy as np
import psi4
from pyscf import scf, gto

class scf_class(object):

    def __init__(self,data,outfile):
        self.level_of_theory = "%s/%s" %(data.method,data.basis)
        self.charge = data.molecular_charge
        self.mult = data.molecular_multiplicity
        self.basis = str(data.basis)
        self.outfile = outfile
        self.keywords = data.keywords

    def psi4_scf(self, geometries):
        """
        Function to run SCF along IRC using PSI4, returns array of energies at each geometry. 
        """
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
            # Code Related to JSON testing #
            #geom_parse = geom.split()
            #del geom_parse[0:2] #Remove Charge and Multiplicity 
            #print(geom_parse)
            geom += "symmetry c1"
            mol = psi4.geometry(geom)
            psi4.set_options(self.keywords)
            #print("pyREX:Single Point Calculation on IRC Point %d" %(count))
            energy, wfn = psi4.energy(self.level_of_theory, return_wfn=True)
            #wfn = psi4.core.Wavefunction.build(mol, self.basis)
            ndocc = wfn.doccpi()[0]
            #print(ndocc)
            eps = np.array(wfn.epsilon_a())
            #print(eps)
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

    def pyscf_scf(self, geometries):
        """
        Function to run SCF along IRC using pySCF, returns array of energies at each geometry.
        """
        if(os.path.isdir('pyscf_output')):
            pass
        else:
            os.makedirs("pyscf_output")
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
            mol = gto.Mole()
            mol.output = 'pyscf_output/irc_%d.dat' %count
            mol.verbose = 0
            mol.atom = geometry
            mol.basis = self.basis
            mol.charge = self.charge
            mol.spin = 0
            mol.build() 
            scf_obj = scf.RHF(mol)
            energy = scf_obj.scf()
            mo_e = scf_obj.mo_energy
            nocc = np.count_nonzero(scf_obj.mo_occ)
            homo_energy = mo_e[nocc - 1]
            lumo_energy = mo_e[nocc]
            energies.append(energy)
            output.write('{:>20} {:>20.4f} {:>20.4f} {:>20.4f}\n'.format(count, energy, homo_energy, lumo_energy))
            count = count+1
            output.close()
        output = open(self.outfile, "a")
        output.write('-------------------------------------------------------------------------------------\n')
        output.close()
    
        return energies

    def opt(self, label, natoms, geom):
        """
        Optimizes individual fragments for strain energy calculations.
        """
        output = open(self.outfile, "a")
        frag = ""
        geom = geom.split('\n')[:(natoms+2)]
        for i in range(natoms+2):
            frag += "%s\n" %geom[i]
        frag += "symmetry c1"
        print("Geometry Optimization on Fragment %s" %label)
        psi4.set_options(self.keywords)
        psi4.geometry(frag)
        psidump = "psi4_output/fragment_%s_opt.out" %label
        psi4.core.set_output_file(psidump, False)
        e = psi4.optimize(self.level_of_theory)
        #psi4.core.clean()
        return e
