import psi4
import numpy as np


class concept_dft(object):

    def potential(self, wavefunctions):
        self.potentials = []
        wfns = wavefunctions
        for wfn in wfns:
            ndocc = wfn.doccpi()[0]
            nelec = ndocc*2.0
            C = np.array(wfn.Ca())
            eps = wfn.epsilon_a()
            eps = np.array([eps.get(x) for x in range(C.shape[0])])
            homo_energy = eps[ndocc - 1]
            lumo_energy = eps[ndocc]
            potential = 0.5*(homo_energy + lumo_energy)
            self.potentials.append(potential)
        return self.potentials
    
    def potential_pyscf(self,frontier_orb_energies):
        self.potentials = []
        for foe in frontier_orb_energies:
            homo_energy = foe[0]
            lumo_energy = foe[1]
            potential = 0.5*(homo_energy + lumo_energy)
            self.potentials.append(potential)
        return self.potentials

    def potential_open_shell(self, wavefunctions):
        """	
        Function for handling open-shell cases for chemical potential
        """
        self.potentials = []
        wfns = wavefunctions
        for wfn in wfns:
            nalpha = wfn.nalpha()
            nbeta = wfn.nbeta()   # Get the number of alpha and beta electrons
            nelec = nalpha + nbeta
            C_a = np.array(wfn.Ca())
            C_b = np.array(wfn.Cb()) # Grab both coefficient matrices
            eps_a = wfn.epsilon_a()
            eps_a = np.array([eps_a.get(x) for x in range(C_a.shape[0])])
            eps_b = wfn.epsilon_b() # Store all orbital energies in vectors
            eps_b = np.array([eps_b.get(x) for x in range(C_b.shape[0])])
            homo_energy = 0.5*(eps_a[nalpha - 1] + eps_b[nbeta - 1])
            lumo_energy = 0.5*(eps_a[nalpha] + eps_b[nbeta]) 
            potential = 0.5*(homo_energy + lumo_energy)
            self.potentials.append(potential)
        return self.potentials
 
    def hardness(self, wavefunctions):
        self.hardness = []
        wfns = wavefunctions
        for wfn in wfns:
            ndocc = wfn.doccpi()[0]
            nelec = ndocc*2.0
            C = np.array(wfn.Ca())
            eps = wfn.epsilon_a()
            eps = np.array([eps.get(x) for x in range(C.shape[0])])
            homo_energy = eps[ndocc - 1]
            lumo_energy = eps[ndocc]
            hardness = lumo_energy - homo_energy
            self.hardness.append(hardness)
        return self.hardness
       
#TODO: Just randomly thought of a better way to do this, just have the SCF class detect the QM Program then let each
#      function calculate the potential, hardness, etc. differently based on that. Rather than having different 
#      functions for every QM program, this is fine for now with just two, but more would be a nightmare. 
    def hardness_pyscf(self, frontier_orb_energies):
        self.hardness = []
        for foe in frontier_orb_energies:
            homo_energy = foe[0]
            lumo_energy = foe[1] 
            hardness = lumo_energy - homo_energy
            self.hardness.append(hardness)
        return self.hardness 

    def electronic_flux(self, step):
        re_flux = -1.0*np.gradient(self.potentials,step)
        return re_flux
