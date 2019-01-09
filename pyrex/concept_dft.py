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
        

    def electronic_flux(self, step):
        re_flux = -1.0*np.gradient(self.potentials,step)
        return re_flux
