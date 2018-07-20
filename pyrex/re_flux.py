import psi4
import numpy as np


class re_flux(object):
    
    def __init__(self, wavefunctions, pol=False):
        self.wavefunctions = wavefunctions
        self.pol = pol

    def potential(self):
        self.potentials = []
        wfns = self.wavefunctions
        for wfn in wfns:
            dimer_wfn = wfn[0]
            ndocc = dimer_wfn.doccpi()[0]
            nelec = ndocc*2.0
            C = np.array(dimer_wfn.Ca())
            eps = dimer_wfn.epsilon_a()
            eps = np.array([eps.get(x) for x in range(C.shape[0])])
            homo_energy = eps[ndocc - 1]
            lumo_energy = eps[ndocc]
            dimer_potential = 0.5*(homo_energy)
            self.potentials.append(dimer_potential)
        return self.potentials

    def pol_transfer(self):
            self.potentials_A = []
            self.potentials_B = []
            wfns = self.wavefunctions
            for wfn in wfns:
            # Chemical Potential for Fragment A
            wfn_A = wfns[1]
            ndocc_A = wfn_A.doccpi()[0]
            nelec_A = ndocc_A*2.0
            C_A = np.array(wfn_A.Ca())
            eps_A = wfn_A.epsilon_a()
            eps_A = np.array([eps_A.get(x) for x in range(C_A.shape[0])])
            homo_energy_A = eps_A[ndocc-1]
            lumo_energy_A = eps_A[ndocc]
            self.potentials_A.append(0.5*(homo_energy_A + lumo_energy_A))
            # Chemical Potential for Fragment B
            wfn_B = wfns[2]
            ndocc_B = wfn_B.doccpi()[0]
            nelec_B = ndocc_B*2.0
            C_B = np.array(wfn_B.Ca())
            eps_B = wfn_B.epsilon_a()
            eps_B = np.array([eps_B.get(x) for x in range(C_B.shape[0])])
            homo_energy_B = eps_A[ndocc-1]
            lumo_energy_B = eps_A[ndocc]
            self.potentials_B.append(0.5*(homo_energy_B + lumo_energy_B))
            return self.potentials_A , self.potentials_B

    def electronic_flux(self, step):
        re_flux = np.gradient(self.potentials,step)
        return re_flux

    def pol_flux(self, step):
        if(self.pol==True):
            re_flux_pol_A = np.gradient(self.potentials_A, step)
            re_flux_pol_B = np.gradient(self.potentials_B, step)
        return re_flux_pol_A , re_flux_pol_B
    

        
