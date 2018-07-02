import psi4
import numpy as np

def potential(wavefunctions, pol=False):
    potentials = []
    for i in range(len(wavefunctions)):
        dimer_wfn = wavefunctions[i][0]
        ndocc = dimer_wfn.doccpi()[0]
        nelec = ndocc*2.0
        C = np.array(dimer_wfn.Ca())
        eps = dimer_wfn.epsilon_a()
        eps = np.array([eps.get(x) for x in range(C.shape[0])])
        homo_energy = eps[ndocc-1]
        lumo_energy = eps[ndocc]
        dimer_potential = 0.5*(homo_energy + lumo_energy)
        dimer_hardness = 0.5*(lumo_energy - homo_energy)
        if(pol==True):
            # Fragment A Chemical Potential
            frag_A_wfn = wavefunctions[i][1]
            ndocc_A = frag_A_wfn.doccpi()[0]
            nelec_A = ndocc_A*2.0
            C_A = np.array(frag_A_wfn.Ca())
            eps_A = frag_A_wfn.epsilon_a()
            eps_A = np.array([eps_A.get(x) for x in range(C_A.shape[0])])
            homo_energy_A = eps_A[ndocc-1]
            lumo_energy_A = eps_A[ndocc]
            frag_A_potential = 0.5*(homo_energy_A + lumo_energy_A)
            frag_A_hardness = 0.5*(lumo_energy_A - homo_energy_A)
            # Fragment B Chemical Potential
            frag_B_wfn = wavefunctions[i][2]
            ndocc_B = frag_B_wfn.doccpi()[0]
            nelec_B = ndocc_B*2.0
            C_B = np.array(frag_B_wfn.Ca())
            eps_B = frag_B_wfn.epsilon_a()
            eps_B = np.array([eps_B.get(x) for x in range(C_B.shape[0])])
            homo_energy_B = eps_B[ndocc-1]
            lumo_energy_B = eps_B[ndocc]
            frag_B_potential = 0.5*(homo_energy_B + lumo_energy_B)
            frag_B_hardness = 0.5*(lumo_energy_B - homo_energy_B)
        else:
            frag_A_potential = 0.0
            frag_B_potential = 0.0
        potentials.append((dimer_potential,frag_A_potential,frag_B_potential))
    return potentials
