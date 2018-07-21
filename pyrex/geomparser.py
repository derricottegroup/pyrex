import numpy as np

class Geomparser(object):

    def __init__(self, natoms, charge, mult, geometries):
        self.charge = charge
        self.mult = mult
        self.geometries = geometries
        self.natoms = natoms

    def geombuilder(self):
        self.mol_inputs = []
        for geometry in self.geometries:
            mol_input = "\n%d %d\n" %(self.charge, self.mult)
            for j in range(self.natoms):
                line = geometry[j]
                mol_input += line
            self.mol_inputs.append(mol_input)
        return self.mol_inputs

    def iso_frag(self, frag_charge, frag_mult, fraglist):
        # This function returns ONLY the isolated fragment
        # DISCLAIMER: SCF convergence can be tough once the framents become significantly distorted
        # from their minimum energy structure. Use with caution.
        self.iso_frags = []
        for geometry in self.geometries:
            iso_frag_geom = "\n%d %d\n" %(frag_charge, frag_mult)
            for j in range(len(fraglist)):
                line = geometry[fraglist[j]]
                iso_frag_geom += line
            self.iso_frags.append(iso_frag_geom)
        return self.iso_frags

    def frag_ghost(self, frag_charge, frag_mult, fraglist):
        # This function returns the fragment with all other atoms as ghost atoms
        # Make array of all atoms and ghost out atoms not in the fragment
        atom_list = np.arange(self.natoms)
        ghost_list = [x for x in atom_list if x not in fraglist]
        self.frag_geoms = []
        for geometry in self.geometries:
            frag_geom = "\n%d %d\n" %(frag_charge, frag_mult)
            for j in range(len(fraglist)):
                line = geometry[fraglist[j]]
                frag_geom += line
            for j in range(len(ghost_list)):
                line = geometry[ghost_list[j]]
                frag_geom += "@"
                frag_geom += line
            self.frag_geoms.append(frag_geom)
        return self.frag_geoms
    
    def sapt_geombuilder(self, charge_A, mult_A, charge_B, mult_B, frag_A, frag_B):
        # Function builds geometries appropriate for the molecule block
        # of a SAPT calculation in PSI4.
        self.sapt_geoms = []
        for geometry in self.geometries:
            sapt_geom = "\n%d %d\n" %(charge_A, mult_A)
            for j in range(len(frag_A)):
                line = geometry[frag_A[j]]
                sapt_geom += line
            sapt_geom += "--\n"
            sapt_geom += "%d %d\n" %(charge_B, mult_B)
            for j in range(len(frag_B)):
                line = geometry[frag_B[j]]
                sapt_geom += line
            self.sapt_geoms.append(sapt_geom)
        return self.sapt_geoms
