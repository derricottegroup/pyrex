import numpy as np
from pyscf import gto

class Geomparser(object):

    def __init__(self, natoms, charge, mult, geometries, coordinates):
        self.charge = charge
        self.mult = mult
        self.geometries = geometries
        self.natoms = natoms
        self.coordinates = coordinates

    def geombuilder(self):
        self.mol_inputs = []
        for geometry in self.geometries:
            mol_input = "\n%d %d\n" %(self.charge, self.mult)
            for j in range(self.natoms):
                line = geometry[j]
                mol_input += line
            self.mol_inputs.append(mol_input)
        return self.mol_inputs

    def pyscf_geombuilder(self):
        self.mol_objs = []
        for geometry in self.geometries:
            mol = ""
            for j in range(self.natoms):
                line = geometry[j]
                mol += line
            self.mol_objs.append(mol)
        return self.mol_objs

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
    
    def atomic_distances(self):
        atomic_distance_csv = open("atomic_distances.csv", "w+")
        bond_matrix = np.zeros((self.natoms,self.natoms))
        atom_labels = []
        unit_x = []
        unit_y = []
        unit_z = []
        geometry_map = {}
        for i in range(len(self.geometries)):
            Distance_Matrix = np.zeros((self.natoms, 3))
            for j in range(self.natoms):
                line = self.geometries[i][j]
                atom_position = line.split()
                atom_labels.append("%s%d" %(atom_position[0],j+1))
                Distance_Matrix[j][0] = atom_position[1]
                Distance_Matrix[j][1] = atom_position[2]
                Distance_Matrix[j][2] = atom_position[3]
            geometry_map[self.coordinates[i]] = Distance_Matrix
        #Make CSV Header
        atomic_distance_csv.write("Coordinate,")
        for i in range(self.natoms):
            for j in range(self.natoms):
                if(i>j):
                    atomic_distance_csv.write("%s-%s," %(atom_labels[i], atom_labels[j]))
                else:
                    pass
        atomic_distance_csv.write("\n")
        for z in range(len(self.geometries)):
            coord = self.coordinates[z]
            Geom = geometry_map[coord]
            atomic_distance_csv.write("%.3f," %coord)
            ex = np.zeros((self.natoms,self.natoms))
            ey = np.zeros((self.natoms,self.natoms))
            ez = np.zeros((self.natoms,self.natoms))
            for i in range(self.natoms):
                for j in range(self.natoms):
                    if(i>j):
                        R = np.sqrt((Geom[i][0]- Geom[j][0])**(2.0) + (Geom[i][1]- Geom[j][1])**(2.0) + (Geom[i][2]- Geom[j][2])**(2.0))
                        ex[i][j] = ex[j][i] = -(Geom[i][0]- Geom[j][0])/R
                        ey[i][j] = ey[j][i] = -(Geom[i][1]- Geom[j][1])/R
                        ez[i][j] = ez[j][i] = -(Geom[i][2]- Geom[j][2])/R
                        bond_matrix[i][j] = R
                        atomic_distance_csv.write("%.7f," %R)
                    else:
                        pass
            unit_x.append(ex)
            unit_y.append(ey)
            unit_z.append(ez)    
            atomic_distance_csv.write("\n")    
        angles_csv = open("angles.csv", "w+")
        angles_csv.write("Coordinate,")
        for i in range(self.natoms):
            for j in range(self.natoms):
                for k in range(self.natoms):
                    if(i>j and j>k and bond_matrix[i][j] < 4.0 and bond_matrix[j][k] < 4.0):
                        angles_csv.write("%s-%s-%s," %(atom_labels[i], atom_labels[j], atom_labels[k]))
        angles_csv.write("\n")
        for z in range(len(self.geometries)):
            unit_vec_x = unit_x[z]
            unit_vec_y = unit_y[z]
            unit_vec_z = unit_z[z]
            for i in range(self.natoms):
                for j in range(self.natoms):
                    for k in range(self.natoms):
                        if(i>j and j>k and bond_matrix[i][j] < 4.0 and bond_matrix[j][k] < 4.0):
                            x = unit_vec_x[j][i]*unit_vec_x[j][k]
                            y = unit_vec_y[j][i]*unit_vec_y[j][k]
                            z = unit_vec_z[j][i]*unit_vec_z[j][k]
                            dot_product = x + y + z
                            angle = np.rad2deg(np.arccos(dot_product))
                            angles_csv.write("%.7f," %angle)
            angles_csv.write("\n") 
