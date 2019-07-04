"""
Class for Atomic Decomposition of the Reaction force
"""

__authors__  = "Wallace D. Derricotte"
__credits__  = ["Wallace D. Derricotte"]

__copyright__  = "(c) 2018-2019, Derricotte Research Group"
__license__    = "MIT License"
__date__       = "2019-01-02"

import numpy as np
import psi4
import os
import sys
import json
from pyscf import gto, scf, dft, grad, solvent, mp
from pyscf.solvent import ddcosmo, ddcosmo_grad

def grad_calc(params,current_geom, mol):
    """
        Uses Psi4/pySCF to calculate the energy gradient and returns the
        gradient and energy. Here any of the keywords the user provides in the
        .json input are used to set the options for the energy calculation.

        Parameters:
        ----------
            params(self) -- contains initialized shared parameters.
            current_geom(np array) -- Matrix of size natoms x 3 containing the geometry.
            mol(psi4.Molecule) -- Psi4 molecule object containing the current molecule.
        Returns:
        -------
            grad_mw(np array) -- Mass weighted gradient matrix of size natoms x 3.
            E(float) -- single-point energy from Psi4 calculation.
    """
    if(params.qm_program=='pyscf'):
        pymol = gto.Mole()
        pymol.verbose = 0
        geom_vec = []
        for i in range(params.natoms):
            atom = [params.symbols[i],]
            atom_coords = []
            for j in range(3):
                atom_coords.append(current_geom[i][j])
            atom_coords = tuple(atom_coords)
            atom.append(atom_coords)
            geom_vec.append(atom)
        #print(geom_vec)
        pymol.atom = geom_vec
        pymol.unit = 'Bohr'
        pymol.basis = params.basis
        pymol.charge = params.molecular_charge
        pymol.spin = params.molecular_multiplicity - 1
        pymol.build()
        if(params.method=="scf"):
            scf_obj = scf.RHF(pymol)
        if(params.method=="dft"):
            scf_obj = dft.RKS(pymol)
            scf_obj.xc = params.xc_functional
        if(params.do_solvent):
            solv_obj = ddcosmo.ddcosmo_for_scf(scf_obj)
            solv_obj.with_solvent.eps = params.eps
            solv_cosmo = solvent.ddCOSMO(solv_obj).run()
            E = solv_cosmo.scf()
            grad = solv_obj.nuc_grad_method().kernel()
        else:
            E = scf_obj.kernel()
            grad = scf_obj.nuc_grad_method().kernel()
    if(params.qm_program=='psi4'):
        mol.set_geometry(psi4.core.Matrix.from_array(current_geom))
        mol.fix_orientation(True)
        mol.fix_com(True)
        mol.reset_point_group('c1')
        grad_method = "%s/%s" %(params.method,params.basis)
        psi4.core.set_output_file("psi4_out.dat", False)
        psi4.set_options(params.keywords)
        psi4.set_num_threads(params.nthreads)
        E, wfn = psi4.energy(grad_method,return_wfn=True)
        psi4.set_num_threads(params.nthreads)
        grad = np.asarray(psi4.gradient(grad_method,ref_wfn=wfn))
        #print(grad)
    return grad

def displacement(geometries):
    displacements = []
    for i in range(len(geometries)):
        if(i<=len(geometries)-2):
            current_mol = psi4.geometry(geometries[i])
            next_mol = psi4.geometry(geometries[i+1])
            current_geom = np.asarray(current_mol.geometry())
            next_geom = np.asarray(next_mol.geometry())
            current_disp = np.subtract(next_geom,current_geom)
            displacements.append(current_disp)
        else:
            pass
    return displacements

def atomic_decomp(params, output_file, geometries):
    disp = displacement(geometries)
    del geometries[-1]
    atom_forces = []
    syms = params.symbols
    count = 0
    output = open(output_file,"a")
    output.write('\n\n--Atomic Contributions to Reaction Work (kcal/mol)--\n')
    for geometry in geometries:
        mol = psi4.geometry(geometry)
        geom_array = np.asarray(mol.geometry())
        grad = grad_calc(params, geom_array, mol)
        force_array = []
        output = open(output_file,"a")
        output.write('\n\nCoordinate %.2f\n' %(params.coordinates[count]))
        output.write('------------------------\n')
        for i in range(params.natoms):
            atom_force = grad[i][0]*disp[count][i][0] + grad[i][1]*disp[count][i][1] + grad[i][2]*disp[count][i][2]
            #print(atom_force)
            output.write('|%s%d|%.5f|   ' %(params.symbols[i], i+1, atom_force*627.509))
            force_array.append(atom_force)
        output.close()
        atom_forces.append(force_array)
        count = count + 1
    atomic_forces_by_atom = []
    for i in range(params.natoms):
        atom_array = []
        for j in range(len(geometries)):
           atom_array.append(atom_forces[j][i])
        atomic_forces_by_atom.append(atom_array)
    return atomic_forces_by_atom
