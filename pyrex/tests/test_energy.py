import pytest
import input_reader
import json
import numpy as np
from scf_class import *


output_filename = "pyrex_output.dat"
# Specifying Input Data for Molecule
inp_data = {
                  "molecule": {
                      "symbols" : ["C", "O", "O", "H", "H"],
                      "molecular_charge" : "0",
                      "molecular_multiplicity" : 1                                           
                      },                               
                  "model": {
                      "method": "scf",
                      "basis": "sto-3g"
                      },
                  "keywords": {
                      "scf_type": "df",
                      "freeze_core" : "True",
                      "reference" : "rhf",
                      "basis": "sto-3g"
                      },
                  "pyrex" : {
                      "irc_filename" : "./ref/co2_h2/full_irc.xyz",
                      "irc_stepsize" : 0.2
                   }
             }

# Store JSON input data and send parameters to input reader
json_input = json.dumps(inp_data)
params = input_reader.Params(json_input)

inp_ref_geom=open("./ref/ref_data_geom.json").read()
ref_data_geom = json.loads(inp_ref_geom)

inp_ref_energy=open("./ref/ref_data_energy.json").read()
ref_data_energy = json.loads(inp_ref_energy)



print("Testing Energy")
def test_psi4_energies():
    params.geoms = ref_data_geom["ref_geoms"]["total_geoms_psi4"]
    params.qm_program = "psi4"
    ref_energies = ref_data_energy["ref_energy"]["psi4_total_energies"]
    scf_instance = scf_class(params, output_filename)
    params.energies, params.wavefunctions = scf_instance.psi4_scf(params.geoms)
    nelec = params.wavefunctions[0].nalpha()
    diff = abs(np.asarray(params.energies) - np.asarray(ref_energies))
    assert max(diff) <= 1.0e-6

def test_pyscf_energies():
    params.geoms = ref_data_geom["ref_geoms"]["total_geoms_pyscf"]
    params.qm_program = "pyscf"
    ref_energies = ref_data_energy["ref_energy"]["pyscf_total_energies"]
    scf_instance = scf_class(params, output_filename)
    params.energies, params.frontier_orb_energies = scf_instance.pyscf_scf(params.geoms)
    diff = abs(np.asarray(params.energies) - np.asarray(ref_energies))
    assert max(diff) <= 1.0e-6

def test_orca_energies():
    params.geoms = ref_data_geom["ref_geoms"]["total_geoms_orca"]
    params.qm_program = "orca"
    params.method = "rhf"
    ref_energies = ref_data_energy["ref_energy"]["orca_total_energies"]
    scf_instance = scf_class(params, output_filename)
    params.energies, params.frontier_orb_energies = scf_instance.orca_scf(params.geoms)
    diff = abs(np.asarray(params.energies) - np.asarray(ref_energies))
    assert max(diff) <= 1.0e-6
