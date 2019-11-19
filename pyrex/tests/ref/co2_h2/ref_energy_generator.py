import input_reader
import json
from scf_class import *

# General Name for Output and Reading in Common Parameters
output_filename = "pyrex_output.dat"
json_input = "input.json"
json_data=open(json_input).read()
params = input_reader.Params(json_data)
json_output = open("ref_data_energy.json", "w+")
data = {}
data["ref_energy"] = {}

inp_ref=open("../ref_data_geom.json").read()
ref_data = json.loads(inp_ref)

scf_instance = scf_class(params, output_filename)

def generate_psi4_energies():
    scf_instance = scf_class(params, output_filename)
    params.geoms = ref_data["ref_geoms"]["total_geoms_psi4"]
    params.energies, params.wavefunctions = scf_instance.psi4_scf(params.geoms)
    nelec = params.wavefunctions[0].nalpha()
    data["ref_energy"]["psi4_total_energies"] = params.energies

def generate_pyscf_energies():
    params.qm_program = "pyscf"
    params.geoms = ref_data["ref_geoms"]["total_geoms_pyscf"]
    scf_instance = scf_class(params, output_filename)
    params.energies, params.frontier_orb_energies = scf_instance.pyscf_scf(params.geoms)
    data["ref_energy"]["pyscf_total_energies"] = params.energies

def generate_orca_energies():
    params.qm_program = "orca"
    params.method = "rhf"
    params.geoms = ref_data["ref_geoms"]["total_geoms_orca"]
    scf_instance = scf_class(params, output_filename)
    params.energies, params.frontier_orb_energies = scf_instance.orca_scf(params.geoms)
    data["ref_energy"]["orca_total_energies"] = params.energies

generate_psi4_energies()
generate_pyscf_energies() 
generate_orca_energies()
json.dump(data, json_output, sort_keys=True, indent=4)
