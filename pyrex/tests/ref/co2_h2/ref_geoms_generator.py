import input_reader
import json
from geomparser import *

# General Name for Output and Reading in Common Parameters 
output_filename = "pyrex_output.dat"
json_input = "input.json"
json_data=open(json_input).read()
params = input_reader.Params(json_data)
json_output = open("ref_data_geom.json", "w+")
data = {}
data["ref_geoms"] = {}
geomparser = Geomparser(params.natoms, params.molecular_charge, params.molecular_multiplicity, params.geometries, params.coordinates)


# Generating total geometries
def generate_psi4_geoms():
    geoms = geomparser.geombuilder()
    data["ref_geoms"]["total_geoms_psi4"] = geoms

def generate_orca_geoms():
    geoms = geomparser.orca_geombuilder()
    data["ref_geoms"]["total_geoms_orca"] = geoms

def generate_pyscf_geoms():
    geoms = geomparser.pyscf_geombuilder()
    data["ref_geoms"]["total_geoms_pyscf"] = geoms

def generate_iso_frags():
    iso_frag_A = geomparser.iso_frag(params.charge_A, params.mult_A, params.frag_A)
    iso_frag_B = geomparser.iso_frag(params.charge_B, params.mult_B, params.frag_B)
    data["ref_geoms"]["iso_frag_A"] = iso_frag_A
    data["ref_geoms"]["iso_frag_B"] = iso_frag_B

def generate_fragghost():
    frag_A_ghost = geomparser.frag_ghost(params.charge_A, params.mult_A, params.frag_A)
    frag_B_ghost = geomparser.frag_ghost(params.charge_B, params.mult_B, params.frag_B)
    data["ref_geoms"]["frag_A_ghost"] = frag_A_ghost
    data["ref_geoms"]["frag_B_ghost"] = frag_B_ghost


# Call all functions to generate reference data
generate_psi4_geoms()
generate_orca_geoms()
generate_pyscf_geoms()
generate_iso_frags()
generate_fragghost()
json.dump(data, json_output, sort_keys=True, indent=4)
