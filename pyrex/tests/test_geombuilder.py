import pytest
import input_reader
import json
from geomparser import *

# Specifying Input Data for Molecule
inp_data = {
          "molecule": {
                "symbols" : ["C", "O", "O", "H", "H"],
                    "molecular_charge" : "0",
                    "molecular_multiplicity" : 1,
                    "fragments" : [[0,1,2],[3,4]],
                    "fragment_charges" : [0,0],
                    "fragment_multiplicities" : [1,1]
                    },
                "pyrex" : {
                    "irc_filename" : "./ref/co2_h2/full_irc.xyz",
                    "irc_stepsize" : 0.2
                    }
                }
# Store JSON input data and send parameters to input reader
json_input = json.dumps(inp_data)
params = input_reader.Params(json_input)
# Read in reference data for testing
json_data=open("./ref/ref_data_geom.json").read()
ref_data = json.loads(json_data)
# Initializing Geomparser class
geomparser = Geomparser(params.natoms, params.molecular_charge, params.molecular_multiplicity, params.geometries, params.coordinates)

# Test PSI4 geometry builder
def test_geombuilder_psi4():
    geoms = geomparser.geombuilder()
    assert geoms==ref_data["ref_geoms"]["total_geoms_psi4"]

# Test ORCA geometry builder
def test_geombuilder_orca():
    geoms = geomparser.orca_geombuilder()
    assert geoms==ref_data["ref_geoms"]["total_geoms_orca"]

# Test pySCF geometry builder
def test_geombuilder_pyscf():
    geoms = geomparser.pyscf_geombuilder()
    assert geoms==ref_data["ref_geoms"]["total_geoms_pyscf"]

# Test iso_frag function for isolating fragments
def test_isofrag_func():
    iso_frag_A = geomparser.iso_frag(params.charge_A, params.mult_A, params.frag_A)
    iso_frag_B = geomparser.iso_frag(params.charge_B, params.mult_B, params.frag_B)
    assert iso_frag_A==ref_data["ref_geoms"]["iso_frag_A"]
    assert iso_frag_B==ref_data["ref_geoms"]["iso_frag_B"]

# Test frag_ghost function for ghosting atoms
def test_fragghost_func():
    frag_A_ghost = geomparser.frag_ghost(params.charge_A, params.mult_A, params.frag_A)
    frag_B_ghost = geomparser.frag_ghost(params.charge_B, params.mult_B, params.frag_B)
    assert frag_A_ghost==ref_data["ref_geoms"]["frag_A_ghost"]
    assert frag_B_ghost==ref_data["ref_geoms"]["frag_B_ghost"]
