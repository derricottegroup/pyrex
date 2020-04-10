from pydantic import BaseModel, Field
from typing import Dict, List, Optional 
import json

class MoleculeOptions(BaseModel):
    symbols: List[str]
    molecular_charge: int = 0
    molecular_multiplicity: int = 1
    geometry: Optional[list] = None 
    fragments: Optional[list] = None
    fragment_charges: Optional[list] = None
    fragment_mults: Optional[list] = None

class ModelOptions(BaseModel):
    method: str = "scf"
    basis: Optional[str] = "sto-3g"

class KeywordOptions(BaseModel):
    keywords: Optional[list] = None

class IrcOptions(BaseModel):
    direction: Optional[str] = None
    step_size: Optional[float] = None
    mode: Optional[int] = None
    normal_mode_file: Optional[str] = None
    grace_period: Optional[int] = None
    e_conv: Optional[float] = 1e-5

class SurfaceScanOptions(BaseModel):
    do_surf_scan: bool = True
    scan_type: Optional[str] = 'relaxed'
    constraint_type: Optional[str] = None
    constrained_atoms: Optional[list] = None
    constrained_values: Optional[list] = None

class FSaptOptions(BaseModel):
    do_fsapt: bool = True
    monomer_A_frags: Optional[list] = None
    monomer_B_frags: Optional[list] = None
    monomer_A_labels: Optional[list] = None
    monomer_B_labels: Optional[list] = None

class PyrexOptions(BaseModel):
    qm_program: str = "psi4"
    do_energy: bool = True
    irc_filename: str
    irc_stepsize: float
    xc_functional: Optional[str] = None
    nthreads: Optional[int] = 1
    do_solvent: Optional[bool] = None
    pcm_solvent: Optional[str] = None
    eps: Optional[float] = None
    do_conceptualdft: Optional[bool] = None
    do_fragility_spec: Optional[bool] = None
    single_connectivity_matrix: Optional[bool] = None
    do_supermolecular: Optional[bool] = None
    energy_read: Optional[str] = None
    force_max: Optional[float] = None
    force_min: Optional[float] = None
    restart: Optional[bool] = None
    do_sapt: Optional[bool] = None
    sapt_method: Optional[str] = None
    do_atomic: Optional[bool] = None
    do_polarization: Optional[bool] = None
    orca_header: Optional[str] = None
    orca_block: Optional[str] = None
    active_site: Optional[list] = None


class PyrexInput(BaseModel): 
    molecule: MoleculeOptions = Field(...)
    model: ModelOptions = Field(...)
    pyrex: PyrexOptions = Field(...)

def clean_empty(d):
    if not isinstance(d, (dict, list)):
        return d
    if isinstance(d, list):
        return [v for v in (clean_empty(v) for v in d) if(v or v == 0)]
    return {k: v for k, v in ((k, clean_empty(v)) for k, v in d.items()) if(v or v == 0)}

def build_inp(PyrexInput):
    d = PyrexInput.dict()
    inp_data = clean_empty(d)
    return inp_data

#TODO: Finish implementing pydantic type recognition! This is pretty cool,
# Merge this idea with the input reader. 

# molecule_options = MoleculeOptions(symbols=["H","O","H"])
# model_options = ModelOptions(method='scf',basis='3-21G')
# pyrex_options = PyrexOptions()

# new = PyrexInput(molecule=molecule_options,model=model_options,pyrex=pyrex_options)
# print(json.dumps(new.dict()))