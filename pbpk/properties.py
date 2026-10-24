import math
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def estimate_parameters_from_smiles(smiles):
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty string.")

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        raise ValueError("Invalid SMILES string: " + str(smiles))

    mw = float(Descriptors.MolWt(mol))
    log_p = float(Crippen.MolLogP(mol))
    tpsa = float(Descriptors.TPSA(mol))

    pat_guanidine = Chem.MolFromSmarts("[NX3,NX2]=C([NX3])[NX3]")
    pat_acid = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    pat_amine = Chem.MolFromSmarts("[NX3;!$(NC=O);!$(N=O);!$(N=C)]")
    pat_aniline = Chem.MolFromSmarts("[NX3;!$(NC=O)][c]")
    pat_phenol = Chem.MolFromSmarts("[OX2H][c]")

    is_guanidine = mol.HasSubstructMatch(pat_guanidine) if pat_guanidine else False
    is_acid = mol.HasSubstructMatch(pat_acid) if pat_acid else False
    is_amine = mol.HasSubstructMatch(pat_amine) if pat_amine else False
    is_aniline = mol.HasSubstructMatch(pat_aniline) if pat_aniline else False
    is_phenol = mol.HasSubstructMatch(pat_phenol) if pat_phenol else False

    pka2 = 0.0

    if is_guanidine:
        compound_type = "monoprotic_base"
        pka1 = 12.4
    elif is_acid and (is_amine or is_guanidine):
        compound_type = "zwitterion"
        pka1 = 9.5
        pka2 = 4.0
    elif is_acid:
        compound_type = "monoprotic_acid"
        pka1 = 4.2
    elif is_amine:
        compound_type = "monoprotic_base"
        pka1 = 9.45
    elif is_aniline:
        compound_type = "monoprotic_base"
        pka1 = 4.6
    elif is_phenol:
        compound_type = "monoprotic_acid"
        pka1 = 9.8
    else:
        compound_type = "neutral"
        pka1 = 0.0

    log_affinity = 0.448 * log_p - 0.478
    fup = 1.0 / (1.0 + 10.0 ** log_affinity)
    fup = max(0.005, min(1.0, float(fup)))

    b2p = 1.0

    log_peff = 0.22 * log_p - 0.008 * tpsa - 0.003 * mw + 0.45
    ka_h = max(0.1, min(3.0, float(2.0 * (10.0 ** log_peff) / 1.75 * 0.36)))

    if log_p < 0.0:
        cl_hep = 0.0
    else:
        cl_int = 10.0 ** (0.45 * log_p + 0.5)
        cl_hep = float((95.0 * fup * cl_int) / (95.0 + fup * cl_int))

    cl_filt = 7.2 * fup
    if is_guanidine or (is_amine and log_p < 0.5):
        cl_ren = float(cl_filt + 24.0 * fup)
    else:
        cl_ren = float(cl_filt)

    has_secretion = bool(is_guanidine or (is_amine and log_p < 0.5))

    return {
        "smiles": smiles.strip(),
        "molecular_weight": mw,
        "log_p": log_p,
        "pka1": pka1,
        "pka2": pka2,
        "compound_type": compound_type,
        "fup": fup,
        "b2p": b2p,
        "ka_h": ka_h,
        "clearance_hep_l_h": cl_hep,
        "clearance_hep_ml_min": float(cl_hep * 1000.0 / 60.0),
        "clearance_ren_l_h": cl_ren,
        "clearance_ren_ml_min": float(cl_ren * 1000.0 / 60.0),
        "renal_mechanism": "filtration_and_tubular_secretion" if has_secretion else "filtration_only",
        "renal_tubular_secretion": has_secretion
    }
