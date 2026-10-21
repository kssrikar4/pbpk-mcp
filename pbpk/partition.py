import math
from pbpk.physiology import get_human_physiology

VALID_COMPOUND_TYPES = {
    "neutral",
    "monoprotic_acid",
    "diprotic_acid",
    "triprotic_acid",
    "monoprotic_base",
    "diprotic_base",
    "triprotic_base",
    "zwitterion"
}

def calculate_rodgers_rowland_kp(compound_type, log_p, pka1=0.0, pka2=0.0, pka3=0.0, fup=1.0, b2p=1.0, species_data=None):
    if not isinstance(compound_type, str) or compound_type.lower() not in VALID_COMPOUND_TYPES:
        raise ValueError("Invalid compound_type: '" + str(compound_type) + "'. Valid types: " + str(sorted(list(VALID_COMPOUND_TYPES))))
    
    ctype = compound_type.lower()
    
    if not isinstance(log_p, (int, float)) or math.isnan(log_p):
        raise ValueError("Invalid log_p: " + str(log_p) + ". Must be a finite number.")
    log_p = float(log_p)
    
    if not isinstance(fup, (int, float)) or fup <= 0.0 or fup > 1.0 or math.isnan(fup):
        raise ValueError("Invalid fup: " + str(fup) + ". Fraction unbound in plasma must be in (0, 1].")
    fup = float(fup)
    
    if not isinstance(b2p, (int, float)) or b2p <= 0.0 or math.isnan(b2p):
        raise ValueError("Invalid b2p: " + str(b2p) + ". Blood to plasma ratio must be strictly positive.")
    b2p = float(b2p)
    
    for pka_name, pka_val in [("pka1", pka1), ("pka2", pka2), ("pka3", pka3)]:
        if not isinstance(pka_val, (int, float)) or math.isnan(pka_val):
            raise ValueError("Invalid " + pka_name + ": " + str(pka_val) + ". Must be a finite number.")
    
    pka1 = float(pka1)
    pka2 = float(pka2)
    pka3 = float(pka3)

    if species_data is None:
        species_data = get_human_physiology()

    hct = species_data["hematocrit"]
    ph_pls = species_data["pH_plasma"]
    ph_rbc = species_data["pH_rbc"]
    ph_iw = species_data["pH_intracellular"]

    f_iw_bc = species_data["blood_cells"]["f_iw"]
    f_nl_bc = species_data["blood_cells"]["f_nl"]
    f_np_bc = species_data["blood_cells"]["f_np"]
    ap_bc = species_data["blood_cells"]["AP_mg_g"]

    f_nl_p = species_data["plasma"]["f_nl"]
    f_np_p = species_data["plasma"]["f_np"]

    is_neutral = False
    is_strong_base = False

    if ctype == "neutral":
        x = 0.0
        y = 0.0
        z = 1.0
        is_neutral = True
    elif ctype == "monoprotic_acid":
        x = 10.0 ** (ph_iw - pka1)
        y = 10.0 ** (ph_pls - pka1)
        z = 1.0
    elif ctype == "diprotic_acid":
        x = 10.0 ** (ph_iw - pka1) + 10.0 ** (2.0 * ph_iw - pka1 - pka2)
        y = 10.0 ** (ph_pls - pka1) + 10.0 ** (2.0 * ph_pls - pka1 - pka2)
        z = 1.0
    elif ctype == "triprotic_acid":
        x = 10.0 ** (ph_iw - pka1) + 10.0 ** (2.0 * ph_iw - pka1 - pka2) + 10.0 ** (3.0 * ph_iw - pka1 - pka2 - pka3)
        y = 10.0 ** (ph_pls - pka1) + 10.0 ** (2.0 * ph_pls - pka1 - pka2) + 10.0 ** (3.0 * ph_pls - pka1 - pka2 - pka3)
        z = 1.0
    elif ctype == "monoprotic_base":
        x = 10.0 ** (pka1 - ph_iw)
        y = 10.0 ** (pka1 - ph_pls)
        z = 10.0 ** (pka1 - ph_rbc)
        is_strong_base = pka1 >= 7.0
    elif ctype == "diprotic_base":
        x = 10.0 ** (pka2 - ph_iw) + 10.0 ** (pka1 + pka2 - 2.0 * ph_iw)
        y = 10.0 ** (pka2 - ph_pls) + 10.0 ** (pka1 + pka2 - 2.0 * ph_pls)
        z = 10.0 ** (pka2 - ph_rbc) + 10.0 ** (pka1 + pka2 - 2.0 * ph_rbc)
        is_strong_base = max(pka1, pka2) >= 7.0
    elif ctype == "triprotic_base":
        x = 10.0 ** (pka3 - ph_iw) + 10.0 ** (pka3 + pka2 - 2.0 * ph_iw) + 10.0 ** (pka1 + pka2 + pka3 - 3.0 * ph_iw)
        y = 10.0 ** (pka3 - ph_pls) + 10.0 ** (pka3 + pka2 - 2.0 * ph_pls) + 10.0 ** (pka1 + pka2 + pka3 - 3.0 * ph_pls)
        z = 10.0 ** (pka3 - ph_rbc) + 10.0 ** (pka3 + pka2 - 2.0 * ph_rbc) + 10.0 ** (pka1 + pka2 + pka3 - 3.0 * ph_rbc)
        is_strong_base = max(pka1, pka2, pka3) >= 7.0
    elif ctype == "zwitterion":
        x = 10.0 ** (pka1 - ph_iw) + 10.0 ** (ph_iw - pka2)
        y = 10.0 ** (pka1 - ph_pls) + 10.0 ** (ph_pls - pka2)
        z = 10.0 ** (pka1 - ph_rbc) + 10.0 ** (ph_rbc - pka2)
        is_strong_base = pka1 >= 7.0

    p_val = 10.0 ** log_p

    if is_strong_base:
        if z <= 0.0:
            raise RuntimeError("RBC ionization factor Z must be positive.")
        kpu_bc = (hct - 1.0 + b2p) / (fup * hct)
        lipid_term_bc = (p_val * f_nl_bc + (0.3 * p_val + 0.7) * f_np_bc) / (1.0 + y)
        term_ap = kpu_bc - (1.0 + z) / (1.0 + y) * f_iw_bc - lipid_term_bc
        ka_ap = max(0.0, term_ap) * (1.0 + y) / (ap_bc * z)
        ka_pr = 0.0
    else:
        ka_ap = 0.0
        lipid_term_p = (p_val * f_nl_p + (0.3 * p_val + 0.7) * f_np_p) / (1.0 + y)
        ka_pr = max(0.0, 1.0 / fup - 1.0 - lipid_term_p)

    kp_dict = {}

    for organ_name, organ_data in species_data["organs"].items():
        few = organ_data["f_ew"]
        fiw = organ_data["f_iw"]
        fnl = organ_data["f_nl"]
        fnp = organ_data["f_np"]
        apt = organ_data["AP_mg_g"]
        alb_ratio = organ_data["albumin_ratio"]
        lp_ratio = organ_data["lipoprotein_ratio"]

        if organ_name == "Fat":
            log_p_ow = 1.115 * log_p - 1.35
            p_eff = 10.0 ** log_p_ow
        else:
            p_eff = p_val

        lipid_term = (p_eff * fnl + (0.3 * p_eff + 0.7) * fnp) / (1.0 + y)

        if is_strong_base:
            binding_term = (ka_ap * apt * x) / (1.0 + y)
        elif is_neutral:
            binding_term = ka_pr * lp_ratio
        else:
            binding_term = ka_pr * alb_ratio

        kp = ((1.0 + x) / (1.0 + y) * fiw + few + lipid_term + binding_term) * fup

        if kp <= 0.0 or math.isnan(kp) or math.isinf(kp):
            raise RuntimeError("Invalid partition coefficient calculated for " + organ_name + ": " + str(kp))

        kp_dict[organ_name] = float(kp)

    return kp_dict
