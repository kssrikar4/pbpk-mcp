import math

HUMAN_ICRP_DATA = {
    "reference_weight_kg": 73.0,
    "reference_height_m": 1.76,
    "hematocrit": 0.47,
    "pH_plasma": 7.4,
    "pH_rbc": 7.22,
    "pH_intracellular": 7.0,
    "blood_cells": {
        "f_iw": 0.63,
        "f_nl": 0.0012,
        "f_np": 0.0033,
        "AP_mg_g": 0.57
    },
    "plasma": {
        "f_nl": 0.0032,
        "f_np": 0.0021
    },
    "vascular": {
        "ArterialBlood": 0.419080314,
        "VenousBlood": 0.963884723,
        "PortalVein": 1.0372237776
    },
    "organs": {
        "Bone": {
            "volume_L": 11.817828072,
            "qspec_L_min_kg": 0.027495,
            "f_ew": 0.1,
            "f_iw": 0.317,
            "f_nl": 0.0174,
            "f_np": 0.0016,
            "AP_mg_g": 0.67,
            "albumin_ratio": 0.1,
            "lipoprotein_ratio": 0.05
        },
        "Brain": {
            "volume_L": 1.5088,
            "qspec_L_min_kg": 0.51675,
            "f_ew": 0.162,
            "f_iw": 0.591,
            "f_nl": 0.0391,
            "f_np": 0.0015,
            "AP_mg_g": 0.4,
            "albumin_ratio": 0.048,
            "lipoprotein_ratio": 0.041
        },
        "Fat": {
            "volume_L": 14.868,
            "qspec_L_min_kg": 0.02184,
            "f_ew": 0.135,
            "f_iw": 0.009,
            "f_nl": 0.853,
            "f_np": 0.0016,
            "AP_mg_g": 0.4,
            "albumin_ratio": 0.049,
            "lipoprotein_ratio": 0.068
        },
        "Gonads": {
            "volume_L": 0.040268154912,
            "qspec_L_min_kg": 0.0806,
            "f_ew": 0.03,
            "f_iw": 0.83,
            "f_nl": 0.0048,
            "f_np": 0.0116,
            "AP_mg_g": 2.45,
            "albumin_ratio": 0.048,
            "lipoprotein_ratio": 0.041
        },
        "Heart": {
            "volume_L": 0.4171997376,
            "qspec_L_min_kg": 0.62335,
            "f_ew": 0.32,
            "f_iw": 0.248,
            "f_nl": 0.0135,
            "f_np": 0.0106,
            "AP_mg_g": 2.25,
            "albumin_ratio": 0.157,
            "lipoprotein_ratio": 0.16
        },
        "Kidney": {
            "volume_L": 0.43771998528,
            "qspec_L_min_kg": 3.02705,
            "f_ew": 0.273,
            "f_iw": 0.399,
            "f_nl": 0.0121,
            "f_np": 0.024,
            "AP_mg_g": 5.03,
            "albumin_ratio": 0.13,
            "lipoprotein_ratio": 0.137
        },
        "LargeIntestine": {
            "volume_L": 0.41217159744,
            "qspec_L_min_kg": 0.6305,
            "f_ew": 0.282,
            "f_iw": 0.456,
            "f_nl": 0.0375,
            "f_np": 0.0124,
            "AP_mg_g": 2.41,
            "albumin_ratio": 0.158,
            "lipoprotein_ratio": 0.141
        },
        "Liver": {
            "volume_L": 2.3578466712,
            "qspec_L_min_kg": 0.1794,
            "f_ew": 0.161,
            "f_iw": 0.481,
            "f_nl": 0.0135,
            "f_np": 0.0238,
            "AP_mg_g": 4.56,
            "albumin_ratio": 0.086,
            "lipoprotein_ratio": 0.161
        },
        "Lung": {
            "volume_L": 1.2943497288,
            "qspec_L_min_kg": 0.0,
            "f_ew": 0.336,
            "f_iw": 0.238,
            "f_nl": 0.0215,
            "f_np": 0.0123,
            "AP_mg_g": 3.91,
            "albumin_ratio": 0.212,
            "lipoprotein_ratio": 0.168
        },
        "Muscle": {
            "volume_L": 32.338641984,
            "qspec_L_min_kg": 0.03419,
            "f_ew": 0.079,
            "f_iw": 0.666,
            "f_nl": 0.022,
            "f_np": 0.0078,
            "AP_mg_g": 2.42,
            "albumin_ratio": 0.064,
            "lipoprotein_ratio": 0.059
        },
        "Pancreas": {
            "volume_L": 0.19026527664,
            "qspec_L_min_kg": 0.3419,
            "f_ew": 0.12,
            "f_iw": 0.521,
            "f_nl": 0.0403,
            "f_np": 0.009,
            "AP_mg_g": 1.67,
            "albumin_ratio": 0.06,
            "lipoprotein_ratio": 0.06
        },
        "Skin": {
            "volume_L": 3.760912944,
            "qspec_L_min_kg": 0.08645,
            "f_ew": 0.382,
            "f_iw": 0.276,
            "f_nl": 0.0603,
            "f_np": 0.0044,
            "AP_mg_g": 1.32,
            "albumin_ratio": 0.277,
            "lipoprotein_ratio": 0.096
        },
        "SmallIntestine": {
            "volume_L": 0.7240974816,
            "qspec_L_min_kg": 0.89765,
            "f_ew": 0.282,
            "f_iw": 0.456,
            "f_nl": 0.0375,
            "f_np": 0.0124,
            "AP_mg_g": 2.41,
            "albumin_ratio": 0.158,
            "lipoprotein_ratio": 0.141
        },
        "Spleen": {
            "volume_L": 0.24341181216,
            "qspec_L_min_kg": 0.801125,
            "f_ew": 0.207,
            "f_iw": 0.355,
            "f_nl": 0.0071,
            "f_np": 0.0107,
            "AP_mg_g": 3.18,
            "albumin_ratio": 0.097,
            "lipoprotein_ratio": 0.207
        },
        "Stomach": {
            "volume_L": 0.16847666928,
            "qspec_L_min_kg": 0.3861,
            "f_ew": 0.282,
            "f_iw": 0.456,
            "f_nl": 0.0375,
            "f_np": 0.0124,
            "AP_mg_g": 2.41,
            "albumin_ratio": 0.158,
            "lipoprotein_ratio": 0.141
        }
    }
}

def get_human_physiology():
    return HUMAN_ICRP_DATA

def scale_physiology(species, weight_kg):
    if not isinstance(species, str) or species.lower() != "human":
        raise ValueError("Unsupported species: '" + str(species) + "'. Only 'Human' is supported.")
    if not isinstance(weight_kg, (int, float)) or weight_kg <= 0.0 or weight_kg > 500.0 or math.isnan(weight_kg):
        raise ValueError("Invalid body weight: " + str(weight_kg) + ". Weight must be positive and within physiological range.")

    ref_weight = HUMAN_ICRP_DATA["reference_weight_kg"]
    weight_ratio = weight_kg / ref_weight

    scaled_volumes = {}
    scaled_flows_l_h = {}

    for organ, data in HUMAN_ICRP_DATA["organs"].items():
        v = data["volume_L"] * weight_ratio
        scaled_volumes[organ] = v
        q_l_min = data["qspec_L_min_kg"] * v
        scaled_flows_l_h[organ] = q_l_min * 60.0

    v_art = HUMAN_ICRP_DATA["vascular"]["ArterialBlood"] * weight_ratio
    v_ven = HUMAN_ICRP_DATA["vascular"]["VenousBlood"] * weight_ratio
    v_pv = HUMAN_ICRP_DATA["vascular"]["PortalVein"] * weight_ratio

    scaled_volumes["ArterialBlood"] = v_art
    scaled_volumes["VenousBlood"] = v_ven
    scaled_volumes["PortalVein"] = v_pv

    q_pv = (scaled_flows_l_h["Pancreas"] +
            scaled_flows_l_h["Spleen"] +
            scaled_flows_l_h["Stomach"] +
            scaled_flows_l_h["SmallIntestine"] +
            scaled_flows_l_h["LargeIntestine"])

    q_ha = scaled_flows_l_h["Liver"]
    q_liver_total = q_ha + q_pv
    scaled_flows_l_h["PortalVein"] = q_pv
    scaled_flows_l_h["HepaticArtery"] = q_ha
    scaled_flows_l_h["LiverTotal"] = q_liver_total

    q_co = (scaled_flows_l_h["Bone"] +
            scaled_flows_l_h["Brain"] +
            scaled_flows_l_h["Fat"] +
            scaled_flows_l_h["Gonads"] +
            scaled_flows_l_h["Heart"] +
            scaled_flows_l_h["Kidney"] +
            scaled_flows_l_h["Muscle"] +
            scaled_flows_l_h["Skin"] +
            q_liver_total)

    scaled_flows_l_h["CardiacOutput"] = q_co
    scaled_flows_l_h["Lung"] = q_co

    return {
        "species": "Human",
        "weight_kg": float(weight_kg),
        "hematocrit": HUMAN_ICRP_DATA["hematocrit"],
        "pH_plasma": HUMAN_ICRP_DATA["pH_plasma"],
        "pH_rbc": HUMAN_ICRP_DATA["pH_rbc"],
        "pH_intracellular": HUMAN_ICRP_DATA["pH_intracellular"],
        "blood_cells": HUMAN_ICRP_DATA["blood_cells"],
        "plasma": HUMAN_ICRP_DATA["plasma"],
        "volumes": scaled_volumes,
        "flows_l_h": scaled_flows_l_h,
        "organs": HUMAN_ICRP_DATA["organs"]
    }
