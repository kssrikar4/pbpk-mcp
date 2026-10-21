import math
import numpy as np
import numba

COMPARTMENT_NAMES = [
    "VenousBlood",
    "Lung",
    "ArterialBlood",
    "Bone",
    "Brain",
    "Fat",
    "Gonads",
    "Heart",
    "Kidney",
    "Muscle",
    "Skin",
    "Pancreas",
    "Spleen",
    "Stomach",
    "SmallIntestine",
    "LargeIntestine",
    "Liver",
    "EliminatedHepatic",
    "EliminatedRenal",
    "GutLumen"
]

@numba.njit(fastmath=True)
def pbpk_rhs_jit(t, y, V, Q, Kp, CL_hep, CL_ren, ka, infusion_rate, infusion_end):
    dy = np.zeros(20, dtype=np.float64)

    c_ven = y[0] / V[0]
    c_art = y[2] / V[2]

    cv_lung = y[1] / (V[1] * Kp[0])
    cv_bone = y[3] / (V[3] * Kp[1])
    cv_brain = y[4] / (V[4] * Kp[2])
    cv_fat = y[5] / (V[5] * Kp[3])
    cv_gonads = y[6] / (V[6] * Kp[4])
    cv_heart = y[7] / (V[7] * Kp[5])
    cv_kidney = y[8] / (V[8] * Kp[6])
    cv_muscle = y[9] / (V[9] * Kp[7])
    cv_skin = y[10] / (V[10] * Kp[8])
    cv_pancreas = y[11] / (V[11] * Kp[9])
    cv_spleen = y[12] / (V[12] * Kp[10])
    cv_stomach = y[13] / (V[13] * Kp[11])
    cv_small_int = y[14] / (V[14] * Kp[12])
    cv_large_int = y[15] / (V[15] * Kp[13])
    cv_liver = y[16] / (V[16] * Kp[14])

    q_bone = Q[0]
    q_brain = Q[1]
    q_fat = Q[2]
    q_gonads = Q[3]
    q_heart = Q[4]
    q_kidney = Q[5]
    q_muscle = Q[6]
    q_skin = Q[7]
    q_pancreas = Q[8]
    q_spleen = Q[9]
    q_stomach = Q[10]
    q_small_int = Q[11]
    q_large_int = Q[12]
    q_ha = Q[13]
    q_pv = Q[14]
    q_liv = Q[15]
    q_co = Q[16]

    c_pv = (q_pancreas * cv_pancreas +
            q_spleen * cv_spleen +
            q_stomach * cv_stomach +
            q_small_int * cv_small_int +
            q_large_int * cv_large_int) / q_pv

    r_hep = CL_hep * cv_liver
    r_ren = CL_ren * c_art

    r_oral = ka * y[19]

    curr_infusion = infusion_rate if t <= infusion_end else 0.0

    r_ven_return = (q_bone * cv_bone +
                    q_brain * cv_brain +
                    q_fat * cv_fat +
                    q_gonads * cv_gonads +
                    q_heart * cv_heart +
                    q_kidney * cv_kidney +
                    q_muscle * cv_muscle +
                    q_skin * cv_skin +
                    q_liv * cv_liver)

    dy[0] = r_ven_return - q_co * c_ven + curr_infusion
    dy[1] = q_co * (c_ven - cv_lung)
    dy[2] = q_co * (cv_lung - c_art)

    dy[3] = q_bone * (c_art - cv_bone)
    dy[4] = q_brain * (c_art - cv_brain)
    dy[5] = q_fat * (c_art - cv_fat)
    dy[6] = q_gonads * (c_art - cv_gonads)
    dy[7] = q_heart * (c_art - cv_heart)
    dy[8] = q_kidney * (c_art - cv_kidney) - r_ren
    dy[9] = q_muscle * (c_art - cv_muscle)
    dy[10] = q_skin * (c_art - cv_skin)

    dy[11] = q_pancreas * (c_art - cv_pancreas)
    dy[12] = q_spleen * (c_art - cv_spleen)
    dy[13] = q_stomach * (c_art - cv_stomach)
    dy[14] = q_small_int * (c_art - cv_small_int) + r_oral
    dy[15] = q_large_int * (c_art - cv_large_int)

    dy[16] = q_ha * c_art + q_pv * c_pv - q_liv * cv_liver - r_hep

    dy[17] = r_hep
    dy[18] = r_ren
    dy[19] = -r_oral

    return dy

class WholeBodyPBPKModel:
    def __init__(self, phys_params, kp_dict):
        self.phys = phys_params
        self.kp = kp_dict

        v = phys_params["volumes"]
        q = phys_params["flows_l_h"]

        self.V = np.array([
            v["VenousBlood"],
            v["Lung"],
            v["ArterialBlood"],
            v["Bone"],
            v["Brain"],
            v["Fat"],
            v["Gonads"],
            v["Heart"],
            v["Kidney"],
            v["Muscle"],
            v["Skin"],
            v["Pancreas"],
            v["Spleen"],
            v["Stomach"],
            v["SmallIntestine"],
            v["LargeIntestine"],
            v["Liver"],
            1.0,
            1.0,
            1.0
        ], dtype=np.float64)

        if np.any(self.V[:17] <= 0.0) or np.any(np.isnan(self.V)):
            raise RuntimeError("Encountered non-positive or NaN organ volume.")

        self.Q = np.array([
            q["Bone"],
            q["Brain"],
            q["Fat"],
            q["Gonads"],
            q["Heart"],
            q["Kidney"],
            q["Muscle"],
            q["Skin"],
            q["Pancreas"],
            q["Spleen"],
            q["Stomach"],
            q["SmallIntestine"],
            q["LargeIntestine"],
            q["HepaticArtery"],
            q["PortalVein"],
            q["LiverTotal"],
            q["CardiacOutput"]
        ], dtype=np.float64)

        if np.any(self.Q <= 0.0) or np.any(np.isnan(self.Q)):
            raise RuntimeError("Encountered non-positive or NaN blood flow rate.")

        self.Kp = np.array([
            kp_dict["Lung"],
            kp_dict["Bone"],
            kp_dict["Brain"],
            kp_dict["Fat"],
            kp_dict["Gonads"],
            kp_dict["Heart"],
            kp_dict["Kidney"],
            kp_dict["Muscle"],
            kp_dict["Skin"],
            kp_dict["Pancreas"],
            kp_dict["Spleen"],
            kp_dict["Stomach"],
            kp_dict["SmallIntestine"],
            kp_dict["LargeIntestine"],
            kp_dict["Liver"]
        ], dtype=np.float64)

        if np.any(self.Kp <= 0.0) or np.any(np.isnan(self.Kp)):
            raise RuntimeError("Encountered non-positive or NaN partition coefficient.")

    def rhs(self, t, y, cl_hep, cl_ren, ka=1.0, infusion_rate=0.0, infusion_end=0.0):
        return pbpk_rhs_jit(t, y, self.V, self.Q, self.Kp, cl_hep, cl_ren, ka, infusion_rate, infusion_end)
