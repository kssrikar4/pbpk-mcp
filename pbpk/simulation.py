import math
import numpy as np
from scipy.integrate import solve_ivp
from concurrent.futures import ProcessPoolExecutor
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp
from pbpk.model import WholeBodyPBPKModel, COMPARTMENT_NAMES

VALID_ROUTES = {"iv_bolus", "iv_infusion", "oral"}

class PBPKSimulation:
    def __init__(self, model, dose_mg, dose_route="iv_bolus", infusion_duration_h=0.0, simulation_time_h=24.0, clearance_hep_l_h=0.0, clearance_ren_l_h=0.0, ka_h=1.0, num_time_points=500):
        if not isinstance(model, WholeBodyPBPKModel):
            raise TypeError("model must be an instance of WholeBodyPBPKModel.")
        if not isinstance(dose_mg, (int, float)) or dose_mg <= 0.0 or math.isnan(dose_mg):
            raise ValueError("dose_mg must be a strictly positive number.")
        if not isinstance(dose_route, str) or dose_route.lower() not in VALID_ROUTES:
            raise ValueError("dose_route must be one of: " + str(sorted(list(VALID_ROUTES))))
        if not isinstance(simulation_time_h, (int, float)) or simulation_time_h <= 0.0 or math.isnan(simulation_time_h):
            raise ValueError("simulation_time_h must be a strictly positive number.")
        if not isinstance(clearance_hep_l_h, (int, float)) or clearance_hep_l_h < 0.0 or math.isnan(clearance_hep_l_h):
            raise ValueError("clearance_hep_l_h must be non-negative.")
        if not isinstance(clearance_ren_l_h, (int, float)) or clearance_ren_l_h < 0.0 or math.isnan(clearance_ren_l_h):
            raise ValueError("clearance_ren_l_h must be non-negative.")
        if not isinstance(infusion_duration_h, (int, float)) or infusion_duration_h < 0.0 or math.isnan(infusion_duration_h):
            raise ValueError("infusion_duration_h must be non-negative.")
        if not isinstance(ka_h, (int, float)) or ka_h <= 0.0 or math.isnan(ka_h):
            raise ValueError("ka_h must be strictly positive.")
        if not isinstance(num_time_points, int) or num_time_points < 10:
            raise ValueError("num_time_points must be an integer >= 10.")

        self.model = model
        self.dose_mg = float(dose_mg)
        self.dose_route = dose_route.lower()
        self.simulation_time_h = float(simulation_time_h)
        self.clearance_hep_l_h = float(clearance_hep_l_h)
        self.clearance_ren_l_h = float(clearance_ren_l_h)
        self.infusion_duration_h = float(infusion_duration_h)
        self.ka_h = float(ka_h)
        self.num_time_points = num_time_points

        if self.dose_route == "iv_infusion" and self.infusion_duration_h <= 0.0:
            raise ValueError("infusion_duration_h must be > 0.0 for iv_infusion.")

    def run(self):
        y0 = np.zeros(20, dtype=np.float64)
        infusion_rate = 0.0
        infusion_end = 0.0

        if self.dose_route == "iv_bolus":
            y0[0] = self.dose_mg
        elif self.dose_route == "iv_infusion":
            infusion_rate = self.dose_mg / self.infusion_duration_h
            infusion_end = self.infusion_duration_h
        elif self.dose_route == "oral":
            y0[19] = self.dose_mg

        t_eval = np.linspace(0.0, self.simulation_time_h, self.num_time_points)

        def odefun(t, y):
            return self.model.rhs(t, y, self.clearance_hep_l_h, self.clearance_ren_l_h, self.ka_h, infusion_rate, infusion_end)

        sol = solve_ivp(
            fun=odefun,
            t_span=(0.0, self.simulation_time_h),
            y0=y0,
            method="Radau",
            t_eval=t_eval,
            rtol=1e-8,
            atol=1e-10
        )

        if not sol.success:
            raise RuntimeError("ODE integration failed: " + str(sol.message))

        times = sol.t
        y = sol.y

        total_mass = np.sum(y, axis=0)

        if self.dose_route == "iv_infusion":
            expected_dose_t = np.minimum(times, self.infusion_duration_h) * infusion_rate
        else:
            expected_dose_t = np.full_like(times, self.dose_mg)

        mass_balance_discrepancy = np.max(np.abs(total_mass - expected_dose_t))

        c_venous_plasma = y[0] / self.model.V[0]
        c_arterial_plasma = y[2] / self.model.V[2]

        organ_concentrations = {}
        organ_names_map = {
            "Lung": 1,
            "Bone": 3,
            "Brain": 4,
            "Fat": 5,
            "Gonads": 6,
            "Heart": 7,
            "Kidney": 8,
            "Muscle": 9,
            "Skin": 10,
            "Pancreas": 11,
            "Spleen": 12,
            "Stomach": 13,
            "SmallIntestine": 14,
            "LargeIntestine": 15,
            "Liver": 16
        }

        for org, idx in organ_names_map.items():
            organ_concentrations[org] = y[idx] / self.model.V[idx]

        c_max = float(np.max(c_venous_plasma))
        t_max = float(times[np.argmax(c_venous_plasma)])
        auc_plasma = float(np.trapezoid(c_venous_plasma, times))

        organ_aucs = {}
        for org, c_arr in organ_concentrations.items():
            organ_aucs[org] = float(np.trapezoid(c_arr, times))

        return {
            "time_h": times.tolist(),
            "venous_plasma_mg_l": c_venous_plasma.tolist(),
            "arterial_plasma_mg_l": c_arterial_plasma.tolist(),
            "organ_concentrations_mg_l": {k: v.tolist() for k, v in organ_concentrations.items()},
            "eliminated_hepatic_mg": y[17].tolist(),
            "eliminated_renal_mg": y[18].tolist(),
            "total_mass_mg": total_mass.tolist(),
            "max_mass_balance_error_mg": float(mass_balance_discrepancy),
            "summary": {
                "c_max_mg_l": c_max,
                "t_max_h": t_max,
                "auc_plasma_mg_h_l": auc_plasma,
                "organ_aucs_mg_h_l": organ_aucs
            }
        }

def _single_subject_simulation_worker(args):
    weight_kg, compound_params, dose_params = args

    scaled_phys = scale_physiology("Human", weight_kg)
    kp_dict = calculate_rodgers_rowland_kp(
        compound_type=compound_params["compound_type"],
        log_p=compound_params["log_p"],
        pka1=compound_params.get("pka1", 0.0),
        pka2=compound_params.get("pka2", 0.0),
        pka3=compound_params.get("pka3", 0.0),
        fup=compound_params.get("fup", 1.0),
        b2p=compound_params.get("b2p", 1.0),
        species_data=scaled_phys
    )

    model = WholeBodyPBPKModel(scaled_phys, kp_dict)

    weight_ratio = weight_kg / 73.0
    cl_hep = dose_params["clearance_hep_l_h"] * (weight_ratio ** 0.75)
    cl_ren = dose_params["clearance_ren_l_h"] * (weight_ratio ** 0.75)

    sim = PBPKSimulation(
        model=model,
        dose_mg=dose_params["dose_mg"],
        dose_route=dose_params.get("dose_route", "iv_bolus"),
        infusion_duration_h=dose_params.get("infusion_duration_h", 0.0),
        simulation_time_h=dose_params.get("simulation_time_h", 24.0),
        clearance_hep_l_h=cl_hep,
        clearance_ren_l_h=cl_ren,
        ka_h=dose_params.get("ka_h", 1.0),
        num_time_points=dose_params.get("num_time_points", 300)
    )

    res = sim.run()
    return {
        "weight_kg": weight_kg,
        "plasma_c": res["venous_plasma_mg_l"],
        "c_max": res["summary"]["c_max_mg_l"],
        "auc": res["summary"]["auc_plasma_mg_h_l"],
        "time_h": res["time_h"]
    }

def run_virtual_population(compound_params, dose_params, num_subjects=10, mean_weight_kg=73.0, weight_cv=0.15, seed=42, n_workers=None):
    if num_subjects <= 0:
        raise ValueError("num_subjects must be a positive integer.")
    if mean_weight_kg <= 0.0:
        raise ValueError("mean_weight_kg must be positive.")

    rng = np.random.default_rng(seed)
    sigma = math.sqrt(math.log(1.0 + weight_cv ** 2))
    mu = math.log(mean_weight_kg) - 0.5 * sigma ** 2
    weights = rng.lognormal(mean=mu, sigma=sigma, size=num_subjects)

    tasks = [(float(w), compound_params, dose_params) for w in weights]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(_single_subject_simulation_worker, tasks))

    all_plasma = np.array([r["plasma_c"] for r in results])
    times = results[0]["time_h"]

    median_profile = np.median(all_plasma, axis=0)
    p05_profile = np.percentile(all_plasma, 5, axis=0)
    p95_profile = np.percentile(all_plasma, 95, axis=0)

    aucs = [r["auc"] for r in results]
    cmaxs = [r["c_max"] for r in results]

    return {
        "num_subjects": num_subjects,
        "time_h": times,
        "plasma_median_mg_l": median_profile.tolist(),
        "plasma_p05_mg_l": p05_profile.tolist(),
        "plasma_p95_mg_l": p95_profile.tolist(),
        "mean_auc": float(np.mean(aucs)),
        "std_auc": float(np.std(aucs)),
        "mean_cmax": float(np.mean(cmaxs)),
        "std_cmax": float(np.std(cmaxs)),
        "individual_weights": [r["weight_kg"] for r in results]
    }
