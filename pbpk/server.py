import json
import os
import sys
import uuid
import tempfile
from fastmcp import FastMCP
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp
from pbpk.model import WholeBodyPBPKModel
from pbpk.simulation import PBPKSimulation, run_virtual_population
from pbpk.visualization import (
    plot_pbpk_curves,
    plot_multicurve_comparison,
    plot_auc_barchart,
    plot_renal_function_vs_exposure
)

mcp = FastMCP("PBPK-Engine")

SIMULATION_STORE = {}
SIM_CACHE_DIR = os.path.join(tempfile.gettempdir(), "pbpk_simulations")
os.makedirs(SIM_CACHE_DIR, exist_ok=True)

def _build_remediation(exc):
    msg = str(exc).lower()
    if "fup" in msg or "fraction unbound" in msg:
        return "Modify 'fup' to be a positive float between 0.0 and 1.0 (exclusive of 0, inclusive of 1, e.g., 0.1 for 10% unbound), then rerun."
    if "b2p" in msg or "blood to plasma" in msg:
        return "Modify 'b2p' to be a strictly positive float (greater than 0.0, e.g., 1.0), then rerun."
    if "weight" in msg:
        return "Modify 'weight_kg' to be a positive human body weight in kg (e.g., 70.0), then rerun."
    if "infusion_duration" in msg:
        return "For 'iv_infusion', set 'infusion_duration_h' to a positive value > 0.0 (e.g., 1.0 or 2.0). For bolus or oral, use 0.0, then rerun."
    if "dose_route" in msg:
        return "Modify 'dose_route' to one of: 'iv_bolus', 'iv_infusion', or 'oral', then rerun."
    if "compound_type" in msg:
        return "Modify 'compound_type' to one of: 'neutral', 'monoprotic_base', 'diprotic_base', 'triprotic_base', 'monoprotic_acid', 'diprotic_acid', 'triprotic_acid', 'zwitterion', then rerun."
    if "clearance" in msg:
        return "Ensure 'clearance_hep_l_h' and 'clearance_ren_l_h' are non-negative floats (>= 0.0), then rerun."
    if "dose_mg" in msg:
        return "Modify 'dose_mg' to be a strictly positive float (> 0.0), then rerun."
    return "Check parameter bounds and data types, adjust the invalid inputs, and rerun."

def _resolve_simulation_data(ref_str):
    s = ref_str.strip()
    if s in SIMULATION_STORE:
        return SIMULATION_STORE[s]
    if os.path.isfile(s):
        with open(s, "r") as f:
            return json.load(f)
    try:
        parsed = json.loads(s)
        if isinstance(parsed, dict):
            sid = parsed.get("simulation_id")
            sfile = parsed.get("simulation_file")
            if sid and sid in SIMULATION_STORE:
                return SIMULATION_STORE[sid]
            if sfile and os.path.isfile(sfile):
                with open(sfile, "r") as f:
                    return json.load(f)
            if "data" in parsed and isinstance(parsed["data"], dict) and "venous_plasma_mg_l" in parsed["data"]:
                return parsed["data"]
            if "venous_plasma_mg_l" in parsed:
                return parsed
    except Exception:
        pass
    return None

@mcp.tool()
def get_physiological_params(species: str, weight_kg: float) -> str:
    try:
        data = scale_physiology(species=species, weight_kg=weight_kg)
        return json.dumps({"status": "success", "data": data}, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

@mcp.tool()
def calculate_partition_coefficients(compound_type: str, log_p: float, pka1: float = 0.0, pka2: float = 0.0, fup: float = 1.0, b2p: float = 1.0) -> str:
    try:
        kp_dict = calculate_rodgers_rowland_kp(
            compound_type=compound_type,
            log_p=log_p,
            pka1=pka1,
            pka2=pka2,
            fup=fup,
            b2p=b2p
        )
        return json.dumps({"status": "success", "data": kp_dict}, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

@mcp.tool()
def run_pbpk_simulation(dose_mg: float, dose_route: str, infusion_duration_h: float, simulation_time_h: float, clearance_hep_l_h: float, clearance_ren_l_h: float, log_p: float, pka1: float, compound_type: str, fup: float, b2p: float, weight_kg: float = 73.0, num_subjects: int = 1) -> str:
    try:
        if num_subjects < 1:
            raise ValueError("num_subjects must be an integer >= 1.")

        compound_params = {
            "compound_type": compound_type,
            "log_p": log_p,
            "pka1": pka1,
            "fup": fup,
            "b2p": b2p
        }

        dose_params = {
            "dose_mg": dose_mg,
            "dose_route": dose_route,
            "infusion_duration_h": infusion_duration_h,
            "simulation_time_h": simulation_time_h,
            "clearance_hep_l_h": clearance_hep_l_h,
            "clearance_ren_l_h": clearance_ren_l_h
        }

        sim_id = "sim_" + uuid.uuid4().hex[:8]
        sim_file = os.path.join(SIM_CACHE_DIR, sim_id + ".json")

        if num_subjects == 1:
            phys = scale_physiology("Human", weight_kg)
            kp_dict = calculate_rodgers_rowland_kp(
                compound_type=compound_type,
                log_p=log_p,
                pka1=pka1,
                fup=fup,
                b2p=b2p,
                species_data=phys
            )
            model = WholeBodyPBPKModel(phys, kp_dict)
            sim = PBPKSimulation(
                model=model,
                dose_mg=dose_mg,
                dose_route=dose_route,
                infusion_duration_h=infusion_duration_h,
                simulation_time_h=simulation_time_h,
                clearance_hep_l_h=clearance_hep_l_h,
                clearance_ren_l_h=clearance_ren_l_h
            )
            res = sim.run()

            with open(sim_file, "w") as f:
                json.dump(res, f)
            SIMULATION_STORE[sim_id] = res

            lightweight_payload = {
                "simulation_id": sim_id,
                "simulation_file": sim_file,
                "summary": res["summary"],
                "max_mass_balance_error_mg": res["max_mass_balance_error_mg"]
            }

            return json.dumps({
                "status": "success",
                "type": "individual",
                "simulation_id": sim_id,
                "simulation_file": sim_file,
                "data": lightweight_payload,
                "summary": res["summary"],
                "max_mass_balance_error_mg": res["max_mass_balance_error_mg"]
            }, indent=2)
        else:
            pop_res = run_virtual_population(
                compound_params=compound_params,
                dose_params=dose_params,
                num_subjects=num_subjects,
                mean_weight_kg=weight_kg
            )

            with open(sim_file, "w") as f:
                json.dump(pop_res, f)
            SIMULATION_STORE[sim_id] = pop_res

            pop_summary = {
                "num_subjects": pop_res["num_subjects"],
                "mean_auc": pop_res["mean_auc"],
                "std_auc": pop_res["std_auc"],
                "mean_cmax": pop_res["mean_cmax"],
                "std_cmax": pop_res["std_cmax"]
            }

            lightweight_payload = {
                "simulation_id": sim_id,
                "simulation_file": sim_file,
                "num_subjects": pop_res["num_subjects"],
                "summary": pop_summary
            }

            return json.dumps({
                "status": "success",
                "type": "population",
                "simulation_id": sim_id,
                "simulation_file": sim_file,
                "data": lightweight_payload,
                "summary": pop_summary
            }, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

@mcp.tool()
def generate_plots(simulation_json: str, output_path: str, target_organ: str = "Liver") -> str:
    try:
        ref_str = simulation_json.strip()
        sim_data = _resolve_simulation_data(ref_str)

        if sim_data is None:
            if "," in ref_str:
                parts = [p.strip() for p in ref_str.split(",") if p.strip()]
                scenarios = []
                for idx, p in enumerate(parts):
                    d = _resolve_simulation_data(p)
                    if d is not None:
                        scenarios.append({"label": "Scenario " + str(idx + 1), "data": d})
                if scenarios:
                    saved_path = plot_multicurve_comparison(scenarios, output_path)
                    return json.dumps({"status": "success", "output_path": saved_path}, indent=2)

            raise ValueError("Could not resolve simulation data from provided simulation_json. Please pass a valid simulation_id, simulation_file, or JSON string.")

        saved_path = plot_pbpk_curves(sim_data, output_path, target_organ=target_organ)
        return json.dumps({"status": "success", "output_path": saved_path}, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

@mcp.tool()
def generate_comparison_plots(simulation_ids: str, labels: str, output_path: str, plot_type: str = "multicurve", renal_metrics: str = "") -> str:
    try:
        id_list = [i.strip() for i in simulation_ids.split(",") if i.strip()]
        label_list = [l.strip() for l in labels.split(",") if l.strip()]

        if len(id_list) != len(label_list):
            raise ValueError("The count of simulation_ids must match the count of labels.")

        resolved_sims = []
        aucs = []

        for sid in id_list:
            d = _resolve_simulation_data(sid)
            if d is None:
                raise ValueError("Could not find simulation for id: " + sid)
            resolved_sims.append(d)
            if "summary" in d and "auc_plasma_mg_h_l" in d["summary"]:
                aucs.append(d["summary"]["auc_plasma_mg_h_l"])
            elif "venous_plasma_mg_l" in d and "time_h" in d:
                import numpy as np
                aucs.append(float(np.trapezoid(d["venous_plasma_mg_l"], d["time_h"])))
            else:
                aucs.append(0.0)

        ptype = plot_type.lower()

        if ptype == "multicurve":
            scenarios = [{"label": label_list[i], "data": resolved_sims[i]} for i in range(len(label_list))]
            saved = plot_multicurve_comparison(scenarios, output_path)
        elif ptype == "auc_bar":
            saved = plot_auc_barchart(label_list, aucs, output_path)
        elif ptype in ("renal_vs_auc", "correlation"):
            if not renal_metrics.strip():
                raise ValueError("renal_metrics (comma-separated numbers) must be provided for renal_vs_auc plot.")
            r_vals = [float(v.strip()) for v in renal_metrics.split(",") if v.strip()]
            if len(r_vals) != len(aucs):
                raise ValueError("The count of renal_metrics must match the count of simulations.")
            saved = plot_renal_function_vs_exposure(r_vals, aucs, output_path)
        else:
            raise ValueError("Unsupported plot_type: " + str(plot_type) + ". Choose from: 'multicurve', 'auc_bar', 'renal_vs_auc'.")

        return json.dumps({"status": "success", "output_path": saved}, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

@mcp.tool()
def compare_renal_impairment(dose_mg: float, log_p: float, pka1: float, compound_type: str, fup: float, b2p: float, healthy_cl_ren_l_h: float, clearance_hep_l_h: float = 0.0, simulation_time_h: float = 48.0, ka_h: float = 0.5, output_dir: str = ".") -> str:
    try:
        os.makedirs(output_dir, exist_ok=True)

        scenarios = [
            ("Healthy (eGFR >= 90)", 1.0, 100.0),
            ("Mild Impairment (eGFR 60-89)", 0.75, 75.0),
            ("Moderate Impairment (eGFR 30-59)", 0.45, 45.0),
            ("Severe Impairment (eGFR 15-29)", 0.20, 20.0)
        ]

        phys = scale_physiology("Human", 73.0)
        kp_dict = calculate_rodgers_rowland_kp(
            compound_type=compound_type,
            log_p=log_p,
            pka1=pka1,
            fup=fup,
            b2p=b2p,
            species_data=phys
        )
        model = WholeBodyPBPKModel(phys, kp_dict)

        results_list = []
        sim_objects = []
        labels = []
        aucs = []
        egfrs = []

        for name, scale_factor, egfr in scenarios:
            cl_ren = healthy_cl_ren_l_h * scale_factor
            sim = PBPKSimulation(
                model=model,
                dose_mg=dose_mg,
                dose_route="oral",
                simulation_time_h=simulation_time_h,
                clearance_hep_l_h=clearance_hep_l_h,
                clearance_ren_l_h=cl_ren,
                ka_h=ka_h
            )
            res = sim.run()

            sid = "sim_renal_" + uuid.uuid4().hex[:6]
            SIMULATION_STORE[sid] = res
            sim_objects.append({"label": name, "data": res})
            labels.append(name)
            aucs.append(res["summary"]["auc_plasma_mg_h_l"])
            egfrs.append(egfr)

            results_list.append({
                "group": name,
                "egfr_ml_min": egfr,
                "renal_clearance_l_h": cl_ren,
                "c_max_mg_l": res["summary"]["c_max_mg_l"],
                "t_max_h": res["summary"]["t_max_h"],
                "auc_plasma_mg_h_l": res["summary"]["auc_plasma_mg_h_l"],
                "half_life_h": res["summary"]["half_life_h"]
            })

        p1 = os.path.join(output_dir, "figure1_plasma_overlay.png")
        p2 = os.path.join(output_dir, "figure2_auc_comparison.png")
        p3 = os.path.join(output_dir, "figure3_renal_vs_exposure.png")

        plot_multicurve_comparison(sim_objects, p1)
        plot_auc_barchart(labels, aucs, p2)
        plot_renal_function_vs_exposure(egfrs, aucs, p3, xlabel="eGFR (mL/min)", ylabel="Plasma AUC (mg·h/L)")

        return json.dumps({
            "status": "success",
            "study": "Renal Impairment Comparison",
            "table": results_list,
            "plots": {
                "figure1_multicurve": p1,
                "figure2_auc_barchart": p2,
                "figure3_renal_vs_exposure": p3
            }
        }, indent=2)
    except Exception as exc:
        return json.dumps({
            "status": "error",
            "error_type": exc.__class__.__name__,
            "error": str(exc),
            "remediation": _build_remediation(exc)
        }, indent=2)

def main():
    mcp.run()

if __name__ == "__main__":
    main()
