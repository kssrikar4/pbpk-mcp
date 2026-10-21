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
from pbpk.visualization import plot_pbpk_curves

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
        sim_data = None

        if ref_str in SIMULATION_STORE:
            sim_data = SIMULATION_STORE[ref_str]
        elif os.path.isfile(ref_str):
            with open(ref_str, "r") as f:
                sim_data = json.load(f)
        else:
            try:
                parsed = json.loads(ref_str)
                if isinstance(parsed, dict):
                    sid = parsed.get("simulation_id")
                    sfile = parsed.get("simulation_file")
                    if sid and sid in SIMULATION_STORE:
                        sim_data = SIMULATION_STORE[sid]
                    elif sfile and os.path.isfile(sfile):
                        with open(sfile, "r") as f:
                            sim_data = json.load(f)
                    elif "data" in parsed and isinstance(parsed["data"], dict) and "venous_plasma_mg_l" in parsed["data"]:
                        sim_data = parsed["data"]
                    elif "venous_plasma_mg_l" in parsed:
                        sim_data = parsed
            except Exception:
                pass

        if sim_data is None:
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

def main():
    mcp.run()

if __name__ == "__main__":
    main()
