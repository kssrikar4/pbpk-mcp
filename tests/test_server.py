import unittest
import json
import os
import tempfile
from pbpk.server import (
    get_physiological_params,
    calculate_partition_coefficients,
    run_pbpk_simulation,
    generate_plots,
    generate_comparison_plots,
    compare_renal_impairment,
    estimate_compound_properties,
    run_smiles_pbpk_simulation,
    convert_clearance_units
)

class TestServerTools(unittest.TestCase):
    def test_get_physiological_params_success(self):
        resp_str = get_physiological_params(species="Human", weight_kg=73.0)
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertIn("data", resp)
        self.assertIn("volumes", resp["data"])
        self.assertIn("flows_l_h", resp["data"])

    def test_get_physiological_params_error(self):
        resp_str = get_physiological_params(species="Rat", weight_kg=0.25)
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "error")
        self.assertIn("error", resp)

    def test_calculate_partition_coefficients_success(self):
        resp_str = calculate_partition_coefficients(
            compound_type="monoprotic_base",
            log_p=2.6,
            pka1=9.45,
            fup=0.13,
            b2p=1.15
        )
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertIn("data", resp)
        self.assertEqual(len(resp["data"]), 15)

    def test_calculate_partition_coefficients_error(self):
        resp_str = calculate_partition_coefficients(
            compound_type="unknown_type",
            log_p=2.0
        )
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "error")
        self.assertIn("error", resp)

    def test_run_pbpk_simulation_individual(self):
        resp_str = run_pbpk_simulation(
            dose_mg=100.0,
            dose_route="iv_bolus",
            infusion_duration_h=0.0,
            simulation_time_h=12.0,
            clearance_hep_l_h=15.0,
            clearance_ren_l_h=5.0,
            log_p=2.6,
            pka1=9.45,
            compound_type="monoprotic_base",
            fup=0.13,
            b2p=1.15,
            weight_kg=73.0,
            num_subjects=1
        )
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertEqual(resp["type"], "individual")
        self.assertIn("data", resp)
        self.assertIn("summary", resp["data"])

    def test_run_pbpk_simulation_population(self):
        resp_str = run_pbpk_simulation(
            dose_mg=100.0,
            dose_route="iv_bolus",
            infusion_duration_h=0.0,
            simulation_time_h=6.0,
            clearance_hep_l_h=15.0,
            clearance_ren_l_h=5.0,
            log_p=2.6,
            pka1=9.45,
            compound_type="monoprotic_base",
            fup=0.13,
            b2p=1.15,
            weight_kg=73.0,
            num_subjects=3
        )
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertEqual(resp["type"], "population")
        self.assertIn("data", resp)
        self.assertEqual(resp["data"]["num_subjects"], 3)

    def test_generate_plots_tool(self):
        sim_resp_str = run_pbpk_simulation(
            dose_mg=50.0,
            dose_route="iv_bolus",
            infusion_duration_h=0.0,
            simulation_time_h=6.0,
            clearance_hep_l_h=10.0,
            clearance_ren_l_h=2.0,
            log_p=2.6,
            pka1=9.45,
            compound_type="monoprotic_base",
            fup=0.13,
            b2p=1.15,
            weight_kg=73.0,
            num_subjects=1
        )
        self.assertLess(len(sim_resp_str), 2000)
        parsed = json.loads(sim_resp_str)
        sim_id = parsed["simulation_id"]
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
            temp_path = tf.name
        try:
            plot_resp_str = generate_plots(
                simulation_json=sim_id,
                output_path=temp_path,
                target_organ="Liver"
            )
            resp = json.loads(plot_resp_str)
            self.assertEqual(resp["status"], "success")
            self.assertEqual(resp["output_path"], temp_path)
            self.assertTrue(os.path.exists(temp_path))
            self.assertGreater(os.path.getsize(temp_path), 500)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_simulation_error_remediation(self):
        err_resp_str = run_pbpk_simulation(
            dose_mg=50.0,
            dose_route="iv_infusion",
            infusion_duration_h=0.0,
            simulation_time_h=6.0,
            clearance_hep_l_h=10.0,
            clearance_ren_l_h=2.0,
            log_p=2.6,
            pka1=9.45,
            compound_type="monoprotic_base",
            fup=0.13,
            b2p=1.15
        )
        resp = json.loads(err_resp_str)
        self.assertEqual(resp["status"], "error")
        self.assertIn("error", resp)
        self.assertIn("remediation", resp)
        self.assertIn("infusion_duration_h", resp["remediation"])

    def test_compare_renal_impairment_tool(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            res_str = compare_renal_impairment(
                dose_mg=500.0,
                log_p=-1.43,
                pka1=12.4,
                compound_type="monoprotic_base",
                fup=1.0,
                b2p=1.0,
                healthy_cl_ren_l_h=30.0,
                simulation_time_h=24.0,
                output_dir=tmpdir
            )
            res = json.loads(res_str)
            self.assertEqual(res["status"], "success")
            self.assertEqual(len(res["table"]), 4)
            self.assertTrue(os.path.exists(res["plots"]["figure1_multicurve"]))
            self.assertTrue(os.path.exists(res["plots"]["figure2_auc_barchart"]))
            self.assertTrue(os.path.exists(res["plots"]["figure3_renal_vs_exposure"]))

    def test_generate_comparison_plots_tool(self):
        s1 = json.loads(run_pbpk_simulation(50.0, "iv_bolus", 0.0, 6.0, 10.0, 2.0, 2.6, 9.45, "monoprotic_base", 0.13, 1.15))["simulation_id"]
        s2 = json.loads(run_pbpk_simulation(50.0, "iv_bolus", 0.0, 6.0, 5.0, 1.0, 2.6, 9.45, "monoprotic_base", 0.13, 1.15))["simulation_id"]
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
            p1 = tf.name
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
            p2 = tf.name
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
            p3 = tf.name
        try:
            r1 = json.loads(generate_comparison_plots(f"{s1},{s2}", "Normal,Impaired", p1, "multicurve"))
            self.assertEqual(r1["status"], "success")
            r2 = json.loads(generate_comparison_plots(f"{s1},{s2}", "Normal,Impaired", p2, "auc_bar"))
            self.assertEqual(r2["status"], "success")
            r3 = json.loads(generate_comparison_plots(f"{s1},{s2}", "Normal,Impaired", p3, "renal_vs_auc", "100,50"))
            self.assertEqual(r3["status"], "success")
        finally:
            for p in [p1, p2, p3]:
                if os.path.exists(p):
                    os.remove(p)

    def test_estimate_compound_properties_tool(self):
        resp_str = estimate_compound_properties("CN(C)C(=N)NC(=N)N")
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertEqual(resp["data"]["compound_type"], "monoprotic_base")
        self.assertEqual(resp["data"]["pka1"], 12.4)

        err_str = estimate_compound_properties("INVALID_SMILES")
        err = json.loads(err_str)
        self.assertEqual(err["status"], "error")
        self.assertIn("remediation", err)

    def test_run_smiles_pbpk_simulation_tool(self):
        resp_str = run_smiles_pbpk_simulation("CN(C)C(=N)NC(=N)N", dose_mg=500.0, dose_route="oral", simulation_time_h=12.0)
        resp = json.loads(resp_str)
        self.assertEqual(resp["status"], "success")
        self.assertEqual(resp["type"], "individual")
        self.assertIn("summary", resp)

    def test_compare_renal_impairment_with_smiles_tool(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            res_str = compare_renal_impairment(
                smiles="CN(C)C(=N)NC(=N)N",
                dose_mg=500.0,
                simulation_time_h=12.0,
                output_dir=tmpdir
            )
            res = json.loads(res_str)
            self.assertEqual(res["status"], "success")
            self.assertEqual(len(res["table"]), 4)
            self.assertTrue(os.path.exists(res["plots"]["figure1_multicurve"]))
            self.assertTrue(os.path.exists(res["plots"]["figure2_auc_barchart"]))
            self.assertTrue(os.path.exists(res["plots"]["figure3_renal_vs_exposure"]))

    def test_convert_clearance_units_tool(self):
        r1 = json.loads(convert_clearance_units(6.0, "L/h", "mL/min"))
        self.assertEqual(r1["status"], "success")
        self.assertEqual(r1["converted_value"], 100.0)

        r2 = json.loads(convert_clearance_units(100.0, "mL/min", "L/h"))
        self.assertEqual(r2["status"], "success")
        self.assertEqual(r2["converted_value"], 6.0)

        r3 = json.loads(convert_clearance_units(100.0, "invalid", "L/h"))
        self.assertEqual(r3["status"], "error")

if __name__ == "__main__":
    unittest.main()
