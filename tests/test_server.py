import unittest
import json
import os
import tempfile
from pbpk.server import (
    get_physiological_params,
    calculate_partition_coefficients,
    run_pbpk_simulation,
    generate_plots
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

if __name__ == "__main__":
    unittest.main()
