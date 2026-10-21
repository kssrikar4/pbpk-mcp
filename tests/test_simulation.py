import unittest
import os
import tempfile
import numpy as np
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp
from pbpk.model import WholeBodyPBPKModel
from pbpk.simulation import PBPKSimulation, run_virtual_population
from pbpk.visualization import plot_pbpk_curves

class TestSimulation(unittest.TestCase):
    def setUp(self):
        self.phys = scale_physiology("Human", 73.0)
        self.kp_dict = calculate_rodgers_rowland_kp(
            compound_type="monoprotic_base",
            log_p=2.6,
            pka1=9.45,
            fup=0.13,
            b2p=1.15,
            species_data=self.phys
        )
        self.model = WholeBodyPBPKModel(self.phys, self.kp_dict)

    def test_individual_simulation_iv_bolus(self):
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=100.0,
            dose_route="iv_bolus",
            simulation_time_h=12.0,
            clearance_hep_l_h=15.0,
            clearance_ren_l_h=5.0
        )
        res = sim.run()
        self.assertIn("summary", res)
        self.assertGreater(res["summary"]["c_max_mg_l"], 0.0)
        self.assertAlmostEqual(res["summary"]["t_max_h"], 0.0, places=4)
        self.assertGreater(res["summary"]["auc_plasma_mg_h_l"], 0.0)
        self.assertEqual(len(res["time_h"]), 500)

    def test_individual_simulation_oral(self):
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=100.0,
            dose_route="oral",
            simulation_time_h=12.0,
            clearance_hep_l_h=15.0,
            clearance_ren_l_h=5.0,
            ka_h=1.5
        )
        res = sim.run()
        self.assertGreater(res["summary"]["t_max_h"], 0.0)
        self.assertGreater(res["summary"]["c_max_mg_l"], 0.0)

    def test_virtual_population_parallel(self):
        compound_params = {
            "compound_type": "monoprotic_base",
            "log_p": 2.6,
            "pka1": 9.45,
            "fup": 0.13,
            "b2p": 1.15
        }
        dose_params = {
            "dose_mg": 100.0,
            "dose_route": "iv_bolus",
            "simulation_time_h": 12.0,
            "clearance_hep_l_h": 15.0,
            "clearance_ren_l_h": 5.0,
            "num_time_points": 200
        }
        pop_res = run_virtual_population(
            compound_params=compound_params,
            dose_params=dose_params,
            num_subjects=4,
            mean_weight_kg=73.0,
            n_workers=2
        )
        self.assertEqual(pop_res["num_subjects"], 4)
        self.assertGreater(pop_res["mean_auc"], 0.0)
        p05 = np.array(pop_res["plasma_p05_mg_l"])
        median = np.array(pop_res["plasma_median_mg_l"])
        p95 = np.array(pop_res["plasma_p95_mg_l"])
        self.assertTrue(np.all(p05 <= median + 1e-9))
        self.assertTrue(np.all(median <= p95 + 1e-9))

    def test_visualization_generation(self):
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=50.0,
            dose_route="iv_bolus",
            simulation_time_h=8.0,
            clearance_hep_l_h=10.0,
            clearance_ren_l_h=2.0
        )
        res = sim.run()
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tf:
            temp_path = tf.name
        try:
            saved = plot_pbpk_curves(res, temp_path, target_organ="Liver")
            self.assertEqual(saved, temp_path)
            self.assertTrue(os.path.exists(temp_path))
            self.assertGreater(os.path.getsize(temp_path), 1000)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

if __name__ == "__main__":
    unittest.main()
