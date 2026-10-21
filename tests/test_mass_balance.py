import unittest
import numpy as np
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp
from pbpk.model import WholeBodyPBPKModel
from pbpk.simulation import PBPKSimulation

class TestMassBalance(unittest.TestCase):
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

    def test_mass_conservation_zero_clearance_iv_bolus(self):
        dose = 100.0
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=dose,
            dose_route="iv_bolus",
            simulation_time_h=24.0,
            clearance_hep_l_h=0.0,
            clearance_ren_l_h=0.0,
            num_time_points=400
        )
        res = sim.run()
        max_err = res["max_mass_balance_error_mg"]
        self.assertLess(max_err, 1e-5)
        totals = np.array(res["total_mass_mg"])
        self.assertTrue(np.allclose(totals, dose, atol=1e-5, rtol=1e-5))

    def test_mass_conservation_with_clearance_iv_bolus(self):
        dose = 150.0
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=dose,
            dose_route="iv_bolus",
            simulation_time_h=24.0,
            clearance_hep_l_h=12.5,
            clearance_ren_l_h=4.2,
            num_time_points=500
        )
        res = sim.run()
        max_err = res["max_mass_balance_error_mg"]
        self.assertLess(max_err, 1e-5)
        totals = np.array(res["total_mass_mg"])
        self.assertTrue(np.allclose(totals, dose, atol=1e-5, rtol=1e-5))

    def test_mass_conservation_iv_infusion(self):
        dose = 200.0
        infusion_dur = 2.0
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=dose,
            dose_route="iv_infusion",
            infusion_duration_h=infusion_dur,
            simulation_time_h=24.0,
            clearance_hep_l_h=10.0,
            clearance_ren_l_h=3.0,
            num_time_points=500
        )
        res = sim.run()
        max_err = res["max_mass_balance_error_mg"]
        self.assertLess(max_err, 1e-5)
        times = np.array(res["time_h"])
        expected = np.minimum(times, infusion_dur) * (dose / infusion_dur)
        totals = np.array(res["total_mass_mg"])
        self.assertTrue(np.allclose(totals, expected, atol=1e-5, rtol=1e-5))

    def test_mass_conservation_oral_route(self):
        dose = 50.0
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=dose,
            dose_route="oral",
            simulation_time_h=24.0,
            clearance_hep_l_h=8.0,
            clearance_ren_l_h=2.0,
            ka_h=1.2,
            num_time_points=500
        )
        res = sim.run()
        max_err = res["max_mass_balance_error_mg"]
        self.assertLess(max_err, 1e-5)
        totals = np.array(res["total_mass_mg"])
        self.assertTrue(np.allclose(totals, dose, atol=1e-5, rtol=1e-5))

if __name__ == "__main__":
    unittest.main()
