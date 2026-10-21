import unittest
import time
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp
from pbpk.model import WholeBodyPBPKModel
from pbpk.simulation import PBPKSimulation, run_virtual_population

class TestPerformance(unittest.TestCase):
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

    def test_rodgers_rowland_throughput(self):
        n_iters = 500
        t0 = time.perf_counter()
        for _ in range(n_iters):
            calculate_rodgers_rowland_kp(
                compound_type="monoprotic_base",
                log_p=2.6,
                pka1=9.45,
                fup=0.13,
                b2p=1.15,
                species_data=self.phys
            )
        t1 = time.perf_counter()
        elapsed = t1 - t0
        time_per_eval_ms = (elapsed / n_iters) * 1000.0
        self.assertLess(time_per_eval_ms, 2.0)

    def test_single_simulation_benchmark(self):
        sim = PBPKSimulation(
            model=self.model,
            dose_mg=100.0,
            dose_route="iv_bolus",
            simulation_time_h=24.0,
            clearance_hep_l_h=15.0,
            clearance_ren_l_h=5.0,
            num_time_points=250
        )
        sim.run()

        n_runs = 5
        t0 = time.perf_counter()
        for _ in range(n_runs):
            res = sim.run()
        t1 = time.perf_counter()
        avg_time_ms = ((t1 - t0) / n_runs) * 1000.0
        self.assertLess(avg_time_ms, 500.0)
        self.assertLess(res["max_mass_balance_error_mg"], 1e-5)

    def test_virtual_population_parallel_benchmark(self):
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
            "simulation_time_h": 24.0,
            "clearance_hep_l_h": 15.0,
            "clearance_ren_l_h": 5.0,
            "num_time_points": 200
        }
        t0 = time.perf_counter()
        pop_res = run_virtual_population(
            compound_params=compound_params,
            dose_params=dose_params,
            num_subjects=6,
            mean_weight_kg=73.0,
            n_workers=2
        )
        t1 = time.perf_counter()
        elapsed_s = t1 - t0
        self.assertEqual(pop_res["num_subjects"], 6)
        self.assertLess(elapsed_s, 10.0)

if __name__ == "__main__":
    unittest.main()
