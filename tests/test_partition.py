import unittest
import math
from pbpk.physiology import scale_physiology
from pbpk.partition import calculate_rodgers_rowland_kp

class TestPartitionCoefficients(unittest.TestCase):
    def setUp(self):
        self.phys = scale_physiology("Human", 73.0)

    def test_neutral_compound(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="neutral",
            log_p=2.82,
            fup=0.015,
            b2p=0.8,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)
            self.assertFalse(math.isnan(val))
        self.assertGreater(kp["Fat"], kp["Muscle"])

    def test_monoprotic_acid(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="monoprotic_acid",
            log_p=2.7,
            pka1=5.05,
            fup=0.01,
            b2p=0.55,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)
            self.assertFalse(math.isnan(val))

    def test_diprotic_acid(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="diprotic_acid",
            log_p=1.5,
            pka1=4.0,
            pka2=6.5,
            fup=0.05,
            b2p=0.7,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)

    def test_monoprotic_base(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="monoprotic_base",
            log_p=2.6,
            pka1=9.45,
            fup=0.13,
            b2p=1.15,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)
        self.assertGreater(kp["Lung"], 1.0)
        self.assertGreater(kp["Liver"], 5.0)

    def test_diprotic_base(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="diprotic_base",
            log_p=2.1,
            pka1=9.1,
            pka2=7.8,
            fup=0.2,
            b2p=1.3,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)

    def test_zwitterion(self):
        kp = calculate_rodgers_rowland_kp(
            compound_type="zwitterion",
            log_p=0.8,
            pka1=9.2,
            pka2=3.5,
            fup=0.4,
            b2p=0.9,
            species_data=self.phys
        )
        self.assertEqual(len(kp), 15)
        for organ, val in kp.items():
            self.assertGreater(val, 0.0)

    def test_invalid_inputs_strict_halt(self):
        with self.assertRaises(ValueError):
            calculate_rodgers_rowland_kp("unknown_type", 2.0)

        with self.assertRaises(ValueError):
            calculate_rodgers_rowland_kp("neutral", 2.0, fup=0.0)

        with self.assertRaises(ValueError):
            calculate_rodgers_rowland_kp("neutral", 2.0, fup=1.5)

        with self.assertRaises(ValueError):
            calculate_rodgers_rowland_kp("neutral", 2.0, b2p=-1.0)

        with self.assertRaises(ValueError):
            calculate_rodgers_rowland_kp("neutral", float("nan"))

        with self.assertRaises(ValueError):
            scale_physiology("Canine", 10.0)

        with self.assertRaises(ValueError):
            scale_physiology("Human", -5.0)

if __name__ == "__main__":
    unittest.main()
