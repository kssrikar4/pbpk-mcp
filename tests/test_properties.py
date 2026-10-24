import unittest
from pbpk.properties import estimate_parameters_from_smiles

class TestPropertiesEstimation(unittest.TestCase):
    def test_metformin_smiles(self):
        props = estimate_parameters_from_smiles("CN(C)C(=N)NC(=N)N")
        self.assertEqual(props["compound_type"], "monoprotic_base")
        self.assertEqual(props["pka1"], 12.4)
        self.assertLess(props["log_p"], 0.0)
        self.assertGreater(props["molecular_weight"], 120.0)
        self.assertLess(props["molecular_weight"], 135.0)
        self.assertEqual(props["clearance_hep_l_h"], 0.0)
        self.assertGreater(props["clearance_ren_l_h"], 10.0)
        self.assertGreater(props["clearance_ren_ml_min"], 100.0)
        self.assertTrue(props["renal_tubular_secretion"])
        self.assertGreater(props["fup"], 0.5)

    def test_propranolol_smiles(self):
        props = estimate_parameters_from_smiles("CC(C)NCC(O)COc1cccc2ccccc12")
        self.assertEqual(props["compound_type"], "monoprotic_base")
        self.assertEqual(props["pka1"], 9.45)
        self.assertGreater(props["log_p"], 1.0)
        self.assertGreater(props["clearance_hep_l_h"], 5.0)
        self.assertGreater(props["clearance_ren_l_h"], 0.0)

    def test_neutral_compound(self):
        props = estimate_parameters_from_smiles("c1ccccc1")
        self.assertEqual(props["compound_type"], "neutral")
        self.assertEqual(props["pka1"], 0.0)

    def test_acid_compound(self):
        props = estimate_parameters_from_smiles("CC(=O)O")
        self.assertEqual(props["compound_type"], "monoprotic_acid")
        self.assertEqual(props["pka1"], 4.2)

    def test_zwitterion_compound(self):
        props = estimate_parameters_from_smiles("NCC(=O)O")
        self.assertEqual(props["compound_type"], "zwitterion")
        self.assertEqual(props["pka1"], 9.5)
        self.assertEqual(props["pka2"], 4.0)

    def test_invalid_smiles(self):
        with self.assertRaises(ValueError):
            estimate_parameters_from_smiles("NOT_A_SMILES")

    def test_empty_smiles(self):
        with self.assertRaises(ValueError):
            estimate_parameters_from_smiles("")

if __name__ == "__main__":
    unittest.main()
