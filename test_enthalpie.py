import unittest

from Enthalpie import ChemicalReaction, element_counts, parse_equation, main


class TestEnthalpie(unittest.TestCase):
    def test_parse(self):
        self.assertEqual(parse_equation("2 H2 + O2 -> 2 H2O"),
                         ({"H2": 2, "O2": 1}, {"H2O": 2}))

    def test_element_counts(self):
        self.assertEqual(element_counts("Ca(OH)2"), {"Ca": 1, "O": 2, "H": 2})
        self.assertEqual(element_counts("Fe2O3(s)"), {"Fe": 2, "O": 3})

    def test_water_formation(self):
        r = ChemicalReaction.from_equation("2 H2 + O2 -> 2 H2O(g)")
        self.assertTrue(r.is_balanced())
        self.assertAlmostEqual(r.calculate_delta_H(), -483.66, places=2)
        self.assertAlmostEqual(r.calculate_delta_S(), -88.83, places=2)
        self.assertLess(r.calculate_delta_G(298.15), 0)

    def test_inversion_temperature(self):
        r = ChemicalReaction.from_equation("CaCO3 -> CaO + CO2")
        self.assertAlmostEqual(r.inversion_temperature(), 1121.5, delta=0.5)
        self.assertGreater(r.calculate_delta_G(1000), 0)
        self.assertLess(r.calculate_delta_G(1200), 0)

    def test_unbalanced(self):
        r = ChemicalReaction.from_equation("H2 + O2 -> H2O(l)")
        self.assertEqual(r.imbalance(), {"O": -1})

    def test_errors(self):
        with self.assertRaises(KeyError):
            ChemicalReaction.from_equation("H2 + Xx2 -> H2O(l)")
        with self.assertRaises(KeyError):
            ChemicalReaction.from_equation("H2 + O2 -> H2O")
        r = ChemicalReaction.from_equation("H2 + Cl2 -> 2 HCl")
        with self.assertRaises(ValueError):
            r.calculate_delta_G(0)

    def test_cli(self):
        self.assertEqual(main(["H2 + Cl2 -> 2 HCl", "-T", "500"]), 0)
        self.assertEqual(main(["inconnu -> H2"]), 1)


if __name__ == "__main__":
    unittest.main()
