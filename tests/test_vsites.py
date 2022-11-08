import martini_openmm.vsites as vs
import unittest


class TestVSiteManager(unittest.TestCase):
    def setUp(self):
        self.vsites = vs.VSiteManager()

    def test_simple_dependency(self):
        site1 = vs.LinearSite({0: 1.0})
        site2 = vs.LinearSite({1: 1.0})
        self.vsites.add(1, site1)
        self.vsites.add(2, site2)

        sites = list(self.vsites.iter())
        self.assertEqual(len(sites), 2)

        ind1, result1 = sites[0]
        ind2, result2 = sites[1]

        # Check that the indices match
        self.assertEqual(ind1, 1)
        self.assertEqual(ind2, 2)

        # Check that the atoms and weights are correct
        self.assertEqual(len(result1.atom_weights), 1)
        self.assertAlmostEqual(result1.atom_weights[0], 1.0)
        self.assertEqual(len(result2.atom_weights), 1)
        self.assertAlmostEqual(result2.atom_weights[0], 1.0)

    def test_two_particle_dependency(self):
        site2 = vs.LinearSite({0: 0.5, 1: 00.5})
        site3 = vs.LinearSite({2: 1.0})
        self.vsites.add(2, site2)
        self.vsites.add(3, site3)

        sites = list(self.vsites.iter())
        self.assertEqual(len(sites), 2)

        ind2, result2 = sites[0]
        ind3, result3 = sites[1]

        # Check that the indices match
        self.assertEqual(ind2, 2)
        self.assertEqual(ind3, 3)

        # Check that the atoms and weights are correct
        self.assertEqual(len(result2.atom_weights), 2)
        self.assertAlmostEqual(result2.atom_weights[0], 0.5)
        self.assertAlmostEqual(result2.atom_weights[1], 0.5)
        self.assertEqual(len(result3.atom_weights), 2)
        self.assertAlmostEqual(result2.atom_weights[0], 0.5)
        self.assertAlmostEqual(result2.atom_weights[1], 0.5)

    def test_complex_dependency(self):
        # particle 0
        # particle 1
        # particle 2
        # particle 3 = 0.5 p0 + 0.5 p1
        # particle 4 = 0.5 p3 + 0.5 p2
        #            = 0.25 p0 + 0.25 p1 + 0.5 p2

        site3 = vs.LinearSite({0: 0.5, 1: 0.5})
        site4 = vs.LinearSite({2: 0.5, 3: 0.5})
        self.vsites.add(3, site3)
        self.vsites.add(4, site4)

        sites = list(self.vsites.iter())
        self.assertEqual(len(sites), 2)

        ind3, result3 = sites[0]
        ind4, result4 = sites[1]

        # Check that the indices match
        self.assertEqual(ind3, 3)
        self.assertEqual(ind4, 4)

        # Check that the atoms and weights are correct
        self.assertEqual(len(result3.atom_weights), 2)
        self.assertAlmostEqual(result3.atom_weights[0], 0.5)
        self.assertAlmostEqual(result3.atom_weights[1], 0.5)
        self.assertEqual(len(result4.atom_weights), 3)
        self.assertAlmostEqual(result4.atom_weights[0], 0.25)
        self.assertAlmostEqual(result4.atom_weights[1], 0.25)
        self.assertAlmostEqual(result4.atom_weights[2], 0.5)
