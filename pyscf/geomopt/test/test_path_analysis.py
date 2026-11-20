# SPDX-License-Identifier: Apache-2.0
import unittest
import numpy as np
from pyscf.geomopt import path_analysis as geometry


class KnownValues(unittest.TestCase):
    def test_distances_and_angles(self):
        path = np.array([[[0, 0, 0], [3, 0, 0], [3, 4, 0], [6, 0, 0]]], dtype=float)
        np.testing.assert_allclose(geometry.distances(path, [(0, 2)]), [[5]])
        np.testing.assert_allclose(geometry.angles(path, [(0, 1, 2), (0, 1, 3)]), [[90, 180]])
        self.assertEqual(geometry.distances(path, []).shape, (1, 0))

    def test_signed_dihedrals_and_degeneracy(self):
        path = np.array([[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 0, 1]],
                         [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 0, -1]],
                         [[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0]],
                         [[-1, 0, 0], [0, 0, 0], [1, 0, 0], [1, 0, 1]]], dtype=float)
        actual = geometry.dihedrals(path, [(0, 1, 2, 3)])[:, 0]
        np.testing.assert_allclose(actual[:3], [90, -90, 180])
        self.assertTrue(np.isnan(actual[3]))
        path[0, 0] = path[0, 1]
        self.assertTrue(np.isnan(geometry.angles(path, [(0, 1, 2)])[0, 0]))

    def test_unwrap_preserves_gaps(self):
        values = np.array([[170], [-170], [np.nan], [-175], [175]])
        expected = [[170], [190], [np.nan], [-175], [-185]]
        np.testing.assert_allclose(geometry.unwrap_dihedrals(values), expected, equal_nan=True)

    def test_rigid_motion_and_reflection(self):
        original = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        path = np.array([original, original @ rotation + [4, -7, 2]])
        np.testing.assert_allclose(geometry.aligned_rmsd(path), 0, atol=1e-14)
        for function, indices in ((geometry.distances, [(0, 3)]),
                                  (geometry.angles, [(1, 0, 2)]),
                                  (geometry.dihedrals, [(1, 0, 2, 3)])):
            result = function(path, indices)
            np.testing.assert_allclose(result[0], result[1], atol=1e-12)
        reflected = original * [-1, 1, 1]
        self.assertGreater(geometry.aligned_rmsd([original, reflected])[1], .1)

    def test_graph_coordinates(self):
        graph = [(0, 1), (1, 2), (2, 3), (1, 0)]
        result = geometry.coordinates_from_bonds(4, graph)
        self.assertEqual(result['distance'], [(0, 1), (1, 2), (2, 3)])
        self.assertEqual(result['angle'], [(0, 1, 2), (1, 2, 3)])
        self.assertEqual(result['dihedral'], [(0, 1, 2, 3)])

    def test_invalid_atom_selections(self):
        path = np.zeros((2, 3, 3))
        for pairs in ([(0, 3)], [(0, 0)], [(0., 1.)]):
            with self.assertRaises(ValueError):
                geometry.distances(path, pairs)
        with self.assertRaises(ValueError):
            geometry.distances(path[0], [(0, 1)])


if __name__ == '__main__':
    unittest.main()
