# Copyright 2026 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import importlib.util
import unittest
import numpy as np
from pyscf import gto, lib
from pyscf.geomopt.geodesic import interpolate


@unittest.skipIf(importlib.util.find_spec('geodesic_interpolate') is None,
                 'geodesic-interpolate library not found')
class KnownValues(unittest.TestCase):
    def setUp(self):
        self.random_state = np.random.get_state()
        np.random.seed(7)
        self.start = gto.M(atom='H 0 0 -1.5; H 0 0 0; H 0 0 .74',
                           basis='sto-3g', spin=1, verbose=0)
        self.end = self.start.copy().set_geom_(
            'H 0 0 -.74; H 0 0 0; H 0 0 1.5')

    def tearDown(self):
        np.random.set_state(self.random_state)

    def test_path_and_endpoints(self):
        before = self.start.atom_coords().copy()
        images = interpolate(self.start, self.end, nimages=5)
        self.assertEqual(len(images), 5)
        # Compare pair distances because the backend removes rigid motion.
        for image, endpoint in ((images[0], self.start), (images[-1], self.end)):
            x = image.atom_coords()
            y = endpoint.atom_coords()
            np.testing.assert_allclose(
                np.linalg.norm(x[:, None] - x, axis=2),
                np.linalg.norm(y[:, None] - y, axis=2), atol=1e-10)
        np.testing.assert_array_equal(self.start.atom_coords(), before)
        # Returned molecules retain electronic settings and own their arrays.
        self.assertEqual(images[2].spin, 1)
        self.assertEqual(images[2].basis, 'sto-3g')
        self.assertFalse(images[2].symmetry)
        images[0].set_geom_('H 0 0 -2; H 0 0 0; H 0 0 2')
        np.testing.assert_array_equal(self.start.atom_coords(), before)

    def test_angstrom_and_bohr(self):
        angstrom = interpolate(self.start, self.end, nimages=3)
        start_bohr = self.start.copy().set_geom_(self.start.atom_coords(), unit='Bohr')
        end_bohr = self.end.copy().set_geom_(self.end.atom_coords(), unit='Bohr')
        np.random.seed(7)
        bohr = interpolate(start_bohr, end_bohr, nimages=3)
        np.testing.assert_allclose([m.atom_coords() for m in angstrom],
                                   [m.atom_coords() for m in bohr], atol=1e-8)

    def test_pyscf_energy_and_gradient(self):
        images = interpolate(self.start, self.end, nimages=3)
        with lib.with_omp_threads(1):
            scanner = self.start.UHF().nuc_grad_method().as_scanner()
            for image in images:
                energy, gradient = scanner(image)
                self.assertTrue(scanner.converged)
                self.assertTrue(np.isfinite(energy))
                self.assertEqual(gradient.shape, (3, 3))
                self.assertTrue(np.isfinite(gradient).all())
                # An isolated molecule has zero net translational force.
                np.testing.assert_allclose(gradient.sum(axis=0), 0, atol=1e-7)

    def test_incompatible_endpoints(self):
        incompatible = gto.M(atom='He 0 0 0; H 0 0 1; H 0 0 2', verbose=0)
        with self.assertRaisesRegex(ValueError, 'labels'):
            interpolate(self.start, incompatible)
        incompatible = self.end.copy()
        incompatible.spin = 3
        with self.assertRaisesRegex(ValueError, 'spin'):
            interpolate(self.start, incompatible)

    def test_invalid_controls(self):
        with self.assertRaises(ValueError):
            interpolate(self.start, self.end, nimages=2)
        with self.assertRaises(ValueError):
            interpolate(self.start, self.end, tol=0)
        with self.assertRaises(ValueError):
            interpolate(self.start, self.end, maxiter=0)


if __name__ == '__main__':
    unittest.main()
