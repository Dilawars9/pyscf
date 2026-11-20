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

'''Molecular path interpolation using the optional geodesic-interpolate package.

This constructs a geometric path, without evaluating electronic energies or
optimizing a minimum-energy path. Atom correspondence must be supplied by the
caller through consistent atom ordering at both endpoints.

Backend: https://github.com/virtualzx-nad/geodesic-interpolate
'''

import operator
import numpy as np
from pyscf import gto

__all__ = ['interpolate']


def interpolate(start, end, nimages=17, tol=2e-3, maxiter=50):
    '''Interpolate between two molecular geometries in internal coordinates.

    Args:
        start, end : pyscf.gto.Mole
            Built molecular endpoints with matching atom labels, atom ordering,
            charge and spin. At least two real atoms are required. Matching
            element lists cannot detect permutations of identical atoms.
            Periodic cells and ghost atoms are unsupported.
        nimages : int
            Number of images including endpoints, at least 3.
        tol : float
            Positive tolerance passed to geodesic smoothing. Initial path
            redistribution uses five times this tolerance.
        maxiter : int
            Positive function-evaluation limit for joint path smoothing.

    Returns:
        list of pyscf.gto.Mole
            Independent molecule copies using the basis, charge, spin and other
            molecular settings of ``start``, with point-group symmetry disabled.
            Coordinates are passed to geodesic-interpolate in Angstrom and the
            returned molecules use Angstrom. Rigid alignment centers and rotates
            images, preserving endpoint internal geometries rather than their
            original Cartesian frames. Neither input molecule is modified.

    Notes:
        Install the optional dependency with ``pip install geodesic-interpolate``.
        The backend uses random initial displacements, so paths can vary between
        calls. This prototype uses joint smoothing, intended for small systems;
        it returns the backend's final iterate, without certifying convergence
        or a minimum-energy path. Inspect the path before scientific use.
    '''
    nimages = operator.index(nimages)
    maxiter = operator.index(maxiter)
    if nimages < 3 or maxiter < 1:
        raise ValueError('nimages must be at least 3 and maxiter must be positive')
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError('tol must be finite and positive')
    for mol in (start, end):
        if not isinstance(mol, gto.Mole) or hasattr(mol, 'pbc_intor'):
            raise TypeError('Only molecular gto.Mole endpoints are supported')
        if not mol._built or mol.natm < 2 or np.any(mol.atom_charges() <= 0):
            raise ValueError('Endpoints must be built with at least two real atoms')
    labels = [start.atom_symbol(i) for i in range(start.natm)]
    if labels != [end.atom_symbol(i) for i in range(end.natm)]:
        raise ValueError('Endpoint atom labels and ordering must match')
    if start.charge != end.charge or start.spin != end.spin:
        raise ValueError('Endpoint charge and spin must match')
    coords = np.array([start.atom_coords(unit='Angstrom'),
                       end.atom_coords(unit='Angstrom')])
    if not np.isfinite(coords).all():
        raise ValueError('Endpoint coordinates must be finite')

    try:
        from geodesic_interpolate.interpolation import redistribute
        from geodesic_interpolate.geodesic import Geodesic
    except ImportError as err:
        raise ImportError('Geodesic interpolation requires the optional package; '
                          'install it with pip install geodesic-interpolate') from err

    symbols = [start.atom_pure_symbol(i) for i in range(start.natm)]
    raw = redistribute(symbols, coords, nimages, tol=tol * 5)
    smoother = Geodesic(symbols, raw, friction=1e-2)
    smoother.smooth(tol=tol, max_iter=maxiter)
    path = np.asarray(smoother.path)
    if path.shape != (nimages, start.natm, 3) or not np.isfinite(path).all():
        raise RuntimeError('Geodesic interpolation returned an invalid path')
    return [start.copy().set_geom_(image.copy(), unit='Angstrom', symmetry=False)
            for image in path]
