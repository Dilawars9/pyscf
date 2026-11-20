# SPDX-License-Identifier: Apache-2.0
'''Geometry measurements for atom-mapped molecular paths.

Coordinates have shape (nframes, natoms, 3). Atom indices are zero based.
Distances retain the coordinate unit; angles and signed dihedrals use degrees.
Undefined angles/dihedrals are NaN, never an invented zero.
'''

from itertools import combinations
import numpy as np

__all__ = ['distances', 'angles', 'dihedrals', 'unwrap_dihedrals',
           'aligned_rmsd', 'coordinates_from_bonds']


def _path(coords):
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 3 or coords.shape[-1] != 3 or not np.isfinite(coords).all():
        raise ValueError('Coordinates must be finite with shape (nframes, natoms, 3)')
    if not coords.shape[0] or not coords.shape[1]:
        raise ValueError('A path must contain frames and atoms')
    return coords


def _indices(indices, width, natoms):
    values = np.asarray(indices)
    if values.size == 0:
        return np.empty((0, width), dtype=int)
    if values.ndim != 2 or values.shape[1] != width or values.dtype.kind not in 'iu':
        raise ValueError('Atom selections must be integer tuples of length %d' % width)
    if (values < 0).any() or (values >= natoms).any():
        raise ValueError('Atom index outside the path')
    if any(len(set(row)) != width for row in values):
        raise ValueError('A measurement must select distinct atoms')
    return values


def distances(coords, pairs):
    '''Return distances for arbitrary atom pairs, shape (nframes, npairs).'''
    coords = _path(coords)
    i, j = _indices(pairs, 2, coords.shape[1]).T
    return np.linalg.norm(coords[:, i] - coords[:, j], axis=-1)


def angles(coords, triples):
    '''Return i-j-k angles in degrees, shape (nframes, ntriples).

    The central atom is j. A zero-length arm gives NaN.
    '''
    coords = _path(coords)
    i, j, k = _indices(triples, 3, coords.shape[1]).T
    a, b = coords[:, i] - coords[:, j], coords[:, k] - coords[:, j]
    denominator = np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)
    cosine = np.full(denominator.shape, np.nan)
    np.divide(np.sum(a * b, axis=-1), denominator, out=cosine, where=denominator > 1e-14)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1)))


def dihedrals(coords, quadruples):
    '''Return signed i-j-k-l torsions in degrees, shape (nframes, nquadruples).

    Look along j to k: the sign follows atan2((axis cross projected i-j)
    dot projected l-k, projected i-j dot projected l-k). Values lie in
    (-180, 180], with numerically exact trans reported as +180.
    A zero central bond or an arm collinear with it gives NaN.
    '''
    coords = _path(coords)
    i, j, k, l = _indices(quadruples, 4, coords.shape[1]).T
    a, axis, b = coords[:, i] - coords[:, j], coords[:, k] - coords[:, j], coords[:, l] - coords[:, k]
    length = np.linalg.norm(axis, axis=-1, keepdims=True)
    unit = np.zeros_like(axis)
    np.divide(axis, length, out=unit, where=length > 1e-14)
    v = a - np.sum(a * unit, axis=-1, keepdims=True) * unit
    w = b - np.sum(b * unit, axis=-1, keepdims=True) * unit
    valid = ((length[..., 0] > 1e-14) &
             (np.linalg.norm(v, axis=-1) > 1e-10 * np.linalg.norm(a, axis=-1)) &
             (np.linalg.norm(w, axis=-1) > 1e-10 * np.linalg.norm(b, axis=-1)))
    result = np.rad2deg(np.arctan2(np.sum(np.cross(unit, v) * w, axis=-1),
                                 np.sum(v * w, axis=-1)))
    result[np.isclose(result, -180, atol=1e-10, rtol=0)] = 180
    return np.where(valid, result, np.nan)


def unwrap_dihedrals(values):
    '''Unwrap degree-valued torsion columns along frames, preserving NaN gaps.

    Each finite segment is unwrapped independently. Assumes successive changes
    smaller than 180 degrees; sparse samples cannot determine a winding number.
    '''
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError('Torsions must have shape (nframes, ncoordinates)')
    result = values.copy()
    for column in range(values.shape[1]):
        finite = np.isfinite(values[:, column])
        edges = np.diff(np.r_[False, finite, False].astype(int))
        for start, stop in zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)):
            result[start:stop, column] = np.rad2deg(np.unwrap(np.deg2rad(values[start:stop, column])))
    return result


def aligned_rmsd(coords, reference=None):
    '''Return atom RMSD after proper rigid alignment, in the coordinate unit.

    All atoms receive equal weights. Reflections are excluded. By default the
    first frame is the reference; otherwise supply an array of shape (natoms, 3).
    '''
    coords = _path(coords)
    ref = coords[0] if reference is None else np.asarray(reference, dtype=float)
    if ref.shape != coords.shape[1:] or not np.isfinite(ref).all():
        raise ValueError('Reference must have shape (natoms, 3) and finite coordinates')
    ref = ref - ref.mean(axis=0)
    values = []
    for frame in coords:
        centered = frame - frame.mean(axis=0)
        u, _, vt = np.linalg.svd(centered.T @ ref)
        if np.linalg.det(u @ vt) < 0:
            u[:, -1] *= -1
        displacement = centered @ u @ vt - ref
        values.append(np.sqrt(np.mean(np.sum(displacement ** 2, axis=-1))))
    return np.asarray(values)


def coordinates_from_bonds(natoms, bonds):
    '''Enumerate unique bond distances, angles and proper torsions from a graph.

    The supplied zero-based graph is fixed; no bond orders or connectivity
    changes are inferred from distances. Reverse-equivalent coordinates appear
    once. Arbitrary contacts can be measured separately with distances().
    '''
    bonds = _indices(bonds, 2, natoms)
    neighbors = [set() for _ in range(natoms)]
    pairs = set()
    for i, j in bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
        pairs.add(tuple(sorted((i, j))))
    triples = {(i, j, k) for j, group in enumerate(neighbors) for i, k in combinations(sorted(group), 2)}
    torsions = set()
    for j, k in sorted(pairs):
        for i in neighbors[j] - {k}:
            for l in neighbors[k] - {j, i}:
                item = (i, j, k, l)
                torsions.add(min(item, item[::-1]))
    return {'distance': sorted(pairs), 'angle': sorted(triples), 'dihedral': sorted(torsions)}
