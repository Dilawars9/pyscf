# SPDX-License-Identifier: Apache-2.0
'''Idealized, atom-mapped butane endpoints, generated without any optimization.

C-C = 1.54 Angstrom, C-H = 1.09 Angstrom, tetrahedral angle arccos(-1/3).
Atom order is four carbons, then H atoms grouped by their parent carbon.
These are explicit model geometries, not experimental or optimized structures.
'''

import numpy as np


def endpoints():
    length = 1.54
    carbons = np.array([[-length / 3, length * np.sqrt(8 / 9), 0],
                        [0, 0, 0], [length, 0, 0],
                        [length * 4 / 3, -length * np.sqrt(8 / 9), 0]])
    xyz = list(carbons.copy())
    bonds = [(0, 1), (1, 2), (2, 3)]
    hydrogens = {}
    for carbon in range(4):
        adjacent = [j for i, j in bonds if i == carbon] + [i for i, j in bonds if j == carbon]
        arms = carbons[adjacent] - carbons[carbon]
        arms /= np.linalg.norm(arms, axis=1, keepdims=True)
        if len(adjacent) == 2:
            normal = np.cross(*arms)
            normal /= np.linalg.norm(normal)
            directions = [-.5 * arms.sum(axis=0) + sign * np.sqrt(2 / 3) * normal for sign in (-1, 1)]
        else:
            axis = arms[0]
            next_carbon = 2 if carbon == 0 else 1
            reference = carbons[next_carbon] - carbons[adjacent[0]]
            u = reference - np.dot(reference, axis) * axis
            u /= np.linalg.norm(u)
            v = np.cross(axis, u)
            directions = [-axis / 3 + np.sqrt(8 / 9) * (np.cos(phi) * u + np.sin(phi) * v)
                          for phi in np.deg2rad([60, 180, 300])]
        hydrogens[carbon] = []
        for direction in directions:
            atom = len(xyz)
            hydrogens[carbon].append(atom)
            xyz.append(carbons[carbon] + 1.09 * direction)
            bonds.append((carbon, atom))
    start = np.array(xyz)
    end = start.copy()
    # Rotate the entire right-hand fragment about the C2-C3 axis by -120 deg.
    group = [3] + hydrogens[2] + hydrogens[3]
    phi = np.deg2rad(-120)
    rotation = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)],
                         [0, np.sin(phi), np.cos(phi)]])
    end[group] = (start[group] - carbons[2]) @ rotation.T + carbons[2]
    return ['C'] * 4 + ['H'] * 10, start, end, bonds
