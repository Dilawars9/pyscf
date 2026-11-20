#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
'''Bend water from 80 to 120 degrees with 0.95 Angstrom endpoint O-H bonds.

This is a molecular deformation example, not a chemical reaction or a geometry
optimization. Interior bond lengths are determined by geodesic interpolation.
'''

import numpy as np
from pyscf import gto
from _common import example_parser, run_example


def water(angle):
    half_angle = np.deg2rad(angle / 2)
    x, z = .95 * np.sin(half_angle), .95 * np.cos(half_angle)
    return gto.M(atom=[('O', (0, 0, 0)), ('H', (-x, 0, z)),
                       ('H', (x, 0, z))], basis='sto-3g', verbose=0)


def main():
    args = example_parser('water_bending').parse_args()
    run_example(water(80), water(120), args, 'Water bending: 80 to 120 degrees')


if __name__ == '__main__':
    main()
