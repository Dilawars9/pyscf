#!/usr/bin/env python
'''Interpolate H + H2 -> H2 + H, then evaluate PySCF energies and gradients.

Requires geodesic-interpolate. This is a small workflow demonstration with a
minimal basis, not a quantitative reaction-barrier calculation. The geometric
interpolator does not use these energies or optimize a minimum-energy path.
'''

from pyscf import gto
from _common import example_parser, run_example


def main():
    args = example_parser('hydrogen_exchange').parse_args()
    reactant = gto.M(atom='H 0 0 -1.5; H 0 0 0; H 0 0 .74',
                     basis='sto-3g', spin=1, verbose=0)
    product = reactant.copy().set_geom_(
        'H 0 0 -.74; H 0 0 0; H 0 0 1.5', unit='Angstrom')
    run_example(reactant, product, args, 'Hydrogen exchange: H + H2 to H2 + H')


if __name__ == '__main__':
    main()
