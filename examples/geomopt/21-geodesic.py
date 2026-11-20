#!/usr/bin/env python
'''Interpolate H + H2 -> H2 + H, then evaluate PySCF energies and gradients.

Requires geodesic-interpolate. This is a small workflow demonstration with a
minimal basis, not a quantitative reaction-barrier calculation. The geometric
interpolator does not use these energies or optimize a minimum-energy path.
'''

from pyscf import gto
from pyscf.geomopt.geodesic import interpolate

reactant = gto.M(atom='H 0 0 -1.5; H 0 0 0; H 0 0 .74',
                 basis='sto-3g', spin=1, verbose=0)
product = reactant.copy().set_geom_(
    'H 0 0 -.74; H 0 0 0; H 0 0 1.5', unit='Angstrom')
images = interpolate(reactant, product, nimages=5)
scanner = reactant.UHF().nuc_grad_method().as_scanner()

print('image   energy / Hartree    max |gradient| / (Hartree/Bohr)')
for i, mol in enumerate(images):
    energy, gradient = scanner(mol)
    if not scanner.converged:
        raise RuntimeError('SCF did not converge for image %d' % i)
    print('%3d    %18.10f    %12.6f' % (i, energy, abs(gradient).max()))
