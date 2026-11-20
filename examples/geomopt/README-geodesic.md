# Geodesic interpolation prototype

This optional adapter connects [geodesic-interpolate](https://github.com/virtualzx-nad/geodesic-interpolate)
to PySCF molecular objects. It generates a geometric path between atom-mapped
endpoints. PySCF can then evaluate an energy and gradient at each image.

Install this branch with the optional dependency using the normal PySCF source
installation procedure:

```sh
python -m pip install -e '.[geodesic]'
python examples/geomopt/21-geodesic.py
```

A source installation requires PySCF's native build dependencies. Installing the
standard PySCF release alone does not provide this prototype module.

```python
from pyscf import gto
from pyscf.geomopt.geodesic import interpolate

reactant = gto.M(atom='H 0 0 -1.5; H 0 0 0; H 0 0 .74',
                 basis='sto-3g', spin=1, verbose=0)
product = reactant.copy().set_geom_(
    'H 0 0 -.74; H 0 0 0; H 0 0 1.5')
images = interpolate(reactant, product, nimages=5)
scanner = reactant.UHF().nuc_grad_method().as_scanner()
energy, gradient = scanner(images[2])  # Hartree; Hartree/Bohr
```

Both endpoints must have matching atom labels, atom order, charge and spin.
The caller must establish atom correspondence, including correspondence between
identical elements. The adapter converts coordinates to Angstrom for the
backend. Returned molecules inherit the first endpoint's basis and molecular
settings, use Angstrom, and disable point-group symmetry. Inputs are unchanged.
Alignment removes rigid translations and rotations; endpoint pair distances
are preserved, while their original Cartesian frames can change.

The current interface uses joint smoothing and is intended initially for small
molecules. It has no energy-based path optimization, NEB, transition-state
search, periodic-cell support or ghost-atom support. The backend uses random
displacements and may produce different paths on repeated calls. It returns
its final iterate even when the iteration limit is reached; the adapter does
not certify geodesic convergence. Inspect paths before using them for research.

The hydrogen-exchange example uses a minimal basis to demonstrate the workflow;
its energy maximum is not a validated activation barrier.

## Local validation

Focused tests cover endpoint pair distances, independent molecule copies,
Angstrom/Bohr consistency, input validation, and PySCF energy/gradient evaluation
including zero net translational force. These are small-system integration
checks, not a full PySCF regression suite or validation of reaction barriers.

Local checks passed with Python 3.9.6, NumPy 2.0.2 and SciPy 1.13.1,
using this PySCF checkout's Python code with native libraries from the PySCF
2.15.0 macOS wheel. The five integration tests and two existing geometry-helper
tests passed against both the published geodesic-interpolate 1.0.0 package and
the selected repository at commit `476aa041f562d37d5fb4121b86817a52d75e4c5d`.
A clean native build of this fork has not been tested.

```sh
python -m pytest -q pyscf/geomopt/test/test_geodesic.py
```

This prototype was developed with OpenAI Codex assistance for implementation,
documentation and testing. Any proposed upstream contribution needs human
review and must follow PySCF's contribution and AI policies.
