# Local validation — 20 November 2026

Nine tests passed: five adapter tests, two existing geometry-helper tests and
two worked-example tests. Both example scripts also ran successfully to generate
the committed results and diagrams. Every saved image had a converged SCF result.

The example tests independently read exported XYZ files, check hydrogen endpoint
separations and water endpoint bond lengths/angles, and recompute selected
electronic energies. Hydrogen endpoint energies agree within 1e-8 Hartree.
Both diagrams were visually inspected for readable labels, units and layout.

| Example | Electronic method | Images | Initial energy / Hartree | Final energy / Hartree |
| --- | --- | --- | --- | --- |
| Hydrogen exchange | UHF/STO-3G | 5 | -1.5735403558 | -1.5735403558 |
| Water bending | RHF/STO-3G | 5 | -74.9393290347 | -74.9480403818 |

Settings and package versions are recorded in each result folder's `run.json`.
The runs used seed 7, Python 3.9.6, PySCF 2.15.0, geodesic-interpolate 1.0.0,
NumPy 2.0.2, SciPy 1.13.1 and Matplotlib 3.9.4 on macOS. The Python source was
loaded from this fork, with native libraries supplied by the PySCF 2.15.0 wheel.

The test command from the repository root, with this fork installed, is:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  pyscf/geomopt/test/test_geodesic.py pyscf/geomopt/test/test_addons.py \
  dilawar7/geodesic/tests
```

The final observed result was `9 passed, 1 warning in 5.33s`. The warning reported
that the wheel had no OpenMP support, so its thread-count context had no effect.
These checks do not cover a clean native build, other platforms, large systems,
minimum-energy paths, transition states, reaction barriers or certification
that the geodesic optimizer reached convergence.
