# Local validation — 20 November 2026

## Geometry-only butane explorer

The new geometry analysis and butane example passed **nine focused tests**:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  pyscf/geomopt/test/test_path_analysis.py \
  dilawar7/geodesic/tests/test_butane_geometry.py
```

The result was `9 passed, 14 warnings in 1.75s`. Warnings were Matplotlib's
Pyparsing deprecations. Tests cover analytic distances/angles/signed torsions,
undefined geometry, torsion unwrapping with gaps, rigid-motion invariance,
reflection exclusion in RMSD, graph enumeration, invalid selections, endpoint
construction and exported-path consistency. A fresh seven-image run is tested
with SCF and electronic-gradient entry points replaced by exceptions, ensuring
this example does not invoke them. That test also rebuilds the standalone HTML.

The committed full run generated 101 images of 14-atom butane in **1.48 seconds**
for interpolation alone on this local environment. Exporting, plotting and
report assembly are outside that timing. It is a single local observation,
not a performance guarantee. All 101 requested images enter the backend's
optimization; there is no extra animation resampling. There are 65 measured
coordinates. Source and backend settings are recorded in
[`results/butane/run.json`](results/butane/run.json).

The run used the sparse backend at commit
`476aa041f562d37d5fb4121b86817a52d75e4c5d`; no tracked modifications were present
in that backend checkout. The prior published-package validation below used
the older dense solver and does not establish equivalent many-image performance.

Browser checks confirmed WebGL molecular rendering, endpoint navigation,
playback, graph-to-frame selection, synchronized markers/readouts, switching
to the Cartesian baseline, coordinate-menu highlighting, custom hydrogen
distances and reverse-order torsions, and picking atoms directly in the model.
The default CSV export produced 102 rows (header plus 101 images) and 185
columns. A copyable-text fallback was inspected for browsers that restrict
downloads. The browser download-event automation stalled, so writing the
export to disk through the browser was not certified.

The report was tested through a temporary loopback HTTP preview. It embeds its
data, scripts and third-party notices, with no external script sources. The
in-app browser blocks `file:` navigation, so direct offline-file opening in
another browser was not tested. The static figure and browser rendering were
visually inspected; the browser console reported no errors in the final check.

No energy-minimized path, reaction barrier, optimizer convergence certificate,
clean native build or cross-platform certification is claimed.

## Earlier energy/gradient workflow examples

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
