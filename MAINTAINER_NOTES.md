# Maintainer notes: geodesic molecular paths

Scope: Dilawar's work on branch `dilawar7/geodesic-interpolation` of this PySCF
fork. These notes describe how to maintain the adapter, geometry analysis,
worked examples, interactive report and notebook. Upstream PySCF policies are
in [CONTRIBUTING.md](CONTRIBUTING.md) and [AI_POLICY.md](AI_POLICY.md).

The [user guide](dilawar7/geodesic/README.md) explains how to use the feature.
The [validation record](dilawar7/geodesic/VALIDATION.md) records checks actually
performed, their environments and their limits. Update that record when new
validation is completed; keep historical results identifiable.

## Code layout and responsibilities

| Component | Location | Responsibility |
| --- | --- | --- |
| PySCF adapter | [`pyscf/geomopt/geodesic.py`](pyscf/geomopt/geodesic.py) | Validate endpoints, convert units, call the optional backend and return molecule copies. |
| Geometry analysis | [`pyscf/geomopt/path_analysis.py`](pyscf/geomopt/path_analysis.py) | Distances, angles, signed dihedrals, unwrapping, aligned RMSD and graph-based coordinate enumeration. |
| Butane construction | [`butane_geometry.py`](dilawar7/geodesic/examples/butane_geometry.py) | Generate idealized endpoints with an explicit atom mapping and bond graph. |
| Path generation and exports | [`03_butane_path.py`](dilawar7/geodesic/examples/03_butane_path.py) | Generate frames, measure both paths, write data/metadata and make static figures. |
| HTML assembly | [`path_report.py`](dilawar7/geodesic/examples/path_report.py) | Embed the template, scripts, data, license notices and viewer checksum. |
| Browser UI | [`report/`](dilawar7/geodesic/report/) | Molecular rendering, frame playback, linked graphs, atom picking and custom-coordinate exports. |
| Colab/Jupyter notebook | [`butane_path_explorer.ipynb`](dilawar7/geodesic/notebooks/butane_path_explorer.ipynb) | Install dependencies, load a pinned feature revision and demonstrate the complete workflow. |
| Generated evidence | [`results/`](dilawar7/geodesic/results/) | Saved paths, tables, figures, reports and run settings. |

The original interpolation algorithm is supplied by the external
[geodesic-interpolate project](https://github.com/virtualzx-nad/geodesic-interpolate).
Our adapter and analysis code remain in `pyscf/geomopt/` so callers can import
them through PySCF. Fork-specific tutorials and results live in `dilawar7/`.

## Contracts to preserve

- Coordinates for analysis have shape `(nframes, natoms, 3)`. Distances and RMSD
  retain the coordinate unit. The adapter passes Angstrom to the backend.
- Python atom indices are **0-based**. Viewer labels, notebook selection forms
  and CSV atom numbers are **1-based**. Keep conversions explicit.
- Angles and signed dihedrals use degrees. The signed-torsion convention is
  documented in `path_analysis.dihedrals()`; exact trans is reported as +180°.
- Undefined angles/dihedrals are NaN in Python, null in JSON and empty fields
  in CSV. Plot lines must break across undefined regions. Unwrap each finite
  segment independently; do not bridge a gap or invent a winding number.
- Aligned RMSD uses equal atom weights and proper rotations, excluding reflections.
- The adapter preserves endpoint internal geometries, while rigid alignment may
  change Cartesian frames. It must leave input molecules unchanged.
- Atom correspondence and reference connectivity are supplied by the caller.
  Matching element lists cannot detect permutations of identical atoms. Display
  bonds do not establish chemical bond orders or a reaction mechanism.
- Keep the butane workflow free of electronic energy and gradient calculations.
  The earlier hydrogen-exchange and water-bending examples intentionally include
  those calculations and are documented separately.

The browser's custom-coordinate formulas are implemented in
[`report.js`](dilawar7/geodesic/report/report.js), allowing the report to work
without a Python server. If a convention or degeneracy rule changes, update
both implementations and compare custom browser readouts against Python,
including torsion signs, wrapping boundaries and undefined selections.

## Dependency and source pins

| Item | Current notebook/report selection | Maintenance implication |
| --- | --- | --- |
| PySCF native package | `pyscf==2.15.0` | The notebook uses its wheel and loads two feature modules from the fork; it does not compile the whole fork. |
| Notebook feature snapshot | `af7df7057c543e3550bf1506f065126346843c8f` | Controls feature modules, example scripts and report assets used by the notebook. |
| Geodesic backend snapshot | `476aa041f562d37d5fb4121b86817a52d75e4c5d` | Includes the sparse joint solver used for the many-image example. |
| Browser renderer | 3Dmol.js `2.5.5` | Vendored in [`report/vendor/`](dilawar7/geodesic/report/vendor/) with project and bundle license notices. |
| Scientific/plotting dependencies | NumPy/SciPy supplied through the package environment; notebook requests Matplotlib `>=3.8` | Record actual resolved versions in each new validation run. These are not a complete environment lock. |

The published geodesic package also identifies itself as version `1.0.0` but
uses the older dense solver. A version number alone does not identify the tested
backend implementation. Preserve the commit pin and notebook's explicit
`--force-reinstall --no-deps` backend installation. The example records backend
provenance and sparse-solver availability in `run.json`.

Notebook setup extends `pyscf.geomopt.__path__` within its kernel to import
`geodesic.py` and `path_analysis.py` from the pinned checkout. Other PySCF modules
and native libraries come from the installed wheel. Restart the kernel after
changing installed code or source pins; package installation does not unload
modules already present in memory.

### Updating the notebook snapshot

1. Make and validate implementation changes on the feature branch. Commit and
   publish the implementation so GitHub can fetch its full commit SHA.
2. Set `FORK_COMMIT` in the notebook's `install` cell to that SHA. A renderer or
   example change also needs a new feature pin: the notebook reads those files
   from the same pinned checkout. Documentation-only edits need no pin change.
3. Change `BACKEND_COMMIT` only when deliberately updating the external backend.
   Check compatibility and sparse-solver behavior rather than relying on its
   package version. Keep the interpolation dependency optional in PySCF.
4. Update current pin descriptions in the [notebook notes](dilawar7/geodesic/notebooks/README.md)
   and any relevant setup instructions. Preserve old pins in historical run
   records; do not relabel an old result as a new validation.
5. Execute the notebook in a fresh kernel, inspect its printed module locations,
   run settings, generated figures, embedded viewer and ZIP contents. Record
   local execution and hosted Colab execution as separate checks.
6. Clear code-cell outputs and execution counts before committing the notebook.
   Preserve cell IDs and Colab form metadata, then check the badge/link still
   opens the published notebook.

When updating 3Dmol.js, retain the selected release's license notices and bundle
provenance. The report builder computes the bundle checksum during assembly.
Rebuild the HTML and inspect WebGL rendering, labels, picking and camera behavior
after changing the renderer; static image tests cannot verify those controls.

## Regenerate examples and reports

Run commands from the repository root in an environment where this fork's
feature modules and PySCF's native libraries are available. The notebook is the
documented wheel-plus-feature setup. A clean native build of the full fork
remains unverified. See the [butane guide](dilawar7/geodesic/BUTANE.md) for the
sparse-backend installation command.

Generate fresh outputs in a separate directory for review:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python dilawar7/geodesic/examples/03_butane_path.py \
  --nimages 101 --seed 7 --maxiter 50 --output /tmp/dilawar-butane-review
```

The script writes coordinates, a Cartesian baseline, 65 predefined coordinate
traces for butane, run metadata, figures and the standalone HTML. It overwrites
files in the selected output directory. Review changes before replacing the
committed [`results/butane/`](dilawar7/geodesic/results/butane/) artifacts. Retain
the generated `run.json` and record the feature commit in `VALIDATION.md` for a
new run. A single local timing is not a performance guarantee.

For report CSS/JavaScript/template changes, rebuild HTML from saved analysis
without rerunning interpolation:

```sh
python dilawar7/geodesic/examples/path_report.py dilawar7/geodesic/results/butane/analysis.json
```

This updates `report.html` beside `analysis.json`; it does not regenerate the
static figures or numerical measurements. Changes to coordinate calculations
require fresh analysis and corresponding result checks.

## Validation workflow

For geometry-analysis or geometry-only example changes, use the focused suite:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q \
  pyscf/geomopt/test/test_path_analysis.py \
  dilawar7/geodesic/tests/test_butane_geometry.py
```

These tests include known geometries, proper rigid-motion invariance, torsion
wrapping/degeneracy, graph enumeration, endpoint mapping and saved-path checks.
The example test replaces electronic calculation entry points with exceptions
while generating a small fresh path. Saved-artifact checks complement a fresh
run; neither alone proves general scientific validity.

For adapter changes, also run `pyscf/geomopt/test/test_geodesic.py`; that file
includes an electronic energy/gradient integration test. For changes to the
earlier energy examples, select their relevant tests in
[`test_worked_examples.py`](dilawar7/geodesic/tests/test_worked_examples.py).
Do not treat historical test counts as the current collection's result.

Notebook validation requires `nbformat`, `nbclient` and `ipykernel` in a chosen
test environment. Register that environment's kernel, then execute into a
scratch directory, retaining the output-bearing notebook as validation evidence:

```sh
python -m pip install nbformat nbclient ipykernel
python -m ipykernel install --user --name dilawar-geodesic --display-name 'Dilawar geodesic'
python - <<'PY'
from pathlib import Path
import tempfile
import nbformat
from nbclient import NotebookClient

source = Path('dilawar7/geodesic/notebooks/butane_path_explorer.ipynb').resolve()
notebook = nbformat.read(source, as_version=4)
nbformat.validate(notebook)
work = Path(tempfile.mkdtemp(prefix='dilawar-notebook-check-'))
NotebookClient(
    notebook, kernel_name='dilawar-geodesic', timeout=300,
    resources={'metadata': {'path': str(work)}},
).execute()
nbformat.write(notebook, work / 'executed.ipynb')
print('Executed notebook and generated results:', work)
PY
```

Notebook setup installs packages into the selected kernel environment and
fetches public pinned sources, so use an isolated environment with internet
access. Keep execution outputs separate from the cleared notebook source.
Check that the viewer cell emits HTML, the plot cells emit PNG images, endpoint
assertions pass, and the ZIP contains the expected generated files.

For viewer changes, inspect rendering, endpoint navigation, playback, linked
graph markers, baseline switching, custom measurements and atom picking in a
browser. Check CSV/XYZ contents and the export-text fallback. Test hosted
Colab's iframe and file-download behavior in Colab itself when claiming support
there. Documentation-only changes need link and formatting checks rather than
fresh numerical runs.

## Current limits and outstanding work

The detailed evidence is in [VALIDATION.md](dilawar7/geodesic/VALIDATION.md).
The original geometry-only suite passed nine tests, and all nine notebook code
cells completed in local Jupyter. The public notebook was also opened in Colab;
loading its source is not execution on a hosted runtime.

- Run a fresh hosted Colab CPU session through installation, generation, iframe
  interaction and the optional ZIP download. Record its actual environment.
- Validate a clean native build of the fork and any additional supported
  platforms. Existing local checks used PySCF wheel-provided native libraries.
- Verify standalone offline-file opening and browser downloads across target
  browsers. Local HTTP rendering and copyable export text have been checked;
  writing exports to disk through the browser was not certified.
- Add explicit optimizer termination/convergence diagnostics before making a
  stable geodesic-convergence claim. More images alone do not establish it.
- Benchmark larger paths/systems before making scaling claims. Additional ring
  or proton-transfer examples need documented, consistently mapped endpoints
  and their own validation.

No minimum-energy path, transition state, reaction barrier or physical time
interpretation has been established. The fixed graph is a visualization and
measurement input, not a bond-order prediction. Endpoint model assumptions and
unit/index conventions must stay visible in examples and reports.

## Upstream contribution and attribution

Preserve attribution to PySCF, the geodesic backend authors and 3Dmol.js, and
retain third-party license notices. Describe our contribution as the adapter,
analysis and workflow around the existing interpolation algorithm.

Follow PySCF's [contribution guidance](CONTRIBUTING.md) and
[AI policy](AI_POLICY.md). Substantial AI assistance must be disclosed. A human
must understand and review generated work and submit an upstream pull request;
the policy does not permit an agent to write code and autonomously submit it.
Changes on this personal feature branch do not establish upstream acceptance.
