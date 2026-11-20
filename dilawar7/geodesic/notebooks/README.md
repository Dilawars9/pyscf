# Colab / Jupyter notebook

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Dilawars9/pyscf/blob/dilawar7/geodesic-interpolation/dilawar7/geodesic/notebooks/butane_path_explorer.ipynb)

Open [`butane_path_explorer.ipynb`](butane_path_explorer.ipynb), select a CPU
runtime, and run the cells from top to bottom. Installation needs internet access.
The default example generates 101 butane images and makes no electronic energy
or gradient calls. There is no GPU requirement or Google Drive mount.

The notebook includes:

- Package installation and a pinned source checkout.
- Controls for image count, seed and interpolation iterations.
- The interactive molecular viewer with linked geometry plots.
- Direct calls to the distance, angle, dihedral and RMSD functions.
- Endpoint consistency checks and static figures.
- A ZIP export, with an optional Colab download cell.

## Installation design

The notebook installs PySCF 2.15.0 from its prebuilt package and adds the pinned
fork's `pyscf/geomopt/` directory to that package's search path within the current
kernel. It imports our two Python feature modules from there, while the native
libraries and other PySCF modules come from the installed wheel. The full fork
is not compiled. The output prints the file locations actually loaded.

Feature code is pinned to `af7df7057c543e3550bf1506f065126346843c8f` and the
sparse geodesic backend to `476aa041f562d37d5fb4121b86817a52d75e4c5d`.
The notebook can evolve independently while these source versions stay fixed
for reproducibility. It is distributed with cleared outputs and no credentials.

The HTML viewer runs inside an isolated iframe and embeds its own library and
data. Rendering depends on the notebook frontend and WebGL support. Static
figures and the downloadable HTML report remain available if the iframe is
restricted. Browser downloads inside an iframe may be restricted; the final
notebook cell provides a separate ZIP download through Colab's file API.

## Code locations

The actual reusable measurement functions are in
[`pyscf/geomopt/path_analysis.py`](../../../pyscf/geomopt/path_analysis.py).
The notebook calls these functions directly; it does not redefine their formulas.
The interpolation adapter is
[`pyscf/geomopt/geodesic.py`](../../../pyscf/geomopt/geodesic.py).

For the distinction between local execution checks and a hosted Colab run,
see the [validation record](../VALIDATION.md).
