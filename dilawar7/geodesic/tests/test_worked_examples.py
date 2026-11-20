# SPDX-License-Identifier: Apache-2.0
'''Check the saved paths against independent endpoint geometry and SCF results.'''

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from pyscf import gto

pytestmark = pytest.mark.skipif(
    any(importlib.util.find_spec(m) is None for m in ('geodesic_interpolate', 'matplotlib')),
    reason='Worked examples require geodesic-interpolate and matplotlib')
EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


def saved_run(script, output):
    subprocess.run([sys.executable, str(EXAMPLES / script), '--output', str(output)],
                   check=True, capture_output=True, text=True)
    lines = (output / 'path.xyz').read_text().splitlines()
    # Both examples have three atoms: count, comment, then three XYZ rows.
    assert len(lines) == 25
    molecules = []
    metadata = json.loads((output / 'run.json').read_text())
    for offset in range(0, len(lines), 5):
        assert lines[offset] == '3'
        molecules.append(gto.M(atom='\n'.join(lines[offset + 2:offset + 5]),
                               unit='Angstrom', basis='sto-3g',
                               spin=metadata['spin'], verbose=0))
    with (output / 'energies.csv').open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == len(molecules) == metadata['nimages'] == 5
    assert [int(row['image']) for row in rows] == list(range(5))
    assert all(row['scf_converged'] == 'True' for row in rows)
    assert all(np.isfinite(float(row['energy_hartree'])) for row in rows)
    assert (output / 'path_and_energy.svg').stat().st_size > 1000
    assert (output / 'path_and_energy.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    return molecules, rows


def test_hydrogen_exchange(tmp_path):
    molecules, rows = saved_run('01_hydrogen_exchange.py', tmp_path)
    for mol, expected in ((molecules[0], [1.5, .74]),
                          (molecules[-1], [.74, 1.5])):
        xyz = mol.atom_coords(unit='Angstrom')
        np.testing.assert_allclose(np.linalg.norm(np.diff(xyz, axis=0), axis=1),
                                   expected, atol=1e-10)
    # Exchanging identical nuclei leaves endpoint energies equal.
    assert abs(float(rows[0]['energy_hartree']) - float(rows[-1]['energy_hartree'])) < 1e-8
    independent = molecules[0].UHF().run(conv_tol=1e-11).e_tot
    assert abs(independent - float(rows[0]['energy_hartree'])) < 1e-8


def test_water_bending(tmp_path):
    molecules, rows = saved_run('02_water_bending.py', tmp_path)
    for mol, angle in ((molecules[0], 80), (molecules[-1], 120)):
        xyz = mol.atom_coords(unit='Angstrom')
        bonds = xyz[1:] - xyz[0]
        lengths = np.linalg.norm(bonds, axis=1)
        np.testing.assert_allclose(lengths, .95, atol=1e-10)
        actual_angle = np.rad2deg(np.arccos(np.dot(*bonds) / np.prod(lengths)))
        assert abs(actual_angle - angle) < 1e-8
    # Independently recompute an interior energy from the exported XYZ data.
    independent = molecules[2].RHF().run(conv_tol=1e-11).e_tot
    assert abs(independent - float(rows[2]['energy_hartree'])) < 1e-8
