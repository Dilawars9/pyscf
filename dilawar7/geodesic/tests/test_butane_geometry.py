# SPDX-License-Identifier: Apache-2.0
'''Validate the idealized endpoint construction and saved 101-image example.'''

import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from pyscf.geomopt import path_analysis as geometry

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('butane_geometry', ROOT / 'examples' / 'butane_geometry.py')
butane = importlib.util.module_from_spec(spec)
spec.loader.exec_module(butane)


def test_endpoint_geometry():
    symbols, start, end, bonds = butane.endpoints()
    path = np.array([start, end])
    assert symbols == ['C'] * 4 + ['H'] * 10
    selections = geometry.coordinates_from_bonds(14, bonds)
    assert [len(selections[key]) for key in ('distance', 'angle', 'dihedral')] == [13, 24, 27]
    expected = [1.54 if max(pair) < 4 else 1.09 for pair in selections['distance']]
    np.testing.assert_allclose(geometry.distances(path, selections['distance']), [expected, expected], atol=1e-12)
    np.testing.assert_allclose(geometry.angles(path, selections['angle']),
                               np.rad2deg(np.arccos(-1 / 3)), atol=1e-12)
    np.testing.assert_allclose(geometry.dihedrals(path, [(0, 1, 2, 3)]), [[180], [60]], atol=1e-12)


def test_saved_path_measurements():
    data = json.loads((ROOT / 'results' / 'butane' / 'analysis.json').read_text())
    coords = np.array(data['coordinates']['geodesic'])
    assert coords.shape == (101, 14, 3)
    assert data['metadata']['electronic_energy_calls'] == 0
    assert data['metadata']['electronic_gradient_calls'] == 0
    assert len(data['measurements']) == 65
    # Independently reconstruct distances directly from the exported frames.
    selected = next(entry for entry in data['measurements'] if entry['id'] == 'distance_1_4')
    independent = np.sqrt(((coords[:, 0] - coords[:, 3]) ** 2).sum(axis=1))
    np.testing.assert_allclose(selected['geodesic'], independent, atol=1e-9)
    # All endpoint pair distances survive the backend's rigid alignment.
    _, start, end, _ = butane.endpoints()
    for actual, expected in ((coords[0], start), (coords[-1], end)):
        np.testing.assert_allclose(np.linalg.norm(actual[:, None] - actual, axis=-1),
                                   np.linalg.norm(expected[:, None] - expected, axis=-1), atol=2e-10)
    # Known endpoints provide an independent angle/torsion reference.
    torsion = next(entry for entry in data['measurements'] if entry['id'] == 'dihedral_1_2_3_4')
    np.testing.assert_allclose([torsion['geodesic'][0], torsion['geodesic'][-1]], [180, 60], atol=1e-8)


def test_example_has_no_electronic_calls(tmp_path, monkeypatch):
    from pyscf.scf import hf, uhf
    from pyscf.grad import rhf as rhf_grad, uhf as uhf_grad

    def forbidden(*args, **kwargs):
        raise AssertionError('Geometry-only examples must not run electronic calculations')

    monkeypatch.setattr(hf, 'kernel', forbidden)
    monkeypatch.setattr(uhf.UHF, 'kernel', forbidden)
    monkeypatch.setattr(rhf_grad, 'grad_elec', forbidden)
    monkeypatch.setattr(uhf_grad, 'grad_elec', forbidden)
    monkeypatch.syspath_prepend(str(ROOT / 'examples'))
    spec = importlib.util.spec_from_file_location('butane_path', ROOT / 'examples' / '03_butane_path.py')
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    monkeypatch.setattr(sys, 'argv', ['03_butane_path.py', '--nimages', '7', '--maxiter', '10',
                                    '--data-only', '--output', str(tmp_path)])
    random_state = np.random.get_state()
    try:
        example.main()
    finally:
        np.random.set_state(random_state)
    data = json.loads((tmp_path / 'analysis.json').read_text())
    assert len(data['coordinates']['geodesic']) == 7
    # Report generation is independent of the expensive interpolation step.
    from path_report import build_report
    build_report(data, tmp_path / 'report.html')
    report = (tmp_path / 'report.html').read_text()
    assert 'All 7 requested images' in report
    assert '__DATA_JSON__' not in report
    assert '<script src=' not in report
    assert '3Dmol' in report
