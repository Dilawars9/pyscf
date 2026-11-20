#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
'''Generate and analyze a geometry-only trans-to-gauche butane path.'''

import argparse
import csv
from importlib.metadata import distribution, version
import json
from pathlib import Path
import platform
import subprocess
import time
from urllib.parse import unquote, urlparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyscf import gto
from pyscf.geomopt.geodesic import interpolate
from pyscf.geomopt import path_analysis as geometry
from butane_geometry import endpoints


def backend_info():
    import geodesic_interpolate.geodesic as backend
    result = {'version': version('geodesic-interpolate'),
              'sparse_solver': hasattr(backend.Geodesic, 'target_deriv_sparse')}
    direct = distribution('geodesic-interpolate').read_text('direct_url.json')
    if direct:
        direct = json.loads(direct)
        result['source_commit'] = direct.get('vcs_info', {}).get('commit_id')
        url = urlparse(direct.get('url', ''))
        if url.scheme == 'file':
            checkout = Path(unquote(url.path))
            commit = subprocess.run(['git', '-C', str(checkout), 'rev-parse', 'HEAD'],
                                    capture_output=True, text=True)
            if commit.returncode == 0:
                result['source_commit'] = commit.stdout.strip()
                dirty = subprocess.run(['git', '-C', str(checkout), 'status', '--porcelain', '--untracked-files=no'],
                                       capture_output=True, text=True, check=True)
                result['source_dirty'] = bool(dirty.stdout.strip())
    return result


def xyz_text(symbols, frames):
    lines = []
    for i, frame in enumerate(frames):
        lines.extend([str(len(symbols)), 'image %d; Angstrom; geometry only' % i])
        lines.extend('%s %.12f %.12f %.12f' % (symbol, *point)
                     for symbol, point in zip(symbols, frame))
    return '\n'.join(lines) + '\n'


def measurements(symbols, bonds, path, baseline):
    selections = geometry.coordinates_from_bonds(len(symbols), bonds)
    selections['distance'].insert(0, (0, 3))  # A contact, not a covalent bond.
    result = []
    for kind, function in [('distance', geometry.distances), ('angle', geometry.angles),
                           ('dihedral', geometry.dihedrals)]:
        indices = selections[kind]
        a, b = function(path, indices), function(baseline, indices)
        wrapped_a, wrapped_b = a.copy(), b.copy()
        if kind == 'dihedral':
            a, b = geometry.unwrap_dihedrals(a), geometry.unwrap_dihedrals(b)
        for column, atoms in enumerate(indices):
            label = '–'.join('%s%d' % (symbols[j], j + 1) for j in atoms)
            entry = dict(id=kind + '_' + '_'.join(str(j + 1) for j in atoms), kind=kind,
                         label=label, atoms=list(atoms), unit='Å' if kind == 'distance' else '°',
                         geodesic=a[:, column].tolist(), cartesian=b[:, column].tolist())
            if kind == 'dihedral':
                entry.update(wrapped_geodesic=wrapped_a[:, column].tolist(),
                             wrapped_cartesian=wrapped_b[:, column].tolist())
            result.append(entry)
    return result


def finite_json(value):
    if isinstance(value, np.generic):
        return finite_json(value.item())
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def static_figure(data, output):
    path = np.array(data['coordinates']['geodesic'])
    # One shared orientation makes the backbone readable in the static views.
    _, _, vt = np.linalg.svd(path[0] - path[0].mean(axis=0))
    if np.linalg.det(vt) < 0:
        vt[-1] *= -1
    path = path @ vt.T
    bonds = data['bonds']
    colors = ['#334155'] * 4 + ['#d6dee9'] * 10
    fig = plt.figure(figsize=(14, 7.5), layout='constrained')
    grid = fig.add_gridspec(2, 3, height_ratios=[1.1, 1])
    fig.suptitle('Butane · trans to gauche', fontsize=23, weight='bold', color='#172944')
    extent = max(np.ptp(path, axis=(0, 1)).max() * .55, 1)
    center = path.mean(axis=(0, 1))
    for panel, frame in enumerate(np.linspace(0, len(path) - 1, 3, dtype=int)):
        ax = fig.add_subplot(grid[0, panel], projection='3d')
        xyz = path[frame]
        for i, j in bonds:
            ax.plot(*xyz[[i, j]].T, color='#64748b' if max(i, j) < 4 else '#b0bdce',
                    lw=5 if max(i, j) < 4 else 2.5, solid_capstyle='round')
        ax.scatter(*xyz.T, c=colors, s=[180] * 4 + [65] * 10, edgecolors='#ffffff',
                   linewidths=.6, depthshade=True)
        ax.set(xlim=(center[0] - extent, center[0] + extent),
               ylim=(center[1] - extent, center[1] + extent),
               zlim=(center[2] - extent, center[2] + extent))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=35, azim=-80)
        ax.set_axis_off()
        ax.set_title('Image %d' % frame, fontsize=12)
    defaults = ['distance_1_4', 'angle_1_2_3', 'dihedral_1_2_3_4']
    for panel, key in enumerate(defaults):
        entry = next(item for item in data['measurements'] if item['id'] == key)
        ax = fig.add_subplot(grid[1, panel])
        ax.plot(entry['geodesic'], color='#1967d2', lw=2.5, label='Geodesic')
        ax.plot(entry['cartesian'], color='#c17c37', ls='--', lw=1.8, label='Cartesian baseline')
        ax.set(title=entry['label'], xlabel='Image index',
               ylabel=entry['kind'].capitalize() + ' (' + entry['unit'] + ')')
        ax.grid(alpha=.15)
        ax.spines[['top', 'right']].set_visible(False)
        if panel == 0:
            ax.legend(frameon=False, fontsize=9)
    fig.supxlabel('%d geodesic images · fixed atom mapping · no energy or gradient calculations' % len(path), fontsize=11)
    fig.savefig(output / 'geometry_summary.png', dpi=170)
    svg = output / 'geometry_summary.svg'
    fig.savefig(svg)
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results' / 'butane')
    parser.add_argument('--nimages', type=int, default=101)
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--maxiter', type=int, default=50)
    parser.add_argument('--data-only', action='store_true', help='Skip report and figure generation')
    args = parser.parse_args()
    symbols, start, end, bonds = endpoints()
    first = gto.M(atom=list(zip(symbols, start)), basis='sto-3g', verbose=0)
    last = first.copy().set_geom_(end, unit='Angstrom')
    np.random.seed(args.seed)
    print('Interpolating %d butane images with backend %s' % (args.nimages, backend_info()), flush=True)
    begin = time.perf_counter()
    images = interpolate(first, last, nimages=args.nimages, maxiter=args.maxiter)
    seconds = time.perf_counter() - begin
    print('Interpolation finished in %.2f s' % seconds, flush=True)
    path = np.array([mol.atom_coords(unit='Angstrom') for mol in images])
    t = np.linspace(0, 1, len(path))[:, None, None]
    baseline = path[0] * (1 - t) + path[-1] * t
    metadata = dict(nimages=len(path), natoms=len(symbols), seed=args.seed,
                    interpolation_seconds=seconds, maxiter=args.maxiter, tolerance=.002,
                    electronic_energy_calls=0, electronic_gradient_calls=0,
                    frames='All requested frames enter geodesic optimization; no display resampling.',
                    endpoint_preparation='Analytic idealized tetrahedral butane, C-C 1.54 A, C-H 1.09 A; '
                                         'trans 180 to gauche +60 deg with fixed atom mapping.',
                    baseline='Cartesian linear interpolation in the geodesic endpoint frames.',
                    connectivity='Fixed reference graph, not inferred bond orders.',
                    backend=backend_info(), python=platform.python_version(),
                    packages={name: version(name) for name in ('pyscf', 'numpy', 'scipy', 'matplotlib')},
                    limitation='No energy-minimized path or certified geodesic convergence.')
    data = finite_json(dict(symbols=symbols, bonds=bonds, metadata=metadata,
                           coordinates={'geodesic': path.round(10).tolist(),
                                        'cartesian': baseline.round(10).tolist()},
                           rmsd=geometry.aligned_rmsd(path).tolist(),
                           measurements=measurements(symbols, bonds, path, baseline)))
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    (output / 'path.xyz').write_text(xyz_text(symbols, path))
    (output / 'cartesian.xyz').write_text(xyz_text(symbols, baseline))
    (output / 'endpoints.xyz').write_text(xyz_text(symbols, [start, end]))
    (output / 'analysis.json').write_text(json.dumps(data, separators=(',', ':'), allow_nan=False) + '\n')
    (output / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    with (output / 'measurements.csv').open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        columns = [(entry, source) for entry in data['measurements'] for source in
                   (('geodesic', 'cartesian', 'wrapped_geodesic', 'wrapped_cartesian')
                    if entry['kind'] == 'dihedral' else ('geodesic', 'cartesian'))]
        writer.writerow(['image', 'aligned_rmsd_A'] +
                        [entry['id'] + '_' + source + ('_A' if entry['kind'] == 'distance' else '_deg')
                         for entry, source in columns])
        for frame in range(len(path)):
            writer.writerow([frame, data['rmsd'][frame]] + [entry[source][frame] for entry, source in columns])
    if not args.data_only:
        from path_report import build_report
        static_figure(data, output)
        build_report(data, output / 'report.html')
    print('Saved %d images and %d measurements in %.2f s to %s' %
          (len(path), len(data['measurements']), seconds, output), flush=True)


if __name__ == '__main__':
    main()
