# SPDX-License-Identifier: Apache-2.0
'''Shared output and plotting for the two worked examples.'''

import argparse
import csv
import json
from importlib.metadata import version
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyscf.geomopt.geodesic import interpolate


def example_parser(name):
    parser = argparse.ArgumentParser(description='Run a small geodesic/PySCF example')
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'results' / name)
    parser.add_argument('--seed', type=int, default=7)
    return parser


def run_example(start, end, args, title):
    # Script-local reproducibility; interpolate() itself does not reseed NumPy.
    np.random.seed(args.seed)
    images = interpolate(start, end, nimages=5)
    method = start.UHF() if start.spin else start.RHF()
    method.conv_tol = 1e-11
    scanner = method.nuc_grad_method().as_scanner()
    rows = []
    for i, mol in enumerate(images):
        energy, gradient = scanner(mol)
        if not scanner.converged:
            raise RuntimeError('SCF did not converge for image %d' % i)
        if not np.isfinite(energy) or not np.isfinite(gradient).all():
            raise RuntimeError('Nonfinite energy or gradient for image %d' % i)
        rows.append([i, float(energy), float(abs(gradient).max()), True])

    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'path.xyz').open('w') as f:
        for i, mol in enumerate(images):
            f.write('%d\nimage %d; coordinates in Angstrom\n' % (mol.natm, i))
            for j, xyz in enumerate(mol.atom_coords(unit='Angstrom')):
                f.write('%s %.12f %.12f %.12f\n' % (mol.atom_pure_symbol(j), *xyz))
    with (output / 'energies.csv').open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['image', 'energy_hartree', 'max_abs_gradient_hartree_per_bohr',
                         'scf_converged'])
        writer.writerows(rows)
    metadata = dict(title=title, seed=args.seed, nimages=len(images), basis=start.basis,
                    charge=start.charge, spin=start.spin, method='UHF' if start.spin else 'RHF',
                    scf_conv_tol=method.conv_tol, coordinate_unit='Angstrom',
                    energy_unit='Hartree', gradient_unit='Hartree/Bohr',
                    python=platform.python_version(),
                    packages={name: version(name) for name in
                              ('pyscf', 'geodesic-interpolate', 'numpy', 'scipy', 'matplotlib')},
                    limitation='Geometric path only; no energy-based path optimization or '
                               'geodesic convergence certification.')
    (output / 'run.json').write_text(json.dumps(metadata, indent=2) + '\n')
    plot_path(images, rows, title, metadata['method'], output)
    print(title)
    print('image   energy / Hartree    max |gradient| / (Hartree/Bohr)')
    for i, energy, gradient, _ in rows:
        print('%3d    %18.10f    %12.6f' % (i, energy, gradient))
    print('Saved XYZ, CSV, metadata and diagrams to %s' % output)


def plot_path(images, rows, title, method, output):
    coords = np.array([m.atom_coords(unit='Angstrom') for m in images])
    colors = ['#dc493a' if images[0].atom_pure_symbol(i) == 'O' else '#3977bd'
              for i in range(images[0].natm)]
    radius = max(np.ptp(coords, axis=(0, 1)).max() * .6, .5)
    center = coords.mean(axis=(0, 1))
    fig = plt.figure(figsize=(12, 6.5), layout='constrained')
    grid = fig.add_gridspec(2, 5, height_ratios=[1.1, 1])
    fig.suptitle(title, fontsize=17, fontweight='bold')
    for i, xyz in enumerate(coords):
        ax = fig.add_subplot(grid[0, i], projection='3d')
        ax.scatter(*xyz.T, c=colors, s=90, depthshade=False)
        for j, point in enumerate(xyz):
            ax.text(*point, '  %s%d' % (images[i].atom_pure_symbol(j), j + 1), fontsize=9)
        ax.set(xlim=(center[0] - radius, center[0] + radius),
               ylim=(center[1] - radius, center[1] + radius),
               zlim=(center[2] - radius, center[2] + radius))
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=18, azim=-65)
        ax.set_axis_off()
        ax.set_title('Image %d' % i, fontsize=11)
    indices = np.array([r[0] for r in rows])
    energy = np.array([r[1] for r in rows])
    ax = fig.add_subplot(grid[1, :3])
    ax.plot(indices, energy - energy[0], 'o-', color='#3977bd')
    ax.set(xlabel='Image index', ylabel='Energy relative to image 0 (Hartree)',
           xticks=indices, title='%s / STO-3G energy scan' % method)
    ax.grid(alpha=.2)
    ax = fig.add_subplot(grid[1, 3:])
    ax.plot(indices, [r[2] for r in rows], 'o-', color='#dc493a')
    ax.set(xlabel='Image index', ylabel='Max |gradient component| (Hartree/Bohr)',
           xticks=indices, title='Electronic energy gradients')
    ax.grid(alpha=.2)
    fig.supxlabel('Geometric interpolation; no energy-based path optimization. '
                  'Atom positions shown at a common scale.', fontsize=10)
    svg = output / 'path_and_energy.svg'
    fig.savefig(svg)
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines()) + '\n')
    fig.savefig(output / 'path_and_energy.png', dpi=160)
    plt.close(fig)
