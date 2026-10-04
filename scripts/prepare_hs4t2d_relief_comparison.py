"""Prepare separate, grid-aligned HS4T2D relief variants; does not run gprMax."""
import argparse
import csv
from decimal import Decimal, ROUND_HALF_UP
import json
from pathlib import Path

import numpy as np

from export_hs4_browser_geometry import Mesh, write_obj
from hs_capsule_identity import read_manifest, verify_file, sha256

ROOT = Path(__file__).resolve().parents[1]
CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_hs4t2d_transect'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists() or args.out.resolve().is_relative_to(CAPSULE.resolve()):
        raise ValueError('use a new directory outside the original capsule')
    records = read_manifest(CAPSULE)
    originals, identities, boxes = [], [], None
    identities.append(verify_file(CAPSULE, records, 'hs4t2d_transect_table.npz'))
    with np.load(CAPSULE / 'hs4t2d_transect_table.npz') as t:
        original_z = t['prof_m'].copy()
    for k in range(1, 122):
        name = f'hs4t2d_t{k:02d}.in'
        identities.append(verify_file(CAPSULE, records, name))
        lines = (CAPSULE / name).read_text('utf-8').splitlines()
        current = [line for line in lines if line.startswith('#box:')]
        if boxes is None:
            boxes = current
        if current != boxes:
            raise ValueError('archived geometry differs across traces')
        originals.append(lines)
    if not np.array_equal([float(line.split()[3]) for line in boxes[1:]], original_z):
        raise ValueError('archive/table mismatch')
    args.out.mkdir(parents=True)
    versions = {}
    for label, factor in [('relief1p5', Decimal('1.5')), ('relief2p0', Decimal('2')),
                          ('flat', Decimal('0'))]:
        # Keep the original extrema midpoint, X footprint and 5 cm grid.
        nominal = [Decimal('9.05') + factor * (Decimal(str(v)) - Decimal('9.05'))
                   for v in original_z]
        changed = np.asarray([float((v / Decimal('.05')).quantize(Decimal('1'), rounding=ROUND_HALF_UP)
                                   * Decimal('.05')) for v in nominal])
        folder = args.out / label
        folder.mkdir()
        for k, lines in enumerate(originals, 1):
            new, bin_index = [], 0
            for line in lines:
                if line.startswith('#box:') and line.split()[-1] == 'cover':
                    words = line.split()
                    words[3] = f'{changed[bin_index]:g}'
                    line = ' '.join(words)
                    bin_index += 1
                # All other input directives, source/receiver positions and precision settings stay unchanged.
                new.append(line)
            if bin_index != 48:
                raise ValueError('unexpected cover bin count')
            (folder / f'hs4t2d_t{k:02d}.in').write_text('\n'.join(new) + '\n', encoding='utf-8')
        rock, cover = Mesh(), Mesh()
        for j, z in enumerate(changed):
            x0, x1 = j * .25, (j + 1) * .25
            rock.quad([(x0,.025,0),(x1,.025,0),(x1,.025,z),(x0,.025,z)])
            cover.quad([(x0,.025,z),(x1,.025,z),(x1,.025,12),(x0,.025,12)])
        write_obj(folder / 'model.obj', [('rock', rock), ('cover', cover)], True)
        (folder / 'hs4_scene.mtl').write_text(
            'newmtl rock\nKd .45 .55 .65\nd 1\nillum 1\n\n'
            'newmtl cover\nKd .75 .58 .32\nd 1\nillum 1\n', encoding='ascii')
        with (folder / 'profile.csv').open('x', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['x0_m','x1_m','original_z_m','nominal_z_m','grid_z_m','cover_depth_m'])
            for j, z in enumerate(changed):
                w.writerow([j*.25,(j+1)*.25,original_z[j],str(nominal[j]),z,12-z])
        versions[label] = {'factor':str(factor),'z_range_m':[float(changed.min()),float(changed.max())],
                           'relief_m':float(np.ptp(changed)),
                           'mean_z_m':float(changed.mean()),
                           'mean_z_change_m':float(changed.mean()-original_z.mean()),
                           'max_adjacent_bin_jump_m':float(np.abs(np.diff(changed)).max()),
                           'z_m':changed.tolist()}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(11,4), constrained_layout=True)
    ax.stairs(12-original_z, np.arange(49)*.25, baseline=None,
              label='Executed baseline (0.40 m relief)', color='gray')
    for label, color in [('relief1p5','tab:blue'),('relief2p0','tab:orange')]:
        ax.stairs(12-np.asarray(versions[label]['z_m']), np.arange(49)*.25, baseline=None,
                  label=f"{label}: {versions[label]['relief_m']:.2f} m relief (NOT RUN)", color=color)
    ax.axhline(2.95, color='black', ls=':', label='Flat control: 2.95 m depth (NOT RUN)')
    ax.set_xlabel('Profile X (m)')
    ax.set_ylabel('Cover depth below Z=12 m surface (m)')
    ax.invert_yaxis()
    ax.legend()
    ax.set_title('Input geometry comparison; common axes; step geometry retained')
    fig.savefig(args.out / 'geometry_comparison.png', dpi=150)
    plt.close(fig)
    for entry in identities:
        verify_file(CAPSULE, records, entry['file'])
    proposal = {'status':'PREPARED_NOT_RUN', 'solver_executed':False,
                'purpose':'Increase interface relief modestly and compare B-scan with executed baseline and flat control',
                'recommended_variant':'relief1p5','reference_z_m':9.05,
                'rounding':'5 cm cells, Decimal ROUND_HALF_UP; nominal and grid coordinates separate',
                'original_z_range_m':[float(original_z.min()),float(original_z.max())],
                'original_mean_z_m':float(original_z.mean()),
                'invariant_fields':'Every input directive except cover-box bottom Z is unchanged; 121 traces, Ey, 600 ns, Debye, 15 m standoff, 5 cm grid',
                'variants':versions,'inputs':identities,'generator_sha256':sha256(__file__),
                'allowed_claims':'Geometry preparation only. Larger relief does not guarantee waveform correspondence.',
                'processing_plan':'Same pinned official 20-170 MHz, 501 points, Hann chain. Shared color scale. Raw, residual and complex flat-reference difference reported separately.',
                'execution_plan':'Freeze selected variant and local solver identity; reproduce archived t61 in FP64; then 121 selected-variant traces and 121 matched flat traces. No reuse of incompatible historical flat BG.'}
    (args.out / 'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n', encoding='utf-8')
    files = [{'file':p.relative_to(args.out).as_posix(),'bytes':p.stat().st_size,'sha256':sha256(p)}
             for p in sorted(args.out.rglob('*')) if p.is_file()]
    (args.out / 'manifest.json').write_text(json.dumps(files,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:{a:b for a,b in v.items() if a!='z_m'} for k,v in versions.items()},indent=2))


if __name__ == '__main__':
    main()
