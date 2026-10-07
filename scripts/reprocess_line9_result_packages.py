"""Reprocess audited Line9 natives into phase-preserving SFCW HDF5; no FDTD.

Private numerical products go to --out; aggregate report and labelled figures
go to --report-out. Both must be fresh directories. Raw packages are read-only.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from review_line9_result_packages import FREQ, sha, save, rel


def process_package(path, audit):
    manifest_path = path / 'package_manifest.json'
    manifest = json.loads(manifest_path.read_text('utf-8'))
    planned = manifest['cases']
    planned_ids = [p['id'] for p in planned]
    rows = audit['records']
    ids = [r['id'] for r in rows]
    if len(set(planned_ids)) != len(planned_ids) or len(set(ids)) != len(ids):
        raise ValueError('Duplicate station identity')
    if not set(ids) <= set(planned_ids) or len(ids) != audit['raw_count']:
        raise ValueError('Completed station inventory mismatch')
    declared = {p['id']: p for p in planned}
    if any(r['chainage_m'] != declared[r['id']]['chainage_m'] for r in rows):
        raise ValueError('Station coordinate mismatch')
    if set(p.parent.name for p in (path/'cases').glob('*/profile.h5')) != set(ids):
        raise ValueError('Unreviewed or missing native trace')
    material = next((path/'geometries').glob('*.json'))
    geometry = next((path/'geometries').glob('*.h5'))
    if sha(material) != audit['materials_sha256'] or sha(geometry) != audit['geometry_sha256']:
        raise ValueError('Model identity changed')
    histories = []
    manifest_x_offsets = set()
    source = receiver = None
    for row in rows:
        raw = path / 'cases' / row['id'] / 'profile.h5'
        if sha(raw) != row['native_sha256']:
            raise ValueError('Original native changed: ' + row['id'])
        card = raw.with_suffix('.in')
        if sha(card) != row['input_sha256']:
            raise ValueError('Reviewed input card changed')
        commands = {s.split(':', 1)[0][1:]: s.split(':', 1)[1].strip().split()
                    for s in card.read_text('utf-8').splitlines() if s.startswith('#') and ':' in s}
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax']) != '4.0.0' or str(h['rxs/rx1/Ez'].dtype) != audit['native_dtype']:
                raise ValueError('Native version or precision mismatch')
            # TMz collapses the nominal half-cell z coordinate to its only plane.
            if int(h.attrs['nx_ny_nz'][2]) != 1:
                raise ValueError('This entrypoint is for the reviewed 2D packages')
            for name, key in [('srcs/src1', 'tx_m'), ('rxs/rx1', 'rx_m')]:
                actual = h[name].attrs['Position']
                card_pos = commands['hertzian_dipole'][1:3] if key == 'tx_m' else commands['rx'][:2]
                np.testing.assert_allclose(actual[:2], np.array(card_pos, float), rtol=0, atol=1e-10)
                # pkg2's manifest uses global x, while its card/native use cropped x.
                # Preserve and report the offset; do not use manifest x to index local H5.
                manifest_x_offsets.add(round(float(declared[row['id']][key][0]-actual[0]), 9))
                if abs(declared[row['id']][key][1]-actual[1]) > 1e-10:
                    raise ValueError('Manifest height differs from native')
                if actual[2] != 0. or declared[row['id']][key][2] != .0125:
                    raise ValueError('Unexpected collapsed 2D plane coordinate')
        src = sf.load_source(raw)
        rx = sf.load_receiver(raw, '/rxs/rx1', 'Ez')
        if src.spatial_scale != .025 or src.source_type.lower().find('hertzian') < 0:
            raise ValueError('Expected reviewed 0.025m Hertzian source')
        if not np.isfinite(src.samples).all() or not np.isfinite(rx.samples).all():
            raise ValueError('Nonfinite native data')
        if source is None:
            source, receiver = src, rx
        else:
            np.testing.assert_array_equal(src.samples, source.samples)
            np.testing.assert_array_equal(src.times, source.times)
            np.testing.assert_array_equal(rx.times, receiver.times)
            if src.spatial_scale != source.spatial_scale:
                raise ValueError('Source spatial normalisation differs')
        histories.append(rx.samples)
    if len(manifest_x_offsets) != 1:
        raise ValueError('Manifest/native coordinate offset is inconsistent')
    # Use the actual source samples and native Yee offsets, never an assumed zero-time impulse.
    receiver = replace(receiver, samples=np.column_stack(histories))
    product = sf.direct_frequency_response(source, receiver, FREQ, tail_taper_fraction=0)
    if not product.source_valid.all():
        raise ValueError('Source spectrum is not supported at every tone')
    product = replace(product, response=product.response/source.spatial_scale)
    # Independent exact off-FFT-grid sum for every trace at seven distributed tones.
    f = FREQ[[0, 83, 167, 250, 333, 417, 500]]
    xf = source.dt * (np.exp(-2j*np.pi*f[:, None]*source.times) @ source.samples)
    yf = receiver.dt * (np.exp(-2j*np.pi*f[:, None]*receiver.times) @ receiver.samples)
    expected = yf/xf[:, None]/source.spatial_scale
    error = rel(product.response[[0, 83, 167, 250, 333, 417, 500]], expected)
    if error > 1e-9:
        raise ValueError('Independent complex DFT mismatch')
    profiles = {}
    inverse_errors = {}
    for window in ['hann', 'gaussian', 'blackman']:
        tr = sf.reconstruct_time_response(product, window=window, zero_pad_factor=8)
        take = np.arange(0, len(tr.time), 97)
        expected = np.exp(2j*np.pi*tr.time[take, None]*FREQ) @ (tr.weights[:, None]*product.response)/len(FREQ)
        e = rel(tr.complex_bandpass[take], expected)
        if e > 1e-9:
            raise ValueError('Independent inverse mismatch')
        profiles[window] = tr
        inverse_errors[window] = e
    return product, profiles, dict(independent_DFT_relative_L2=error,
        independent_inverse_relative_L2=inverse_errors, planned=planned,
        manifest_sha256=sha(manifest_path),
        manifest_x_minus_native_x_m=manifest_x_offsets.pop())


def main(root, review, out, report_out):
    if out.resolve() == report_out.resolve() or out.exists() or report_out.exists():
        raise ValueError('Use two distinct fresh output directories')
    audits = sorted(review.glob('*_audit.json'))
    if len(audits) != 3:
        raise ValueError('Expected all three reviewed packages')
    out.mkdir(parents=True)
    report_out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    report = dict(calls_solver=False, calls_training=False, script_sha256=sha(__file__),
        official_processing_sha256=sha(sf.__file__), frequency_MHz=[20., 170.], step_MHz=.3,
        frequency_count=501, primary_window='hann', diagnostic_windows=['gaussian', 'blackman'],
        no_background_removal=True, no_AGC=True, no_trace_normalisation=True,
        no_time_partition=True, no_tail_taper=True, no_layer_guided_signal_extraction=True,
        no_extra_time_shift=True, packages=[])
    for ap in audits:
        audit = json.loads(ap.read_text('utf-8'))
        product, profiles, evidence = process_package(root/audit['package'], audit)
        destination = out/(audit['package']+'_corrected_sfcw.h5')
        sf.write_sfcw_output(destination, product, profiles['hann'])
        rows = audit['records']
        with h5py.File(destination, 'a') as h:
            h.attrs['NativeDtype'] = audit['native_dtype']
            h.attrs['DerivedDtypeDoesNotUpgradeNativePrecision'] = True
            h.attrs['SourceSpatialNormalisationApplied'] = True
            h.attrs['ResponseUnits'] = '(V/m)/(A*m); 2D field proxy, not port S21'
            h.attrs['NoBackgroundOrGainOrTimePartition'] = True
            h.attrs['PlotReference'] = 'Gaussian full-record peak, per-package, no independent panel normalisation'
            h.attrs['AuditSHA256'] = sha(ap)
            h.attrs['ManifestXMinusNativeXM'] = evidence['manifest_x_minus_native_x_m']
            h.attrs['ScriptSHA256'] = sha(__file__)
            h.attrs['OfficialProcessingSHA256'] = sha(sf.__file__)
            text_type = h5py.string_dtype('utf-8')
            h.create_dataset('station_id', data=[r['id'] for r in rows], dtype=text_type)
            h.create_dataset('native_sha256', data=[r['native_sha256'] for r in rows], dtype=text_type)
            h.create_dataset('chainage_m', data=[r['chainage_m'] for r in rows])
            h.create_dataset('planned_station_id', data=[p['id'] for p in evidence['planned']], dtype=text_type)
            h.create_dataset('planned_chainage_m', data=[p['chainage_m'] for p in evidence['planned']])
            h.create_dataset('planned_completed_mask', data=[p['id'] in {r['id'] for r in rows} for p in evidence['planned']])
            for window in ['gaussian', 'blackman']:
                tr = profiles[window]
                g = h.create_group('diagnostic_'+window)
                g.attrs['Window'] = window
                for name in ['weights', 'complex_envelope', 'complex_bandpass', 'real_bandpass']:
                    g.create_dataset(name, data=getattr(tr, name), compression='gzip')
        # Write-read verification of the actual deliverable, not only in-memory transforms.
        with h5py.File(destination) as h:
            np.testing.assert_array_equal(h['frequency'][:], FREQ)
            np.testing.assert_array_equal(h['response'][:], product.response)
            np.testing.assert_array_equal(h['time_response/real_bandpass'][:], 2*h['time_response/complex_bandpass'][:].real)
            np.testing.assert_allclose(abs(h['time_response/complex_envelope'][:]), abs(h['time_response/complex_bandpass'][:]), rtol=1e-14, atol=0)
            assert int(h['planned_completed_mask'][:].sum()) == len(rows)
        fig, axes = plt.subplots(2, 3, figsize=(17, 10), layout='constrained')
        positions = np.array([r['chainage_m'] for r in rows])
        order = np.argsort(positions)
        xe = np.r_[positions[order][0]-.1, (positions[order][:-1]+positions[order][1:])/2, positions[order][-1]+.1]
        reference = float(np.max(abs(profiles['gaussian'].complex_bandpass)))
        for col, name in enumerate(['gaussian', 'hann', 'blackman']):
            tr = profiles[name]
            for row, (lo, hi) in enumerate([(0, 600), (120, 450)]):
                keep = (tr.time*1e9 >= lo) & (tr.time*1e9 <= hi)
                tt = tr.time[keep]*1e9
                te = np.r_[tt[0]-(tt[1]-tt[0])/2, (tt[:-1]+tt[1:])/2, tt[-1]+(tt[1]-tt[0])/2]
                db = 20*np.log10(np.maximum(abs(tr.complex_bandpass[keep][:, order])/reference, 1e-12))
                pic = axes[row, col].pcolormesh(xe, te, db, cmap='gray_r', vmin=-110, vmax=-35, rasterized=True)
                axes[row, col].set(xlabel='原测线里程 / m（仅完成区）', ylabel='时间 / ns', ylim=(hi, lo),
                    title=f'{name}总场；'+('0–600ns' if row == 0 else '120–450ns局部；相同数值/色标'))
                axes[row, col].invert_xaxis()
                fig.colorbar(pic, ax=axes[row, col], label='复幅度dB / 同一Gaussian总场峰值')
        fig.suptitle(audit['package']+f'：{len(rows)}/{len(evidence["planned"])}道；native '+audit['native_dtype']+
            '\n修正SFCW交付：20–170MHz/0.3MHz/501点；Hann为当前诊断基线，另两窗对照\n无AGC/去背景/包络相减/叠道/时间分割；未算区域不补道；本图无模型参考线')
        pic_path = report_out/(audit['package']+'_corrected_sfcw.png')
        fig.savefig(pic_path, dpi=125)
        plt.close(fig)
        evidence.pop('planned')
        item = dict(package=audit['package'], native_count=len(rows), planned_count=audit['planned_count'],
            native_dtype=audit['native_dtype'], native_hashes_verified=True, audit_sha256=sha(ap),
            private_hdf5_path=str(destination.resolve()), private_hdf5_sha256=sha(destination),
            plot_sha256=sha(pic_path), global_gaussian_reference=reference, **evidence)
        report['packages'].append(item)
        print(json.dumps(item, ensure_ascii=False), flush=True)
    report['limits'] = 'Preserves complex numerical data, not a physical validation or target recovery certificate. Hann/Blackman are explicit diagnostic processing choices, not verified instrument settings. Private originals unchanged; no FP32-to-FP64 solver upgrade.'
    save(report_out/'corrected_processing_report.json', report)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'out', 'report-out']:
        p.add_argument('--'+name, type=Path, required=True)
    args = p.parse_args()
    main(args.root, args.review, args.out, args.report_out)
