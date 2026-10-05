"""Read-only Line9 comparison audit; no field fitting, solver, or training.

Aggregates/hashes are reviewable evidence. Field-derived curves/images are
written only to the explicitly separate local, gitignored plot directory.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.signal import hilbert
from field_profile_metrics import late_peak_metrics
from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]
CHECKS = ROOT/'artifacts/research_checks'


def old_metric(env, time):
    mask = time >= 300
    peaks = time[mask][np.argmax(env[mask], axis=0)]
    ratio = env[mask].max(axis=0)/env[time < 150].max(axis=0)
    return {'median_dB': float(np.median(20*np.log10(ratio))),
            'fraction_at_last_sample': float(np.mean(peaks == time[-1])),
            'fraction_in_last50ns': float(np.mean(peaks > time[-1]-50))}


def field(csv):
    with csv.open(encoding='utf-8-sig') as stream:
        head = [stream.readline() for _ in range(4)]
    numbers = [float(s.split('=')[1].split(',')[0]) for s in head]
    ns, span, nt, spacing = int(numbers[0]), numbers[1], int(numbers[2]), numbers[3]
    data = np.loadtxt(csv, delimiter=',', skiprows=4)
    if data.shape != (ns*nt, 5) or not np.isfinite(data).all():
        raise ValueError('finite stacked five-column payload required')
    y = data[:, 3].reshape(nt, ns).T.copy()
    height = data[:, 4].reshape(nt, ns)
    info = {'file_name': csv.name, 'sha256': sha256(csv), 'shape_sample_trace': [ns, nt],
            'time_window_ns': span, 'nominal_trace_spacing_m': spacing,
            'nominal_station_span_m': (nt-1)*spacing,
            'declared_height_median_m': float(np.median(height[:, 0])),
            'declared_height_minmax_m': [float(height.min()), float(height.max())],
            'traces_with_within_trace_height_change': int(np.count_nonzero(np.ptp(height, axis=1))),
            'time_axis_scope': 'Inclusive endpoint is a provisional display convention; original time vector unavailable.',
            'upstream_processing': 'UNKNOWN; Preprocessed directory name is not evidence of gain/background subtraction.'}
    del data, height
    time = np.linspace(0, span, ns)
    old_env = abs(hilbert(y, axis=0))
    metrics = {'historical_unbounded_hilbert': old_metric(old_env, time),
               'hilbert_guard50ns': late_peak_metrics(old_env, time, late_stop_ns=span-50),
               'hilbert_guard100ns': late_peak_metrics(old_env, time, late_stop_ns=span-100),
               'raw_abs_guard50ns': late_peak_metrics(abs(y), time, late_stop_ns=span-50)}
    chosen = None
    for pad in (100, ns, 2*ns):
        env = abs(hilbert(np.pad(y, ((pad, pad), (0, 0)), mode='reflect'), axis=0))[pad:pad+ns]
        metrics[f'reflect_pad{pad}_guard50ns'] = late_peak_metrics(env, time, late_stop_ns=span-50)
        if pad == 2*ns:
            chosen = env
    exclusive = np.arange(ns)*span/ns
    metrics['raw_abs_exclusive_axis_guard50ns'] = late_peak_metrics(abs(y), exclusive, late_stop_ns=span-50)
    event = (time >= 400)&(time <= 500)
    anchors = {}
    for label, env in [('old_unpadded', old_env), ('reflect_pad1002', chosen)]:
        mean = env.mean(axis=1); db = 20*np.log10(mean/mean.max())
        anchors[label] = {'whole_profile_mean_envelope_peak_ns': float(time[np.argmax(mean)]),
                         'mean_log_mean_envelope_in400_500ns_dB': float(db[event].mean()),
                         'definition': 'Mean across times of log(mean across traces); not a per-trace target peak or clean interface response.'}
    info['peak_metrics'] = metrics; info['interface_band_anchors'] = anchors
    return info, time, y, old_env, chosen


def simulation_products():
    from analyze_hs4_height_wavefield import response
    from check_hs4_v4_factor_evidence import direct_response
    from sfcw_official_loader_v0_2 import verify_official_runtime
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    verify_official_runtime()
    capsule = CHECKS/'2026-10-05_hs4_debye_scan_r1'
    cp = capsule/'execution_contract.json'
    contract = json.loads(cp.read_text('utf-8'))
    verified = json.loads((capsule/'completed_verification.json').read_text('utf-8'))
    if verified['status'] != 'PASS' or verified['contract_sha256'] != sha256(cp):
        raise ValueError('complete native simulation identity required')
    rows, products, hashes, dft = {}, {}, {}, {}
    for tag in ('de7.878_tau6.4567_s0.003', 'de0.5_tau6.4567'):
        for role in ('rough', 'halfspace'):
            group = next(g for g in contract['groups'] if g['id'] == tag+'_'+role)
            folder = capsule/group['id']
            audit = json.loads((folder/'audit.json').read_text('utf-8'))
            for row in audit:
                raw = folder/row['file']; digest = sha256(raw)
                if digest != row['sha256']:
                    raise ValueError('native raw differs from archived audit')
                r = response(raw); key = (tag, role, row['station_index']); rows[key] = r
                take = np.array([0, 83, 167, 250, 333, 417, 500])
                exact = direct_response(raw, r.frequency[take])
                error = float(np.linalg.norm(exact-r.response[take])/np.linalg.norm(exact))
                if error > 1e-9:
                    raise ValueError('independent actual-source direct DFT mismatch')
                dft['/'.join((tag, role, row['file']))] = error
                hashes['/'.join((tag, role, row['file']))] = digest
        for si in range(1, 122, 10):
            r = rows[tag, 'rough', si]; h0 = rows[tag, 'halfspace', si]
            products[tag, si, 'total'] = reconstruct_time_response(r, window='hann', zero_pad_factor=8)
            products[tag, si, 'contrast'] = reconstruct_time_response(replace(r, response=r.response-h0.response), window='hann', zero_pad_factor=8)
    return {'contract_sha256': sha256(cp), 'raw_identities': hashes,
            'independent_direct_DFT_max_relative_L2': max(dft.values()),
            'material_changes': 'Candidate changes BOTH sigma .003->.001 and Delta7.878->.5; no unique fitted field material.',
            'geometry': '8m altitude; flat ground; approx3m cover;0.8m basement relief;13 stations over6m;2D line source.',
            'observable': 'Ey/current ((V/m)/A), not measured port S21; actual-source50120-170MHz,200ns tail,Hann/8x.',
            'scope': 'No newly simulated18.5m depth. No measured gain/port calibration/antenna validation.'}, products


def aperture_common_station_check():
    from analyze_hs4_height_wavefield import response
    from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
    aperture = CHECKS/'2026-10-05_hs4_aperture_scan_r1'
    material = CHECKS/'2026-10-05_hs4_material_scan_r1'
    results = {}
    tag = 'eps18.017_s0.003_debye'
    for i in range(13):
        selected = []
        for role in ('rough', 'halfspace'):
            old_folder = material/(tag+'_'+role)
            old_row = json.loads((old_folder/'audit.json').read_text('utf-8'))[i]
            station = i+25
            seg = 1 if station <= 31 else 2
            new_folder = aperture/f'{tag}_{role}_seg{seg}'
            new_row = next(r for r in json.loads((new_folder/'audit.json').read_text('utf-8')) if r['station_index'] == station)
            pair = []
            for folder, row in [(old_folder, old_row), (new_folder, new_row)]:
                p = folder/row['file']
                if sha256(p) != row['sha256']:
                    raise ValueError('common-station raw identity differs')
                pair.append(response(p))
            if not np.allclose(np.array(new_row['source_m'])-np.array(old_row['source_m']), [36, 0, 0], atol=1e-12, rtol=0):
                raise ValueError('stations not physically matched after translation')
            selected.append(pair)
        old = selected[0][0]; new = selected[0][1]
        a = reconstruct_time_response(replace(old, response=old.response-selected[1][0].response), window='hann', zero_pad_factor=8)
        b = reconstruct_time_response(replace(new, response=new.response-selected[1][1].response), window='hann', zero_pad_factor=8)
        time = a.time*1e9; mask = (time >= 80)&(time <= 220)
        results[str(i)] = float(np.linalg.norm(b.real_bandpass[mask]-a.real_bandpass[mask])/np.linalg.norm(a.real_bandpass[mask]))
    return {'old13_stations_in_new61': True, 'matched_interface_window_ns': [80, 220],
            'signed_relative_L2_by_station': results, 'maximum': max(results.values()),
            'interpretation': 'Wider sampling adds stations; it does not physically change the common un-imaged traces. Full-span correlation includes added flat flanks.'}


def plot(folder, time, y, old, reflected, products):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    folder.mkdir(parents=True)
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    ref = old[time < 150].max(axis=0)
    for name, env in [('未延拓Hilbert包络', old), ('镜像延拓后的包络（敏感性对照）', reflected), ('原始绝对幅值（未做Hilbert）', abs(y))]:
        db = np.median(20*np.log10(np.maximum(env/ref, 1e-12)), axis=1)
        ax.plot(time, db, label=name)
    ax.axvspan(650, 700, color='red', alpha=.12, label='末端50ns：排除峰值统计')
    ax.axvspan(400, 500, color='grey', alpha=.1, label='原报告解释条带窗；未经地质标定')
    ax.set(xlim=(280, 700), ylim=(-90, -15), xlabel='CSV显示时间（ns，端点约定待核验）', ylabel='逐道相对早峰的dB中位曲线', title='Line9末端包络伪影：旧−26.6dB来自记录末端')
    ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.savefig(folder/'endpoint_artifact.png', dpi=140); plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), layout='constrained')
    field_db = 20*np.log10(np.maximum(reflected/reflected[time < 150].max(), 1e-12))
    axes[0, 0].imshow(field_db, aspect='auto', cmap='gray', vmin=-80, vmax=0, extent=[0, 2377*.09093, 700, 0])
    axes[0, 0].set_title('实测Line9：2378道，约216m；导出实数A-scan')
    axes[1, 0].imshow(field_db, aspect='auto', cmap='gray', vmin=-80, vmax=-30, extent=[0, 2377*.09093, 700, 0])
    axes[1, 0].set_title('同一实测数据；仅收窄显示色标，未去背景')
    axes[1, 0].set_ylim(550, 250)
    labels = [('de7.878_tau6.4567_s0.003', '归档材料：Δε7.878 / σ0.003'), ('de0.5_tau6.4567', '诊断候选：Δε0.5 / σ0.001')]
    for col, (tag, title) in enumerate(labels, 1):
        total = np.stack([2*abs(products[tag, si, 'total'].complex_envelope) for si in range(1, 122, 10)], axis=1)
        contrast = np.stack([2*abs(products[tag, si, 'contrast'].complex_envelope) for si in range(1, 122, 10)], axis=1)
        st = products[tag, 1, 'total'].time*1e9; take = st <= 600
        for row, data in [(0, total), (1, contrast)]:
            db = 20*np.log10(np.maximum(data/total.max(), 1e-12))
            axes[row, col].imshow(db[take], aspect='auto', cmap='gray', vmin=-80, vmax=0 if row == 0 else -30, extent=[0, 6, st[take][-1], 0])
            axes[row, col].set_title(title+'\n'+('模型总场' if row == 0 else '模型起伏−全覆盖层：理想参考对比'))
        axes[1, col].set_ylim(220, 80)
    for ax in axes.flat:
        ax.set_xlabel('各自测段的名义距离（m）'); ax.set_ylabel('各自接收时间（ns）')
    fig.suptitle('观察口径核查：每列相对自身总场早峰；无绝对标定，时间/场景不能逐像素对应\n实测显示时间未核验；模型8m航高/约3m覆层/6m测段，未复现营山设备及真实地层')
    fig.savefig(folder/'field_vs_model_scope.png', dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--csv', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    ap.add_argument('--plots-dir', required=True, type=Path)
    a = ap.parse_args()
    if a.out.exists() or a.plots_dir.exists():
        raise ValueError('new outputs required')
    if not a.plots_dir.resolve().is_relative_to((ROOT/'artifacts/local_checks').resolve()):
        raise ValueError('field-derived plots must remain under ignored local_checks')
    info, time, y, old, reflected = field(a.csv)
    simulation, products = simulation_products()
    aperture = aperture_common_station_check()
    from analyze_hs4_material_scan import real_anchors
    fixed = real_anchors(a.csv)
    if abs(fixed['late_peak_over_early_peak_db_median']-info['peak_metrics']['hilbert_guard50ns']['median_dB']) > 1e-12:
        raise ValueError('live anchor helper correction not reflected')
    result = {'status': 'COMPLETED_COMPARISON_AUDIT_NOT_PHYSICAL_CALIBRATION',
              'solver_called': False, 'training_called': False, 'field_parameters_fitted': False,
              'code_sha256': sha256(__file__), 'field_csv': info, 'simulation': simulation,
              'common_station_aperture_check': aperture,
              'live_anchor_helper_late_peak_db': fixed['late_peak_over_early_peak_db_median'],
              'retracted_claim': 'Measured late peak -26.6dB and the 30-45dB sim/field gap inferred from it: dominated by FFT Hilbert endpoint artifact.',
              'limits': ['No unique material inversion from amplitude anchors.', 'Field total stripe mean and simulated ideal contrast peak are different estimators.',
                         'CSV is time_real; no second IFFT or frequency reconstruction applied.', 'No hardware-calibrated correspondence or known site geometry.']}
    a.out.mkdir(parents=True)
    (a.out/'summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    plot(a.plots_dir, time, y, old, reflected, products)
    print(json.dumps({'status': result['status'], 'field': info, 'aperture': aperture,
                      'independent_DFT_max': simulation['independent_direct_DFT_max_relative_L2']}, indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
