"""Bounded official 4.0.1 free-space test; no Line9 model changes.

Prepare once; explicitly run selected cases; analyse completed pairs.
This tests exterior propagation only, not Wang2026's incident-surface solver.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import re
import shutil
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / 'artifacts/research_checks/2026-10-10_air_domain_ksir_r1'
PRIVATE = ROOT / 'artifacts/local_checks/2026-10-10_air_domain_ksir_r1'
CASES = {f'{domain}_{mesh}': (domain, dl)
         for mesh, dl in [('d10', .1), ('d05', .05)]
         for domain in ['full', 'compact']}
BASE_CASES = dict(CASES)
CASES.update({name + '_t480': values for name, values in BASE_CASES.items()})


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def prepare():
    if (PUBLIC / 'contract.json').exists():
        raise RuntimeError('Refusing to overwrite frozen study')
    entries = []
    for name, (domain, dl) in BASE_CASES.items():
        z = 2 + dl / 2
        length = 12 if domain == 'full' else 4
        lines = [f'#title: Air exterior diagnostic {name}',
                 f'#domain: {length} 4 4', f'#dx_dy_dz: {dl} {dl} {dl}',
                 '#time_window: 120e-9', f'#pml_cells: {round(1 / dl)}',
                 '#waveform: impulse 1 1 impulse',
                 '#hertzian_dipole: z 2 2 2 impulse',
                 '#ntff_surface: 1.4 1.4 1.4 2.6 2.6 2.6 surface',
                 f'#ksir_time_rx: 6 2 {z} surface exterior4 Ez simulation',
                 f'#ksir_time_rx: 10 2 {z} surface exterior8 Ez simulation']
        if domain == 'full':
            lines += ['#rx: 6 2 2 direct4 Ez', '#rx: 10 2 2 direct8 Ez']
        path = PUBLIC / name / 'model.in'
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        entries.append(dict(id=name, input=str(path.relative_to(ROOT)), input_sha256=sha(path),
                            dl_m=dl, cells=round(length / dl)*round(4 / dl)**2,
                            domain_m=[length, 4, 4], native_Ez_position_m=[2, 2, z]))
    save(PUBLIC / 'contract.json', dict(
        authorization='User 2026-10-10: 本地仿真研究一轮，空气传播省算研究',
        scope='Four small homogeneous-air 3D diagnostic cases; no ground, no new production route',
        solver='official gprMax 4.0.1 serial CPU double, OMP_NUM_THREADS=4',
        waveform='official impulse 1 1 impulse; source start 0',
        pml_physical_thickness_m=1, surface_bounds_m=[[1.4]*3, [2.6]*3],
        observation_distances_m=[4, 8], time_window_s=120e-9,
        sfcw=dict(method='official direct', start_hz=20e6, stop_hz=170e6,
                  step_hz=.3e6, count=501, tail_taper_fraction=0,
                  reconstruction_window='gaussian', gaussian_sigma=.2,
                  time_shift_s=0, fitted_amplitude_or_delay=False),
        comparison='Same-mesh direct Yee Ez vs physical-position KSIR Ez; same-mesh full vs compact KSIR',
        caveats=['No underground scattering or incident-surface injection tested',
                 '2D Line9 incompatible with official 3D-only KSIR',
                 'Closed surface must enclose all scattering in homogeneous lossless exterior',
                 'Impulse is broadband; evaluate only specified band and check mesh sensitivity'],
        cases=entries))
    print('Prepared four cases; no solver launched')


def prepare_extension():
    target = PUBLIC / 'time_extension_contract.json'
    if target.exists():
        raise RuntimeError('Refusing to overwrite time extension')
    original = json.loads((PUBLIC / 'contract.json').read_text(encoding='utf-8'))
    entries = []
    for case in original['cases']:
        old = ROOT / case['input']
        assert sha(old) == case['input_sha256']
        path = PUBLIC / (case['id'] + '_t480') / 'model.in'
        path.parent.mkdir(parents=True, exist_ok=True)
        content = old.read_text(encoding='utf-8').replace('#time_window: 120e-9', '#time_window: 480e-9')
        path.write_text(content, encoding='utf-8')
        entries.append(dict(case, id=case['id'] + '_t480', input=str(path.relative_to(ROOT)),
                            input_sha256=sha(path), original_input_sha256=sha(old)))
    save(target, dict(original, cases=entries, time_window_s=480e-9,
                      scope='Four additional time-window diagnostics; other case inputs unchanged',
                      reason='120ns terminal decay failed and fine-grid response did not converge; retain failures'))
    print('Prepared four separate 480ns cases; original inputs retained')


def run_case(name):
    contract_path = 'time_extension_contract.json' if name.endswith('_t480') else 'contract.json'
    contract = json.loads((PUBLIC / contract_path).read_text(encoding='utf-8'))
    case = next(c for c in contract['cases'] if c['id'] == name)
    inp = ROOT / case['input']
    assert sha(inp) == case['input_sha256']
    directory = PRIVATE / name
    if directory.exists():
        raise RuntimeError('Refusing to repeat an existing attempt')
    directory.mkdir(parents=True)
    output = directory / 'model'
    command = [sys.executable, '-m', 'gprMax', str(inp), '-cpu_precision', 'double', '-o', str(output)]
    env = dict(os.environ, OMP_NUM_THREADS='4', MPLBACKEND='Agg')
    record = dict(case=name, command=command, input_sha256=sha(inp),
                  started_utc=datetime.now(timezone.utc).isoformat(), python=sys.executable)
    save(directory / 'execution.json', record)
    start = time.perf_counter()
    with (directory / 'solver.log').open('w', encoding='utf-8') as log:
        result = subprocess.run(command, cwd=ROOT, env=env, stdout=log,
                                stderr=subprocess.STDOUT, timeout=900)
    record.update(exit_code=result.returncode, wall_seconds=time.perf_counter()-start,
                  finished_utc=datetime.now(timezone.utc).isoformat())
    record['solver_log_sha256'] = sha(directory / 'solver.log')
    files = list(directory.glob('*.h5'))
    record['outputs'] = [dict(path=str(p.relative_to(ROOT)), sha256=sha(p), bytes=p.stat().st_size) for p in files]
    save(directory / 'execution.json', record)
    print(json.dumps(record, ensure_ascii=False))
    if result.returncode or len(files) != 1:
        raise RuntimeError('Solver failed or ambiguous output; inspect log')


def analyse():
    import h5py
    import numpy as np
    from gprMax.toolboxes.SFCW import processing as sf
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    frequencies = 20e6 + .3e6 * np.arange(501)
    spectra, rows, runs, direct_audits = {}, [], [], []

    def independent_audit(fr):
        selected = np.array([0, 1, 100, 250, 400, 500])
        def dft(signal):
            return np.exp(-2j*np.pi*frequencies[selected, None]*signal.times) @ signal.samples * signal.dt
        calculated = dft(fr.receiver) / dft(fr.source)
        error = np.linalg.norm(calculated-fr.response[selected])/np.linalg.norm(calculated)
        assert error < 1e-11
        direct_audits.append(dict(path=fr.receiver.path, relative_error=float(error)))
    for name in CASES:
        directory = PRIVATE / name
        if not (directory / 'execution.json').exists():
            continue
        execution = json.loads((directory / 'execution.json').read_text(encoding='utf-8'))
        if execution.get('exit_code') != 0:
            continue
        raw = ROOT / execution['outputs'][0]['path']
        assert sha(raw) == execution['outputs'][0]['sha256']
        assert sha(execution['command'][3]) == execution['input_sha256']
        archive = PUBLIC / 'native' / (name + '.h5')
        archive.parent.mkdir(parents=True, exist_ok=True)
        if archive.exists():
            assert sha(archive) == sha(raw)
        else:
            shutil.copyfile(raw, archive)
        execution['public_native'] = str(archive.relative_to(ROOT))
        source = sf.load_source(raw)
        log = (directory / 'solver.log').read_text(encoding='utf-8', errors='replace')
        execution['solver_log_sha256'] = sha(directory / 'solver.log')
        execution['solver_duration_text'] = re.findall(r'Time taken: (.+)', log)[-1]
        execution['estimated_memory_MB'] = float(re.findall(r'Memory used \(estimated\): ~([\d.]+) MB', log)[-1])
        runs.append(execution)
        with h5py.File(raw) as h:
            assert str(h.attrs['gprMax']) == '4.0.1'
            expected_dl = CASES[name][1]
            assert np.allclose(h.attrs['dx_dy_dz'], expected_dl, rtol=0, atol=1e-14)
            assert source.samples.dtype == np.float64
            assert h['srcs/src1/excitation/samples'].dtype == np.float64
            assert np.count_nonzero(source.samples) == 1 and source.samples[0] == 1
            assert source.time_offset == .5*float(h.attrs['dt'])
            execution['native_source_time_offset_s'] = source.time_offset
            execution['native_iterations'] = int(h.attrs['Iterations'])
            execution['native_cells'] = int(np.prod(h.attrs['nx_ny_nz']))
            paths = []
            h.visititems(lambda n, o: paths.append(n) if isinstance(o, h5py.Group)
                         and 'fully_supported_lengths' in o else None)
            for path in paths:
                group = h[path]
                label = str(group.attrs.get('ID', path.rsplit('/', 1)[-1]))
                distance = 4 if '4' in label else 8
                fully_supported = int(group['fully_supported_lengths'][0])
                # Compare the same absolute physical observation interval as a native Rx.
                # Retarded KSIR buffers extend beyond the FDTD iteration count.
                length = min(fully_supported, int(h.attrs['Iterations']))
                data = group['fields/Ez'][0, :length]
                assert data.dtype == np.float64 and np.isfinite(data).all()
                times = group['times'][:length] + group['time_origins'][0]
                signal = sf.SampledSignal(path=path, samples=data, dt=float(times[1]-times[0]),
                                          time_offset=float(times[0]), quantity='Ez', units='V/m')
                fr = sf.direct_frequency_response(source, signal, frequencies, tail_taper_fraction=0)
                assert fr.source_valid.all()
                independent_audit(fr)
                spectra[name, distance, 'ksir'] = fr
                rows.append(dict(case=name, distance_m=distance, kind='ksir', raw_sha256=sha(raw),
                                 dtype=str(data.dtype), shape=list(data.shape), dt_s=signal.dt,
                                 time_offset_s=signal.time_offset, fully_supported_length=fully_supported,
                                 common_record_length=length,
                                 tail_db=fr.receiver_tail_relative_db,
                                 terminal_decay_ok=bool(group['terminal_decay_ok'][0])))
            if name.startswith('full'):
                for distance, rx in [(4, 'rx1'), (8, 'rx2')]:
                    assert np.allclose(h['rxs/' + rx].attrs['Position'], [2+distance, 2, 2])
                    signal = sf.load_receiver(raw, '/rxs/' + rx, 'Ez')
                    assert h['rxs/' + rx + '/Ez'].dtype == np.float64
                    fr = sf.direct_frequency_response(source, signal, frequencies, tail_taper_fraction=0)
                    independent_audit(fr)
                    spectra[name, distance, 'direct'] = fr
                    rows.append(dict(case=name, distance_m=distance, kind='direct', dtype='float64',
                                     shape=list(signal.samples.shape), dt_s=signal.dt,
                                     time_offset_s=signal.time_offset, tail_db=fr.receiver_tail_relative_db))
    comparisons, convergence = [], []
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for col, distance in enumerate([4, 8]):
        for mesh, style in [('d10', ':'), ('d05', '--'), ('d10_t480', '-.'), ('d05_t480', '-')]:
            label_mesh = f"网格{10 if mesh.startswith('d10') else 5}cm，{480 if 't480' in mesh else 120}ns"
            key = ('full_' + mesh, distance, 'direct')
            if key not in spectra or ('compact_' + mesh, distance, 'ksir') not in spectra:
                continue
            reference = spectra[key]
            for domain in ['full', 'compact']:
                fr = spectra[domain + '_' + mesh, distance, 'ksir']
                ratio = fr.response / reference.response
                band_errors = {label: float(np.linalg.norm((fr.response-reference.response)[sel])/np.linalg.norm(reference.response[sel]))
                               for label, sel in [('20_40MHz', frequencies<=40e6), ('40_170MHz', frequencies>40e6)]}
                comparisons.append(dict(mesh=mesh, distance_m=distance, comparison=domain+'_ksir_vs_full_direct',
                    complex_relative_l2=float(np.linalg.norm(fr.response-reference.response)/np.linalg.norm(reference.response)),
                    band_complex_relative_l2=band_errors,
                    amplitude_ratio_min=float(np.min(abs(ratio))), amplitude_ratio_max=float(np.max(abs(ratio))),
                    phase_error_max_deg=float(np.max(abs(np.angle(ratio, deg=True))))))
            compact = spectra['compact_' + mesh, distance, 'ksir']
            fullksir = spectra['full_' + mesh, distance, 'ksir']
            comparisons.append(dict(mesh=mesh, distance_m=distance, comparison='compact_vs_full_ksir',
                complex_relative_l2=float(np.linalg.norm(compact.response-fullksir.response)/np.linalg.norm(fullksir.response))))
            for fr, label in [(reference, '完整空气域直接接收'), (compact, '缩小空气域KSIR重建')]:
                tr = sf.reconstruct_time_response(fr, window='gaussian', gaussian_sigma=.2, zero_pad_factor=32)
                sel = tr.time < 100e-9
                # Per unit current moment removes the declared edge-length change across meshes.
                if 't480' in mesh:
                    axes[0, col].plot(tr.time[sel]*1e9, abs(tr.complex_envelope[sel])/fr.source.spatial_scale,
                                      style, label=f'{label}，{label_mesh}')
            axes[1, col].plot(frequencies/1e6, np.angle(compact.response/reference.response, deg=True),
                              style, label=label_mesh)
        axes[0, col].axvline(distance/299792458*1e9, color='grey', alpha=.6, label='几何空气到时 R/c')
        axes[0, col].set(title=f'距源{distance}米：官方direct SFCW包络，保留物理时延',
                         xlabel='从发射起算的时间（ns）', ylabel='包络 / 源电流矩 [(V/m)/(A·m)]')
        axes[1, col].set(title=f'距源{distance}米：缩域重建相对直接接收的相位误差',
                         xlabel='频率（MHz）', ylabel='相位差（度）')
        for ax in axes[:, col]:
            ax.legend(fontsize=8); ax.grid(alpha=.2)
        for mesh in ['d10', 'd05']:
            for domain, kind in [('full', 'direct'), ('full', 'ksir'), ('compact', 'ksir')]:
                short = spectra.get((domain+'_'+mesh, distance, kind))
                long = spectra.get((domain+'_'+mesh+'_t480', distance, kind))
                if short is not None and long is not None:
                    convergence.append(dict(mesh=mesh, domain=domain, kind=kind, distance_m=distance,
                        time_window_120_vs_480_complex_relative_l2=float(np.linalg.norm(short.response-long.response)/np.linalg.norm(long.response))))
    fig.suptitle('空气传播省算诊断：均匀空气三维偶极；两个接收点，不是测线B-scan\n未包含地表、地下反射、天线实体或论文虚拟入射面')
    fig.tight_layout()
    fig.savefig(PUBLIC / 'air_domain_comparison.png', dpi=150)
    plt.close(fig)
    report = dict(audit=rows, comparisons=comparisons, time_window_sensitivity=convergence,
                  independent_six_tone_direct_DFT_audits=direct_audits, execution=runs,
                  completed_cases=len(runs), official_processing_sha256=sha(Path(sf.__file__)),
                  analysis_script_sha256=sha(Path(__file__)),
                  processing_settings=dict(method='official_direct', frequency_count=501,
                      start_hz=20e6, step_hz=.3e6, stop_hz=170e6,
                      tail_taper_fraction=0, common_absolute_record=True,
                      reconstruction_window='gaussian', gaussian_sigma=.2,
                      zero_pad_factor=32, normalise_window=True, time_shift_s=0),
                  analytic_air_delays_ns={str(r): r/299792458*1e9 for r in [4, 8]},
                  claims='Free-space exterior reconstruction only; not underground model validation')
    save(PUBLIC / 'results.json', report)
    arrays = {'frequency_hz': frequencies}
    arrays.update({f'{name}_r{distance}_{kind}': fr.response for (name, distance, kind), fr in spectra.items()})
    np.savez_compressed(PUBLIC / 'spectra.npz', **arrays)
    print(json.dumps(dict(completed_cases=len(runs), comparisons=comparisons), ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['prepare', 'prepare-extension', 'run', 'analyse'])
    parser.add_argument('--case', choices=CASES)
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare()
    elif args.action == 'prepare-extension':
        prepare_extension()
    elif args.action == 'run':
        if not args.case:
            parser.error('--case is required')
        run_case(args.case)
    else:
        analyse()
