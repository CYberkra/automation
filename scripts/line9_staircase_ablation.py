"""Approved pure-staircase isolation ablation: 3 surfaces x 3 grids, homogeneous cover; never retried.

Design v0.2 (user-approved 2026-10-11, "批准，开跑"): isolate the surface voxelisation
(staircase) error from physical roughness and from non-staircase grid effects.
- Single air/cover interface; cover is the production Debye cover copied verbatim from
  line9_materials_v0_2 (sha256 fdf02830...); no layering, no landslide body (H0-equivalent).
- Domain 40 x 18 m, time window 900 ns (bands [20,120] and [150,450] ns only; the deep band
  is not evaluated because the homogeneous half-space has no deep structure).
- Surfaces: flat (control, zero staircase), 5-degree planar slope (pure staircase),
  fBm micro-relief (same statistical family as the S5 batch: spectral synthesis H=0.5,
  seed 20261011, RMS target 0.1 m, clip +/-0.2 m, defined continuous at 0.625 cm then
  block-averaged to coarser columns; a new realisation at this domain length, disclosed).
- Grids 2.5 / 1.25 / 0.625 cm; the 0.625 cm solve is the proxy smooth reference (not an
  analytic truth). Same continuous surface, same quantisation rule (round to cell), identical
  absolute source/rx coordinates (multiples of 0.625 cm), identical physical PML thickness
  (2 m -> 80/160/320 cells).
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
C0 = 299792458.
DOMAIN = (40.0, 18.0)
SURFACE_MEAN = 8.0
ANT_Y = 13.0
SRC_X, RX_X = 19.35, 20.65
TIME_WINDOW = 900e-9
PML_PHYS = 2.0
SLOPE_DEG = 5.0
SEED, HURST, RMS_M, CLIP_M = 20261011, 0.5, 0.1, 0.2
GRIDS = {'d25': 0.025, 'd125': 0.0125, 'd625': 0.00625}
DB_ID = 'staircase_ablation_materials_v0_1'
PROD_DB_SHA256 = 'fdf02830012ea2fe4dfe8820096491a2bafd0e8729e0abb5db63cfb85b98d43d'
APPROVAL = ('2026-10-11 用户显性批准（“批准，开跑”）：纯阶梯误差隔离消融 v0.2。均匀覆盖层半空间'
            '（生产 Debye cover 逐字复制）、40x18m 域、900ns、平/5度坡/fBm 微起伏三种地表、'
            '2.5/1.25/0.625cm 三网格共 9 解，H0-only，无滑体无分层；0.625cm 解为代理光滑参考。'
            '目的：分离地表体素化误差与真实起伏物理效应及非阶梯网格效应，服务真实地形导入的网格决策。')


def dt_of(dl):
    return dl / (C0 * np.sqrt(2))


def fractal_profile(n, hurst=HURST, seed=SEED):
    """1D fBm height profile via spectral synthesis, zero mean, unit std (same family as S5)."""
    rng = np.random.default_rng(seed)
    freqs = np.fft.rfftfreq(n)
    amp = np.zeros_like(freqs, dtype=complex)
    amp[1:] = freqs[1:] ** -(hurst + 0.5)
    amp *= np.exp(1j * rng.uniform(0, 2 * np.pi, len(freqs)))
    h = np.fft.irfft(amp, n)
    h -= h.mean()
    return h / h.std()


def surface_height(surface, x):
    """Continuous surface elevation y_s(x) in metres."""
    if surface == 'flat':
        return np.full_like(x, SURFACE_MEAN)
    if surface == 'slope':
        return SURFACE_MEAN + (x - 20.0) * np.tan(np.deg2rad(SLOPE_DEG))
    if surface == 'rough':
        raise AssertionError('rough profile is defined on the fine grid, use rough_fine()')
    raise ValueError(surface)


def rough_fine():
    dl_fine = GRIDS['d625']
    n = int(round(DOMAIN[0] / dl_fine))
    x = (np.arange(n) + 0.5) * dl_fine
    h = np.clip(fractal_profile(n) * RMS_M, -CLIP_M, CLIP_M)
    return x, SURFACE_MEAN + h


def material_database():
    return {
        'schema': 'gprMax-material-database', 'schema_version': 1,
        'database': {'id': DB_ID, 'name': 'Staircase ablation minimal database v0.1',
                     'version': '0.1.0',
                     'description': 'air + production Debye cover copied verbatim from '
                                    f'line9_materials_v0_2 (sha256 {PROD_DB_SHA256}); not site calibration'},
        'materials': {
            'material_000_air': {
                'name': 'air', 'model': 'constant',
                'base': {'relative_permittivity': 1.0, 'electric_conductivity_s_per_m': 0.0,
                         'relative_permeability': 1.0, 'magnetic_conductivity_s_per_m': 0.0},
                'metadata': {'original_id': 'air',
                             'source': f'copied verbatim from line9_materials_v0_2 sha256 {PROD_DB_SHA256}',
                             'geological_role': 'air'}},
            'material_001_cover': {
                'name': 'cover', 'model': 'debye',
                'base': {'relative_permittivity': 11.0, 'electric_conductivity_s_per_m': 0.001,
                         'relative_permeability': 1.0, 'magnetic_conductivity_s_per_m': 0.0},
                'metadata': {'original_id': 'cover',
                             'source': f'copied verbatim from line9_materials_v0_2 sha256 {PROD_DB_SHA256}',
                             'geological_role': 'silty clay',
                             'parameter_status': 'FIXED_RESEARCH_NOT_SITE_CALIBRATION'},
                'poles': [{'relative_permittivity_difference': 0.5, 'relaxation_time_s': 6.4567e-09}],
                'averagable': True}}}


def build_geometry(surface, dl, out):
    nx = int(round(DOMAIN[0] / dl))
    ny = int(round(DOMAIN[1] / dl))
    if surface == 'rough':
        xf, yf = rough_fine()
        factor = int(round(dl / GRIDS['d625']))
        y_cont = yf.reshape(nx, factor).mean(axis=1)  # block mean == sample at coarse centre
    else:
        y_cont = surface_height(surface, (np.arange(nx) + 0.5) * dl)
    k = np.round(y_cont / dl).astype(int)
    assert (k > 0).all() and (k < ny).all()
    data = np.zeros((nx, ny, 1), dtype=np.int16)
    for i in range(nx):
        data[i, :k[i], 0] = 1
    with h5py.File(out, 'w') as h:
        h.create_dataset('data', data=data, compression='gzip')
        h.create_dataset('material_keys', data=[b'material_000_air', b'material_001_cover'])
        h.attrs['MaterialDatabase'] = DB_ID
        h.attrs['dx_dy_dz'] = [dl, dl, dl]
        h.attrs['origin_xyz'] = [0., 0., 0.]
        h.attrs['shape_nxyz'] = [nx, ny, 1]
        h.attrs['batch_note'] = ('staircase ablation geometry; surface quantised by round(y/dl); '
                                 '0=air 1=cover; not a site model')
    stats = dict(surface=surface, dl_m=dl, nx=nx, ny=ny,
                 surface_min_m=float(k.min() * dl), surface_max_m=float(k.max() * dl),
                 surface_std_m=float((k * dl - SURFACE_MEAN).std()))
    return stats, (nx, ny)


def build(a):
    assert not a.out.exists()
    a.out.mkdir(parents=True)
    geo_dir = a.out / 'geometries'
    geo_dir.mkdir()
    db_path = geo_dir / f'{DB_ID}.json'
    db_path.write_text(json.dumps(material_database(), ensure_ascii=False, indent=2) + '\n',
                       encoding='utf-8')
    groups, stats = [], {}
    for surface in ['flat', 'slope', 'rough']:
        for label, dl in GRIDS.items():
            cid = f'stair_{surface}_{label}'
            geo = geo_dir / f'{surface}_{label}.h5'
            st, (nx, ny) = build_geometry(surface, dl, geo)
            stats[cid] = st
            card = a.out / 'cases' / cid / 'profile.in'
            card.parent.mkdir(parents=True)
            pml = int(round(PML_PHYS / dl))
            card.write_text(
                f'#title: {cid}; approved pure-staircase isolation ablation v0.2\n'
                f'#domain_mode: TM\n'
                f'#domain: {DOMAIN[0]} {DOMAIN[1]} inf\n'
                f'#dx_dy_dz: {dl} {dl} {dl}\n'
                f'#time_window: {TIME_WINDOW!r}\n'
                f'#omp_threads: 2\n'
                f'#pml_cells: {pml} {pml} 0 {pml} {pml} 0\n'
                f'#pml_formulation: HORIPML\n'
                f'#dispersive_averaging: y\n'
                f'#waveform: impulse 1 1 impulse\n'
                f'#hertzian_dipole: z {SRC_X} {ANT_Y} {dl / 2!r} impulse\n'
                f'#rx: {RX_X} {ANT_Y} {dl / 2!r} staircase_rx Ez\n'
                f'#geometry_objects_read: 0 0 0 ../../geometries/{surface}_{label}.h5 {DB_ID} y\n',
                encoding='utf-8')
            dt = dt_of(dl)
            g = dict(id=cid, surface=surface, grid=label,
                     input=str(card.relative_to(a.out)).replace('\\', '/'),
                     geometry=str(geo.relative_to(a.out)).replace('\\', '/'),
                     material=str(db_path.relative_to(a.out)).replace('\\', '/'),
                     # official TM mode records the invariant-axis position as 0.0
                     tx_m=[SRC_X, ANT_Y, 0.], rx_m=[RX_X, ANT_Y, 0.],
                     native_shape=[nx, ny, 1], spacing_m=[dl, dl, dl], dt_s=float(dt),
                     expected_samples=int(np.ceil(TIME_WINDOW / dt)) + 1,
                     pml_cells=pml,
                     source_type='impulse', source_frequency_Hz=1., source_amplitude_A=1.)
            for key in ['input', 'geometry', 'material']:
                g[key + '_sha256'] = hashlib.sha256((a.out / g[key]).read_bytes()).hexdigest()
            groups.append(g)
    m = dict(status='PREPARED_APPROVED_STAIRCASE_ABLATION', approval_basis=APPROVAL,
             execution_authorized=True, dt_s=float(dt_of(GRIDS['d25'])), groups=groups,
             surface_stats=stats, no_retry=True,
             design=dict(domain_m=list(DOMAIN), surface_mean_m=SURFACE_MEAN, antenna_y_m=ANT_Y,
                         src_x_m=SRC_X, rx_x_m=RX_X, time_window_s=TIME_WINDOW,
                         pml_physical_m=PML_PHYS, slope_deg=SLOPE_DEG,
                         fbm=dict(seed=SEED, hurst=HURST, rms_target_m=RMS_M, clip_m=CLIP_M,
                                  fine_grid_m=GRIDS['d625'],
                                  note='same statistical family as S5 batch; new realisation at 40 m length'),
                         quantisation_rule='k = round(y_s(x_centre)/dl); cells below are cover'),
             SFCW=dict(method='direct', frequency_start_Hz=20e6, frequency_stop_Hz=170e6,
                       frequency_step_Hz=300e3, tones=501, windows=['hann', 'blackman'],
                       zero_pad_factor=8, tail_taper_fraction=0, time_shift_s=0),
             analysis_prereg=dict(
                 primary='E_s(d) = RMS_501tones|S_{s,d}(f) - S_{s,d625}(f)| / RMS_501tones|S_{s,d625}(f)|, '
                         's in {flat,slope,rough}, d in {d25,d125}; complex direct DFT per tone',
                 convergence='ratio E(d25)/E(d125) per surface',
                 staircase_attribution='excess of E_slope and E_rough over E_flat (flat carries '
                                       'dispersion/source/PML grid effects only)',
                 secondary='band rel_delta_rms on Hann envelope, [20,120] and [150,450] ns, plus peak-time shifts',
                 anchors='S5 physical roughness effect 0.42-0.73 in-band; thresholds stay null, magnitudes only'),
             limits=('Ablation only: homogeneous cover half-space, no layering, no landslide body; '
                     '0.625 cm is a proxy smooth reference, not analytic truth; deep band [450,900] ns '
                     'not evaluated; cross-grid comparison valid, absolute waveforms not comparable to '
                     'production batches (domain/window/material truncation); not a site claim.'))
    (a.out / 'manifest.json').write_text(json.dumps(m, ensure_ascii=False, indent=2) + '\n',
                                         encoding='utf-8')
    print(json.dumps(dict(status=m['status'], groups=len(groups))))


def freeze(a):
    import line9_v401_version_controls as native
    m = json.loads((a.package / 'manifest.json').read_text('utf-8'))
    assert m['approval_basis'] == APPROVAL and m['execution_authorized']
    assert len(m['groups']) == 9
    a.action = 'prepare'
    native.prepare(a)
    p = a.out / 'execution_contract.json'
    c = json.loads(p.read_text('utf-8'))
    c['approval_basis'] = APPROVAL
    c['study_manifest'] = m
    c['limits'] = m['limits']
    c['code_identities'][str(Path(__file__).resolve())] = native.sha(__file__)
    dt_by_id = {g['id']: g['dt_s'] for g in m['groups']}
    for g in c['groups']:
        g['dt_s'] = dt_by_id[g['id']]  # per-grid dt; prepare() applies one manifest dt to all
    c['max_group_wall_s'] = 3600
    c['max_batch_wall_s'] = 14400
    native.save(p, c)
    native.save(a.out / 'preflight_verification.json', native.audit(p))
    print(json.dumps(dict(contract_sha256=native.sha(p), groups=9, version=c['expected_version'])))


def preflight(a):
    """geometry_only build of all 9 prepared inputs; validates the fresh H5 format, no solve."""
    import os
    import line9_v401_version_controls as native
    p = a.out / 'execution_contract.json'
    native.audit(p)
    import gprMax
    groups = json.loads(p.read_text('utf-8'))['groups']
    scratch = (a.out / 'geometry_preflight').resolve()
    scratch.mkdir(exist_ok=True)
    cwd = os.getcwd()
    results = []
    for g in groups:
        card = Path(g['input'])
        os.chdir(card.parent)  # ../../geometries resolves relative to the case dir
        try:
            gprMax.run(inputfile=card.name, outputfile=str(scratch / g['id']),
                       geometry_only=True, hide_progress_bars=True)
        finally:
            os.chdir(cwd)
        results.append(g['id'])
    print(json.dumps(dict(preflight='GEOMETRY_BUILD_OK', groups=results)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['build', 'freeze', 'preflight', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--package', type=Path)
    a = p.parse_args()
    a.out = a.out.resolve()
    if a.package:
        a.package = a.package.resolve()
    if a.action == 'build':
        build(a)
    else:
        import line9_v401_version_controls as native
        if a.action == 'freeze':
            freeze(a)
        elif a.action == 'preflight':
            preflight(a)
        elif a.action == 'run':
            native.supervisor.audit = native.audit
            native.supervisor.run(a.out / 'execution_contract.json')
        else:
            print(json.dumps(native.audit(a.out / 'execution_contract.json', True)))


if __name__ == '__main__':
    main()
