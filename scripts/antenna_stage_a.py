"""Stage A bowtie antenna characterization batch (approved 2026-10-12).

Approved scope (user, 2026-10-12, before sleep): realistic antenna model for the
folded-bowtie UAV-GPR device (equal-length centre-fed bowtie approximation:
total length 3.35 m, 60 deg flare, 5 cm-class feed, no resistive loading,
labelled "unfocused original-antenna approximation" = conservative lower bound),
three solves: A1 free space (arms along y), A1b free space (arms along x),
A2 flat cover half-space with Tx+passive Rx bowties at 1.3 m baseline.
Outputs: official NTFF far-field patterns at 20/95/170 MHz (A1/A1b),
feed-port voltage/current (ringing), near-field probes, Rx gap field (A2).

Deviations from the pre-approval sketch, forced by the 16 GB VRAM budget
(user pre-authorized domain resizing so ROG can run it): PML 60 cells
(1.5 m) instead of 80; domains 7.0^3 m (A1/A1b) and 6.0 x 7.75 x 9.75 m (A2).
Pattern extraction uses the official V4 NTFF commands instead of mid-field
probe rings as primary product; probe rings/lines are kept as cross-check.
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import sys

import h5py
import numpy as np

import gprMax
import line9_v401_version_controls as native
import hs4_station_grid_controls as supervisor
from hs_capsule_identity import sha256 as sha

ROOT = Path(__file__).resolve().parents[1]
DB_SOURCE = ROOT/'model_current/geometries/line9_materials_v0_2.json'
DB_ID = 'line9_materials_v0_2'

APPROVAL = ('2026-10-12 用户批准「做真实天线吧」及阶段A方案（自由空间领结两取向+平地覆盖层双天线），'
            '并明确允许为适配 ROG 算力调整域大小等配套。领结为设备方幻灯片的等长中馈近似：'
            '全长3.35m、60°张角、无电阻加载，标注为未聚焦原天线近似（保守下限）；'
            '设备为偶极（折叠领结/等效磁回路，2L+d≈λ_low/4，λ_low=20MHz），无更详细结构资料。')

DL = 0.025
C0 = 299792458.0
DT = DL/(C0*math.sqrt(3.0))
WINDOW_S = 6e-7
STEPS = int(math.ceil(WINDOW_S/DT))+1
PML = 60
ARM = 1.65          # each arm length beyond the 2*DL feed gap -> total 3.35 m
HALFW = 0.95        # end half-width -> 2*atan(0.95/1.65) = 59.9 deg flare
STRIP = 0.15        # feed-strip length per arm (guarantees feed-edge connectivity)
TONES_HZ = [20e6, 95e6, 170e6]
FEED_R_OHM = 50.0

A1_SHAPE = (280, 280, 280)              # 7.0^3 m
A2_SHAPE = (240, 310, 390)              # 6.0 x 7.75 x 9.75 m
A2_COVER_TOP_M = 2.5                    # surface z; bottom PML 1.5 m + 1.0 m inner cover
A2_FEED_Z = 7.5                         # 5.0 m AGL
A2_TX_X, A2_RX_X, A2_CY = 2.35, 3.65, 3.875   # 1.3 m baseline along x, arms along y


def fmt(values):
    return ' '.join(f'{v:.12g}' for v in values)


def make_geometry(path, shape, cover_top_m=None):
    nx, ny, nz = shape
    data = np.zeros(shape, dtype=np.int16)
    if cover_top_m is not None:
        top = int(round(cover_top_m/DL))
        assert abs(top*DL-cover_top_m) < 1e-12
        data[:, :, :top] = 1
    with h5py.File(path, 'x') as h:
        h.create_dataset('data', data=data, chunks=(min(nx, 8), min(ny, 64), min(nz, 64)),
                         compression='gzip', compression_opts=1)
        h.create_dataset('material_keys', data=np.asarray(['material_000_air', 'material_001_cover'], dtype='S'))
        h.attrs['dx_dy_dz'] = [DL]*3
        h.attrs['origin_xyz'] = [0., 0., 0.]
        h.attrs['shape_nxyz'] = list(shape)
        h.attrs['MaterialDatabase'] = DB_ID
        h.attrs['MaterialDatabaseSchemaVersion'] = 1
        h.attrs['ModelContent'] = 'antenna_stage_a; air%s' % ('' if cover_top_m is None else f'; flat cover z<{cover_top_m} m')
    with h5py.File(path, 'r') as h:
        assert np.array_equal(h['data'][:], data)
    return dict(sha256=sha(path), cover_top_m=cover_top_m)


def bowtie_lines(feed, axis, source_waveform, passive_id=None, window_s=WINDOW_S):
    """Two triangular PEC arms along `axis` ('x' or 'y') in the horizontal plane z=feed[2].

    Feed connectivity: a point-apex #triangle rasterizes to zero cells within ~2 cells of
    the apex (confirmed 2026-10-12 with feedcell.vtkhdf), which left the feed edge floating
    in the r1 batches. Each arm therefore gets a rectangular feed strip (width 2*DL, the
    approved 5 cm-class feed width; STRIP long, well past the sliver gap) that overlaps the
    triangle and touches the voltage-source edge endpoints.
    """
    fx, fy, fz = feed
    lines = []
    if axis == 'y':
        apex_p, apex_m = (fx, fy+DL, fz), (fx, fy-DL, fz)
        c1p, c2p = (fx-HALFW, fy+DL+ARM, fz), (fx+HALFW, fy+DL+ARM, fz)
        c1m, c2m = (fx-HALFW, fy-DL-ARM, fz), (fx+HALFW, fy-DL-ARM, fz)
        # source y-edge spans [fy, fy+DL]; arm+ must touch node fy+DL, arm- node fy.
        # Triangles rasterize into z-layer [fz, fz+DL]; strips must use the same layer.
        strips = [(fx-DL, fy+DL, fz, fx+DL, fy+DL+STRIP, fz+DL),
                  (fx-DL, fy-STRIP, fz, fx+DL, fy, fz+DL)]
    else:
        apex_p, apex_m = (fx+DL, fy, fz), (fx-DL, fy, fz)
        c1p, c2p = (fx+DL+ARM, fy-HALFW, fz), (fx+DL+ARM, fy+HALFW, fz)
        c1m, c2m = (fx-DL-ARM, fy-HALFW, fz), (fx-DL-ARM, fy+HALFW, fz)
        strips = [(fx+DL, fy-DL, fz, fx+DL+STRIP, fy+DL, fz+DL),
                  (fx-STRIP, fy-DL, fz, fx, fy+DL, fz+DL)]
    for apex, c1, c2 in [(apex_p, c1p, c2p), (apex_m, c1m, c2m)]:
        lines.append(f'#triangle: {fmt(apex)} {fmt(c1)} {fmt(c2)} {DL:.12g} pec')
    for x1, y1, z1, x2, y2, z2 in strips:
        lines.append(f'#box: {fmt((x1, y1, z1))} {fmt((x2, y2, z2))} pec')
    if passive_id is None:
        lines.append(f'#voltage_source: {axis} {fmt(feed)} {FEED_R_OHM:.12g} {source_waveform} 0 {window_s:.12g} feed nyquist')
    else:
        comp = 'E'+axis
        lines.append(f'#rx: {fmt(feed)} {passive_id} {comp}')
    return lines


def ring_rx_lines(center, plane, radius, step_deg=15, prefix='pr'):
    cx, cy, cz = center
    lines = []
    for i, deg in enumerate(range(0, 360, step_deg)):
        a = math.radians(deg)
        u, v = radius*math.cos(a), radius*math.sin(a)
        if plane == 'yz':
            p = (cx, cy+u, cz+v)
        else:
            p = (cx+u, cy, cz+v)
        lines.append(f'#rx: {fmt(p)} {prefix}{plane}{i:02d} Ex Ey Ez')
    return lines


def ntff_lines(feed, arm_span, axis):
    fx, fy, fz = feed
    lo = [fx-0.4, fy-0.4, fz-0.4]
    hi = [fx+0.4, fy+0.4, fz+0.4]
    i = 1 if axis == 'y' else 0
    lo[i] = round((feed[i]-arm_span-0.1)/DL)*DL
    hi[i] = round((feed[i]+arm_span+0.1)/DL)*DL
    tones = ' '.join(f'{f:.12g}' for f in TONES_HZ)
    return [f'#ntff_surface: {fmt(lo)} {fmt(hi)} stage_a_surface',
            f'#ntff_frequency: stage_a_surface stage_a_pattern {tones} rectangular',
            '#ntff_antenna_ports: stage_a_pattern feed',
            '#ntff_far_field_array: 0 180 5 0 355 5 stage_a_pattern full_sphere Etheta Ephi gain_dbi directivity_dbi']


def case_text(case_id, shape, feed, axis, geometry_file, extra_lines, ntff=False):
    nx, ny, nz = shape
    lines = [f'#title: {case_id}; Stage A bowtie antenna characterization; unfocused approximation',
             f'#domain: {fmt([nx*DL, ny*DL, nz*DL])}',
             f'#dx_dy_dz: {fmt([DL]*3)}',
             f'#time_window: {WINDOW_S:.12g}',
             '#omp_threads: 2',
             f'#pml_cells: {PML} {PML} {PML} {PML} {PML} {PML}',
             '#pml_formulation: HORIPML',
             '#dispersive_averaging: y',
             '#waveform: impulse 1 1 impulse',
             # geometry read must precede every in-file geometry command: a later
             # #geometry_objects_read overwrites the PEC antenna cells (verified 2026-10-12,
             # owlast/owfirst rasterization test: 0 vs 62 pec cells).
             f'#geometry_objects_read: 0 0 0 ../../geometries/{geometry_file} {DB_ID} n',
             *bowtie_lines(feed, axis, 'impulse'),
             *extra_lines, '']
    if ntff:
        lines[-1:-1] = ntff_lines(feed, ARM+DL, axis)
    return '\n'.join(lines)


def grid_check(feed, shape, label):
    for i, v in enumerate(feed):
        assert abs(v/DL-round(v/DL)) < 1e-9, f'{label} feed off grid'
        assert PML*DL < v < shape[i]*DL-PML*DL, f'{label} feed inside PML'


def prepare(a):
    assert not a.out.exists()
    (a.out/'cases').mkdir(parents=True)
    (a.out/'geometries').mkdir()
    db = a.out/'geometries'/DB_SOURCE.name
    shutil.copyfile(DB_SOURCE, db)

    geos = {
        'A1': make_geometry(a.out/'geometries/A1.h5', A1_SHAPE),
        'A2': make_geometry(a.out/'geometries/A2.h5', A2_SHAPE, A2_COVER_TOP_M),
    }
    specs = []
    feed1 = (3.5, 3.5, 3.5)
    grid_check(feed1, A1_SHAPE, 'A1')
    specs.append(dict(id='A1_free_y', shape=A1_SHAPE, geo='A1.h5', feed=feed1, axis='y', ntff=True,
                      extra=ring_rx_lines(feed1, 'yz', 1.5)+ring_rx_lines(feed1, 'xz', 1.5)))
    specs.append(dict(id='A1b_free_x', shape=A1_SHAPE, geo='A1.h5', feed=feed1, axis='x', ntff=True,
                      extra=ring_rx_lines(feed1, 'yz', 1.5)+ring_rx_lines(feed1, 'xz', 1.5)))
    tx = (A2_TX_X, A2_CY, A2_FEED_Z)
    rx = (A2_RX_X, A2_CY, A2_FEED_Z)
    grid_check(tx, A2_SHAPE, 'A2tx')
    grid_check(rx, A2_SHAPE, 'A2rx')
    probes = [f'#rx: {fmt((A2_TX_X, A2_CY, z))} dz{i:02d} Ex Ey Ez'
              for i, z in enumerate(np.arange(2.75, 7.26, 0.25))]
    probes += [f'#rx: {fmt((x, A2_CY, A2_FEED_Z))} hx{i:02d} Ex Ey Ez'
               for i, x in enumerate(np.arange(1.6, 4.41, 0.2))]
    specs.append(dict(id='A2_ground_pair', shape=A2_SHAPE, geo='A2.h5', feed=tx, axis='y', ntff=False,
                      extra=[*bowtie_lines(rx, 'y', None, passive_id='rxgap'), *probes], rx_gap=rx))

    groups = []
    for s in specs:
        text = case_text(s['id'], s['shape'], s['feed'], s['axis'], s['geo'], s['extra'], s['ntff'])
        case = a.out/'cases'/s['id']
        case.mkdir()
        p = case/'profile.in'
        p.write_text(text, encoding='utf-8')
        g = dict(id=s['id'], input=f'cases/{s["id"]}/profile.in', geometry=f'geometries/{s["geo"]}',
                 material='geometries/'+DB_SOURCE.name,
                 input_sha256=sha(p), geometry_sha256=geos[s['geo'][:2]]['sha256'], material_sha256=sha(db),
                 native_shape=list(s['shape']), spacing_m=[DL]*3, expected_samples=STEPS,
                 tx_m=list(s['feed']), rx_m=[list(s['feed'])], antenna_axis=s['axis'],
                 ntff=bool(s['ntff']), tones_hz=TONES_HZ, source_type='impulse',
                 source_frequency_Hz=1., source_amplitude_A=1.)
        if 'rx_gap' in s:
            g['rx_gap_m'] = list(s['rx_gap'])
        groups.append(g)

    manifest = dict(status='PREPARED_APPROVED_ANTENNA_STAGE_A', approval_basis=APPROVAL,
                    execution_authorized=True, groups=groups, no_retry=True, dt_s=DT,
                    tones_hz=TONES_HZ, antenna_model='unfocused_centre_fed_bowtie_approximation',
                    bowtie=dict(total_length_m=2*(ARM+DL), flare_deg=2*math.degrees(math.atan(HALFW/ARM)),
                                gap_m=2*DL, thickness_m=DL, feed_resistance_ohm=FEED_R_OHM,
                                loading='none; conservative lower bound vs device folded/focused antenna'),
                    generator_sha256=sha(Path(__file__)),
                    limits='Stage A antenna characterization only; not a survey-line simulation. '
                           'Bowtie is an equal-length centre-fed approximation of the folded device antenna '
                           '(no balun/folding/focusing detail available). PML 60 cells (1.5 m). '
                           'A2 flat cover uses the line9_research_materials cover; flat ground only.')
    native.save(a.out/'manifest.json', manifest)
    print(json.dumps(dict(groups=len(groups), expected_samples=STEPS, dt_s=DT)))


def audit(path, completed=False):
    c = json.loads(Path(path).read_text('utf-8'))
    assert gprMax.__version__ == c['expected_version']
    assert Path(sys.executable).resolve() == Path(c['python']).resolve()
    native_dir = Path(gprMax.__file__).parent
    for category in ['code_identities', 'file_identities', 'binary_identities']:
        for p, digest in c[category].items():
            assert sha(p) == digest, (category, p)
    for p, digest in c['source_identities'].items():
        assert sha(native_dir/p) == digest, p
    rows = []
    for g in c['groups']:
        row = dict(id=g['id'], input_sha256=sha(g['input']))
        assert row['input_sha256'] == g['input_sha256'], g['id']
        if completed:
            raw = Path(g['input']).with_suffix('.h5')
            with h5py.File(raw) as h:
                assert str(h.attrs['gprMax']) == c['expected_version']
                np.testing.assert_array_equal(h.attrs['nx_ny_nz'], g['native_shape'])
                np.testing.assert_array_equal(h.attrs['dx_dy_dz'], g['spacing_m'])
                assert abs(float(h.attrs['dt'])/g['dt_s']-1) < 1e-13
                assert h.attrs['Iterations'] == g['expected_samples']
                rxs = h['rxs']
                assert len(rxs) >= 1
                for name in rxs:
                    for comp in rxs[name]:
                        v = rxs[name][comp][:]
                        assert v.dtype == np.float64 and v.shape == (g['expected_samples'],) and np.isfinite(v).all()
                if 'srcs' in h:
                    exc = h['srcs/src1/excitation']
                    assert exc.attrs['WaveformType'] == 'impulse'
                    samples = exc['samples'][:]
                    np.testing.assert_array_equal(np.flatnonzero(samples), [0])
                rx_count = len(rxs)
            extras = sorted(p.name for p in raw.parent.iterdir() if p.name != raw.name and p.suffix not in ('.in', '.log', '.jsonl'))
            row.update(native_sha256=sha(raw), rx_count=rx_count, extra_outputs=extras)
        rows.append(row)
    return dict(status='PASS_NATIVE_IDENTITY_NOT_PHYSICAL_CERTIFICATION', completed=completed,
                contract_sha256=sha(path), version=c['expected_version'], groups=rows)


def freeze(a):
    manifest = json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert manifest['approval_basis'] == APPROVAL and manifest['execution_authorized']
    assert len(manifest['groups']) == 3
    a.action = 'prepare'
    native.prepare(a)
    p = a.out/'execution_contract.json'
    c = json.loads(p.read_text('utf-8'))
    c.update(approval_basis=APPROVAL, study_manifest=manifest, limits=manifest['limits'],
             max_batch_wall_s=10800, max_group_wall_s=2400)
    c['code_identities'][str(Path(__file__).resolve())] = sha(Path(__file__))
    native.save(p, c)
    native.save(a.out/'preflight_verification.json', audit(p))
    print(json.dumps(dict(contract_sha256=sha(p), groups=3, max_batch_wall_s=10800)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=['prepare', 'freeze', 'run', 'verify'])
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--package', type=Path)
    a = p.parse_args()
    a.out = a.out.resolve()
    if a.action == 'prepare':
        prepare(a)
    elif a.action == 'freeze':
        a.package = a.package.resolve()
        freeze(a)
    elif a.action == 'run':
        supervisor.audit = audit
        supervisor.run(a.out/'execution_contract.json')
    else:
        print(json.dumps(audit(a.out/'execution_contract.json', True)))
