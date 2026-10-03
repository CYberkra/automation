"""CPU-only HS4 3D morphology audit; diagnostics are not physical truth.

Checks 38 raw traces, actual coordinates, full input geometry and source
identity, then compares several declared picks/processing variants against
local-depth approximations. Never invokes a solver or accesses C5/C8.
"""
import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import brentq
from scipy.signal import find_peaks

from hs_capsule_identity import read_manifest, verify_file, sha256
from sfcw_official_loader_v0_2 import (
    verify_official_runtime, trace_time_response, FREQ, OFFICIAL_PROCESSING_SHA256)
from gprMax.toolboxes.SFCW.processing import (
    load_source, load_receiver, direct_frequency_response, reconstruct_time_response)

ROOT = Path(__file__).resolve().parents[1]
HS = ROOT / 'artifacts/research_checks/2026-10-02_halfspace_standard_hs'
WINDOWS = [(160., 220.), (170., 215.), (180., 210.)]
LAGS = np.linspace(-20., 20., 401)  # 0.1 ns interpolation, not physical resolution
C0 = 299792458.
EPS0 = 8.8541878128e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def correlation(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if np.std(a) == 0 or np.std(b) == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def cover_index(f):
    omega = 2 * np.pi * np.asarray(f)
    return np.sqrt(18.017 + 7.878 / (1 + 1j * omega * 6.4567e-9)
                   - 1j * .003 / (EPS0 * omega))


def local_depth_times(depths, offset, anchor):
    """Single-frequency planar-layer ray approximation, NOT rough-interface truth."""
    frequency = 95e6
    n = float(cover_index(frequency).real)
    dn = float((cover_index(frequency + 1e3).real - cover_index(frequency - 1e3).real) / 2e3)
    group_n = n + frequency * dn
    def time(depth):
        p = brentq(lambda p: 2 * (15 * p / np.sqrt(1 - p*p)
                                  + depth * p / np.sqrt(n*n - p*p)) - offset,
                   0., 1. - 1e-12)
        return 2 * (15 / np.sqrt(1 - p*p)
                    + depth * group_n / np.sqrt(1 - (p/n)**2)) / C0 * 1e9
    baseline = time(3.)
    return np.array([anchor + time(float(d)) - baseline for d in depths]), {
        'role': 'local planar-ray approximation only; not a validation reference',
        'frequency_Hz': frequency, 'cover_phase_index': n, 'cover_group_index': group_n,
        'air_height_m': 15., 'offset_m': offset, 'flat_anchor_ns': anchor,
        'flat_3m_group_ray_ns': baseline,
        'limitations': 'single frequency, local planar depth; no rough-interface lateral scattering, full footprint, or pulse distortion'}


def full_field_ray_proxy(source, receiver, table, approximation):
    """Patch-centre optical-path minimum, without scattering/amplitude model.

    Select by phase optical path at 95 MHz; report corresponding group time.
    This proxy is NOT arrival truth; vertical faces and multiple paths are
    omitted, and the earliest path need not dominate the measured envelope.
    """
    n, ng = approximation['cover_phase_index'], approximation['cover_group_index']
    coords = (np.arange(48)+.5)*.25
    xx, yy = np.meshgrid(coords, coords, indexing='ij')
    points = np.stack((xx.ravel(), yy.ravel()), axis=1)
    depth = (12-table).ravel()[None,:]
    def half_path(position):
        lateral = np.linalg.norm(position[:,None,:2]-points[None,:,:], axis=2)
        lo, hi = np.zeros_like(lateral), np.full_like(lateral,1.-1e-12)
        for _ in range(60):
            p=(lo+hi)/2
            distance=15*p/np.sqrt(1-p*p)+depth*p/np.sqrt(n*n-p*p)
            lo=np.where(distance<lateral,p,lo); hi=np.where(distance>=lateral,p,hi)
        p=(lo+hi)/2
        air=15/np.sqrt(1-p*p); cover=depth/np.sqrt(1-(p/n)**2)
        return (air+cover*n)/C0*1e9, (air+cover*ng)/C0*1e9
    phase_a,group_a=half_path(source); phase_b,group_b=half_path(receiver)
    optical=phase_a+phase_b; group=group_a+group_b
    ix=optical.argmin(axis=1)
    time=group[np.arange(len(source)),ix]-approximation['flat_3m_group_ray_ns']+approximation['flat_anchor_ns']
    point=points[ix]; midpoint=(source+receiver)/2
    return {'role':'full-field patch optical-path minimum proxy only; NOT arrival truth',
            'group_time_proxy_ns':time.tolist(),'selected_patch_xy_m':point.tolist(),
            'selected_patch_cover_depth_m':depth.ravel()[ix].tolist(),
            'lateral_distance_from_midpoint_m':np.linalg.norm(point-midpoint[:,:2],axis=1).tolist(),
            'limitations':'patch-centre discretization; no reflection amplitudes, full-wave interference, vertical step faces or validation of source/PML; earliest optical path need not dominate the envelope'}


def remove_rank(matrix, rank):
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    u, singular, vh = np.linalg.svd(centered, full_matrices=False)
    return centered - (u[:, :rank] * singular[:rank]) @ vh[:rank], singular


def lag_fit(t, matrix, ref_t, template, window):
    mask = (t >= window[0]) & (t <= window[1])
    target = matrix[mask] - matrix[mask].mean(axis=0, keepdims=True)
    require(np.all(np.linalg.norm(target, axis=0) > 0), 'zero-energy correlation input')
    model = np.array([np.interp(t[mask] - lag, ref_t, template) for lag in LAGS])
    model -= model.mean(axis=1, keepdims=True)
    score = model @ target / (np.linalg.norm(model, axis=1)[:, None]
                              * np.linalg.norm(target, axis=0)[None, :])
    ix = np.argmax(score, axis=0)
    alternatives=[]
    for j,best in enumerate(ix):
        candidates=find_peaks(score[:,j])[0]
        candidates=candidates[candidates!=best]
        if len(candidates):
            second=int(candidates[np.argmax(score[candidates,j])])
            alternatives.append({'lag_ns':float(LAGS[second]),'correlation':float(score[second,j]),
                                 'best_minus_alternative_correlation':float(score[best,j]-score[second,j])})
        else:
            alternatives.append(None)
    return {'lag_ns': LAGS[ix].tolist(), 'signed_normalized_correlation': score[ix, np.arange(matrix.shape[1])].tolist(),
            'best_alternative_local_maximum':alternatives,
            'search_boundary_count': int(np.count_nonzero((ix == 0) | (ix == len(LAGS)-1))),
            'search_ns': [-20., 20.], 'interpolation_step_ns': .1,
            'convention': 'positive lag = target later; real signed correlation; de-mean each window; linear interpolation; no periodic roll'}


def sanity_checks():
    t = np.arange(0., 300., .2)
    def pulse(t):
        return np.exp(-((t-187)/6)**2) * np.cos(2*np.pi*.095*(t-187))
    target = pulse(t-4.3)[:, None]
    fit = lag_fit(t, target, t, pulse(t), (160., 220.))
    require(abs(fit['lag_ns'][0]-4.3) <= .10001, 'known delay recovery failed')
    _, params = local_depth_times([2.8, 3., 3.2], 1.3, 187.)
    times, _ = local_depth_times([2.8, 3., 3.2], 1.3, 187.)
    require(np.all(np.diff(times) > 0) and abs(times[1]-187.) < 1e-10, 'local-depth sign/anchor failed')
    require(params['cover_group_index'] > 0, 'invalid group index')
    source=np.array([[6.,5.35,27.]])
    receiver=np.array([[6.,6.65,27.]])
    proxy=full_field_ray_proxy(source,receiver,np.full((48,48),9.),params)
    require(abs(proxy['group_time_proxy_ns'][0]-187.) < .01,'flat full-field proxy consistency failed')
    return {'known_analytic_pulse_delay_4p3ns': fit, 'local_depth_monotonic_anchor_check': True,
            'flat_full_field_proxy_anchor_check':True,
            'evidence_level': 'array/mathematical sanity only'}


def geometry_lines(path):
    return [line for line in path.read_text(encoding='utf-8').splitlines() if line.startswith('#box:')]


def numerics_lines(path):
    return [line for line in path.read_text(encoding='utf-8').splitlines()
            if line.startswith('#') and not line.startswith(('#title:', '#box:', '#hertzian_dipole:', '#rx:', '#geometry_view:'))]


def load_group(prefix, count, saved_name, table):
    raw, sources, receivers, input_positions = [], [], [], []
    first_src = load_source(HS / f'{prefix}_t01.h5')
    first_rx = load_receiver(HS / f'{prefix}_t01.h5', receiver_path='/rxs/rx1', component='Ex')
    boxes = geometry_lines(HS / 'hs4_rough_halfspace.in')
    numerics = numerics_lines(HS / 'hs4_rough_halfspace.in')
    require(len(boxes) == 2305, 'unexpected full 3D geometry')
    # Check all cover patches independently of the generator.
    parsed_table = np.full((48, 48), np.nan)
    for line in boxes[1:]:
        fields = line.split()[1:]
        x0, y0, z0, x1, y1, z1 = map(float, fields[:6])
        bi, bj = int(round(x0/.25)), int(round(y0/.25))
        require(fields[6] == 'cover' and abs(x1-x0-.25) < 1e-12
                and abs(y1-y0-.25) < 1e-12 and z1 == 12, 'unexpected cover patch')
        require(np.isnan(parsed_table[bi,bj]), 'duplicate cover patch')
        parsed_table[bi,bj] = z0
    require(np.array_equal(parsed_table, table), 'input cover boxes differ from saved table')
    for k in range(1, count+1):
        path = HS / f'{prefix}_t{k:02d}.h5'
        src, rx = load_source(path), load_receiver(path, receiver_path='/rxs/rx1', component='Ex')
        require(src.dt == first_src.dt and src.time_offset == first_src.time_offset
                and np.array_equal(src.samples, first_src.samples), 'source identity differs')
        require(rx.dt == first_rx.dt and rx.time_offset == first_rx.time_offset
                and rx.samples.shape == first_rx.samples.shape, 'receiver axes differ')
        require(geometry_lines(path.with_suffix('.in')) == boxes, 'geometry differs between traces')
        require(numerics_lines(path.with_suffix('.in')) == numerics, 'numerics/material commands differ')
        text = path.with_suffix('.in').read_text(encoding='utf-8').splitlines()
        tx = next(v for v in text if v.startswith('#hertzian_dipole:')).split()
        r = next(v for v in text if v.startswith('#rx:')).split()
        require(tx[1] == 'x' and r[-1] == 'Ex', 'polarization/component differs')
        input_positions.append({'source': list(map(float, tx[2:5])), 'receiver': list(map(float, r[1:4]))})
        with h5py.File(path, 'r') as h:
            require(h['rxs/rx1/Ex'].dtype == np.float64, 'raw receiver not FP64')
            require(np.array_equal(h.attrs['nx_ny_nz'],[240,240,660])
                    and np.array_equal(h.attrs['dx_dy_dz'],[.05,.05,.05]), 'H5 grid differs from 3D input')
            require(h.attrs['gprMax']=='4.0.0', 'H5 solver version differs')
            sources.append(np.asarray(h['srcs/src1'].attrs['Position']).tolist())
            receivers.append(np.asarray(h['rxs/rx1'].attrs['Position']).tolist())
        raw.append(rx.samples)
    raw = np.stack(raw, axis=1)
    n = min(len(raw), int(np.floor(1200e-9 / first_rx.dt)) + 1)
    response = direct_frequency_response(first_src, replace(first_rx, samples=raw[:n]), FREQ,
        tail_taper_fraction=(round(200e-9/first_rx.dt)-.25)/n)
    reconstructed = reconstruct_time_response(response, zero_pad_factor=8, window='hann')
    t, signed, env_c = reconstructed.time*1e9, reconstructed.real_bandpass, reconstructed.complex_envelope
    with np.load(HS/saved_name) as saved:
        require(np.array_equal(t, saved['t']) and np.array_equal(signed, saved['S']), 'saved reconstruction differs')
    return t, signed, env_c, np.asarray(sources), np.asarray(receivers), input_positions


def audit_group(prefix, count, saved, table, ref_t, template, anchor, ref_positions):
    t, signed, env_c, source, receiver, nominal = load_group(prefix, count, saved, table)
    midpoint = (source + receiver)/2
    xb = np.floor(midpoint[:,0]/.25).astype(int)
    geometry = {}
    for role, points in [('Tx',source),('midpoint',midpoint),('Rx',receiver)]:
        z = table[xb, np.floor(points[:,1]/.25).astype(int)]
        geometry[role] = {'z0_m': z.tolist(), 'cover_depth_m': (12-z).tolist(),
                          'depth_range_m': [float((12-z).min()), float((12-z).max())]}
    expected, approx = local_depth_times(geometry['midpoint']['cover_depth_m'], 1.3, anchor)
    full_field=full_field_ray_proxy(source,receiver,table,approx)
    whole = (t >= 0) & (t <= 250)
    tm, sm = t[whole], signed[whole]
    r1, singular = remove_rank(sm, 1)
    r2, _ = remove_rank(sm, 2)
    crop = (t >= 140) & (t <= 240)
    crop_r1, _ = remove_rank(signed[crop], 1)
    centered = sm-sm.mean(axis=0,keepdims=True)
    demean_space = centered - centered.mean(axis=1, keepdims=True)
    variants = [('raw',tm,sm), ('spatial_mean_removed',tm,demean_space),
                ('rank1_full0_250',tm,r1), ('rank2_full0_250',tm,r2),
                ('rank1_crop140_240',t[crop],crop_r1)]
    picks = {}
    for name, axis, matrix in variants:
        picks[name] = {}
        for lo,hi in WINDOWS:
            mask = (axis >= lo) & (axis <= hi)
            ix = np.argmax(np.abs(matrix[mask]), axis=0)
            ridge = axis[mask][ix]
            fit = lag_fit(axis, matrix, ref_t, template, (lo,hi))
            picks[name][f'{lo:g}_{hi:g}'] = {
                'signed_abs_peak_ns': ridge.tolist(), 'peak_window_edge_count': int(np.count_nonzero((ix == 0)|(ix == mask.sum()-1))),
                'peak_vs_local_depth_correlation': correlation(ridge,expected), 'peak_span_ns': float(np.ptp(ridge)),
                'template_fit': fit, 'fit_vs_local_depth_correlation': correlation(fit['lag_ns'],expected)}
            picks[name][f'{lo:g}_{hi:g}']['fit_vs_full_field_proxy_correlation']=correlation(fit['lag_ns'],full_field['group_time_proxy_ns'])
    official = 2*np.abs(env_c)
    mask = (t >= 160) & (t <= 220)
    idx = np.argmax(official[mask], axis=0)
    env_ridge = t[mask][idx]
    ground = (t >= 85) & (t <= 120)
    ground_ridge = t[ground][np.argmax(official[ground], axis=0)]
    anchor_receiver = np.asarray(ref_positions['receiver'])
    anchor_source = np.asarray(ref_positions['source'])
    matched = np.all(source == anchor_source,axis=1) & np.all(receiver == anchor_receiver,axis=1)
    results = {'traces': count, 'component':'Ex', 'raw_dtype':'float64', 'saved_rebuild_max_abs_error':0.,
        'H5_solver_version':'4.0.0','H5_grid_cells':[240,240,660],'H5_cell_size_m':[.05,.05,.05],
        'all_input_geometries_match_table':True, 'source_signals_identical':True,
        'all_numerics_and_material_commands_match':True,
        'actual_sources':source.tolist(), 'actual_receivers':receiver.tolist(), 'nominal_input_positions':nominal,
        'profile_axis':'actual Tx/Rx midpoint y (m); transverse x from raw H5',
        'geometry':geometry, 'local_depth_approx_ns':expected.tolist(), 'local_depth_approx':approx,
        'full_field_ray_proxy':full_field,
        'matched_HS1_HS2_reference_trace_numbers': (np.flatnonzero(matched)+1).tolist(),
        'matched_reference_note':'HS1/HS2 are single central traces; only exact-coordinate matches have matched acquisition. Template fitting elsewhere is diagnostic, not clean truth.',
        'rank1_removed_total_centered_energy_fraction':float(1-np.sum(r1*r1)/np.sum(centered*centered)),
        'rank1_residual_total_energy_fraction':float(np.sum(r1*r1)/np.sum(centered*centered)),
        'rank1_vs_spatial_mean_residual_relative_norm':float(np.linalg.norm(r1-demean_space)/np.linalg.norm(r1)),
        'first_three_singular_values':singular[:3].tolist(),
        'official_interface_window_envelope_peak_ns':env_ridge.tolist(),
        'official_envelope_window_edge_count':int(np.count_nonzero((idx == 0)|(idx == mask.sum()-1))),
        'official_envelope_peak_vs_local_depth_correlation':correlation(env_ridge,expected),
        'official_ground_window_peak_ns':ground_ridge.tolist(), 'picks':picks,
        'raw_bandlimited_edge_center_difference_over_direct':{},
        'nearest_x_PML_inner_edge_distance_m':float(np.min(np.minimum(source[:,0]-1,11-source[:,0]))),
        'full_x_domain_m':12., 'formal_morphology_validation':'NOT_ESTABLISHED'}
    direct_peak = float(np.max(official[(t>=0)&(t<=10)]))
    difference = signed[:,[0,-1]]-signed[:,[count//2]]
    for name,(lo,hi) in {'direct':(0,20),'pre_ground':(20,85),'ground':(85,120),'interface':(160,220)}.items():
        m=(t>=lo)&(t<=hi)
        results['raw_bandlimited_edge_center_difference_over_direct'][name]=float(np.max(np.abs(difference[m]))/direct_peak)
    return results, {'t':t,'signed':signed,'complex_envelope':env_c,'tm':tm,'r1':r1,'r2':r2,
                     'crop_t':t[crop],'crop_r1':crop_r1,'midpoint_y':midpoint[:,1],'expected':expected}


def plot(out, name, r, arrays):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    x=arrays['midpoint_y']
    for role in ('Tx','midpoint','Rx'):
        ax[0,0].plot(x,r['geometry'][role]['cover_depth_m'],'.-',label=role+' sample')
    ax[0,0].plot(x,r['full_field_ray_proxy']['selected_patch_cover_depth_m'],'k:',label='full-field proxy patch')
    ax[0,0].set_ylabel('Cover depth below surface (m)'); ax[0,0].invert_yaxis()
    ax[0,0].set_title('Actual model depth; ray approximation is not truth'); ax[0,0].legend()
    rawmask=(arrays['t']>=160)&(arrays['t']<=220)
    env=2*np.abs(arrays['complex_envelope'][rawmask])
    im=ax[0,1].pcolormesh(x,arrays['t'][rawmask],env,cmap='viridis',shading='nearest')
    fig.colorbar(im,ax=ax[0,1],label='Official 2|complex envelope|')
    ax[0,1].plot(x,r['official_interface_window_envelope_peak_ns'],'w.-',label='window max')
    ax[0,1].plot(x,arrays['expected'],'r--',label='local-depth approximation')
    ax[0,1].plot(x,r['full_field_ray_proxy']['group_time_proxy_ns'],'k:',label='full-field ray proxy')
    ax[0,1].set_ylim(220,160); ax[0,1].set_ylabel('Time (ns)'); ax[0,1].legend()
    ax[0,1].set_title('Raw official envelope; no background/SVD')
    ref=r['local_depth_approx']['flat_anchor_ns']
    ax[0,2].plot(x,arrays['expected']-ref,'k--',label='local-depth approximation')
    ax[0,2].plot(x,np.asarray(r['full_field_ray_proxy']['group_time_proxy_ns'])-ref,'k:',label='full-field ray proxy')
    for variant in ['raw','spatial_mean_removed','rank1_full0_250','rank2_full0_250','rank1_crop140_240']:
        p=r['picks'][variant]['160_220']['template_fit']
        ax[0,2].plot(x,p['lag_ns'],'.-',label=variant)
    ax[0,2].set_ylabel('Template lag (ns); positive=later'); ax[0,2].set_title('Diagnostic picks depend on processing'); ax[0,2].legend(fontsize=7)
    bound=float(np.quantile(np.abs(arrays['r1']),.999))
    for col,(axis,residual,title) in enumerate([(arrays['tm'],arrays['r1'],'Rank 1 removed (0-250 ns fit)'),
                                                (arrays['tm'],arrays['r2'],'Rank 2 removed (0-250 ns fit)'),
                                                (arrays['crop_t'],arrays['crop_r1'],'Rank 1 removed (140-240 ns fit)')]):
        m=(axis>=160)&(axis<=220)
        im=ax[1,col].pcolormesh(x,axis[m],residual[m],cmap='gray',shading='nearest',vmin=-bound,vmax=bound)
        ax[1,col].plot(x,arrays['expected'],'r--',label='local-depth approximation')
        ax[1,col].plot(x,r['full_field_ray_proxy']['group_time_proxy_ns'],'b:',label='full-field ray proxy')
        ax[1,col].set_ylim(220,160); ax[1,col].set_ylabel('Time (ns)'); ax[1,col].set_title(title)
        fig.colorbar(im,ax=ax[1,col],label='Signed residual; common fixed scale')
    for a in ax.flat:
        a.set_xlabel('Actual Tx/Rx midpoint y (m)'); a.grid(alpha=.15)
    fig.suptitle(name+': 3D morphology NOT certified; template/picks are diagnostics only')
    fig.savefig(out/f'{name}_shape_audit.png',dpi=140,bbox_inches='tight'); plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,help='NEW directory outside the input capsule')
    args=ap.parse_args(); out=Path(args.out)
    require(not out.exists() and not out.resolve().is_relative_to(HS.resolve()),'unsafe/existing output')
    verify_official_runtime(); sanity=sanity_checks()
    manifest_digest=sha256(HS/'manifest.json'); records=read_manifest(HS)
    identities=[verify_file(HS,records,k) for k in records]
    with np.load(HS/'hs4_interface_binned_table.npz') as z:
        table=z['z0_m'].copy()
    flat=trace_time_response(HS/'hs1_flat_halfspace.h5')
    bg=trace_time_response(HS/'hs2_coveronly_halfspace.h5')
    ref_t=flat.time*1e9
    require(np.array_equal(ref_t,bg.time*1e9),'reference axes differ')
    envelope=2*np.abs(flat.complex_envelope-bg.complex_envelope)
    mask=(ref_t>=160)&(ref_t<=220); anchor=float(ref_t[mask][np.argmax(envelope[mask])])
    template=flat.real_bandpass-bg.real_bandpass
    with h5py.File(HS/'hs1_flat_halfspace.h5','r') as h:
        positions={'source':np.asarray(h['srcs/src1'].attrs['Position']).tolist(),
                   'receiver':np.asarray(h['rxs/rx1'].attrs['Position']).tolist()}
    with h5py.File(HS/'hs2_coveronly_halfspace.h5','r') as h:
        require(np.array_equal(h['srcs/src1'].attrs['Position'],positions['source'])
                and np.array_equal(h['rxs/rx1'].attrs['Position'],positions['receiver']),
                'HS1/HS2 positions differ')
    source_flat=load_source(HS/'hs1_flat_halfspace.h5')
    source_bg=load_source(HS/'hs2_coveronly_halfspace.h5')
    require(source_flat.dt==source_bg.dt and source_flat.time_offset==source_bg.time_offset
            and np.array_equal(source_flat.samples,source_bg.samples),'HS1/HS2 source differs')
    groups={}; data={}
    for name,prefix,count,npz in [('col6','hs4_col6',25,'hs4_col6_bscan_official.npz'),
                                  ('co13','hs4_co13',13,'hs4_co13_bscan_official.npz')]:
        groups[name],data[name]=audit_group(prefix,count,npz,table,ref_t,template,anchor,positions)
    for k in records:
        verify_file(HS,records,k)
    require(sha256(HS/'manifest.json')==manifest_digest,'manifest changed')
    out.mkdir(parents=True)
    matched_contrast={}
    for name in groups:
        matches=groups[name]['matched_HS1_HS2_reference_trace_numbers']
        for trace_number in matches:
            k=trace_number-1
            c=data[name]['complex_envelope'][:,k]-bg.complex_envelope
            signed=data[name]['signed'][:,k]-bg.real_bandpass
            env=2*np.abs(c); envelope_peak=float(ref_t[mask][np.argmax(env[mask])])
            matched_contrast[f'{name}_t{trace_number:02d}']={
                'role':'exact-position matched media contrast only; not whole-line clean truth',
                'rough_contrast_env_peak_ns':envelope_peak,'flat_contrast_env_peak_ns':anchor,
                'envelope_peak_difference_ns':envelope_peak-anchor,
                'rough_contrast_env_peak':float(env[mask].max()),
                'flat_contrast_env_peak':float(envelope[mask].max()),
                'signed_template_fit':lag_fit(ref_t,signed[:,None],ref_t,template,(160.,220.)),
                'local_depth_approx_ns':groups[name]['local_depth_approx_ns'][k]}
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
            for label,sig,e in [('HS1-HS2 flat contrast',template,envelope),('HS4-HS2 rough contrast',signed,env)]:
                axes[0].plot(ref_t[mask],sig[mask],label=label)
                axes[1].plot(ref_t[mask],e[mask],label=label)
            for a in axes:
                a.set_xlabel('Time (ns)');a.legend(fontsize=8);a.grid(alpha=.2)
            axes[0].set_ylabel('Signed contrast');axes[1].set_ylabel('Official contrast envelope 2|c-c_BG|')
            fig.suptitle('Only CO13 t07 has this exact matched acquisition; waveform peaks are not a morphology gate')
            fig.savefig(out/f'{name}_t{trace_number:02d}_matched_contrast.png',dpi=140,bbox_inches='tight');plt.close(fig)
    for name in groups:
        plot(out,name,groups[name],data[name])
        np.savez_compressed(out/f'{name}_audit_arrays.npz',**data[name])
        r=groups[name]
        with (out/f'{name}_per_trace.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.writer(stream); writer.writerow(['trace','midpoint_y_m','depth_midpoint_m','local_approx_ns','full_field_ray_proxy_ns','official_env_window_peak_ns','rank1_abs_peak_ns','rank1_template_lag_ns'])
            for k,x in enumerate(data[name]['midpoint_y']):
                pick=r['picks']['rank1_full0_250']['160_220']
                writer.writerow([k+1,x,r['geometry']['midpoint']['cover_depth_m'][k],data[name]['expected'][k],r['full_field_ray_proxy']['group_time_proxy_ns'][k],r['official_interface_window_envelope_peak_ns'][k],pick['signed_abs_peak_ns'][k],pick['template_fit']['lag_ns'][k]])
    cells=C0/(FREQ*np.real(cover_index(FREQ))*.05)
    metrics={'schema':'hs4-3d-shape-audit/0.1','solver_executed':False,'input_capsule_unchanged':True,
        'input_manifest_sha256':manifest_digest,'input_identities':identities,'script_sha256':sha256(Path(__file__)),
        'identity_helper_sha256':sha256(Path(__file__).with_name('hs_capsule_identity.py')),
        'loader_sha256':sha256(Path(__file__).with_name('sfcw_official_loader_v0_2.py')),
        'official_processing_sha256':OFFICIAL_PROCESSING_SHA256,'sanity_checks':sanity,'reference_positions':positions,
        'groups':groups,'matched_position_contrasts':matched_contrast,
        'validation':'3D morphology not established; no matched-position COL6 planar line or near-x-edge domain convergence control',
        'spatial_resolution_diagnostic':{'grid_m':.05,'cover_cells_per_phase_wavelength_at_170MHz':float(cells[-1]),
            'frequencies_Hz_with_fewer_than_10_cells_per_wavelength':FREQ[cells<10].tolist(),
            'role':'resolution screening only; 10 cells is official guidance, not a convergence certificate. Hann reduces high-end weights but does not prove fidelity.'},
        'external_source':'https://docs.gprmax.com/en/latest/inc_SFCW.html (4.0.1 current docs; runtime is pinned local 4.0.0 code)'}
    (out/'metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=1,allow_nan=False),encoding='utf-8')
    (out/'manifest.json').write_text(json.dumps([{'file':p.name,'sha256':sha256(p),'bytes':p.stat().st_size}
        for p in sorted(out.iterdir()) if p.is_file()],indent=1),encoding='utf-8')
    print(json.dumps({k:{'env_edge_picks':v['official_envelope_window_edge_count'],
        'residual_energy_fraction':v['rank1_residual_total_energy_fraction'],
        'rank1_peak_depth_correlation':v['picks']['rank1_full0_250']['160_220']['peak_vs_local_depth_correlation'],
        'rank1_fit_depth_correlation':v['picks']['rank1_full0_250']['160_220']['fit_vs_local_depth_correlation'],
        'matched_reference_traces':v['matched_HS1_HS2_reference_trace_numbers']} for k,v in groups.items()},indent=1))


if __name__=='__main__':
    main()
