"""Audit/comparison of a completed local 8 m replay (CPU only)."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from hs_capsule_identity import sha256
from analyze_hs4_height_wavefield import response
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response
from analyze_hs4_height8m_wavefield import interface_relief

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT/'artifacts/research_checks/2026-10-04_hs4_height8m_wavefield_c'
BANDS = {'interface': (8.5, 9.6), 'cover': (9.6, 12),
         'air_low': (12, 19), 'antenna': (19, 21)}


def relative_l2(a, b):
    den = np.linalg.norm(b)
    return None if den == 0 else float(np.linalg.norm(a-b)/den)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--capsule', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new analysis directory required')
    c = json.loads((args.capsule/'execution_contract.json').read_text('utf-8'))
    v = json.loads((args.capsule/'completed_verification.json').read_text('utf-8'))
    old_v = json.loads((REFERENCE/'completed_verification.json').read_text('utf-8'))
    if v['status'] != 'PASS' or not v['completed'] or v['contract_sha256'] != sha256(args.capsule/'execution_contract.json'):
        raise ValueError('completed verified local capsule required')
    events = [json.loads(s) for s in (args.capsule/'execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1]['status'] != 'COMPLETED':
        raise ValueError('completed solver execution required')
    verify_official_runtime()
    spectra, rows = {}, {}
    for group, old in zip(v['groups'], old_v['groups']):
        if group['id'] != old['id']:
            raise ValueError('case pairing differs')
        name = group['id']
        local, ref = args.capsule/name/'profile.h5', REFERENCE/name/'profile.h5'
        if sha256(local) != group['raw_sha256'] or sha256(ref) != old['raw_sha256']:
            raise ValueError('raw receiver identity differs')
        with h5py.File(local) as h, h5py.File(ref) as r:
            if h.attrs['dt'] != r.attrs['dt']:
                raise ValueError('native clocks differ')
            a, b = h['rxs/rx1/Ey'][:], r['rxs/rx1/Ey'][:]
            if a.dtype != np.float64 or b.dtype != np.float64:
                raise ValueError('native double required')
        new, previous = response(local), response(ref)
        spectra[name] = new
        rows[name] = {'native_Ey_bit_equal': a.tobytes()==b.tobytes(),
                      'native_Ey_max_abs_difference': float(np.max(np.abs(a-b))),
                      'native_Ey_relative_L2': relative_l2(a,b),
                      'official501_response_relative_L2': relative_l2(new.response, previous.response),
                      'local_raw_sha256': group['raw_sha256'], 'reference_raw_sha256': old['raw_sha256']}
    a = spectra['mid8_rough']; b = spectra['mid8_halfspace']
    contrast = a.response-b.response
    old_contrast = response(REFERENCE/'mid8_rough/profile.h5').response-response(REFERENCE/'mid8_halfspace/profile.h5').response
    new_time = reconstruct_time_response(replace(a,response=contrast),window='hann',zero_pad_factor=8)
    old_time = reconstruct_time_response(replace(a,response=old_contrast),window='hann',zero_pad_factor=8)
    time_ns = new_time.time*1e9
    # Predeclared 15m window shifted by the vertical air-path proxy, not event picking.
    lo, hi = np.array([160.,220.])-2*7/299792458*1e9
    mask = (time_ns>=lo)&(time_ns<=hi)
    energy = {name: [] for name in BANDS}
    z = 8+.15*(np.arange(134)+.5)
    iters = np.asarray(c['snapshot_iterations'])
    case_map = {g['id']:g for g in v['groups']}
    snapshot_byte_equal = {'mid8_rough':0, 'mid8_halfspace':0}
    frame_times=iters*c['dt_s']*1e9
    shown={int(np.argmin(abs(frame_times-t))):None for t in (87,98,131,140)}
    global_diff_max=0.
    for index, iteration in enumerate(iters):
        fields=[]
        for name in snapshot_byte_equal:
            entry=case_map[name]['snapshots'][index]
            old=next(g for g in old_v['groups'] if g['id']==name)['snapshots'][index]
            path=args.capsule/name/entry['file']
            actual=sha256(path)
            if actual!=entry['sha256']:
                raise ValueError('local frame identity differs')
            snapshot_byte_equal[name]+=int(actual==old['sha256'])
            with h5py.File(path) as h:
                if h.attrs['iteration']!=iteration:
                    raise ValueError('frame timeline differs')
                fields.append(h['Ey'][:,0,:])
        diff=fields[0]-fields[1]
        global_diff_max=max(global_diff_max,float(abs(diff).max()))
        if index in shown:
            shown[index]=diff.copy()
        for name,(z0,z1) in BANDS.items():
            energy[name].append(float(np.sum(diff[:,(z>=z0)&(z<z1)]**2)))
    peaks={name:float(iters[np.argmax(values)]*c['dt_s']*1e9) for name,values in energy.items()}
    result={'status':'COMPLETED_LOCAL_REPLAY_AND_WAVEFIELD_DIAGNOSTIC',
            'code_sha256':sha256(__file__), 'contract_sha256':sha256(args.capsule/'execution_contract.json'),
            'changed_binary_count':len(c['changed_binary_files_from_reference']),
            'receiver_comparison':rows, 'snapshot_file_byte_equal_count':snapshot_byte_equal,
            'snapshot_count_per_case':len(iters), 'spatial_band_Ey_square_peak_ns':peaks,
            'contrast_official501_relative_L2':relative_l2(contrast,old_contrast),
            'contrast_complex_envelope_relative_L2_in_fixed_window':
                relative_l2(new_time.complex_envelope[mask],old_time.complex_envelope[mask]),
            'comparison_window_ns':[float(lo),float(hi)],
            'passive_receiver_bit_identical':v['passive_receiver_bit_identical'],
            'source_all501_valid':True, 'physical_attribution_certified':False,
            'limitations':['NativeRicker spatialEy squared is not directional energy flux.',
                           'One centre station at8m; no new height/grid/3D acceptance.',
                           'Identical field files certify this replay, not continuum or field validity.']}
    args.out.mkdir(parents=True)
    (args.out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    np.savez_compressed(args.out/'receiver_arrays.npz',time_ns=time_ns,
                        contrast_complex_envelope=new_time.complex_envelope,
                        reference_contrast_complex_envelope=old_time.complex_envelope,
                        snapshot_time_ns=iters*c['dt_s']*1e9, **{name+'_Ey_square':np.asarray(values) for name,values in energy.items()})
    fig,axes=plt.subplots(1,4,figsize=(14,5),sharex=True,sharey=True,layout='constrained')
    xx=13.5+.15*(np.arange(60)+.5)
    iface_x,iface_z=interface_relief()
    norm=SymLogNorm(linthresh=global_diff_max*.01,vmin=-global_diff_max,vmax=global_diff_max)
    for ax,(index,field) in zip(axes,shown.items()):
        im=ax.pcolormesh(xx,z,field.T,cmap='RdBu_r',norm=norm,shading='nearest')
        ax.plot(iface_x,iface_z,'k-',lw=.8); ax.axhline(12,color='k',ls='--',lw=.7)
        ax.scatter([17.6,18.9],[20,20],marker='x',c='green')
        ax.set(xlim=(13.5,22.5),ylim=(8,28.1),xlabel='x (m)',title=f'{frame_times[index]:.1f} ns')
    axes[0].set_ylabel('z (m)')
    fig.colorbar(im,ax=axes,label='rough - full-cover Ey; shared SymLog')
    fig.suptitle('Local V4/double, 8 m: interface-contrast field; raw Ricker95, not an SFCW image')
    fig.savefig(args.out/'local_wavefield_atlas.png',dpi=130); plt.close(fig)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
