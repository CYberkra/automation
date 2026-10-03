"""Third-grid centre contraction, plus explicitly limited spatial ray diagnostics."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

from analyze_hs4t2d_joint_grid import compare, completed_study
from analyze_hs4t2d_cause_controls import load_response, verify_manifest, weighted_comparison
from hs4t2d_path_diagnostic import ray_timings
from hs_capsule_identity import sha256
from sfcw_official_loader_v0_2 import verify_official_runtime, OFFICIAL_PROCESSING_SHA256

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--third', type=Path, required=True)
    p.add_argument('--scan-results', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    if args.out.exists():
        raise ValueError('new results required')
    verify_official_runtime()
    c = completed_study(args.third, 'centre')
    if c.get('resolution_level') != 'third_1over60m':
        raise ValueError('wrong third-resolution study')
    verify_manifest(args.scan_results)
    report_before = json.loads((args.scan_results/'results.json').read_text('utf-8'))
    if report_before['station_count'] != 13 or report_before['studies'][0]['contract_sha256'] != c['prerequisite']['contract_sha256']:
        raise ValueError('wrong complete2.5cm scan')
    with np.load(args.scan_results/'joint_grid_arrays.npz') as d:
        time, x = d['time_ns'].copy(), d['midpoint_m'].copy()
        centre = int(np.flatnonzero(d['old_station_indices']==61)[0])
        f = d['frequency_Hz'].copy()
        coarse = d['coarse_frequency_complex'][:,centre:centre+1].copy()
        fine = d['fine_frequency_complex'][:,centre:centre+1].copy()
        full_fine_signed = d['fine_signed'].copy()
    responses, source_audit = {}, {}
    for g in c['groups']:
        r,rows = load_response(args.third/g['id'])
        if not r.source_valid.all() or not np.array_equal(r.frequency,f) or len(rows)!=1:
            raise ValueError('invalid third response')
        if not np.allclose(rows[0]['source_m'],[17.6,0,27],rtol=0,atol=1e-12) or not np.allclose(rows[0]['receiver_m'],[18.9,0,27],rtol=0,atol=1e-12):
            raise ValueError('third station shifted')
        expected_source = next(iter(report_before['source_audit'].values()))
        if r.source.spatial_scale != expected_source['spatial_scale'] or r.source.quantity != expected_source['source_quantity'] or r.source.units != expected_source['source_units']:
            raise ValueError('source scaling incompatible')
        responses['halfspace' if g['halfspace'] else 'rough'] = r
        source_audit[g['id']] = {'quantity':r.source.quantity,'units':r.source.units,
            'spatial_scale':r.source.spatial_scale,'dt_s':r.source.dt,'all501source_bins_valid':True}
    third = responses['rough'].response-responses['halfspace'].response
    spectra = {'coarse5cm':coarse,'fine2p5cm':fine,'third1over60m':third}
    products = {name:reconstruct_time_response(replace(responses['rough'],response=s),window='hann',zero_pad_factor=8)
                for name,s in spectra.items()}
    if any(not np.array_equal(p.time*1e9,time) for p in products.values()):
        raise ValueError('time axes differ')
    mask = (time>=160)&(time<=220)
    metrics = {}
    for name,a,b in [('5_to_2p5','coarse5cm','fine2p5cm'),('2p5_to_third','fine2p5cm','third1over60m')]:
        metrics[name] = {'hann_complex_spectrum':weighted_comparison(spectra[a],spectra[b]),
                         'signed_band':compare(products[a].real_bandpass[mask],products[b].real_bandpass[mask]),
                         'phase_aware_complex_envelope':compare(products[a].complex_envelope[mask],products[b].complex_envelope[mask]),
                         'envelope_magnitude':compare(2*np.abs(products[a].complex_envelope[mask]),2*np.abs(products[b].complex_envelope[mask]))}
    profile_path = ROOT/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv'
    profile = np.genfromtxt(profile_path,delimiter=',',names=True)
    paths = {str(height):ray_timings(profile['grid_z_m'],x,height) for height in (15,2)}
    path_stats = {height:{'minimum_ray_group_time_span_ns':float(np.ptp(r['minimum_ray_group_time_ns'])),
                           'vertical_group_time_span_ns':float(np.ptp(r['vertical_group_time_ns'])),
                           'minimum_ray_reflection_x_range_m':[min(r['minimum_ray_reflection_x_m']),max(r['minimum_ray_reflection_x_m'])]}
                  for height,r in paths.items()}
    # Exploratory track of one named positive lobe, separate from envelope maxima.
    # The fixed interval defines an observable only, not an identified interface.
    lobe_mask = (time>=160)&(time<=180)
    positive_lobe = time[lobe_mask][np.argmax(full_fine_signed[lobe_mask],axis=0)]
    args.out.mkdir(parents=True)
    payload = {'time_ns':time,'midpoint_m':x,'frequency_Hz':f,'fine_scan_positive_lobe_ns':positive_lobe}
    peaks = {}
    for name,p in products.items():
        payload[name+'_signed'] = p.real_bandpass
        payload[name+'_complex_envelope'] = p.complex_envelope
        payload[name+'_frequency_complex'] = spectra[name]
        peaks[name] = float(time[mask][np.argmax(np.abs(p.complex_envelope[mask,0]))])
    np.savez_compressed(args.out/'third_grid_arrays.npz',**payload)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    fig,axes = plt.subplots(2,1,figsize=(10,7),constrained_layout=True)
    for name,p in products.items():
        axes[0].plot(time[mask],p.real_bandpass[mask,0],label=name)
        axes[1].plot(time[mask],2*np.abs(p.complex_envelope[mask,0]),label=name)
    for ax in axes:
        ax.legend(); ax.set_xlim(160,220); ax.set_xlabel('时间（ns）'); ax.set_ylabel('场/电流响应，无归一化')
    axes[0].set_title('保留符号和绝对时轴'); axes[1].set_title('包络模值仅作补充')
    fig.suptitle('中心道三档联合X/Z网格 / 同一实际15m条件 / 相同物理PML厚度\n20–170MHz / Hann / 匹配全覆盖层复谱差分')
    fig.savefig(args.out/'third_grid_center.png',dpi=145); plt.close(fig)
    fig,axes = plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    q = np.r_[profile['x0_m'],profile['x1_m'][-1]]
    axes[0].stairs(12-profile['grid_z_m'],q,baseline=None,label='实际界面深度')
    axes[0].set_ylim(3.4,2.4); axes[0].set_ylabel('距地表深度（m）'); axes[0].set_title('候选最早路径落点')
    for height,r in paths.items():
        reflect = r['minimum_ray_reflection_x_m']
        depths = 12-profile['grid_z_m'][np.clip(np.floor(np.array(reflect)/.25).astype(int),0,47)]
        axes[0].scatter(reflect,depths,label=f'{height}m高度，13候选落点')
        vertical=np.array(r['vertical_group_time_ns']); earliest=np.array(r['minimum_ray_group_time_ns'])
        axes[1].plot(x,vertical-vertical[centre],label=f'{height}m，局部垂直代理')
        axes[1].plot(x,earliest-earliest[centre],linestyle='--',label=f'{height}m，最早候选路径')
    axes[1].set_ylabel('相对中心站时延（ns）'); axes[1].set_title('同一高度内仅减去中心常量')
    axes[2].plot(x,positive_lobe,label='2.5cm图中160–180ns正波瓣最大值')
    axes[2].plot(x,paths['15']['minimum_ray_group_time_ns'],label='15m最早路径群时延代理')
    axes[2].plot(x,paths['15']['vertical_group_time_ns'],label='15m局部垂直群时延代理')
    axes[2].set_ylabel('时间（ns）'); axes[2].set_title('不拟合逐道时延，曲线不是事件真值')
    for ax in axes:
        ax.set_xlabel('原模型X（m）'); ax.legend(fontsize=8)
    fig.suptitle('空间路径诊断：95MHz色散折射率、平地表折射；不含散射振幅/反射相位\n最早路径≠最强回波；曲线接近不能证明条带来自某段界面')
    fig.savefig(args.out/'spatial_path_diagnostic.png',dpi=145); plt.close(fig)
    report = {'status':'COMPLETED','new_third_grid_traces':2,'scope':'Centre contraction and limited ray-geometry diagnostics; not full-space convergence or unique event attribution.',
              'contract_sha256':sha256(args.third/'execution_contract.json'),
              'scan_manifest_sha256':sha256(args.scan_results/'manifest.json'),
              'profile_sha256':sha256(profile_path),'official_processing_sha256':OFFICIAL_PROCESSING_SHA256,
              'code_identities':{name:sha256(ROOT/'scripts'/name) for name in ['analyze_hs4t2d_joint_grid_third.py','analyze_hs4t2d_joint_grid.py','hs4t2d_path_diagnostic.py']},
              'source_audit':source_audit,'comparisons':metrics,'window_envelope_peak_ns':peaks,
              'peak_scope':'160–220ns sampled maxima, may switch lobes and limited by reconstruction time grid; not certified event arrival.',
              'path_diagnostics_95MHz':paths,'path_span_diagnostics':path_stats,
              'exploratory_positive_lobe':{'window_ns':[160,180],'times_ns':positive_lobe.tolist(),'span_ns':float(np.ptp(positive_lobe)),
                                         'scope':'post-scan observable, no physical event label or fitting/selection qualification'},
              'limitations':['Third grid covers centre only. Adjacent differences are not absolute error bounds.',
                             'Ray minima omit reflection strength, interference, finite-band phase and strong multiple scattering; do not identify the dominant B-scan event.',
                             '2D invariant line source is not a finite3D antenna. No3D/field/training validation.']}
    (args.out/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)} for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'comparisons':metrics,'window_envelope_peak_ns':peaks,'path_span_diagnostics':path_stats},indent=2))


if __name__ == '__main__':
    main()
