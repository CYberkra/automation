"""Existing-data height comparison; CPU only, no production geometry correction."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from audit_line9_postprocessing import FREQ, inverse, weights, sha, save
from review_line9_result_packages import C0, indices
from diagnose_line9_layer_kinematics import primary_paths
from trace_line9_time_origin import correlation


def phase_advance(response, delay_s):
    return response * np.exp(2j*np.pi*FREQ[:, None]*delay_s)


def summarize(low, high, masks):
    result = {}
    coefficients = {}
    for key, mask in masks.items():
        per_trace = []; ratios = []; fitted = []
        for j in range(low.shape[1]):
            a, b = low[mask[:, j], j], high[mask[:, j], j]
            per_trace.append(correlation(a, b))
            ratios.append(20*np.log10(np.linalg.norm(b)/np.linalg.norm(a)))
            fitted.append(np.vdot(a, b)/np.vdot(a, a))
        coefficients[key] = np.array(fitted)
        result[key] = dict(global_complex_correlation=correlation(low[mask], high[mask]),
            per_trace_complex_correlation_p10_50_90=np.percentile(per_trace, [10,50,90]).tolist(),
            per_trace_complex_correlation=per_trace,
            gate_total_field_norm_ratio_high_over_low_db_p10_50_90=np.percentile(ratios,[10,50,90]).tolist(),
            gate_total_field_norm_ratio_high_over_low_db=ratios,
            fitted_coefficient_real_imag=np.column_stack([np.real(fitted),np.imag(fitted)]).tolist())
    # Fit exclusively in the cover-base gate, evaluate the separate basal gate.
    mask = masks['basal_sand']
    predicted = low * coefficients['cover_base'][None, :]
    result['cover_fit_tested_on_basal'] = dict(
        global_complex_correlation=correlation(predicted[mask], high[mask]),
        relative_L2=float(np.linalg.norm(predicted[mask]-high[mask])/np.linalg.norm(high[mask])),
        phase_basal_over_cover_coefficient_deg_p10_50_90=np.percentile(
            np.angle(coefficients['basal_sand']/coefficients['cover_base'],deg=True),[10,50,90]).tolist(),
        scope='Diagnostic using the low-height simulation as reference, not a single-scan deployable calibration')
    return result


def checks():
    # Independent delayed synthetic primary. Inverse phase advance must preserve amplitude.
    source = np.exp(-2j*np.pi*FREQ[:,None]*np.array([[160.,185.]])*1e-9)
    delays = np.array([[61.,64.]])*1e-9
    delayed = source*np.exp(-2j*np.pi*FREQ[:,None]*delays)
    reconstructed = phase_advance(delayed,delays)
    error = float(np.linalg.norm(source-reconstructed)/np.linalg.norm(source))
    amplitude_error = float(np.max(abs(abs(delayed)-abs(reconstructed))))
    assert error < 1e-12 and amplitude_error < 1e-12
    zero_error = float(np.max(abs(phase_advance(source,np.zeros((1,2)))-source)))
    assert zero_error == 0
    rng=np.random.default_rng(108)
    low=rng.normal(size=(60,3))+1j*rng.normal(size=(60,3))
    masks={k:np.broadcast_to(((np.arange(60)>=start)&(np.arange(60)<start+10))[:,None],low.shape)
           for k,start in [('cover_base',0),('first_sand',20),('basal_sand',40)]}
    high=low*np.array([.5+.3j,.3-.4j,.4+.2j])[None,:]
    constant=summarize(low,high,masks)['cover_fit_tested_on_basal']['relative_L2']
    assert constant < 1e-12
    high[masks['basal_sand']] *= 1j
    depth_dependent=summarize(low,high,masks)['cover_fit_tested_on_basal']['relative_L2']
    assert abs(depth_dependent-np.sqrt(2)) < 1e-12
    return dict(known_delay_complex_relative_L2=error,phase_only_magnitude_max_error=amplitude_error,
                zero_delay_max_error=zero_error,common_complex_gain_transfer_relative_L2=constant,
                depth_dependent_phase_reject_relative_L2=depth_dependent)


def main(root, review, cache, ray_path, out):
    if out.exists():
        raise ValueError('Use a fresh output directory')
    audits = [json.loads((review/(name+'_audit.json')).read_text('utf-8')) for name in
              ['line9_pkg2_x160-110_agl35_142st','line9_pkg3_x160-120_agl1_70st']]
    a,b = audits; count = len(b['records']); assert count == 70
    rays = json.loads(ray_path.read_text('utf-8'))
    rp = {p['package']: p for p in rays['packages']}
    fields = []; source_info = []; geometries = []
    for audit in audits:
        package = root/audit['package']; cp = cache/(audit['package']+'.npz')
        gf = next((package/'geometries').glob('*.h5'))
        mf = next((package/'geometries').glob('*.json'))
        assert sha(gf) == audit['geometry_sha256']
        assert sha(mf) == audit['materials_sha256']
        assert sha(cp) == rp[audit['package']]['private_cache_sha256']
        with np.load(cp) as z:
            np.testing.assert_array_equal(z['frequency_Hz'],FREQ)
            np.testing.assert_array_equal(z['ids'],[r['id'] for r in audit['records']])
            np.testing.assert_allclose(z['chainage_m'],[r['chainage_m'] for r in audit['records']],atol=1e-10,rtol=0)
            fields.append(z['response'][:,:count].copy())
        positions = []
        for row in audit['records'][:count]:
            raw = package/'cases'/row['id']/'profile.h5'
            assert sha(raw) == row['native_sha256']
            with h5py.File(raw) as h:
                positions.append([h['srcs/src1'].attrs['Position'],h['rxs/rx1'].attrs['Position']])
        with h5py.File(gf) as h:
            geometries.append(h['data'][:,:,0])
        source_info.append(dict(package=audit['package'],cache_sha256=sha(cp),
                               geometry_sha256=sha(gf),materials_sha256=sha(mf),positions=np.array(positions)))
    assert a['materials_sha256'] == b['materials_sha256']
    np.testing.assert_allclose(source_info[0]['positions'][:,:,0],source_info[1]['positions'][:,:,0],atol=1e-10,rtol=0)
    for ra,rb in zip(a['records'][:count],b['records']):
        assert abs(ra['chainage_m']-rb['chainage_m']) < 1e-10
        np.testing.assert_allclose([[v[k] for k in ['y_m','above','below']] for v in ra['geometry']['boundaries']],
                                   [[v[k] for k in ['y_m','above','below']] for v in rb['geometry']['boundaries']],atol=1e-10,rtol=0)
    # Same origin and entire common voxel region; right domain edge still differs by 10 m.
    nx = min(g.shape[0] for g in geometries)
    assert np.array_equal(geometries[0][:nx],geometries[1][:nx])
    high, low = fields
    materials_file=next((root/b['package']/'geometries').glob('*.json'))
    materials=json.loads(materials_file.read_text('utf-8'))['materials']
    n=np.array([indices(materials,f) for f in FREQ])
    air_delay=np.array([ra['geometry']['midpoint_agl_m']-rb['geometry']['midpoint_agl_m']
                        for ra,rb in zip(a['records'][:count],b['records'])])[None,:]*2/C0
    path_arrays=[]
    for audit in audits:
        records=rp[audit['package']]['records']
        assert [r['station_index'] for r in records] == list(range(count))
        assert all(r['geometry_check']=='NO_OUTSIDE_BAND_MISMATCH' for r in records)
        path_arrays.append(np.array([r['path_length_by_material_m'] for r in records]))
    ray_delay=n.real@(path_arrays[0]-path_arrays[1]).T/C0
    variants={'air_delay_only':phase_advance(high,air_delay),
              'basal_ray_phase_only_proxy':phase_advance(high,ray_delay)}
    primary=[primary_paths(r['geometry'],materials) for r in b['records']]
    result=dict(calls_solver=False,calls_training=False,script_sha256=sha(__file__),self_checks=checks(),
                source_git_base='d2fbf2896cda132b1fabb825d2ce2462e267a350',
                ray_json_sha256=sha(ray_path),sources=[{k:v for k,v in d.items() if k!='positions'} for d in source_info],
                paired_station_count=count,common_voxel_region_equal=True,
                right_domain_boundary_difference_m=10.,height_is_only_changed_condition=False,
                chainage_m=[r['chainage_m'] for r in b['records']],
                air_advance_ns=(air_delay.ravel()*1e9).tolist(),
                basal_ray_minus_air_delay95_ns=((ray_delay[250]-air_delay.ravel())*1e9).tolist(),
                variants={})
    processed={}; centers_by_window={}
    for window in ['hann','blackman']:
        vl,t=inverse(low,FREQ,weights(window,501)); centers={}; masks={}
        for key in ['cover_base','first_sand','basal_sand']:
            spectra=np.column_stack([p[r['geometry']['boundaries'].index(r['geometry'][key])]
                                     for p,r in zip(primary,b['records'])])
            model,_=inverse(spectra,FREQ,weights(window,501))
            centers[key]=t[np.argmax(abs(model),axis=0)]*1e9
            masks[key]=abs(t[:,None]*1e9-centers[key][None,:])<=12
        centers_by_window[window]=centers
        for name,field in variants.items():
            vh,_=inverse(field,FREQ,weights(window,501))
            metrics=summarize(vl,vh,masks)
            if name!='air_delay_only':
                # A basal path is not an appropriate shallow alignment/calibration.
                metrics={'basal_sand':metrics['basal_sand']}
            metrics['low_model_primary_peak_ns']={k:v.tolist() for k,v in centers.items()}
            result['variants'][name+'_'+window]=metrics
            if window=='hann':processed[name]=vh
        if window=='hann':processed['low']=vl
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    x=np.array(result['chainage_m']);order=np.argsort(x);keep=(t*1e9>=135)&(t*1e9<=315)
    # One fixed physical display reference shared across the two height runs.
    reference=float(np.max(abs(processed['low'])))
    titles=['1m：Hann原始总场','高航高：仅扣空气双程时延','高航高：底砂二维路径相位对齐（诊断）']
    fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained')
    for column,(name,title) in enumerate(zip(['low','air_delay_only','basal_ray_phase_only_proxy'],titles)):
        db=20*np.log10(np.maximum(abs(processed[name][keep][:,order])/reference,1e-12))
        for row in range(2):
            ax=axs[row,column]
            pic=ax.pcolormesh(x[order],t[keep]*1e9,db,cmap='gray_r',vmin=-100,vmax=-35,shading='nearest',rasterized=True)
            if row:
                for key,label,color in [('cover_base','覆盖层底','#009e73'),('first_sand','首层砂岩','#e69f00'),('basal_sand','底砂','#d55e00')]:
                    ax.plot(x[order],centers_by_window['hann'][key][order],color=color,lw=1,label=label)
                ax.legend(fontsize=8,loc='lower right')
            ax.set_title(title+('；无参考线' if row==0 else '；叠低航高近似界面峰位'),fontsize=11)
            ax.set_ylim(315,135);ax.invert_xaxis();ax.set_xlabel('剖面横坐标 X / m（采集方向向左）');ax.set_ylabel('相对低航高的时间 / ns')
    fig.colorbar(pic,ax=axs,label='复包络幅度 / dB；共同1m全时峰值参考；无AGC/逐道归一化')
    fig.suptitle('70共同站位：全频501点Hann；总场，不是配对差分；右边界不同10m，不能作纯航高因果实验',fontsize=13)
    fig.savefig(out/'height_pair_common_scale.png',dpi=160);plt.close(fig)
    metrics=result['variants']['air_delay_only_hann']
    fig,axs=plt.subplots(3,1,figsize=(12,10),layout='constrained')
    for key,label in [('cover_base','覆盖层底窗'),('first_sand','首层砂岩窗'),('basal_sand','底砂窗')]:
        axs[0].plot(x,metrics[key]['per_trace_complex_correlation'],label=label)
        axs[1].plot(x,metrics[key]['gate_total_field_norm_ratio_high_over_low_db'],label=label)
    c=np.array(metrics['cover_base']['fitted_coefficient_real_imag']);cb=c[:,0]+1j*c[:,1]
    c=np.array(metrics['basal_sand']['fitted_coefficient_real_imag']);cd=c[:,0]+1j*c[:,1]
    axs[2].plot(x,np.angle(cd/cb,deg=True),label='底砂/覆盖层底复系数相位差')
    axs[0].set_ylabel('单道复数相关（模）');axs[0].set_ylim(0,1.02)
    axs[1].set_ylabel('高/低航高窗总场范数比 / dB')
    axs[2].set_ylabel('相位差 / 度');axs[2].set_xlabel('剖面横坐标 X / m')
    for ax in axs:ax.legend();ax.invert_xaxis();ax.grid(alpha=.2)
    fig.suptitle('只对齐空气时延：单道相关会忽略每道相位/幅度差；图中各窗由模型预先给定，不是盲检')
    fig.savefig(out/'height_pair_phase_amplitude.png',dpi=160);plt.close(fig)
    result['display_shared_reference']=reference
    result['limits']=('Same common geology, but domains differ by 10m: not height-only attribution. '
        'All gates are model-informed +-12ns; total-field norms are neither target amplitude nor SNR. '
        'Per-trace correlation discards per-trace scalar phase/amplitude; global correlation retains their variation. '
        'Basal ray phase advance excludes absorption compensation and uses a target-specific approximate branch; '
        'not a production whole-trace correction or migration, and not an independent target identification. '
        'Native FP32 precision unchanged. No new FDTD, training, model-guided deletion or deep-delay fitting.')
    save(out/'height_pair.json',result)
    print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk in ['basal_sand','cover_fit_tested_on_basal']}
                      for k,v in result['variants'].items()},ensure_ascii=False)[:1000])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','review','cache','ray','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.root,a.review,a.cache,a.ray,a.out)
