"""Audit native identities and render the approved line with official direct SFCW."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import shutil

import gprMax
from gprMax.toolboxes.SFCW import processing as sf
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import audit_source, save

FREQ = 20e6+np.arange(501)*300000.


def relative(x, y):
    return float(np.linalg.norm(x-y)/np.linalg.norm(y))


def inspect(path, g):
    with h5py.File(path) as h:
        assert h.attrs['gprMax'] == '4.0.1' and h.attrs['Iterations'] == 40703
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'], [8400,1700,1])
        np.testing.assert_array_equal(h.attrs['dx_dy_dz'], [.025]*3)
        assert abs(float(h.attrs['dt'])/g['dt_s']-1) < 1e-13
        x = h['rxs/rx1/Ez'][:]; src = h['srcs/src1/excitation/samples'][:]
        assert x.dtype == src.dtype == np.float64 and x.shape == src.shape == (40703,)
        assert np.isfinite(x).all() and np.isfinite(src).all() and np.any(x)
        audit_source(h['srcs/src1/excitation'], g, float(h.attrs['dt']))
        for key, position in [('srcs/src1',g['tx_m']), ('rxs/rx1',g['rx_m'])]:
            np.testing.assert_allclose(h[key].attrs['Position'], position, rtol=0, atol=1e-12)
    return x


def transfer(path, full_check):
    s = sf.load_source(path); r = sf.load_receiver(path,'/rxs/rx1','Ez')
    assert s.spatial_scale == .025 and s.time_offset == .5*s.dt and r.time_offset == 0
    q = sf.direct_frequency_response(s, r, FREQ, tail_taper_fraction=0)
    assert q.source_valid.all()
    tones = np.arange(501) if full_check else np.unique(np.r_[np.arange(0,501,31),500])
    independent = np.empty(len(tones), complex)
    for k in range(0,len(tones),16):
        f = FREQ[tones[k:k+16],None]
        independent[k:k+16] = (r.dt*(np.exp(-2j*np.pi*f*r.times)@r.samples))/(
            s.dt*(np.exp(-2j*np.pi*f*s.times)@s.samples))/.025
    q = replace(q,response=q.response/.025)
    error = relative(q.response[tones], independent); assert error < 1e-9
    return q, dict(independent_DFT_tones=len(tones), relative_L2=error,
                   receiver_final_5pct_peak_relative_db=q.receiver_tail_relative_db)


def main(a):
    assert not a.out.exists() and not a.numerical.exists() and gprMax.__version__ == '4.0.1'
    c = json.loads((a.raw/'execution_contract.json').read_text('utf-8'))
    m = c['study_manifest']; xs = np.array([s['profile_x_m'] for s in m['stations']])
    assert len(xs) == 41 and np.array_equal(xs, np.arange(190.,169.99,-.5))
    assert sha(sf.__file__) == c['source_identities']['toolboxes/SFCW/processing.py']
    retrieval = json.loads((a.raw/'retrieval_manifest.json').read_text('utf-8-sig'))
    for record in retrieval:
        p = a.raw/record['file']
        assert sha(p) == record['sha256'] and p.stat().st_size == record['bytes']
    events = [json.loads(line) for line in (a.raw/'execution.jsonl').read_text('utf-8-sig').splitlines()]
    receipts = {e['group']: e for e in events if e['status'] == 'COMPLETED' and 'group' in e}
    if not a.allow_partial:
        v = json.loads((a.raw/'completed_verification.json').read_text('utf-8'))
        assert v['completed'] and v['contract_sha256'] == sha(a.raw/'execution_contract.json')
        assert len(v['groups']) == len(receipts) == 42
        for g in v['groups']:
            assert g['native_sha256'] == receipts[g['id']]['raw_sha256']
    rows = []; all_q = {}; paths = {}
    for g in c['groups']+m['reused_solver_outputs']:
        reused = 'native_path' in g
        if not reused and g['id'] not in receipts:
            continue
        path = Path(g['native_path']) if reused else a.raw/(g['id']+'.h5')
        expected = g['native_sha256'] if reused else receipts[g['id']]['raw_sha256']
        assert sha(path) == expected
        inspect(path,g)
        x = 190. if reused else g['profile_x_m']; role = g['role']
        assert (x,role) not in all_q
        q, check = transfer(path, x in [190.,180.,170.])
        all_q[x,role] = q; paths[x,role] = path
        rows.append(dict(id=g['id'],profile_x_m=x,role=role,reused=reused,
                         native_sha256=sha(path),native_dtype='float64',**check))
    base = all_q[190.,'H1']; z = np.full((501,41), np.nan+1j*np.nan)
    for k,x in enumerate(xs):
        if (x,'H1') in all_q:
            z[:,k] = all_q[x,'H1'].response
    valid = np.isfinite(z).all(axis=0)
    if not a.allow_partial: assert valid.all()
    anchors = [x for x in [190.,180.,170.] if (x,'H1') in all_q and (x,'H0') in all_q]
    delta = np.column_stack([all_q[x,'H1'].response-all_q[x,'H0'].response for x in anchors])
    # Do not pass missing traces through transforms: reconstruct only native completed columns.
    a.out.mkdir(parents=True); a.numerical.parent.mkdir(parents=True,exist_ok=True)
    profiles = {}; metrics = {}; inverse_errors = {}
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ); h.create_dataset('profile_x_m',data=xs)
        h.create_dataset('H1_response',data=z); h.create_dataset('H1_computed_mask',data=valid)
        h.create_dataset('anchor_x_m',data=anchors); h.create_dataset('anchor_delta_response',data=delta)
        for window in ['hann','blackman']:
            p = sf.reconstruct_time_response(replace(base,response=z[:,valid]),window=window,zero_pad_factor=8,time_shift=0)
            t = p.time*1e9; full = np.full((len(t),41),np.nan+1j*np.nan); full[:,valid] = p.complex_bandpass
            d = sf.reconstruct_time_response(replace(base,response=delta),window=window,zero_pad_factor=8,time_shift=0)
            k = np.arange(7,len(t),101)
            direct = np.exp(2j*np.pi*p.time[k,None]*FREQ)@(z[:,valid]*p.weights[:,None])/p.weights.sum()
            inverse_errors[window] = relative(p.complex_bandpass[k],direct)
            assert inverse_errors[window] < 1e-9
            h.create_dataset(f'{window}/H1_complex_bandpass',data=full)
            real = np.full(full.shape,np.nan); real[:,valid] = p.real_bandpass
            h.create_dataset(f'{window}/H1_real_bandpass',data=real)
            h.create_dataset(f'{window}/anchor_delta_complex_bandpass',data=d.complex_bandpass)
            if window == 'hann': h.create_dataset('time_s',data=p.time)
            profiles[window] = (p,full,d)
            deep = (t >= 250)&(t <= 550)
            line_peaks = []
            for j,x in enumerate(xs):
                if not valid[j]: continue
                segment = abs(full[deep,j]); peak = int(np.argmax(segment))
                line_peaks.append(dict(profile_x_m=float(x),deep_peak_abs=float(segment[peak]),deep_peak_ns=float(t[deep][peak])))
            paired = []
            for j,x in enumerate(anchors):
                h0 = sf.reconstruct_time_response(all_q[x,'H0'],window=window,zero_pad_factor=8,time_shift=0)
                target = abs(d.complex_bandpass[deep,j]); peak = int(np.argmax(target))
                paired.append(dict(profile_x_m=x,delta_peak_abs=float(target[peak]),delta_peak_ns=float(t[deep][peak]),
                    H0_peak_abs=float(abs(h0.complex_bandpass[deep]).max()),
                    delta_over_H0_peak=float(target.max()/abs(h0.complex_bandpass[deep]).max()),
                    vertical_95MHz_guide_ns=m['stations'][list(xs).index(x)]['guide_ns']['basal_sandstone']))
            metrics[window] = dict(predeclared_deep_window_ns=[250,550], H1_peaks=line_peaks, matched_anchors=paired)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    gray = plt.get_cmap('gray').copy(); gray.set_bad('white')
    count = int(valid.sum())
    for window,(p,full,d) in profiles.items():
        t = p.time*1e9
        fig,axs = plt.subplots(3,1,figsize=(13,13),layout='constrained',gridspec_kw={'height_ratios':[1,1, .8]})
        for ax,bounds,title in [(axs[0],[0,600],'H1总场：早段及深部'),(axs[1],[250,550],'H1总场：深部放大')]:
            mask = (t>=bounds[0])&(t<=bounds[1]); image_data = full[mask].real
            limit = float(np.nanmax(abs(image_data)))
            im = ax.imshow(image_data,cmap=gray,aspect='auto',vmin=-limit,vmax=limit,
                           extent=[190.25,169.75,t[mask][-1],t[mask][0]],interpolation='nearest')
            ax.set(xlabel='剖面 x / m（沿190→170方向）',ylabel='SFCW 时间 / ns',
                   title=title+'；本面板全部站位共用绝对灰度，未计算列留白')
            fig.colorbar(im,ax=ax,label='Re(复带通) / (V/m)/(A·m)')
        for label,key,color in [('覆盖层底：95MHz垂直群时延参考','cover_base','tab:blue'),
                                 ('底砂顶：95MHz垂直群时延参考','basal_sandstone','tab:green')]:
            target_ax = axs[0] if key == 'cover_base' else axs[1]
            target_ax.plot(xs,[s['guide_ns'][key] for s in m['stations']],ls=':',color=color,lw=1,label=label)
        axs[0].set_ylim(600,0); axs[1].set_ylim(550,250)
        axs[0].legend(fontsize=8)
        axs[1].legend(fontsize=8)
        ground = np.array([s['interfaces_y_m']['1'] for s in m['stations']])
        cover = np.array([s['interfaces_y_m']['2'] for s in m['stations']])
        rock = np.array([s['interfaces_y_m']['3'] for s in m['stations']])
        for low,high,label,color in [(rock-3,rock,'砂岩','#c9b784'),(rock,cover,'泥岩','#dab692'),
                                     (cover,ground,'粉质粘土','#decad4'),(ground,ground+2,'空气','#f4f4f4')]:
            axs[2].fill_between(xs,low,high,label=label,color=color)
        for y in [ground,cover,rock]: axs[2].plot(xs,y,color='black',lw=.7)
        axs[2].set(xlim=[190.25,169.75],ylim=[float(rock.min()-2),float(ground.max()+1)],
                   xlabel='剖面 x / m',ylabel='模型内部 y / m（非绝对高程）',title='同一H5模型的局部地质界面；不是雷达深度转换')
        axs[2].legend(ncol=4,fontsize=8)
        fig.suptitle(f'190→170m 短测线：已计算{count}/41道H1，0.5m间距，约8m离地\n'
            f'泥岩 εr18、σ0.006 S/m；4.0.1原生FP64；官方impulse/direct，20–170MHz、501点；{window}，补零8；无增益/拟合')
        fig.savefig(a.out/f'bscan_{window}_gray.png',dpi=140); plt.close(fig)
        fig,axs = plt.subplots(2,len(anchors),figsize=(5*len(anchors),10),layout='constrained',squeeze=False)
        for row,bounds in enumerate([[120,550],[250,550]]):
            display = (t>=bounds[0])&(t<=bounds[1])
            for j,x in enumerate(anchors):
                ax = axs[row,j]
                for role,color in [('H1','black'),('H0','tab:blue')]:
                    profile = sf.reconstruct_time_response(all_q[x,role],window=window,zero_pad_factor=8,time_shift=0)
                    ax.plot(t[display],abs(profile.complex_bandpass[display]),label=role+'总场包络',color=color,lw=.9)
                ax.plot(t[display],abs(d.complex_bandpass[display,j]),label='同站H1−H0复差场包络',color='tab:red',lw=.9)
                ax.axvline(m['stations'][list(xs).index(x)]['guide_ns']['basal_sandstone'],color='tab:green',ls=':',label='底砂垂直参考')
                ax.set(title=f'剖面{x:g}m，同站配对 / '+('结构窗' if row==0 else '深部放大'),
                       xlabel='SFCW 时间 / ns',ylabel='幅度 / (V/m)/(A·m)')
                ax.legend(fontsize=7)
            common = max(ax.get_ylim()[1] for ax in axs[row])
            for ax in axs[row]: ax.set_ylim(0,common)
        fig.suptitle(f'已计算{len(anchors)}处对照站独立展示，各行共用幅度标尺 / {window}\n'
                    '上行120–550ns结构显示，下行250–550ns深部放大；数值无增益，H0包含地质响应，各站独立配对')
        fig.savefig(a.out/f'anchor_pairs_{window}.png',dpi=140); plt.close(fig)
    if a.cli:
        with h5py.File(a.cli) as cli, h5py.File(a.numerical) as h:
            q = all_q[170.,'H1']
            np.testing.assert_array_equal(cli['frequency'][:], FREQ)
            assert relative(cli['response'][:]/.025,q.response) < 1e-12
            cli_error = relative(cli['time_response/complex_bandpass'][:]/.025,
                                 h['hann/H1_complex_bandpass'][:,-1])
            assert cli_error < 1e-12
    else:
        cli_error = None
    result = dict(status='COMPLETE_41_STATION_LINE' if valid.all() else 'PARTIAL_COMPLETED_STATIONS_ONLY',
        script_sha256=sha(__file__),official_processing_sha256=sha(sf.__file__),contract_sha256=sha(a.raw/'execution_contract.json'),
        retrieved_files_checked=len(retrieval),completed_new_runs=len(receipts),reused_runs=2,H1_station_count=count,
        matched_anchor_count=len(anchors),native_audits=rows,metrics=metrics,inverse_relative_L2=inverse_errors,
        official_CLI_complex_profile_relative_L2=cli_error,numerical_sha256=sha(a.numerical),
        anchor_display_window_ns=[120,550],
        method='official_direct',SFCW=m['SFCW'],limits=m['limits'])
    save(a.out/'analysis.json',result)
    for record in retrieval:
        src = a.raw/record['file']
        if src.suffix == '.h5': continue
        dst = a.out/record['file']; dst.parent.mkdir(parents=True,exist_ok=True)
        assert not dst.exists(); shutil.copyfile(src,dst)
    shutil.copyfile(a.raw/'retrieval_manifest.json',a.out/'retrieval_manifest.json')
    print(json.dumps(dict(status=result['status'],completed_new_runs=len(receipts),H1_stations=count,
                         matched_anchors=len(anchors),inverse_relative_L2=inverse_errors,metrics=metrics),ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['raw','out','numerical']:
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--cli',type=Path)
    p.add_argument('--allow-partial',action='store_true')
    main(p.parse_args())
