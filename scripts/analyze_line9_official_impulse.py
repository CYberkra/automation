"""Official direct SFCW and causal-convolution audit of the two impulse cases."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.signal import fftconvolve
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import audit_source

FREQ = 20e6 + np.arange(501)*300000.


def relative(a, b):
    return float(np.linalg.norm(a-b)/np.linalg.norm(b))


def response(p):
    s = sf.load_source(p)
    r = sf.load_receiver(p, '/rxs/rx1', 'Ez')
    q = sf.direct_frequency_response(s, r, FREQ, tail_taper_fraction=0)
    assert q.source_valid.all() and s.spatial_scale == .025
    exact = np.empty(501, complex)
    for k in range(0, 501, 16):
        f = FREQ[k:k+16, None]
        exact[k:k+16] = (r.dt*(np.exp(-2j*np.pi*f*r.times)@r.samples))/(
            s.dt*(np.exp(-2j*np.pi*f*s.times)@s.samples))/.025
    q = replace(q, response=q.response/.025)
    error = relative(q.response, exact)
    assert error < 1e-9
    return q, error


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    c = json.loads((a.new/'execution_contract.json').read_text('utf-8'))
    v = json.loads((a.new/'completed_verification.json').read_text('utf-8'))
    oldc = json.loads((a.old/'execution_contract.json').read_text('utf-8'))
    oldv = json.loads((a.old/'snapshot_manifest.json').read_text('utf-8'))
    assert v['completed'] and v['contract_sha256'] == sha(a.new/'execution_contract.json')
    assert oldv['contract_sha256'] == sha(a.old/'execution_contract.json')
    old_groups = {g['id']: g for g in oldc['groups']}
    old_records = {r['id']: r for r in oldv['records']}
    rows = []; qs = [[], []]; raw = []
    for g, record in zip(c['groups'], v['groups']):
        name = g['id']; oldg = old_groups[name]
        for key in ['geometry_sha256', 'material_sha256', 'tx_m', 'rx_m', 'native_shape', 'dt_s', 'expected_samples']:
            assert g[key] == oldg[key], (name, key)
        pn = a.new/(name+'.h5'); po = a.old/(name+'.h5')
        assert sha(pn) == record['native_sha256'] and sha(po) == old_records[name]['native_sha256']
        with h5py.File(pn) as n, h5py.File(po) as o:
            for key in ['gprMax', 'nx_ny_nz', 'dx_dy_dz', 'dt', 'Iterations']:
                np.testing.assert_array_equal(n.attrs[key], o.attrs[key])
            assert str(n.attrs['gprMax']) == '4.0.1'
            xn = n['rxs/rx1/Ez'][:]; xo = o['rxs/rx1/Ez'][:]
            sn = n['srcs/src1/excitation/samples'][:]; so = o['srcs/src1/excitation/samples'][:]
            assert xn.dtype == xo.dtype == sn.dtype == np.float64 and xn.shape == (20352,)
            assert np.isfinite(xn).all()
            audit_source(n['srcs/src1/excitation'], g, float(n.attrs['dt']))
            audit_source(o['srcs/src1/excitation'], oldg, float(o.attrs['dt']))
            for key in ['srcs/src1', 'rxs/rx1']:
                np.testing.assert_array_equal(n[key].attrs['Position'], o[key].attrs['Position'])
            # Full causal convolution BEFORE cropping to the Ricker observation window.
            full_convolution = fftconvolve(xn/sn[0], so, mode='full')
            reconstructed = full_convolution[:len(xo)]
            lti = relative(reconstructed, xo)
            peak_error = float(abs(reconstructed-xo).max()/abs(xo).max())
            assert lti < 1e-8 and peak_error < 1e-8
            raw.append((xo, reconstructed, float(n.attrs['dt'])))
        errors = []; tails = []
        for j, p in enumerate([po, pn]):
            q, error = response(p); qs[j].append(q)
            errors.append(error); tails.append(q.receiver_tail_relative_db)
        # Identify the frequency discrepancy due to cropping that SAME known
        # causal convolution. This does not recover unknown impulse arrivals
        # after the simulated time window or certify that window as sufficient.
        source = sf.load_source(po)
        omitted = np.empty(501, complex)
        times = np.arange(len(xo), len(full_convolution))*source.dt
        for k in range(0, 501, 16):
            f = FREQ[k:k+16, None]
            omitted[k:k+16] = source.dt*(np.exp(-2j*np.pi*f*times)@full_convolution[len(xo):])
        omitted /= qs[0][-1].source_spectrum*.025
        crop_identity_error = float(np.linalg.norm(qs[1][-1].response-qs[0][-1].response-omitted)/np.linalg.norm(qs[1][-1].response))
        assert crop_identity_error < 1e-8
        rows.append(dict(id=name, impulse_native_sha256=sha(pn), ricker_native_sha256=sha(po),
            geometry_sha256=g['geometry_sha256'], material_sha256=g['material_sha256'],
            causal_convolution_relative_L2=lti, causal_convolution_peak_relative_error=peak_error,
            independent_DFT_relative_L2=errors, receiver_final_5pct_peak_relative_db=tails,
            final_5pct_below_official_minus60dB=[bool(x < -60) for x in tails],
            known_convolution_crop_identity_relative_error=crop_identity_error,
            full_501tone_relative_source_change=relative(qs[1][-1].response,qs[0][-1].response)))
    a.out.mkdir(parents=True)
    metrics = {}; profiles = {}
    for window in ['hann', 'blackman']:
        profiles[window] = []
        for j in range(2):
            z = np.column_stack([q.response for q in qs[j]])
            z = np.column_stack([z, z[:, 1]-z[:, 0]])
            profile = sf.reconstruct_time_response(replace(qs[j][0], response=z), window=window, zero_pad_factor=8, time_shift=0)
            profiles[window].append(profile)
            k = np.arange(7, len(profile.time), 101)
            independent = np.exp(2j*np.pi*profile.time[k, None]*FREQ)@(z*profile.weights[:, None])/profile.weights.sum()
            assert relative(profile.complex_bandpass[k], independent) < 1e-9
        old, new = profiles[window]; t = new.time*1e9; gates = {}
        for gate, bounds in [('early',[0,120]), ('deep',[300,450]), ('basal', c['groups'][0]['basal_gate_ns'])]:
            mask = (t >= bounds[0]) & (t <= bounds[1])
            gates[gate] = dict(bounds_ns=bounds, configurations={})
            for j, name in enumerate(['H0','H1','H1_minus_H0']):
                x = old.complex_bandpass[mask, j]; y = new.complex_bandpass[mask, j]
                gates[gate]['configurations'][name] = dict(relative_complex_change=relative(y,x),
                    old_peak_ns=float(t[mask][np.argmax(abs(x))]), new_peak_ns=float(t[mask][np.argmax(abs(y))]),
                    new_over_old_peak=float(abs(y).max()/abs(x).max()))
        metrics[window] = gates
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    for window, (old, new) in profiles.items():
        t = new.time*1e9
        fig, axs = plt.subplots(2,2,figsize=(13,9),layout='constrained')
        for ax, bounds, title in [(axs[0,0],[0,450],'总场：早段与深部'), (axs[0,1],[300,450],'总场：深部时间窗')]:
            mask = (t>=bounds[0]) & (t<=bounds[1])
            data = np.column_stack([old.complex_bandpass[mask,:2], new.complex_bandpass[mask,:2]]).real
            lim = float(abs(data).max())
            im = ax.imshow(data,cmap='gray',aspect='auto',vmin=-lim,vmax=lim,extent=[-.5,3.5,t[mask][-1],t[mask][0]])
            ax.set_xticks(range(4),['旧Ricker H0','旧Ricker H1','新impulse H0','新impulse H1'],rotation=12,fontsize=8)
            ax.set(title=title+'；四列共用绝对灰度',ylabel='SFCW 时间 / ns')
            fig.colorbar(im,ax=ax,label='Re(复带通)，(V/m)/(A·m)')
        for j, label in [(0,'H0：无连通底砂'),(1,'H1：含连通底砂'),(2,'底砂配对差场 H1−H0')]:
            ax = axs[1,0] if j<2 else axs[1,1]
            mask = (t>=300)&(t<=450)
            for p, kind, ls in [(old,'旧Ricker','-'),(new,'新impulse','--')]:
                ax.plot(t[mask],p.complex_bandpass[mask,j].real,ls,label=kind+' '+label,lw=1)
            ax.set(xlabel='SFCW 时间 / ns',ylabel='Re(复带通)，(V/m)/(A·m)')
        axs[1,0].set_title('总场波形：源归一化后直接比较，无拟合')
        axs[1,1].set_title('配对复差场单独标尺；不代表全图干净真值')
        for ax in axs[1]:
            ax.axvspan(*c['groups'][0]['basal_gate_ns'],color='gray',alpha=.12,label='既有底砂时间窗')
            ax.legend(fontsize=7)
        fig.suptitle(f'190 m同站位 / 约8 m离地 / gprMax4.0.1 FP64 / 官方direct SFCW\n20–170 MHz，0.3 MHz，501点；{window}，补零8倍，无尾窗/时移/AGC；配置列不是连续空间B-scan')
        fig.savefig(a.out/f'official_impulse_{window}_gray.png',dpi=140);plt.close(fig)
    fig, axs = plt.subplots(2,1,figsize=(12,7),layout='constrained')
    for ax, (xo, recon, dt), g in zip(axs,raw,c['groups']):
        t=np.arange(len(xo))*dt*1e9; mask=t<=450
        ax.plot(t[mask],xo[mask],label='直接Ricker正演，40 A',lw=1)
        ax.plot(t[mask],recon[mask],'--',label='新impulse响应与实际Ricker源做因果卷积',lw=1)
        ax.set(title=g['id']+' 原生时域线性校核',xlabel='FDTD 时间 / ns',ylabel='Ez / V/m');ax.legend(fontsize=8)
    fig.suptitle('完整卷积后裁到相同20352样本；未移动、拟合或归一化接收波形')
    fig.savefig(a.out/'official_impulse_convolution.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ)
        for j,name in enumerate(['ricker','impulse']):
            z=np.column_stack([q.response for q in qs[j]])
            h.create_dataset(name+'/response',data=z)
            for w,pp in profiles.items():
                h.create_dataset(name+'/'+w+'/complex_bandpass',data=pp[j].complex_bandpass)
                h.create_dataset(name+'/'+w+'/real_bandpass',data=pp[j].real_bandpass)
        h.create_dataset('time_s',data=profiles['hann'][0].time)
    result=dict(status='PASS_NATIVE_OFFICIAL_SFCW_AND_CAUSAL_CONVOLUTION_NOT_FIELD_VALIDATION',
        script_sha256=sha(__file__), execution_contract_sha256=sha(a.new/'execution_contract.json'),
        completed_verification_sha256=sha(a.new/'completed_verification.json'),
        old_snapshot_manifest_sha256=sha(a.old/'snapshot_manifest.json'),rows=rows,metrics=metrics,
        numerical_sha256=sha(a.numerical),method='official_direct',frequency_Hz=[20e6,170e6],frequency_step_Hz=300000,
        tones=501,normalization='actual source spectrum and 0.025 m current moment',
        display='Re(complex_bandpass); official real_bandpass is twice this and saved separately',
        tail_taper=False,time_shift=False,AGC=False,fit=False,
        limits='Two configurations at one station only, not a spatial scan. Impulse end-of-record tail is not fully decayed by the official -60 dB advisory; crop algebra only explains the known finite convolution, not unknown arrivals after 1200 ns. No unique-path, mesh-convergence, finite-3D antenna, field permittivity or field-performance certification.')
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'rows':rows,'metrics':metrics},ensure_ascii=False))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['new','old','out','numerical']:parser.add_argument('--'+name,type=Path,required=True)
    main(parser.parse_args())
