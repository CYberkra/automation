"""Audit an immutable observer solve and analyze collocated native TM total fields.

These are point observations of one excitation, not a moving-antenna B-scan.
Flux signs and waveform correlations do not certify a unique ray or bounce order.
"""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha

C = 299792458.
BANDS = ['air_surface','cover_upper','cover_middle','cover_lower','mud_below','air_34','air_38']
LABELS = ['地表上方空气','覆盖层上部','覆盖层中部','覆盖层下部','底界面下方泥岩','空气 y=34.025m','空气 y=38.025m']


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def collocate(e, hx, hx_below, hy, hy_left):
    """Native H is half a step behind E; neighboring H brackets E spatially."""
    arrays = [np.asarray(x) for x in [e,hx,hx_below,hy,hy_left]]
    if len({x.shape for x in arrays})!=1 or arrays[0].ndim!=1 or arrays[0].size<2:
        raise ValueError('Equal one-dimensional native histories with at least two samples required')
    if not all(np.isfinite(x).all() for x in arrays):
        raise ValueError('Nonfinite field history')
    hh=.5*(arrays[1]+arrays[2]);hv=.5*(arrays[3]+arrays[4])
    return arrays[0][:-1], .5*(hh[:-1]+hh[1:]), .5*(hv[:-1]+hv[1:])


def main(a):
    from audit_line9_native_probe_inputs import check
    from analyze_line9_v401_version_controls import response, inverse, FREQ
    from gprMax.toolboxes.SFCW import processing as sf
    assert not a.out.exists() and not a.arrays.exists()
    check(a.package)
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    cert=json.loads((a.source/'completed_verification.json').read_text('utf-8'))
    assert cert['completed'] and len(cert['groups'])==1
    assert sha(a.source/'execution_contract.json')==cert['contract_sha256']=='8b538b98c839ccc0b74af7a00dc89b5e372fec48ce1439a9418dd0fc11612820'
    raw=a.source/'profile.h5';baseline=a.package/'baseline/native_H0.h5'
    assert sha(raw)==cert['groups'][0]['native_sha256']
    assert sha(baseline)==m['baseline_sha256']['native_H0.h5']
    ee=[];xx=[];yy=[]
    with h5py.File(raw) as h,h5py.File(baseline) as old:
        assert str(h.attrs['gprMax'])=='4.0.1' and h.attrs['nrx']==862 and len(h['rxs'])==862
        dt=float(h.attrs['dt']);assert dt==m['dt_s']
        for key in ['dt','Iterations','dx_dy_dz','nx_ny_nz']:
            np.testing.assert_array_equal(h.attrs[key],old.attrs[key])
        for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
            assert h[key][:].tobytes()==old[key][:].tobytes()
        np.testing.assert_array_equal(h['rxs/rx1'].attrs['Position'],old['rxs/rx1'].attrs['Position'])
        rx=h['rxs/rx1/Ez'][:-1]
        for p in m['probes']:
            rr=[]
            for anchor in p['anchors']:
                r=h[f"rxs/rx{anchor['receiver_index']}"]
                assert r.attrs['Name']==anchor['name'] and set(r.keys())==set(anchor['outputs'])
                np.testing.assert_array_equal(r.attrs['GridPosition'],anchor['coord'])
                np.testing.assert_allclose(r.attrs['Position'],anchor['position_m'],atol=1e-12,rtol=0)
                row={}
                for field in anchor['outputs']:
                    d=r[field];v=d[:]
                    assert v.dtype==np.float64 and v.shape==(20352,) and np.isfinite(v).all()
                    assert d.attrs['Quantity']==field and d.attrs['SampleInterval']==dt
                    assert d.attrs['TimeSampleOffset']==(0 if field=='Ez' else -.5*dt)
                    row[field]=v
                rr.append(row)
            e,hx,hy=collocate(rr[0]['Ez'],rr[0]['Hx'],rr[1]['Hx'],rr[0]['Hy'],rr[2]['Hy'])
            ee.append(e);xx.append(hx);yy.append(hy)
    e=np.array(ee);hx=np.array(xx);hy=np.array(yy);sx=-e*hy;sy=e*hx
    t=np.arange(e.shape[1])*dt;ns=t*1e9;index={p['id']:i for i,p in enumerate(m['probes'])}
    gate_rows=[]
    for k,p in enumerate(m['probes']):
        gates={}
        for name,bounds in m['gates_native_ns'].items():
            mask=(ns>=bounds[0])&(ns<=bounds[1]);ix=np.flatnonzero(mask)
            integ=[float(np.sum(sx[k,mask])*dt),float(np.sum(sy[k,mask])*dt)]
            pos=ix[np.argmax(sy[k,mask])];neg=ix[np.argmin(sy[k,mask])]
            gates[name]={'integrated_Sx_J_m2':integ[0],'integrated_Sy_J_m2':integ[1],
                         'direction_deg':float(np.degrees(np.arctan2(integ[1],integ[0]))),
                         'max_Sy_W_m2':float(sy[k,pos]),'max_Sy_time_ns':float(ns[pos]),
                         'min_Sy_W_m2':float(sy[k,neg]),'min_Sy_time_ns':float(ns[neg])}
        gate_rows.append({'id':p['id'],'band':p['band'],'x_m':p['x_m'],'y_m':p['y_m'],'gates':gates})
    correlations=[];mask=(ns>=300)&(ns<=370);b=rx[mask]
    for p in m['probes']:
        if p['band']!='air_surface':continue
        k=index[p['id']];delay=np.hypot(170.65-p['x_m'],39.025-p['y_m'])/C
        z=np.interp(t[mask]-delay,t,e[k])
        corr=float(z@b/np.linalg.norm(z)/np.linalg.norm(b))
        correlations.append({'id':p['id'],'x_m':p['x_m'],'delay_ns':float(delay*1e9),'signed_cosine':corr})
    # Explicit post-observation exploratory gates: describe extrema, not certify pulses.
    windows=[('first_down',[50,110],'min'),('first_up',[120,175],'max'),('next_down',[190,240],'min'),('late_up_candidate',[250,290],'max')]
    events=[]
    for p in m['probes']:
        if p['band']!='cover_middle':continue
        k=index[p['id']];q={}
        for name,(lo,hi),mode in windows:
            ix=np.flatnonzero((ns>=lo)&(ns<=hi));j=ix[getattr(np,'arg'+mode)(sy[k,ix])]
            near=(ns>=ns[j]-4)&(ns<=ns[j]+4)
            q[name]={'time_ns':float(ns[j]),'Sy_W_m2':float(sy[k,j]),'local8ns_signed_J_m2':float(np.sum(sy[k,near])*dt)}
        events.append({'id':p['id'],'x_m':p['x_m'],'events':q})
    z0,d0=response(baseline,.025);z1,d1=response(raw,.025)
    assert np.array_equal(z0,z1)
    profiles={};inverse_errors={}
    for window in ['hann','blackman']:
        w=sf.spectral_window(window,501);z,st=inverse(z1[:,None],w)
        pick=np.arange(7,len(st),101)
        direct=np.exp(2j*np.pi*st[pick,None]*FREQ)@(z1*w)/501
        err=float(np.linalg.norm(z[pick,0]-direct)/np.linalg.norm(direct));assert err<1e-9
        inverse_errors[window]=err;profiles[window]=z[:,0]
    a.out.mkdir(parents=True);a.arrays.parent.mkdir(parents=True,exist_ok=True)
    np.savez(a.arrays,time_s=t,Ez=e,Hx=hx,Hy=hy,Sx=sx,Sy=sy,rx=rx,sfcw_time_s=st,**profiles)
    summary={'status':'PASS_NATIVE_OBSERVER_AND_TOTAL_FIELD_DIAGNOSTICS','script_sha256':sha(__file__),
             'contract_sha256':cert['contract_sha256'],'manifest_sha256':sha(a.package/'manifest.json'),
             'native_sha256':sha(raw),'baseline_sha256':sha(baseline),'native_source_rx_bitwise_equal':True,
             'raw_receiver_count':862,'collocated_points':287,'collocated_samples':20351,'dt_s':dt,
             'independent_DFT_relative_L2':[d0,d1],'inverse_errors':inverse_errors,'SFCW_complex_response_equal':True,
             'gates':gate_rows,'surface_delay_correlations':correlations,'exploratory_cover_middle_extrema':events,
             'exploratory_windows':windows,'arrays_sha256':sha(a.arrays),
             'limits':'Native total E/H flux; no separated-ray or bounce count, no flux conservation certificate; correlations are 41 searched positions with geometric delays and no fitted lag/amplitude; not site/material/finite3D/whole-line certification.'}
    save(a.out/'analysis.json',summary)
    plot(a.out,m,index,ns,e,sy,rx,correlations,events,st,profiles)
    print(json.dumps({'status':summary['status'],'best_surface_correlation':max(correlations,key=lambda r:r['signed_cosine']),
                      'SFCW_complex_response_equal':True,'events_at160m':next(r for r in events if r['x_m']==160)},ensure_ascii=False))


def plot(out,m,index,ns,e,sy,rx,correlations,events,st,profiles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    def image(ax,band,lo,hi,lim):
        kk=[index[p['id']] for p in m['probes'] if p['band']==band];mask=(ns>=lo)&(ns<=hi)
        data=sy[kk][:,mask].T
        im=ax.imshow(data,cmap='RdBu_r',vmin=-lim,vmax=lim,aspect='auto',extent=[153.75,174.25,hi,lo])
        ax.set(xlabel='模型 x / m（剖面里程=x+20m）',ylabel='原生时间 / ns')
        clipped=np.mean(abs(data)>lim)*100
        ax.set_title(LABELS[BANDS.index(band)]+f'；超过色标截色 {clipped:.2f}%')
        return im
    fig,axs=plt.subplots(3,1,figsize=(12,12),layout='constrained')
    for ax,band,lo,hi,lim in zip(axs,['cover_middle','cover_middle','air_surface'],[50,240,260],[240,330,350],[1.,1e-4,1e-5]):
        im=image(ax,band,lo,hi,lim);fig.colorbar(im,ax=ax,label='总场 Sy / (W/m²)；红向上，蓝向下')
    fig.suptitle('190m高损耗H0：覆盖层内明确下行→上行→再次下行；晚返回另用弱场尺度\n固定一次激发的41列探针时空图；不是移动天线B-scan，不能按总场符号数反射阶次')
    fig.savefig(out/'native_cover_transport.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(4,2,figsize=(14,14),layout='constrained')
    for ax,band in zip(axs.flat,BANDS):
        im=image(ax,band,230,380,3e-5);fig.colorbar(im,ax=ax,label='总场 Sy / (W/m²)')
    axs.flat[-1].axis('off');axs.flat[-1].text(.03,.75,'七个面板统一 ±3e-5 W/m²\n红：向上；蓝：向下\n覆盖层残留强波有截色，见各标题\n同位四邻点空间/原生相邻时间步对齐\n未经分波：局部正峰不等于独立多次波',fontsize=12)
    fig.suptitle('地表两侧、覆盖层内及空气的晚波：相同物理色标，保留符号\n原始总场；无AGC、逐道归一化、去背景或时延拟合')
    fig.savefig(out/'native_seven_bands_late.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    for ax,x in zip(axs.flat,[158,160,164,170]):
        p=next(p for p in m['probes'] if p['band']=='cover_middle' and p['x_m']==x);k=index[p['id']]
        mask=(ns>=115)&(ns<=240);ax.plot(ns[mask],sy[k,mask],lw=1,label='覆盖层中部总场 Sy')
        ax.axhline(0,color='gray',lw=.7);ax.set(xlabel='原生时间 / ns',ylabel='Sy / (W/m²)',title=f'模型 x={x}m / y={p["y_m"]:.3f}m');ax.set_ylim(-.1,.8)
    fig.suptitle('四个示例点的首次上行与再次下行；四图同尺度\n位置是事后探索示例，完整41列和冻结宽窗数字另存；不能把每个峰数作一次反射')
    fig.savefig(out/'native_midcover_signed_examples.png',dpi=145);plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(13,10),layout='constrained')
    axs[0,0].plot([r['x_m'] for r in correlations],[r['signed_cosine'] for r in correlations],'.-')
    axs[0,0].axvspan(162.14,162.72,color='gray',alpha=.2,label='此前方向外推候选范围')
    axs[0,0].set(xlabel='地表空气探针模型 x / m',ylabel='带符号余弦相关',title='完整41位置：仅补几何空气传播时延');axs[0,0].legend(fontsize=9)
    mask=(ns>=300)&(ns<=370);axs[0,1].plot(ns[mask],rx[mask],color='black',lw=1,label='主Rx 原幅度')
    ax2=axs[0,1].twinx()
    for x in [160,162,164]:
        row=next(r for r in correlations if r['x_m']==x);k=index[row['id']]
        ax2.plot(ns[mask],np.interp(ns[mask]-row['delay_ns'],ns,e[k]),lw=.9,label=f'地表x={x}m（几何延时）')
    axs[0,1].set(xlabel='主Rx原生时间 / ns',ylabel='Rx Ez / (V/m)',title='晚波波形比较：两个坐标轴，未拟合振幅');ax2.set_ylabel('地表探针 Ez / (V/m)');ax2.legend(fontsize=8)
    smask=(st*1e9>=150)&(st*1e9<=450);data=np.column_stack([profiles['hann'][smask].real]*2+[np.zeros(smask.sum())]);lim=float(abs(data).max())
    im=axs[1,0].imshow(data,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[-.5,2.5,450,150])
    axs[1,0].set_xticks(range(3),['原H0','探针H0','新−原']);axs[1,0].set(ylabel='SFCW时间 / ns',title='同站配置列：20–170MHz精确501点 Hann');fig.colorbar(im,ax=axs[1,0],label='(V/m)/(A·m)')
    for window,z in profiles.items():axs[1,1].plot(st[smask]*1e9,z[smask].real,lw=1,label=window)
    axs[1,1].set(xlabel='SFCW时间 / ns',ylabel='双极 / (V/m)/(A·m)',title='新增观察器没有改变原有SFCW响应');axs[1,1].legend()
    fig.suptitle('候选非局部地表返回与主Rx的联系；相关平台不定位唯一出射点\n4.0.1 原生FP64 / 190m高损耗H0（无底部连通砂岩）；单站配置列不是连续空间B-scan')
    fig.savefig(out/'native_surface_rx_and_sfcw.png',dpi=145);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','package','out','arrays']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
