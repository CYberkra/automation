"""Read-only native review: source/observable, target-specific windows and truncation.

No solver, no field fitting. Historical cross-machine contracts are preserved;
local files are resolved by capsule/group ID and checked against archived hashes.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256
from green_halfspace_hed_v0_1 import C_SI,epsr_cover

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/research_checks'
F=np.linspace(20e6,170e6,501)
DEPTHS={'depth3':3.,'depth10':10.,'depth18_5':18.5}
PROCESSING_SHA='adad556f09140956f0ee19d3038430e06a9ae8be6a826dcad723096d99624a3b'


def relative(a,b):
    denominator=np.linalg.norm(b)
    if denominator==0:
        raise ValueError('undefined relative norm')
    return float(np.linalg.norm(a-b)/denominator)


def response(source,receiver,tail_ns):
    fraction=max(0,round(tail_ns*1e-9/receiver.dt)-.25)/len(receiver.samples)
    product=sf.direct_frequency_response(source,receiver,F,tail_taper_fraction=fraction)
    if not product.source_valid.all():
        raise ValueError('invalid source frequency')
    # Length is from the native polarized current element, not x/z mesh spacing.
    return replace(product,response=product.response/source.spatial_scale)


def load_capsule(name):
    study=BASE/name
    cp=study/'execution_contract.json'
    c=json.loads(cp.read_text('utf-8'))
    audit=json.loads((study/'completed_verification.json').read_text('utf-8'))
    if audit['status']!='PASS' or audit['contract_sha256']!=sha256(cp):
        raise ValueError('completed contract identity required')
    verified={r['id']:r['raw_sha256'] for r in audit['groups']}
    records={};identities=[]
    for g in c['groups']:
        raw=study/g['id']/'profile.h5';inp=raw.with_suffix('.in')
        if sha256(raw)!=verified[g['id']] or sha256(inp)!=g['input_sha256']:
            raise ValueError('archived input/native identity mismatch')
        pol=g.get('polarisation','y')
        source=sf.load_source(raw);receiver=sf.load_receiver(raw,component='E'+pol)
        with h5py.File(raw) as h:
            if str(h.attrs['gprMax'])!='4.0.0' or receiver.samples.dtype!=np.float64:
                raise ValueError('native V4 double required')
            if not np.array_equal(h['srcs/src1'].attrs['GridPosition'],
                   np.rint(np.array(g['tx_m'])/h.attrs['dx_dy_dz']).astype(int)):
                raise ValueError('native source position mismatch')
            receiver_position=np.array(g['receivers_m'][0])
            if not np.array_equal(h['rxs/rx1'].attrs['GridPosition'],
                   np.rint(receiver_position/h.attrs['dx_dy_dz']).astype(int)):
                raise ValueError('native receiver position mismatch')
            native_spacing=np.array(h.attrs['dx_dy_dz'])
        if (not np.isfinite(receiver.samples).all() or not np.isfinite(source.samples).all()
                or source.spatial_scale!=native_spacing['xyz'.index(pol)]):
            raise ValueError('native source length/finite values mismatch')
        product=response(source,receiver,20)
        take=np.array([0,83,167,250,333,417,500])
        n=round(20e-9/receiver.dt)
        y=receiver.samples.copy();y[-n:]*=.5*(1+np.cos(np.linspace(0,np.pi,n)))
        xf=np.exp(-2j*np.pi*F[take,None]*source.times[None,:])@source.samples
        yf=np.exp(-2j*np.pi*F[take,None]*receiver.times[None,:])@y
        difference=relative(product.response[take],yf/xf/source.spatial_scale)
        if difference>1e-9:
            raise ValueError('independent actual-clock DFT mismatch')
        records[g['id']]=dict(group=g,source=source,receiver=receiver,product=product)
        identities.append(dict(id=g['id'],input_sha256=sha256(inp),raw_sha256=sha256(raw),
             independent_DFT_relative_L2=difference,source_length_m=source.spatial_scale,
             historical_transfer_divisor_m=g['spacing_m'],
             historical_response_over_native_length_response=source.spatial_scale/g['spacing_m'],
             source_abs_peak_time_ns=float(source.times[np.argmax(abs(source.samples))]*1e9)))
    return c,records,dict(capsule=name,contract_sha256=sha256(cp),records=identities,
                        path_resolution='Local capsule/group/profile files; archived absolute paths untouched')


def raw_tail(record,background,depth):
    r=record['receiver'];b=background['receiver'];v=r.samples-b.samples
    peak=np.max(abs(v));time=r.times*1e9
    mask=time>=time[-1]*.9
    lower=(16+2*depth*np.sqrt(18.017))/C_SI*1e9
    early=time<lower-20
    return dict(native_end_ns=float(time[-1]),front_path_lower_bound_ns=float(lower),
        last_sample_over_isolated_peak=float(abs(v[-1])/peak),
        last_tenth_rms_over_isolated_peak=float(np.sqrt(np.mean(v[mask]**2))/peak),
        recorded_isolated_energy_fraction_in_last_tenth=float(np.sum(v[mask]**2)/np.sum(v**2)),
        pre_front_minus20ns_peak_over_isolated_peak=float(np.max(abs(v[early]))/peak),
        no_claim_of_complete_tail=True)


def main(out):
    if out.exists() or sha256(sf.__file__)!=PROCESSING_SHA:
        raise ValueError('new output and reviewed official processing required')
    out.mkdir(parents=True)
    capsules={};identities=[]
    for key,name in [('layout','2026-10-05_layout_pol_3d_r1'),
                     ('depth','2026-10-05_depth_factor_2d_r1'),
                     ('bscan','2026-10-05_depth_bscan_2d_r1'),
                     ('halfspace','2026-10-05_halfspace_point_3d_r2')]:
        c,r,identity=load_capsule(name);capsules[key]=(c,r);identities.append(identity)
        print('Verified native capsule',name,len(r),flush=True)
    # Layout: compare target contrasts as well as total fields. Ricker source
    # centre delay is measured from native excitation, not fitted to the echo.
    _,layout=capsules['layout'];layout_rows=[];window_rows=[]
    for scenario in ('flat','rough'):
        for pol in ('x','y'):
            ia,ca=f'{scenario}_inline_{pol}',f'{scenario}_cross_{pol}'
            ib,cb=f'fullcover_inline_{pol}',f'fullcover_cross_{pol}'
            t=layout[ia]['receiver'].times*1e9
            wa=layout[ia]['receiver'].samples;wc=layout[ca]['receiver'].samples
            a=wa-layout[ib]['receiver'].samples;b=wc-layout[cb]['receiver'].samples
            m=t>=46;mi=(t>=120)&(t<200)
            sa=layout[ia]['product'].response-layout[ib]['product'].response
            sb=layout[ca]['product'].response-layout[cb]['product'].response
            dena=np.linalg.norm(layout[ia]['product'].response)
            denb=np.linalg.norm(layout[ca]['product'].response)
            layout_rows.append(dict(scenario=scenario,polarisation=pol,
                total_raw_post46_relative_L2=relative(wa[m],wc[m]),
                target_raw_120_200_relative_L2=relative(a[mi],b[mi]),
                target_complex_band_relative_L2=relative(sa,sb),
                target_band_amplitude_inline_over_cross=float(np.linalg.norm(sa)/np.linalg.norm(sb)),
                target_over_total_band_inline=float(np.linalg.norm(sa)/dena),
                target_over_total_band_cross=float(np.linalg.norm(sb)/denb)))
    for orientation in ('inline','cross'):
        x=layout['fullcover_'+orientation+'_x'];y=layout['fullcover_'+orientation+'_y']
        t=x['receiver'].times*1e9
        source_delay=float(x['source'].times[np.argmax(abs(x['source'].samples))]*1e9)
        for name,start,end in [('direct',3,9),('surface',46,62)]:
            old=(t>=start)&(t<end);shift=(t>=start+source_delay)&(t<end+source_delay)
            xv=x['receiver'].samples;yv=y['receiver'].samples
            window_rows.append(dict(layout=orientation,event=name,
                old_window_ns=[start,end],source_peak_delay_ns=source_delay,
                source_delay_included_window_ns=[start+source_delay,end+source_delay],
                raw_absolute_peak_x_over_y_old=float(np.max(abs(xv[old]))/np.max(abs(yv[old]))),
                raw_absolute_peak_x_over_y_shifted=float(np.max(abs(xv[shift]))/np.max(abs(yv[shift])))))
    # Depth: matched complex contrasts, tail-window sensitivity and corrected
    # propagation phase. Approximate cylindrical spreading is still a hypothesis.
    _,depth=capsules['depth'];depth_rows=[];k=2*np.pi*F/C_SI*np.sqrt(epsr_cover(F))
    iso={case:depth[case]['product'].response-depth['fullcover']['product'].response for case in DEPTHS}
    no_taper={name:response(r['source'],r['receiver'],0) for name,r in depth.items()}
    for case,d in DEPTHS.items():
        contrast0=no_taper[case].response-no_taper['fullcover'].response
        delta=d-3
        n_eff=np.sqrt(epsr_cover(F).real)
        spread=np.sqrt((8/n_eff+6)/(8/n_eff+2*d))
        old=iso['depth3']*spread*np.exp(2*k.imag*delta)
        corrected=iso['depth3']*spread*np.exp(-2j*k*delta)
        tails={}
        for tail in (0,10,20,50):
            a=response(depth[case]['source'],depth[case]['receiver'],tail)
            b=response(depth['fullcover']['source'],depth['fullcover']['receiver'],tail)
            value=a.response-b.response
            tails[str(tail)]=dict(complex_relative_L2_to_no_taper=relative(value,contrast0),
                                  band_L2_amplitude_ratio_to_no_taper=float(np.linalg.norm(value)/np.linalg.norm(contrast0)))
        depth_rows.append(dict(case=case,depth_m=d,tail=raw_tail(depth[case],depth['fullcover'],d),
            tail_taper_sensitivity_ns=tails,
            old_attenuation_only_complex_extrapolation_error=relative(old,iso[case]),
            propagation_phase_included_complex_extrapolation_error=relative(corrected,iso[case]),
            propagation_phase_model='exp(-2j*k(f)*(depth-3)); retains the original approximate spreading, not an exact layered Green function'))
    # 52-trace B-scan: retain both raw Ricker and genuine source-normalized
    # SFCW reconstructions. Every branch uses its own declared shared background
    # early-peak reference because their source/units differ.
    c,bscan=capsules['bscan'];stations=c['geometry']['stations_tx_x_m'];count=len(stations)
    cases=('fullcover',*DEPTHS)
    matrices={name:[] for name in ('raw','hann','rectangular')};products={};native=None;time_sf=None
    for case in cases:
        raw=[];han=[];rect=[]
        for index in range(count):
            rec=bscan[f'{case}_s{index:02d}'];native=rec['receiver'].times*1e9
            raw.append(rec['receiver'].samples)
            tr=sf.reconstruct_time_response(rec['product'],window='hann',zero_pad_factor=8)
            rr=sf.reconstruct_time_response(rec['product'],window='rectangular',zero_pad_factor=8)
            time_sf=tr.time*1e9;han.append(tr.real_bandpass);rect.append(rr.real_bandpass)
        products['raw',case]=np.stack(raw,axis=1)
        products['hann',case]=np.stack(han,axis=1)
        products['rectangular',case]=np.stack(rect,axis=1)
    axes_time={'raw':native,'hann':time_sf,'rectangular':time_sf}
    bscan_rows=[];centre_match=[]
    for case,d in DEPTHS.items():
        for branch in axes_time:
            t=axes_time[branch];all_t=(t>=0)&(t<=650);early=(t>=0)&(t<100)
            target=products[branch,case]-products[branch,'fullcover']
            denominator=np.max(abs(products[branch,'fullcover'][early]),axis=0)
            numerator=np.max(abs(target[all_t]),axis=0)
            ratio=np.linalg.norm(target[all_t],axis=0)/np.linalg.norm(products[branch,'fullcover'][all_t],axis=0)
            bscan_rows.append(dict(case=case,branch=branch,
                target_peak_over_background_early_peak_median_dB=float(np.median(20*np.log10(numerator/denominator))),
                target_over_background_record_L2_median=float(np.median(ratio)),
                ratio_definition='Declared matched target contrast / background, not field clean/SNR',
                native_tail_unresolved=(case=='depth18_5')))
        centre_match.append(dict(case=case,centre_native_relative_L2_to_single_depth=relative(
            products['raw',case][:,count//2],depth[case]['receiver'].samples)))
    # Re-evaluate half-space gate on the isolated reflection, which is not part
    # of the prior PASS predicate. Existing independent-reference result is kept.
    old_half=BASE/'2026-10-05_halfspace_point_3d_r2/analysis/summary.json'
    hsummary=json.loads(old_half.read_text('utf-8'))
    cross=hsummary['reflection_isolation']
    half_gate=dict(total_field_metrics=hsummary['metrics'],reflection=cross,
         isolated_reflection_passes_same_5percent_5degree_limits=bool(
             cross['fdtd_isolated_reflection_vs_sommerfeld_relative_L2']<=.05 and
             cross['fdtd_isolated_reflection_vs_sommerfeld_max_phase_deg']<=5),
         old_analysis_sha256=sha256(old_half),
         scope='Existing reference result verified by native hashes; reflection gate was omitted from historical status predicate')
    # Scientific plots only; no field image/data exported.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei'];plt.rcParams['axes.unicode_minus']=False
    labels={'fullcover':'全覆盖背景','depth3':'平界面3m','depth10':'平界面10m','depth18_5':'平界面18.5m'}
    fig,axes=plt.subplots(3,4,figsize=(15,9),layout='constrained')
    for r,branch in enumerate(axes_time):
        t=axes_time[branch];mask=(t>=0)&(t<=650);early=(t>=0)&(t<100)
        scale=np.max(abs(products[branch,'fullcover'][early]))
        for col,case in enumerate(cases):
            ax=axes[r,col]
            image=ax.imshow(products[branch,case][mask]/scale,aspect='auto',cmap='RdBu_r',vmin=-.005,vmax=.005,
                extent=[stations[0]-.25,stations[-1]+.25,t[mask][-1],t[mask][0]])
            ax.set_title(labels[case]+': '+{'raw':'原始Ricker总场','hann':'SFCW Hann总场','rectangular':'SFCW矩形窗总场'}[branch],fontsize=9)
            if col==0:ax.set_ylabel('时间(ns)')
            if r==2:ax.set_xlabel('发射x(m)，接收=x+1.3m')
            fig.colorbar(image,ax=ax,shrink=.8,label='本行背景早峰的相对幅值')
    fig.suptitle('同一52道native：原始Ricker 与 20–170MHz/501点 SFCW\n每行共同背景早峰作固定参考；20ns尾渐消；窗是机制假设，非已确认设备链；18.5m尾端未闭合')
    fig.savefig(out/'depth_observation_bscans.png',dpi=135);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    for case,d in DEPTHS.items():
        r=depth[case]['receiver'];b=depth['fullcover']['receiver'];v=r.samples-b.samples
        axes[0].plot(r.times*1e9,v/np.max(abs(v)),label=labels[case])
    axes[0].axvspan(630,650,color='orange',alpha=.2,label='20ns渐消区')
    axes[0].set(title='单道隔离响应：各自峰值参考，仅检查尾端',xlabel='时间(ns)',ylabel='相对各自隔离峰',xlim=(550,650))
    x=np.arange(4)
    for row in depth_rows:
        axes[1].plot(x,[row['tail_taper_sensitivity_ns'][str(n)]['complex_relative_L2_to_no_taper'] for n in (0,10,20,50)],'o-',label=labels[row['case']])
    axes[1].set(title='尾渐消对目标复谱的影响',xticks=x,xticklabels=['0ns','10ns','20ns','50ns'],ylabel='相对未渐消复谱L2差')
    vals=[r['target_complex_band_relative_L2'] for r in layout_rows]
    axes[2].bar(np.arange(4),vals)
    axes[2].set(title='三维界面对比：沿线vs横线',xticks=np.arange(4),
        xticklabels=[r['scenario']+'-'+r['polarisation'] for r in layout_rows],ylabel='目标复谱相对L2差')
    for ax in axes:ax.legend(fontsize=7) if ax!=axes[2] else None;ax.grid(alpha=.2)
    fig.suptitle('目标响应必须单独验收；总场指标与DFT一致性不能给弱回波完整性/误差预算')
    fig.savefig(out/'target_specific_diagnostics.png',dpi=135);plt.close(fig)
    np.savez_compressed(out/'reconstructed_profiles.npz',native_time_ns=native,sfcw_time_ns=time_sf,
         stations_tx_x_m=stations,**{branch+'_'+case:value for (branch,case),value in products.items()})
    summary=dict(status='COMPLETED_READ_ONLY_OBSERVATION_AUDIT_NOT_PHYSICAL_ATTRIBUTION',
        source_git_base='fd0af55',calls_solver=False,calls_training=False,field_fitted=False,
        code_sha256=sha256(__file__),official_processing_sha256=sha256(sf.__file__),
        native_identities=identities,layout_target_comparisons=layout_rows,source_delayed_windows=window_rows,
        depth_tail_and_extrapolation=depth_rows,bscan_observation_comparisons=bscan_rows,
        single_vs_bscan_centre_native_comparisons=centre_match,halfspace_gate_review=half_gate,
        sfcw_assumptions=dict(frequency_Hz=[20e6,170e6],tones=501,zero_pad_factor=8,
          window_branches=['hann','rectangular'],receiver_tail_taper_ns=20,
          domain='source-current-element-normalised field, not port S21',
          historical_instrument_window_unknown=True,native_record_end_ns=650,missing_tail_not_recovered=True),
        numerical_floor_note='DFT agreement and zero pre-arrival differences are not an FDTD target error budget')
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(layout=layout_rows,depth=depth_rows,bscan=bscan_rows,
                         windows=window_rows,half_reflection_pass=half_gate['isolated_reflection_passes_same_5percent_5degree_limits']),indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
