"""Exploratory surface re-reflection sequences from immutable native observers."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha

EPS0=8.8541878128e-12
C0=299792458.


def group_index(material,f=95e6):
    base=material['base'];sigma=base['electric_conductivity_s_per_m']
    er=base['relative_permittivity']-1j*sigma/(2*np.pi*f*EPS0)
    derivative=1j*sigma/(2*np.pi*EPS0*f*f)
    for pole in material.get('poles',[]):
        d=pole['relative_permittivity_difference'];tau=pole['relaxation_time_s'];den=1+2j*np.pi*f*tau
        er+=d/den;derivative-=2j*np.pi*tau*d/den**2
    n=np.sqrt(er)
    return float((n+f*derivative/(2*n)).real)


def main(a):
    assert not a.out.exists()
    read=lambda p:json.loads(p.read_text('utf-8'))
    m=read(a.package/'manifest.json');s=read(a.parent/'analysis.json');audit=read(a.parent/'independent_audit.json')
    assert audit['analysis_sha256']==sha(a.parent/'analysis.json') and audit['status'].startswith('PASS')
    assert sha(a.arrays)==s['arrays_sha256'] and sha(a.raw)==s['native_sha256']
    db=read(a.package/'baseline/materials.json')['materials'];ng=group_index(db['material_001_cover'])
    with np.load(a.arrays) as h:
        t=h['time_s']*1e9;sy=h['Sy'];sx=h['Sx']
    lookup={(p['band'],p['x_m']):i for i,p in enumerate(m['probes'])};rows=[]
    def extremum(k,lo,hi,sign):
        ids=np.flatnonzero((t>=lo)&(t<=hi));j=ids[np.argmax(sy[k,ids]) if sign=='up' else np.argmin(sy[k,ids])]
        win=(t>=t[j]-2)&(t<=t[j]+2)
        return {'index':int(j),'time_ns':float(t[j]),'Sy_W_m2':float(sy[k,j]),
                'local4ns_signed_J_m2':float(np.sum(sy[k,win])*s['dt_s']),
                'local4ns_flux_direction_deg':float(np.degrees(np.arctan2(np.sum(sy[k,win]),np.sum(sx[k,win]))))}
    for x in sorted({p['x_m'] for p in m['probes']}):
        ku=lookup['cover_upper',x];km=lookup['cover_middle',x];ka=lookup['air_surface',x]
        up=extremum(ku,160,210,'up');down=extremum(ku,up['time_ns']+2,up['time_ns']+15,'down')
        midup=extremum(km,120,175,'up');middown=extremum(km,up['time_ns']+15,up['time_ns']+70,'down')
        airup=extremum(ka,up['time_ns'],up['time_ns']+15,'up')
        p=m['probes'][ku];q=m['probes'][km]
        rows.append({'x_m':x,'upper_y_m':p['y_m'],'middle_y_m':q['y_m'],
                     'upper_up':up,'upper_later_down':down,'middle_first_up':midup,'middle_later_down':middown,'air_first_up':airup,
                     'upper_up_to_down_ns':down['time_ns']-up['time_ns'],
                     'upper_stencil_surface_roundtrip_proxy_ns':2*(p['surface_m']-p['y_m'])*ng/C0*1e9,
                     'middle_up_to_upper_up_ns':up['time_ns']-midup['time_ns'],
                     'upper_down_to_middle_down_ns':middown['time_ns']-down['time_ns'],
                     'vertical_upper_middle_group95_proxy_ns':abs(p['y_m']-q['y_m'])*ng/C0*1e9})
    a.out.mkdir(parents=True)
    summary={'status':'EXPLORATORY_NATIVE_SURFACE_REREFLECTION_SEQUENCE','script_sha256':sha(__file__),
             'parent_analysis_sha256':sha(a.parent/'analysis.json'),'native_sha256':sha(a.raw),'arrays_sha256':sha(a.arrays),
             'cover_group_index95':ng,'rows':rows,'rows_with_expected_local4ns_signs':sum(r['upper_up']['local4ns_signed_J_m2']>0 and r['upper_later_down']['local4ns_signed_J_m2']<0 and r['middle_first_up']['local4ns_signed_J_m2']>0 and r['middle_later_down']['local4ns_signed_J_m2']<0 for r in rows),
             'exploratory_gates':{'upper_up':[160,210],'upper_down_relative_to_up':[2,15],'middle_first_up':[120,175],'middle_down_relative_to_upper_up':[15,70],'air_up_relative_to_upper_up':[0,15],'flux_integral_relative_to_peak':[-2,2]},
             'limits':'Post-observation extrema. Near-interface superposition changes apparent peak time and net flux. Vertical95MHz group proxies are not ray solutions; no unique bounce order or reflected-energy fraction.'}
    (a.out/'analysis.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    for ax,x in zip(axs.flat,[158,160,164,170]):
        r=next(r for r in rows if r['x_m']==x);mask=(t>=135)&(t<=235)
        for band,label in [('cover_upper','覆盖层上部'),('cover_middle','覆盖层中部'),('air_surface','地表上方空气')]:
            ax.plot(t[mask],sy[lookup[band,x],mask],label=label,lw=1)
        ax.axhline(0,color='gray',lw=.7);ax.set_ylim(-.13,.8)
        ax.set(title=f'模型 x={x}m；三个真实探针同尺度',xlabel='原生时间 / ns',ylabel='总场 Sy / (W/m²)');ax.legend(fontsize=9)
    fig.suptitle('首次上行到达地表附近后，覆盖层上部先转为下行，中部随后出现下行波包\n同一次激发、总场原幅值；四位置为探索示例，非移动天线B-scan，峰值不独立数阶次')
    fig.savefig(a.out/'surface_rereflection_signed_traces.png',dpi=140);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(12,9),layout='constrained');xs=[r['x_m'] for r in rows]
    for key,label in [('middle_first_up','中部首次上行'),('upper_up','上部上行'),('air_first_up','空气上行'),('upper_later_down','上部后续下行'),('middle_later_down','中部后续下行')]:
        axs[0].plot(xs,[r[key]['time_ns'] for r in rows],'.-',label=label,lw=1)
    axs[0].set(xlabel='模型 x / m',ylabel='原生极值时间 / ns',title='全部41列：事件窗已声明为事后探索');axs[0].legend(ncol=2,fontsize=9)
    for key,label in [('upper_up_to_down_ns','上部上行→后续下行'),('upper_stencil_surface_roundtrip_proxy_ns','上部点到地表往返：垂直95MHz群速近似')]:
        axs[1].plot(xs,[r[key] for r in rows],label=label,lw=1)
    axs[1].set(xlabel='模型 x / m',ylabel='时间差 / ns',title='近界面峰偏移含入射/反射叠加，差值不能单独认证射线');axs[1].legend(fontsize=9)
    fig.suptitle('地表再次反射的原生时序相容性；不是唯一往返次数或能量预算')
    fig.savefig(a.out/'surface_rereflection_all_positions.png',dpi=140);plt.close(fig)
    print(json.dumps({'status':summary['status'],'sign_consistent_rows':summary['rows_with_expected_local4ns_signs'],'example160':next(r for r in rows if r['x_m']==160)}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['package','parent','arrays','raw','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
