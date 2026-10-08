"""Read-only 194-trace v5 path-time overlay, no event identification or solver."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from review_line9_result_packages import indices,ray_time,C0
from analyze_line9_v401_version_controls import inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window


def main(a):
    assert not a.out.exists()
    review=json.loads(a.review.read_text('utf-8'));records=review['records']
    assert review['package']=='line9_pkg1_v5_d02m_194st' and len(records)==194
    material=a.package/'geometries/line9_research_materials_v1_smoothed.json'
    assert sha(material)==review['materials_sha256'];materials=json.loads(material.read_text('utf-8'))['materials']
    with h5py.File(a.corrected) as h:
        np.testing.assert_array_equal(h['frequency'][:],FREQ)
        assert h.attrs['NoBackgroundOrGainOrTimePartition'] and h.attrs['TailTaperFraction']==0 and h.attrs['NativeDtype']=='float32'
        z=h['response'][:];positions=h['chainage_m'][:];ids=[x.decode() for x in h['station_id'][:]]
        assert ids==[r['id'] for r in records]
        native_hash=[x.decode() for x in h['native_sha256'][:]]
    selected=[0,83,167,250,333,417,500];errors=[];height=[];cover=[];mud=[];dx=[]
    for j,row in enumerate(records):
        p=a.package/'cases'/row['id']/'profile.h5';assert sha(p)==row['native_sha256']==native_hash[j]
        with h5py.File(p) as h:
            x=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:];dt=float(h.attrs['dt'])
            assert x.dtype==s.dtype==np.float32 and h.attrs['gprMax']=='4.0.0'
            tr=np.arange(len(x))*dt+h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
            ts=np.arange(len(s))*dt+h['srcs/src1/excitation'].attrs['TimeSampleOffset']
            f=FREQ[selected,None]
            direct=(np.exp(-2j*np.pi*f*tr)@x)/(np.exp(-2j*np.pi*f*ts)@s)/.025
            err=float(np.linalg.norm(direct-z[selected,j])/np.linalg.norm(direct));assert err<1e-9;errors.append(err)
            tx=h['srcs/src1'].attrs['Position'];rx=h['rxs/rx1'].attrs['Position']
        g=row['geometry'];assert len(g['boundaries'])==3
        height.append(g['midpoint_agl_m']);cover.append(g['cover_base']['depth_m']);mud.append(g['basal_sand']['depth_m']-cover[-1]);dx.append(abs(tx[0]-rx[0]))
    height=np.array(height);cover=np.array(cover);mud=np.array(mud)
    n=np.array([indices(materials,f) for f in FREQ]);kc=2*np.pi*FREQ[:,None]*n/C0
    kg=2/.025*np.arcsin(n*.025/(C0*dt)*np.sin(np.pi*FREQ*dt)[:,None])
    r01=(n[:,0]-n[:,1])/(n[:,0]+n[:,1]);r12=(n[:,1]-n[:,2])/(n[:,1]+n[:,2]);r23=(n[:,2]-n[:,3])/(n[:,2]+n[:,3])
    predictions={};curves={};profiles={};stats={}
    for model,k in [('continuum',kc),('bulk_Yee',kg)]:
        primary=[];multiple=[]
        for h,d1,d2,sep in zip(height,cover,mud,dx):
            # Local horizontal symmetric-ray correction is fixed from geometry.
            bpath=np.array([h,d1,d2,0]);mpath=np.array([h,2*d1,0,0])
            bcorr=ray_time(bpath,n[250],sep)*1e-9-2*np.dot(bpath,n[250].real)/C0
            mcorr=ray_time(mpath,n[250],sep)*1e-9-2*np.dot(mpath,n[250].real)/C0
            primary.append((1-r01*r01)*(1-r12*r12)*r23*np.exp(-2j*(k@bpath))*np.exp(-2j*np.pi*FREQ*bcorr))
            multiple.append(-(1-r01*r01)*r01*r12*r12*np.exp(-2j*(k@mpath))*np.exp(-2j*np.pi*FREQ*mcorr))
        predictions[model]=np.column_stack([*primary,*multiple])
    for window in ['hann','blackman']:
        w=spectral_window(window,501);actual,t=inverse(z,w);profiles[window]=actual;curves[window]={}
        for model,pred in predictions.items():
            x,_=inverse(pred,w);peaks=t[np.argmax(abs(x),axis=0)]*1e9
            b=peaks[:194];m=peaks[194:];sep=abs(b-m)
            curves[window][model]=dict(basal_primary_ns=b.tolist(),cover_second_trip_ns=m.tolist())
            stats[window+'_'+model]=dict(absolute_separation_ns_percentiles_0_50_100=np.percentile(sep,[0,50,100]).tolist(),
                stations_separation_le_reciprocal_bandwidth=int(np.sum(sep<=1e9/150e6)),
                stations_separation_le_twice_reciprocal_bandwidth=int(np.sum(sep<=2e9/150e6)),
                note='Descriptive delay comparison, not a validated resolution criterion or measured event classification.')
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    order=np.argsort(positions);xp=positions[order];xe=np.r_[xp[0]-(xp[1]-xp[0])/2,(xp[:-1]+xp[1:])/2,xp[-1]+(xp[-1]-xp[-2])/2]
    keep=(t*1e9>=250)&(t*1e9<=430);tt=t[keep]*1e9;te=np.r_[tt[0]-(tt[1]-tt[0])/2,(tt[:-1]+tt[1:])/2,tt[-1]+(tt[-1]-tt[-2])/2]
    limit=max(float(abs(profiles[w][keep].real).max()) for w in profiles)
    fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
    for ax,window in zip(axes,['hann','blackman']):
        im=ax.pcolormesh(xe,te,profiles[window][keep][:,order].real,cmap='gray',vmin=-limit,vmax=limit,shading='flat',rasterized=True)
        for model,style in [('continuum',':'),('bulk_Yee','--')]:
            q=curves[window][model]
            for key,color,label in [('basal_primary_ns','limegreen','底砂一次'),('cover_second_trip_ns','orange','覆盖层第二次往返')]:
                ax.plot(xp,np.array(q[key])[order],style,color=color,lw=1.1,label=label+('连续预测' if model=='continuum' else '体网格传播预测'))
        ax.set(ylim=(430,250),xlabel='剖面里程 / m（仅原194道完成区）',ylabel='SFCW时间 / ns',title=window+' 原始总场，共同绝对灰度');ax.invert_xaxis();ax.legend(fontsize=8)
        fig.colorbar(im,ax=ax,label='共同参考 (V/m)/(A·m)')
    fig.suptitle('原v5 194道 / AGL8m / 原生4.0.0 FP32 / 精确501点，无AGC/去背景\n曲线为几何驱动局部平层候选路径，非观测回波身份认证；未计算的其余测线不补齐')
    fig.savefig(a.out/'v5_original_194_path_time_overlay.png',dpi=140);plt.close(fig)
    result=dict(status='OFFLINE_MODEL_INFORMED_PATH_TIME_OVERLAY_NOT_EVENT_ATTRIBUTION',script_sha256=sha(__file__),
        corrected_sha256=sha(a.corrected),review_sha256=sha(a.review),material_sha256=sha(material),native_count=194,
        original_dtype='float32',solver_version='4.0.0',new_solves=0,calls_training=False,
        independent_seven_frequency_DFT_all_native_relative_L2_max=max(errors),chainage_m=positions.tolist(),
        curves=curves,statistics=stats,AGC=False,fit=False,
        limits='Local normal-incidence Fresnel paths with geometric bistatic correction; optional normal bulk Yee propagation. No angular spreading,slopes,lateral scattering,other multiples,discrete-interface errors or event identity. Seven-frequency audit of archived corrected spectra, not new all-frequency FDTD certification. Native FP32 remains FP32 evidence.')
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(stats))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','review','corrected','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
