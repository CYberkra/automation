"""Same planar constitutive stack: 2D TE line versus full3D horizontal point source."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from scipy.special import hankel2
from line9_planar_point_layer_green import point_integral,layer_parts
from diagnose_line9_v5_planar_green import integral,MU0
from diagnose_line9_direct_dimensionality import point_dipole
from review_line9_result_packages import indices,C0
from analyze_line9_v401_version_controls import FREQ,inverse
from hs_capsule_identity import sha256 as sha


def main(a):
    assert not a.out.exists() and not a.arrays.exists()
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    geometry=a.package/'baseline/geometry.h5';database=a.package/'baseline/materials.json'
    assert sha(geometry)==m['baseline_sha256']['geometry.h5'] and sha(database)==m['baseline_sha256']['materials.json']
    db=json.loads(database.read_text('utf-8'))['materials']
    tx=np.array([169.35,38.95]);rx=np.array([170.65,39.025])
    with h5py.File(geometry) as h:
        col=h['data'][6800,:,0];surface=(np.flatnonzero(col!=0)[-1]+1)*.025;bottom=(np.flatnonzero(col==2)[-1]+1)*.025
    hs=tx[1]-surface;hr=rx[1]-surface;offset=float(abs(rx[0]-tx[0]));ng=np.array([indices(db,f) for f in FREQ]);eps1=ng[:,1]**2;eps2=ng[:,2]**2;k0=2*np.pi*FREQ/C0
    thicknesses=[5.6,float(surface-bottom),7.3]
    contract={'status':'FROZEN_CPU_PLANAR_MECHANISM_CONTRACT_NOT_FDTD_RUN','code_sha256':sha(__file__),
              'helper_sha256':{n:sha(Path(__file__).with_name(n)) for n in ['line9_planar_point_layer_green.py','diagnose_line9_v5_planar_green.py','diagnose_line9_direct_dimensionality.py','review_line9_result_packages.py','analyze_line9_v401_version_controls.py']},
              'material_sha256':sha(database),'geometry_sha256':sha(geometry),'midpoint_surface_m':surface,'midpoint_bottom_m':bottom,
              'air_source_rx_heights_m':[hs,hr],'offset_m':offset,'cover_thicknesses_m':thicknesses,
              'factor':'2D uniform cross-axis TE line current versus3D transverse horizontal point current element; infinite planar interfaces both cases',
              'frequency_Hz':FREQ.tolist(),'gates_ns':{'first_cover':[150,250],'late':[300,450]},'windows':['hann','blackman'],
              'quadrature_orders':[256,512],'quadrature_relative_L2_tolerance':1e-6,'tail_exponent':36.,
              'point_bearing_from_dipole_rad':float(np.pi/2),'units':{'line':'(V/m)/A','point':'(V/m)/(A*m)'},
              'metric':'Each source type: exact cover-twice norm / own cover-once norm; do not compare absolute source powers or field amplitudes across dimensions',
              'limitations':'Planar surrogate only, not actual nonflat volume/finite geology/antenna/field validation; no new gprMax or stopped-attempt restart.'}
    a.out.mkdir(parents=True);(a.out/'contract.json').write_text(json.dumps(contract,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    responses=[];convergence=[]
    for thickness in thicknesses:
        kernel=lambda q,ky:layer_parts(q,ky,k0,eps1,eps2,thickness)
        values=[]
        for order in contract['quadrature_orders']:
            line=-FREQ[:,None]*MU0*integral(FREQ,hs+hr,offset,lambda q,ky:kernel(q,ky)[0],order=order)
            point=point_integral(FREQ,hs+hr,offset,np.pi/2,kernel,order=order)
            values.append(np.stack([line,point],axis=1))
        err=np.linalg.norm(values[1]-values[0],axis=0)/np.linalg.norm(values[1],axis=0)
        convergence.append(err.tolist())
        assert np.max(err)<contract['quadrature_relative_L2_tolerance'],err
        responses.append(values[1])
    response=np.stack(responses,axis=1) # frequency,thickness,dimension,component
    direct_line=-2*np.pi*FREQ*MU0/4*hankel2(0,k0*np.hypot(offset,hs-hr))
    direct_point=point_dipole(FREQ,np.hypot(offset,hs-hr),np.zeros(501))
    direct=np.column_stack([direct_line,direct_point]);profiles={};metrics={};inverse_errors={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();z,t=inverse(response.reshape(501,-1),w);z=z.reshape(len(t),3,2,4)
        d,_=inverse(direct,w);total=z.sum(axis=-1)+d[:,None,:]
        profiles[window]=(z,total);rows=[]
        first=(t*1e9>=150)&(t*1e9<=250);late=(t*1e9>=300)&(t*1e9<=450)
        for j,thickness in enumerate(thicknesses):
            rr={}
            for dim,label in enumerate(['2D_TE_line','3D_transverse_point']):
                once=z[:,j,dim,1];twice=z[:,j,dim,2];higher=z[:,j,dim,3]
                rr[label]={'cover_once_first_window_norm':float(np.linalg.norm(once[first])),
                           'cover_twice_late_window_norm':float(np.linalg.norm(twice[late])),
                           'twice_late_over_once_first':float(np.linalg.norm(twice[late])/np.linalg.norm(once[first])),
                           'twice_full_over_once_full':float(np.linalg.norm(twice)/np.linalg.norm(once)),
                           'higher_late_over_twice_late':float(np.linalg.norm(higher[late])/np.linalg.norm(twice[late])),
                           'total_late_over_once_first':float(np.linalg.norm(total[late,j,dim])/np.linalg.norm(once[first])),
                           'once_peak_ns':float(t[np.argmax(abs(once))]*1e9),
                           'twice_peak_ns':float(t[np.argmax(abs(twice))]*1e9)}
            r=rr['3D_transverse_point']['twice_late_over_once_first']/rr['2D_TE_line']['twice_late_over_once_first']
            rows.append({'cover_thickness_m':thickness,'dimensions':rr,'3D_over2D_relative_return_ratio':r,'relative_return_change_dB':float(20*np.log10(r))})
        metrics[window]=rows
        pick=np.arange(7,len(t),101)
        direct_sum=np.exp(2j*np.pi*t[pick,None]*FREQ)@(response.reshape(501,-1)*w[:,None])/501
        er=float(np.linalg.norm(z[pick].reshape(len(pick),-1)-direct_sum)/np.linalg.norm(direct_sum));assert er<1e-9
        inverse_errors[window]=er
    a.arrays.parent.mkdir(parents=True,exist_ok=True);np.savez(a.arrays,frequency_Hz=FREQ,response=response,direct=direct,time_s=t,
        **{w+'_components':q[0] for w,q in profiles.items()},**{w+'_total':q[1] for w,q in profiles.items()})
    summary={'status':'PLANAR_CPU_DIMENSIONAL_MECHANISM_DIAGNOSTIC_NOT_NONFLAT3D_CERTIFICATION',
             'script_sha256':sha(__file__),'contract_sha256':sha(a.out/'contract.json'),'arrays_sha256':sha(a.arrays),
             'quadrature256_512_relative_L2_by_thickness_dimension_component':convergence,'inverse_direct_sum_relative_L2':inverse_errors,
             'metrics':metrics,'gprMax_new_solves':0,'source_amplitude_or_delay_fit':False,'site_parameters_changed':False,
             'limits':contract['limitations']}
    (a.out/'analysis.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained');tn=t*1e9
    labels=['地表项','覆盖层底一次返回','覆盖层两次返回项','更高阶剩余'];z,total=profiles['hann'];show=(tn>=100)&(tn<=450)
    for dim,label,unit in [(0,'二维 TE线源','(V/m)/A'),(1,'三维横向点偶极','(V/m)/(A·m)')]:
        ax=axs[0,dim]
        for k,l in enumerate(labels):ax.semilogy(tn[show],abs(z[show,1,dim,k]),lw=1,label=l)
        ax.set(xlabel='SFCW时间 / ns',ylabel='分量包络 / '+unit,title=f'同一平层6.575m覆盖层；{label}，独立物理单位');ax.legend(fontsize=9)
    for window,style in [('hann','-'),('blackman','--')]:
        for dim,label in [('2D_TE_line','二维 TE线源'),('3D_transverse_point','三维横向点偶极')]:
            axs[1,0].plot(thicknesses,[r['dimensions'][dim]['twice_late_over_once_first'] for r in metrics[window]],style+'o',label=label+' / '+window)
        axs[1,1].plot(thicknesses,[r['relative_return_change_dB'] for r in metrics[window]],style+'o',label=window)
    axs[1,0].set(xlabel='平层覆盖层厚度 / m',ylabel='两次返回晚窗范数 / 本源一次返回早窗范数',title='源强度消去的相对返回量；不是能量占比');axs[1,0].legend(fontsize=8)
    axs[1,1].set(xlabel='平层覆盖层厚度 / m',ylabel='三维/二维 相对返回变化 / dB',title='只比较各自“晚/早”比值；未比较源功率');axs[1,1].axhline(0,color='gray',lw=.7);axs[1,1].legend()
    fig.suptitle('CPU角谱参考：三维点源是否令覆盖层重复返回消失？\n无限平界面、相同材料/收发高度/精确501频点；不代表实际起伏二维模型的三维校核，无AGC或拟合')
    fig.savefig(a.out/'planar_dimensional_return_comparison.png',dpi=140);plt.close(fig)
    print(json.dumps({'status':summary['status'],'quadrature_max_error':float(np.max(convergence)),'metrics':metrics}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['package','out','arrays']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
