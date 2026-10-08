"""Exact SFCW factorial comparison; no field fitting or production selection."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window

IDS=['span3_dc003','span3_dc0003','span03_dc003','span03_dc0003']
NAMES=['跨度3 / DC3mS/m','跨度3 / DC0.3mS/m','跨度0.3 / DC3mS/m','跨度0.3 / DC0.3mS/m']


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    c=json.loads((a.source/'execution_contract.json').read_text('utf-8'))
    v=json.loads((a.source/'completed_verification.json').read_text('utf-8'))
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert c['study_manifest']==m and v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    paths={r['id']:(a.package/r['path'],r['native_sha256']) for r in m['reused']}
    for g,r in zip(c['groups'],v['groups']):
        assert g['id']==r['id'];paths[g['id']]=(a.source/(g['id']+'.h5'),r['native_sha256'])
    vals=[];native=[]
    for name in IDS:
        pair=[]
        for role in ['H0','H1']:
            label=name+'_'+role;p,digest=paths[label];assert sha(p)==digest
            with h5py.File(p) as h:
                assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64
            z,error=response(p,.025);pair.append(z)
            native.append(dict(id=label,native_sha256=digest,independent_DFT_relative_L2=error))
        vals.extend([pair[0],pair[1],pair[1]-pair[0]])
    z=np.column_stack(vals);metrics={};profiles={}
    for window in ['hann','blackman']:
        x,t=inverse(z,spectral_window(window,501));profiles[window]=x;metrics[window]={}
        for gate in ['basal','deep']:
            lo,hi=m[gate+'_gate_ns'];k=(t*1e9>=lo)&(t*1e9<=hi)
            rows={}
            for j,name in enumerate(IDS):
                h0,h1,d=x[k,3*j:3*j+3].T
                norm=np.linalg.norm
                rows[name]=dict(H0_norm=float(norm(h0)),H1_norm=float(norm(h1)),delta_norm=float(norm(d)),
                    H0_over_delta=float(norm(h0)/norm(d)),
                    H1_vs_delta_correlation=float(abs(np.vdot(h1,d))/(norm(h1)*norm(d))),
                    delta_vs_original_correlation=float(abs(np.vdot(x[k,2],d))/(norm(x[k,2])*norm(d))),
                    delta_peak_time_ns=float(t[k][np.argmax(abs(d))]*1e9))
            d=[rows[n]['delta_norm'] for n in IDS]
            metrics[window][gate]=dict(cases=rows,effects=dict(
                DC_reduction_at_span3=d[1]/d[0],DC_reduction_at_span03=d[3]/d[2],
                span_reduction_at_DC003=d[2]/d[0],span_reduction_at_DC0003=d[3]/d[1],
                both_reductions_over_original=d[3]/d[0],
                multiplicative_interaction=d[3]*d[0]/(d[2]*d[1])))
    # Constitutive budgets for actual frozen mudstone assignments. No path/template fit.
    parent=json.loads((a.parent/'baseline_H1/geometries/line9_research_materials_v1_smoothed.json').read_text('utf-8'))
    byid={g['id']:g for g in m['groups']};budgets=[];material_curves=[]
    for name in IDS:
        if name.startswith('span3_'):
            db=parent;dc=.003 if name=='span3_dc003' else .0003
        else:
            db=json.loads((a.package/byid[name+'_H1']['material']).read_text('utf-8'))
            dc=db['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m']
        mud=db['materials']['material_002_mudstone'];pole=mud['poles'][0]
        de=pole['relative_permittivity_difference'];tau=pole['relaxation_time_s']
        er=mud['base']['relative_permittivity']+de/(1+2j*np.pi*FREQ*tau)-1j*dc/(2*np.pi*FREQ*8.8541878128e-12)
        alpha=-(2*np.pi*FREQ/299792458*np.sqrt(er)).imag
        material_curves.append((er,-20*14*alpha/np.log(10)))
        at95=int(np.flatnonzero(FREQ==95e6)[0])
        budgets.append(dict(id=name,epsilon_real_20_95_170=er.real[[0,at95,-1]].tolist(),
            epsilon_imag_20_95_170=(-er.imag[[0,at95,-1]]).tolist(),
            sigma_DC_S_m=dc,epsilon_infinity=mud['base']['relative_permittivity'],delta_epsilon=de,tau_s=tau,
            absorption_7m_two_way_dB_20_95_170=(-20*14*alpha[[0,at95,-1]]/np.log(10)).tolist()))
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,x in profiles.items():
        k=(t*1e9>=300)&(t*1e9<=450)
        fig,axes=plt.subplots(2,2,figsize=(16,9),layout='constrained')
        # Common limit across all twelve raw/difference columns; no independent normalization.
        lim=float(abs(x[k].real).max())
        for ax,offset,title in [(axes[0,0],0,'H0：仅替换底连通砂岩后的地质响应'),(axes[0,1],1,'H1：保留原底砂的总场'),(axes[1,0],2,'H1−H0：底砂替换差场')]:
            im=ax.imshow(x[k,offset::3].real,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',interpolation='nearest',extent=[-.5,3.5,450,300])
            ax.set_xticks(range(4),NAMES,rotation=12,fontsize=8)
            ax.set(title=title,ylabel='SFCW时间 / ns')
            for y in m['basal_gate_ns']:ax.axhline(y,color='green',ls=':',lw=1)
            fig.colorbar(im,ax=ax,label='共同绝对幅度 (V/m)/(A·m)')
        for j,label in enumerate(NAMES):axes[1,1].plot(t[k]*1e9,abs(x[k,3*j+2]),label=label,lw=1)
        for y in m['basal_gate_ns']:axes[1,1].axvline(y,color='green',ls=':',lw=1)
        axes[1,1].set(title='四组底砂差场包络，共同幅度',xlabel='SFCW时间 / ns',ylabel='(V/m)/(A·m)');axes[1,1].legend(fontsize=9)
        fig.suptitle('原非平v5 / X198.6m / AGL8m / 4.0.1 FP64 / '+window+'\n泥岩极化与DC 2×2；同站配置列，不是空间B-scan；三灰度面板同标尺；无AGC/拟合；正式材料未换')
        fig.savefig(a.out/f'polarization_dc_{window}_gray.png',dpi=140);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    for j,((er,attenuation),label) in enumerate(zip(material_curves,NAMES)):
        if j in [0,2]:axes[0].plot(FREQ/1e6,er.real,label='跨度3' if j==0 else '跨度0.3')
        axes[1].plot(FREQ/1e6,-er.imag,label=label)
        axes[2].plot(FREQ/1e6,attenuation,label=label)
    for ax,title,ylabel in zip(axes,['中心实部相同，整段实部不同','总虚部 = 极化 + DC','仅泥岩7m往返的连续吸收预算'],['ε′','ε″','场幅吸收 / dB']):
        ax.set(title=title,xlabel='频率 / MHz',ylabel=ylabel);ax.legend(fontsize=8);ax.axvline(95,color='gray',ls=':',lw=.8)
    fig.suptitle('四组泥岩研究假设；完整复介电谱必须联合检查\n解析预算不含扩散/界面/侧向散射/仪器，不是现场SNR；正式材料未换')
    fig.savefig(a.out/'polarization_dc_material_budget.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=z)
        h.attrs['column_ids']=json.dumps([n+'_'+r for n in IDS for r in ['H0','H1','delta']])
    result=dict(status='COMPLETED_NONFLAT_POLARIZATION_DC_FACTORIAL_NOT_SITE_CALIBRATION',script_sha256=sha(__file__),
        contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
        numerical_sha256=sha(a.numerical),native=native,metrics=metrics,constitutive_budget=budgets,
        gate_ns={g:m[g+'_gate_ns'] for g in ['basal','deep']},new_solves=4,reused_native=4,
        AGC=False,fit=False,production_materials_changed=False,limits=m['limits'])
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(metrics=metrics,constitutive_budget=budgets)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','parent','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
