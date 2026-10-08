"""Single-factor mudstone DC-conductivity counterfactual, not field calibration."""
import argparse
import copy
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from analyze_line9_v401_version_controls import response,inverse,FREQ
from gprMax.toolboxes.SFCW.processing import spectral_window


def prepare(a):
    assert not a.out.exists()
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    assert m['spacing_m']==.025 and m['width_m']==80
    shutil.copytree(a.package,a.out)
    for g in m['groups']:
        for key in ['input','geometry','material']:assert sha(a.package/g[key])==g[key+'_sha256']
        path=a.out/g['material'];db=json.loads(path.read_text('utf-8'));before=copy.deepcopy(db)
        mud=db['materials']['material_002_mudstone'];assert mud['base']['electric_conductivity_s_per_m']==.003
        mud['base']['electric_conductivity_s_per_m']=.0003
        check=copy.deepcopy(db);check['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m']=.003
        assert check==before, 'Only DC conductivity may change'
        metadata=mud['metadata'];metadata['parent_parameter_design_sha256']=metadata.pop('parameter_design_sha256')
        metadata['parameter_status']='DC_CONDUCTIVITY_COUNTERFACTUAL_NOT_SITE_CALIBRATION'
        db['database']['description']='Diagnostic only: mudstone sigmaDC .003 -> .0003 S/m; real dielectric spectrum/poles and all other materials unchanged'
        path.write_text(json.dumps(db,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');g['material_sha256']=sha(path)
        assert sha(a.out/g['input'])==g['input_sha256'] and sha(a.out/g['geometry'])==g['geometry_sha256']
    m.update(approval_basis='User authorizes configured V4.0.1 and continued autonomous research; two bounded single-factor material diagnostic solves.',
             generator_sha256=sha(__file__),parent_manifest_sha256=sha(a.package/'manifest.json'),
             diagnostic='Mudstone sigmaDC .003 -> .0003 S/m only; full real permittivity spectrum, Debye pole strengths/times, other materials, source and voxels unchanged.',
             predeclared_question='Does the assumed DC loss materially suppress bottom response even after separating numerical propagation error? Report H0/delta and unfitted total profiles in original fixed gate; do not use field data or adjust parameters after seeing results.',
             limits='Counterfactual tests sensitivity of an unmeasured research assumption; .0003 is not a measured/selected production value or field fit. Planar single station, not nonflat/full-line/finite-3D validation.')
    (a.out/'manifest.json').write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(groups=2,only_physical_change='mudstone conductivity .003 -> .0003 S/m')))


def analyze(a):
    assert not a.out.exists() and not a.numerical.exists()
    c=json.loads((a.source/'execution_contract.json').read_text('utf-8'));v=json.loads((a.source/'completed_verification.json').read_text('utf-8'))
    bc=json.loads((a.baseline/'execution_contract.json').read_text('utf-8'))
    assert v['contract_sha256']==sha(a.source/'execution_contract.json') and v['completed']
    responses=[];audit=[]
    for g,record,old in zip(c['groups'],v['groups'],bc['groups']):
        assert g['input_sha256']==old['input_sha256'] and g['geometry_sha256']==old['geometry_sha256']
        p=a.source/(g['id']+'.h5');assert sha(p)==record['native_sha256']
        with h5py.File(p) as h:assert h.attrs['gprMax']=='4.0.1' and h['rxs/rx1/Ez'].dtype==np.float64
        z,err=response(p,reference_scale=.025);responses.append(z);audit.append(err)
    changed=np.column_stack([*responses,responses[1]-responses[0]])
    with h5py.File(a.baseline/'sfcw.h5') as h:np.testing.assert_array_equal(h['frequency_Hz'][:],FREQ);baseline=h['response'][:]
    m=c['study_manifest'];metrics={};profiles={}
    for window in ['hann','blackman']:
        w=spectral_window(window,501);old,t=inverse(baseline,w);new,_=inverse(changed,w);profiles[window]=(old,new)
        mask=(t*1e9>=m['basal_gate_ns'][0])&(t*1e9<=m['basal_gate_ns'][1]);norm=lambda z:np.linalg.norm(z[mask])
        metrics[window]=dict(baseline_H0_over_delta=float(norm(old[:,0])/norm(old[:,2])),
                             changed_H0_over_delta=float(norm(new[:,0])/norm(new[:,2])),
                             changed_delta_over_baseline_delta=float(norm(new[:,2])/norm(old[:,2])),
                             changed_H0_over_baseline_H0=float(norm(new[:,0])/norm(old[:,0])),
                             changed_H1_vs_delta_correlation=float(abs(np.vdot(new[mask,2],new[mask,1]))/(norm(new[:,2])*norm(new[:,1]))),
                             baseline_H1_vs_delta_correlation=float(abs(np.vdot(old[mask,2],old[mask,1]))/(norm(old[:,2])*norm(old[:,1]))))
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,(old,new) in profiles.items():
        fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained');mask=(t*1e9>=340)&(t*1e9<=405)
        for ax,indices,title in [(axes[0],[0,1],'原始H0/H1总场，共同灰度'),(axes[1],[2],'配对底砂差场，共同灰度')]:
            values=np.column_stack([old[mask][:,indices],new[mask][:,indices]]).real;lim=float(abs(values).max())
            im=ax.imshow(values,cmap='gray',vmin=-lim,vmax=lim,aspect='auto',extent=[-.5,values.shape[1]-.5,405,340]);ax.set_xticks(range(values.shape[1]),['σ.003 H0','σ.003 H1','σ.0003 H0','σ.0003 H1'] if len(indices)==2 else ['σ.003 差场','σ.0003 差场'],rotation=18,fontsize=8);ax.set(title=title,ylabel='SFCW时间 / ns');ax.axhline(371.756487,color='green',ls=':',lw=1);fig.colorbar(im,ax=ax,label='共同参考 (V/m)/(A·m)')
        for z,label,style in [(old[:,1],'σ.003 总场','-'),(old[:,2],'σ.003 底砂差场','--'),(new[:,1],'σ.0003 总场','-'),(new[:,2],'σ.0003 底砂差场','--')]:axes[2].plot(t[mask]*1e9,abs(z[mask]),style,lw=1,label=label)
        axes[2].set(title='共同绝对幅度，无AGC或拟合',xlabel='SFCW时间 / ns',ylabel='包络');axes[2].legend(fontsize=8)
        fig.suptitle('泥岩DC电导率单因素敏感性 / '+window+' / 4.0.1 FP64 / 80m平层\nε′完整谱、Debye极点、几何和航高不变；同站配置列，非空间测线；不是实测拟合')
        fig.savefig(a.out/f'v401_conductivity_{window}_gray.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=changed)
    result=dict(status='COMPLETED_SINGLE_FACTOR_CONDUCTIVITY_SENSITIVITY_NOT_FIELD_CALIBRATION',script_sha256=sha(__file__),
                contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
                numerical_sha256=sha(a.numerical),baseline_numerical_sha256=sha(a.baseline/'sfcw.h5'),
                metrics=metrics,independent_DFT_relative_L2=audit,limits=m['limits'],AGC=False,fit=False)
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print(json.dumps(metrics))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','analyze'])
    for name in ['package','out','source','baseline','numerical']:p.add_argument('--'+name,type=Path,required=name=='out')
    a=p.parse_args();prepare(a) if a.action=='prepare' else analyze(a)
