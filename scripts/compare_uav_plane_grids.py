"""Geometry-preserving grid comparison after an independent transverse-width control."""
import argparse
import json
from pathlib import Path
import numpy as np
from analyze_uav_local_plane import C,F,OMEGA,reference,spectrum
from hs_capsule_identity import sha256

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/research_checks'
PREFIX='2026-10-05_local_uav_'


def main(out):
    if out.exists():
        raise ValueError('new output required')
    responses={};summaries={};identities={}
    for name,result in [('plane_r5','plane_results_r5'),('plane_domain_r7','plane_results_domain_r7'),
                        ('plane_fine_r6','plane_results_fine_r6')]:
        study=BASE/(PREFIX+name)
        cp=study/'execution_contract.json';c=json.loads(cp.read_text('utf-8'))
        verification=json.loads((study/'completed_verification.json').read_text('utf-8'))
        sp=BASE/(PREFIX+result)/'summary.json';summary=json.loads(sp.read_text('utf-8'))
        if summary['contract_sha256']!=sha256(cp) or verification['contract_sha256']!=sha256(cp):
            raise ValueError('contract identity differs')
        summaries[name]=summary;identities[name]={'contract_sha256':sha256(cp),'analysis_sha256':sha256(sp)}
        fields={}
        for g,v in zip(c['groups'],verification['groups']):
            raw=Path(g['input']).with_suffix('.h5')
            if sha256(raw)!=v['raw_sha256'] or sha256(raw)!=summary['raw_identities'][g['id']]:
                raise ValueError('native identity differs')
            for i in (1,3):
                fields[g['id'],i]=spectrum(raw,i)[0]
        for case in ('halfspace','slab','conductive','debye'):
            responses[name,case,'reflection']=(fields[case,1]-fields['air',1])/fields['air',1]
            responses[name,case,'transmission']=fields[case,3]/fields['air',3]
    domain={}
    for case in ('halfspace','slab','conductive','debye'):
        for kind in ('reflection','transmission'):
            a,b=responses['plane_r5',case,kind],responses['plane_domain_r7',case,kind]
            domain[case+'_'+kind]=float(np.linalg.norm(a-b)/np.linalg.norm(a))
    if max(domain.values())>1e-10:
        raise ValueError('transverse width change is not negligible; cannot assign refinement to grid alone')
    coarse=summaries['plane_domain_r7']['metrics'];fine=summaries['plane_fine_r6']['metrics']
    if not all(b['complex_relative_L2']<a['complex_relative_L2']*.5 for a,b in zip(coarse,fine)):
        raise ValueError('declared grid-error decrease not achieved')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,3,figsize=(13,4.1),layout='constrained')
    prediction_residual={}
    k0=OMEGA/C
    for name,dl,label in [('plane_domain_r7',.05,'5cm'),('plane_fine_r6',.0125,'1.25cm')]:
        dt=dl/(C*np.sqrt(2));ka=2*np.arcsin(np.sqrt(2)*np.sin(OMEGA*dt/2))/dl
        k1=2*np.arcsin(3*np.sqrt(2)*np.sin(OMEGA*dt/2))/dl
        predicted=[(-2*(ka-k0)*4+k0*dl)*180/np.pi,
                   (-(k1-ka-2*k0)*4-k0*dl)*180/np.pi]
        for panel,kind,pred in zip(axes[:2],('reflection','transmission'),predicted):
            exact=reference('halfspace')[0 if kind=='reflection' else 1]
            actual=np.angle(responses[name,'halfspace',kind]/exact)*180/np.pi
            line=panel.plot(F/1e6,actual,label=label+' FDTD')[0]
            panel.plot(F/1e6,pred,'--',color=line.get_color(),label=label+' 无拟合预测')
            prediction_residual[label+'_'+kind]=float(abs(actual-pred).max())
    axes[0].set(title='无损半空间：反射相位误差',xlabel='频率(MHz)',ylabel='相对连续解析解(度)')
    axes[1].set(title='无损半空间：透射相位误差',xlabel='频率(MHz)',ylabel='相对连续解析解(度)')
    labels=['半空间反射','半空间透射','平层反射','平层透射','电导反射','电导透射','Debye反射','Debye透射']
    x=np.arange(8)
    axes[2].bar(x-.18,[r['complex_relative_L2']*100 for r in coarse],.36,label='5cm')
    axes[2].bar(x+.18,[r['complex_relative_L2']*100 for r in fine],.36,label='1.25cm')
    axes[2].axhline(5,color='gray',ls='--',label='预声明5%诊断门限')
    axes[2].set(title='全部反射/透射复频响误差',xticks=x,xticklabels=labels,ylabel='相对L2(%)')
    axes[2].tick_params(axis='x',labelrotation=55)
    for ax in axes:
        ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('平层网格对照：宽域/窄域同网格差≤1e−14；细化降低误差，未拟合材料或时延')
    fig.savefig(out/'grid_comparison.png',dpi=140);plt.close(fig)
    result=dict(status='PASS_LOCAL_PLANE_GRID_DIAGNOSTIC',domain_change_complex_relative_L2=domain,
                coarse_metrics=coarse,fine_metrics=fine,prediction_maximum_phase_residual_deg=prediction_residual,
                identities=identities,code_sha256=sha256(__file__),
                prediction='1D lossless Yee dispersion plus half-cell dielectric-step position; parameters from grid only',
                scope='Normal-incidence plane check only; no historical UAV B-scan root-cause or 3D convergence certification')
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({'status':result['status'],'maximum_domain_effect':max(domain.values()),
                      'prediction_phase_residual_deg':prediction_residual},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
