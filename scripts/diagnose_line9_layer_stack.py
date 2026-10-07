"""Local 1D complete stack reflection proxy on existing Line9 models; no FDTD."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_line9_postprocessing import FREQ, inverse, weights, sha, save
from review_line9_result_packages import C0, indices
from diagnose_line9_layer_kinematics import primary_paths
from trace_line9_time_origin import correlation


def stack_response(n, thickness, frequency):
    """Return full and single-reflection spectra at the top of the first layer.

    n: [frequency, layer], including upper/lower half spaces. thickness[j]
    is the thickness of layer j; first entry is source height, last unused.
    exp(-i*k*d) with Im(n)<=0 is the passive engineering convention.
    """
    n=np.asarray(n,dtype=np.complex128)
    r=(n[:,:-1]-n[:,1:])/(n[:,:-1]+n[:,1:])
    full=r[:,-1].copy()
    for j in range(n.shape[1]-3,-1,-1):
        roundtrip=np.exp(-4j*np.pi*frequency*n[:,j+1]*thickness[j+1]/C0)
        full=(r[:,j]+full*roundtrip)/(1+r[:,j]*full*roundtrip)
    air=np.exp(-4j*np.pi*frequency*n[:,0]*thickness[0]/C0)
    full*=air
    single=np.zeros_like(full); transmission=np.ones_like(full);distance=np.zeros_like(full)
    for j in range(n.shape[1]-1):
        distance+=n[:,j]*thickness[j]
        single+=transmission*r[:,j]*np.exp(-4j*np.pi*frequency*distance/C0)
        transmission*=1-r[:,j]**2
    return full,single


def admittance_response(n,thickness,frequency):
    """Independent transmission-line input-admittance formula."""
    y=n[:,-1].copy()
    for j in range(n.shape[1]-2,0,-1):
        tangent=np.tan(2*np.pi*frequency*n[:,j]*thickness[j]/C0)
        y=n[:,j]*(y+1j*n[:,j]*tangent)/(n[:,j]+1j*y*tangent)
    return (n[:,0]-y)/(n[:,0]+y)*np.exp(-4j*np.pi*frequency*n[:,0]*thickness[0]/C0)


def checks():
    n=np.broadcast_to(np.array([1.,3.,2.]),(len(FREQ),3))
    thickness=np.array([2.,7.,0.]);full,single=stack_response(n,thickness,FREQ)
    z=np.exp(-4j*np.pi*FREQ*3*7/C0)
    air=np.exp(-4j*np.pi*FREQ*2/C0)
    # Independent slab expression and its single-reflection approximation.
    expected=(-.5+.2*z)/(1-.1*z)*air
    expected_single=(-.5+.75*.2*z)*air
    err=float(np.max(abs(full-expected)));singleerr=float(np.max(abs(single-expected_single)))
    assert err<1e-12 and singleerr<1e-12
    assert max(abs(full))<=1
    zero,_=stack_response(np.ones((len(FREQ),3)),thickness,FREQ)
    assert np.max(abs(zero))==0
    half,_=stack_response(n[:,:2],np.array([2.,0.]),FREQ)
    halferr=float(np.max(abs(half+.5*air)))
    assert halferr<1e-12
    return dict(slab_complex_max_error=err,single_reflection_max_error=singleerr,
                equal_media_max_reflection=float(np.max(abs(zero))),halfspace_max_error=halferr)


def main(root,review,cache,out):
    if out.exists():raise ValueError('Fresh directory required')
    result=dict(calls_solver=False,calls_training=False,script_sha256=sha(__file__),self_checks=checks(),packages=[])
    for ap in sorted(review.glob('*_audit.json')):
        audit=json.loads(ap.read_text('utf-8'));rows=audit['records']
        mf=next((root/audit['package']/'geometries').glob('*.json'))
        assert sha(mf)==audit['materials_sha256']
        materials=json.loads(mf.read_text('utf-8'))['materials']
        n=np.array([indices(materials,f) for f in FREQ])
        full=[];single=[];paths=[];independent_errors=[]
        for row in rows:
            g=row['geometry'];bounds=g['boundaries']
            media=[b['above'] for b in bounds]+[bounds[-1]['below']]
            thickness=np.r_[g['midpoint_agl_m'],np.diff([b['depth_m'] for b in bounds]),0.]
            hf,hs=stack_response(n[:,media],thickness,FREQ)
            hi=admittance_response(n[:,media],thickness,FREQ)
            independent_errors.append(float(np.linalg.norm(hi-hf)/np.linalg.norm(hf)))
            # One common bistatic surface phase correction; both models use it.
            correction=g['surface']['time95_ns']*1e-9-2*thickness[0]/C0
            phase=np.exp(-2j*np.pi*FREQ*correction)
            full.append(hf*phase);single.append(hs*phase)
            paths.append(primary_paths(g,materials))
        full=np.column_stack(full);single=np.column_stack(single)
        cp=cache/(audit['package']+'.npz')
        with np.load(cp) as z:
            np.testing.assert_array_equal(z['frequency_Hz'],FREQ)
            np.testing.assert_array_equal(z['ids'],[r['id'] for r in rows])
            actual=z['response'].copy()
        assert max(independent_errors)<1e-10
        item=dict(package=audit['package'],audit_sha256=sha(ap),cache_sha256=sha(cp),interfaces={},
                  independent_admittance_relative_L2_max=max(independent_errors))
        for window in ['hann','blackman']:
            w=weights(window,501);vf,t=inverse(full,FREQ,w);vs,_=inverse(single,FREQ,w)
            va,_=inverse(actual,FREQ,w)
            for key in ['cover_base','first_sand','basal_sand']:
                hp=np.column_stack([p[r['geometry']['boundaries'].index(r['geometry'][key])]
                                     for p,r in zip(paths,rows)])
                vp,_=inverse(hp,FREQ,w)
                centers=t[np.argmax(abs(vp),axis=0)]*1e9
                mask=abs(t[:,None]*1e9-centers[None,:])<=12
                item['interfaces'][key+'_'+window]=dict(
                    single_target_primary_complex_correlation=correlation(vp[mask],va[mask]),
                    sum_all_primaries_complex_correlation=correlation(vs[mask],va[mask]),
                    full_stack_complex_correlation=correlation(vf[mask],va[mask]),
                    multiple_only_proxy_complex_correlation=correlation((vf-vs)[mask],va[mask]),
                    multiple_over_full_proxy_norm=float(np.linalg.norm((vf-vs)[mask])/np.linalg.norm(vf[mask])),
                    gate='Existing primary-predicted +-12ns, not fitted to observations')
                if key=='basal_sand':
                    split=len(rows)//2
                    train=mask.copy();train[:,split:]=False
                    test=mask.copy();test[:,:split]=False
                    item['interfaces'][key+'_'+window]['first_half_fit_second_half_test']={}
                    for name,model in [('sum_all_primaries',vs),('full_stack',vf)]:
                        coefficient=np.vdot(model[train],va[train])/np.vdot(model[train],model[train])
                        prediction=model*coefficient
                        item['interfaces'][key+'_'+window]['first_half_fit_second_half_test'][name]=dict(
                            global_coefficient=[float(coefficient.real),float(coefficient.imag)],
                            complex_correlation=correlation(prediction[test],va[test]),
                            relative_L2=float(np.linalg.norm(prediction[test]-va[test])/np.linalg.norm(va[test])))
                    if window=='hann':
                        plot_data=(va,vs,vf,t,centers,train)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
        plt.rcParams['axes.unicode_minus']=False
        va,vs,vf,t,centers,train=plot_data
        x=np.array([r['chainage_m'] for r in rows]);order=np.argsort(x)
        keep=(t*1e9>=min(centers)-40)&(t*1e9<=max(centers)+40)
        display=[va]
        for model in [vs,vf]:
            coef=np.vdot(model[train],va[train])/np.vdot(model[train],model[train])
            display.append(model*coef)
        ref=float(np.max(abs(va)))
        fig,axs=plt.subplots(1,3,figsize=(18,6),layout='constrained')
        for ax,value,title in zip(axs,display,['真实Hann总场：未改动','全部一次反射近似：前半估1个复系数','完整层间多次近似：前半估1个复系数']):
            db=20*np.log10(np.maximum(abs(value[keep][:,order])/ref,1e-12))
            pic=ax.pcolormesh(x[order],t[keep]*1e9,db,cmap='gray_r',vmin=-110,vmax=-45,shading='nearest',rasterized=True)
            ax.plot(x[order],centers[order],color='#d55e00',lw=1,label='底砂一次反射近似峰位')
            ax.axvline((x[len(rows)//2-1]+x[len(rows)//2])/2,color='#0072b2',ls='--',label='拟合/未拟合站位分界')
            ax.set_title(title,fontsize=11);ax.set_xlabel('剖面横坐标 X / m');ax.set_ylabel('时间 / ns')
            ax.set_ylim(float(t[keep][-1]*1e9),float(t[keep][0]*1e9));ax.invert_xaxis();ax.legend(fontsize=8)
        fig.colorbar(pic,ax=axs,label='共同实际Hann全时峰值参考 / dB；无逐道归一化')
        fig.suptitle(audit['package']+'：右两列是1D理论对照，不是修复后的数据；缺二维扩散/绕射/边界',fontsize=12)
        out.mkdir(parents=True,exist_ok=True)
        fig.savefig(out/(audit['package']+'_stack_comparison.png'),dpi=160);plt.close(fig)
        result['packages'].append(item)
        print(json.dumps(item,ensure_ascii=False),flush=True)
    result['limits']=('All internal reflections and transmission factors for a local 1D normal-incidence stack, '
        'not the 2D Green function. No lateral scattering, slope, diffraction, spreading, air/PML echoes, '
        'numerical precision or target-ablation certificate. A model-only norm ratio is not an observed energy fraction. '
        'No removal of any observed component, model-guided filtering, imaging or additional FDTD.')
    save(out/'layer_stack.json',result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','review','cache','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();main(a.root,a.review,a.cache,a.out)
