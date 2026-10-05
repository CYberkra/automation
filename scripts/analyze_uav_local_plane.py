"""Independent continuous Fresnel/slab diagnostic, without amplitude/phase fitting."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256
from gprMax.toolboxes.SFCW.processing import load_receiver, engineering_dft

C = 299792458.
EPS = 8.8541878188e-12
F = np.linspace(20e6,170e6,501)
OMEGA = 2*np.pi*F


def reference(case):
    er = np.full(501,9.,dtype=complex)
    if case in ('conductive','debye'):
        er -= 1j*.001/(OMEGA*EPS)
    if case == 'debye':
        er += 1/(1+1j*OMEGA*6e-9)
    n1 = np.sqrt(er); n2 = 2.
    k0 = OMEGA/C; k1 = k0*n1; k2 = k0*n2
    r01 = (1-n1)/(1+n1); r12 = (n1-n2)/(n1+n2)
    t01 = 2/(1+n1); t12 = 2*n1/(n1+n2)
    if case == 'halfspace':
        reflection = r01
        transmission = t01*np.exp(-1j*(k1-k0)*4)
    else:
        propagation = np.exp(-2j*k1*3)
        reflection = (r01+r12*propagation)/(1+r01*r12*propagation)
        transmission = t01*t12*np.exp(-1j*(k1*3+k2-k0*4))/(1+r01*r12*propagation)
    return reflection*np.exp(-2j*k0*4), transmission


def spectrum(path, i):
    s = load_receiver(path, receiver_path=f'/rxs/rx{i}',component='Ey')
    values = s.samples.copy()
    n = round(20e-9/s.dt)
    values[-n:] *= .5*(1+np.cos(np.linspace(0,np.pi,n)))
    result = engineering_dft(values,s.dt,F,time_offset=s.time_offset)
    take = np.array([0,83,167,250,333,417,500])
    independent = s.dt*np.exp(-2j*np.pi*F[take,None]*s.times[None,:])@values
    difference = np.linalg.norm(result[take]-independent)/np.linalg.norm(independent)
    if difference > 1e-9:
        raise ValueError('independent actual-clock DFT mismatch')
    return result,s,float(difference)


def main(study,out):
    if out.exists():
        raise ValueError('fresh output required')
    cp = study/'execution_contract.json'
    c = json.loads(cp.read_text('utf-8'))
    completed = json.loads((study/'completed_verification.json').read_text('utf-8'))
    if completed['status'] != 'PASS' or completed['contract_sha256'] != sha256(cp):
        raise ValueError('completed identity required')
    spectra={}; waveforms={}; dfts={}; transverse={}; raw_identities={}
    for g,verified in zip(c['groups'],completed['groups']):
        raw = Path(g['input']).with_suffix('.h5')
        if sha256(raw) != verified['raw_sha256']:
            raise ValueError('raw changed')
        raw_identities[g['id']] = sha256(raw)
        for i in (1,2,3):
            spectra[g['id'],i],waveforms[g['id'],i],dfts[g['id'],i] = spectrum(raw,i)
        a,b = spectra[g['id'],1],spectra[g['id'],2]
        transverse[g['id']] = float(np.linalg.norm(a-b)/np.linalg.norm(a))
        with h5py.File(raw) as h:
            unwanted = max(np.max(abs(h['rxs/rx1/'+component][:])) for component in ('Ex','Ez'))
            if unwanted > 1e-8*np.max(abs(h['rxs/rx1/Ey'][:])):
                raise ValueError('unexpected polarisation')
    metrics=[]; products={}
    for case in ('halfspace','slab','conductive','debye'):
        reflection = (spectra[case,1]-spectra['air',1])/spectra['air',1]
        transmission = spectra[case,3]/spectra['air',3]
        for name,value,exact in zip(('reflection','transmission'),(reflection,transmission),reference(case)):
            err = float(np.linalg.norm(value-exact)/np.linalg.norm(exact))
            phase = float(abs(np.angle(value/exact)*180/np.pi).max())
            metrics.append(dict(case=case,observation=name,complex_relative_L2=err,
                                maximum_phase_error_deg=phase,
                                passes_declared_diagnostic=bool(err<=.05 and phase<=5)))
            products[case,name]=(value,exact)
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei']
    plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    names=dict(halfspace='无损半空间',slab='无损三米平层',conductive='平层加电导',debye='平层加电导与Debye')
    for case in names:
        r,ex=products[case,'reflection']; t,tx=products[case,'transmission']
        line=ax[0,0].plot(F/1e6,abs(r),label=names[case])[0]
        ax[0,0].plot(F/1e6,abs(ex),'--',color=line.get_color())
        ax[0,1].plot(F/1e6,np.angle(t/tx)*180/np.pi,label=names[case])
        ax[1,0].plot(F/1e6,20*np.log10(abs(t)),label=names[case])
        ax[1,0].plot(F/1e6,20*np.log10(abs(tx)),'--',color=line.get_color())
        s=waveforms[case,1]; air=waveforms['air',1]
        ax[1,1].plot(s.times*1e9,s.samples-air.samples,label=names[case])
    ax[0,0].set(title='上方反射/空气入射：实线FDTD，虚线解析',xlabel='频率(MHz)',ylabel='幅度比')
    ax[0,1].set(title='下方透射相位误差；不拟合时延',xlabel='频率(MHz)',ylabel='度')
    ax[1,0].set(title='下方透射/同点空气：实线FDTD，虚线解析',xlabel='频率(MHz)',ylabel='幅度比(dB)')
    ax[1,1].set(title='上方散射场 A-scan；单站，非B-scan',xlabel='接收时间(ns)',ylabel='Ey(V/m)',xlim=(0,180))
    for panel in ax.flat:
        panel.legend(fontsize=7); panel.grid(alpha=.2)
    fig.suptitle('第二轮辅助：正入射平层解析自检；不能认证悬空点源或端口S21')
    fig.savefig(out/'plane_comparison.png',dpi=140);plt.close(fig)
    result=dict(status='PASS' if all(r['passes_declared_diagnostic'] for r in metrics) else 'DIAGNOSTIC_FAILED',
                metrics=metrics,transverse_uniformity_relative_L2=transverse,
                maximum_independent_DFT_error=max(dfts.values()),raw_identities=raw_identities,
                contract_sha256=sha256(cp),analysis_code_sha256=sha256(__file__),
                scope='Normal incidence only; not suspended source, cross-track layout, antenna or field validation',
                fit_parameters=False)
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();main(a.study,a.out)
