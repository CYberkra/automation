"""Positive-Debye endpoint-span loss bound; no solver or field fitting."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

EPS0=8.8541878128e-12
C0=299792458.


def main(a):
    assert not a.out.exists()
    design=json.loads(a.design.read_text('utf-8'))['materials']['mudstone']
    span=design['target_band_endpoint_difference'];centre=design['target_center_epsilon_real']
    # Use ns for polynomial conditioning; x = tau_ns**2.
    w=2*np.pi*np.array([20,95,170])*1e-3
    wa,wc,wb=w;aa,cc,bb=w*w;A=aa+bb;B=aa*bb
    coefficients=[B*cc,3*B-A*cc,A-3*cc,-1]
    roots=np.roots(coefficients)
    positive=[float(r.real) for r in roots if abs(r.imag)<1e-10 and r.real>0]
    assert len(positive)==1
    def ratio(t):
        return wc*(1+A*t*t+B*t**4)/((bb-aa)*t*(1+cc*t*t))
    tau=np.sqrt(positive[0]);loss=span*ratio(tau)
    scalar=minimize_scalar(lambda logt:span*ratio(np.exp(logt)),method='bounded',
                           bounds=(np.log(.001),np.log(1e6)),options={'xatol':1e-13})
    assert abs(scalar.fun-loss)<1e-10 and abs(np.exp(scalar.x)-tau)<1e-5
    def parameters(t):
        support=1/(1+(w*t)**2)
        delta=span/(support[0]-support[2]);infinity=centre-delta*support[1]
        return infinity,delta
    infinity,delta=parameters(tau)
    assert infinity>=1 and delta>0
    samples=[]
    for label,t in [('fixed_tau',design['debye'][0]['tau_s']*1e9),('minimum_positive_Debye',tau)]:
        ei,de=parameters(t)
        for sigma in [.003,.0003,0]:
            pol=span*ratio(t);dc=sigma/(2*np.pi*95e6*EPS0);total=pol+dc
            alpha=-np.imag(2*np.pi*95e6/C0*np.sqrt(centre-1j*total))
            # independent real-valued square-root expression
            alt=2*np.pi*95e6/C0*np.sqrt((np.hypot(centre,total)-centre)/2)
            np.testing.assert_allclose(alpha,alt,rtol=1e-12,atol=0)
            samples.append(dict(label=label,tau_ns=float(t),epsilon_infinity=float(ei),delta_epsilon=float(de),
                                sigma_DC_S_m=sigma,epsilon_real95=centre,epsilon_imag_polarization95=float(pol),
                                epsilon_imag_DC95=float(dc),sigma_polarization95_S_m=float(pol*2*np.pi*95e6*EPS0),
                                sigma_effective95_S_m=float(total*2*np.pi*95e6*EPS0),
                                two_way_7m_field_absorption_dB=float(-20*alpha*14/np.log(10))))
    # Random positive mixtures: each pole's loss/span >= minimum ratio.
    rng=np.random.default_rng(20261008)
    min_slack=np.inf
    for _ in range(200):
        ts=np.exp(rng.uniform(np.log(.01),np.log(1e4),20))
        drops=rng.dirichlet(np.ones(20))*span
        value=np.sum(drops*np.array([ratio(t) for t in ts]))
        min_slack=min(min_slack,value-loss)
        assert value>=loss-1e-12
    result=dict(status='PASS_POSITIVE_DEBYE_SPAN_LOSS_BOUND_NOT_SITE_MEASUREMENT',
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                design_sha256=hashlib.sha256(a.design.read_bytes()).hexdigest(),
                endpoint_MHz=[20,170],centre_MHz=95,real_endpoint_drop=span,centre_epsilon_real=centre,
                cubic_coefficients_descending=coefficients,positive_tau_squared_ns2=positive,
                minimizing_tau_ns=float(tau),minimum_polarization_epsilon_imag95=float(loss),
                independent_scalar_minimum=float(scalar.fun),mixture_checks=200,minimum_mixture_slack=float(min_slack),
                cases=samples,
                proof='For each positive Debye pole, loss95=drop20to170*r(tau); all drops are nonnegative. Therefore sum(loss95)>=r_min*sum(drop)=r_min*3. r diverges at tau->0 and tau->infinity; its derivative numerator is the stated cubic in tau^2, which has exactly one positive root. Constant epsilon_infinity cancels from drop and loss; sigmaDC>=0 only adds loss. A positive single pole with epsilon_infinity>=1 attains this lower bound at epsilon_real95=12.',
                limits='Assumes a sum of ordinary Debye poles with positive strengths/times, nonnegative DC and mu_r=1. Not a bound on arbitrary causal materials, field detectability or discrete FDTD error. 7m is this controlled station, not whole-line thickness.')
    a.out.mkdir(parents=True)
    (a.out/'analysis.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
    ts=np.geomspace(.6,100,1500);values=span*np.array([ratio(t) for t in ts])
    axes[0].semilogx(ts,values,label='保持ε′95=12，ε′20−ε′170=3')
    axes[0].plot([tau],[loss],'ro',label=f'最低点τ={tau:.3f} ns')
    axes[0].axvline(8,color='gray',ls='--',label='当前τ=8 ns')
    axes[0].set(xlabel='单极弛豫时间 / ns',ylabel='95 MHz极化损耗 ε″',title='只换τ，无法任意降低损耗');axes[0].legend(fontsize=9)
    fixed=samples[:3];labels=['σDC=.003','σDC=.0003','σDC=0']
    axes[1].bar(labels,[q['epsilon_imag_polarization95'] for q in fixed],label='Debye极化损耗')
    axes[1].bar(labels,[q['epsilon_imag_DC95'] for q in fixed],bottom=[q['epsilon_imag_polarization95'] for q in fixed],label='DC损耗')
    for j,q in enumerate(fixed):axes[1].text(j,q['epsilon_imag_polarization95']+q['epsilon_imag_DC95']+.04,f"7m往返 {q['two_way_7m_field_absorption_dB']:.1f} dB",ha='center',fontsize=9)
    axes[1].set(ylabel='95 MHz总损耗 ε″',ylim=(0,2.15),title='当前配方：降DC后仍有极化损耗');axes[1].legend()
    fig.suptitle('泥岩研究假设的损耗约束 / 连续介质解析，无FDTD、无实测拟合\n正参数Debye和固定实部跨度下成立；吸收预算不等于实际SNR')
    fig.savefig(a.out/'mudstone_positive_debye_loss_bound.png',dpi=150);plt.close(fig)
    print(json.dumps(dict(minimum_epsilon_imag95=loss,minimizing_tau_ns=tau,cases=samples)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--design',type=Path,default=Path('configs/research/mudstone_design_v0_1.json'))
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
