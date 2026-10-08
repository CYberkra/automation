"""Independent scalar Debye loss minimization and exact rational Sturm root count."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path


def trim(p):
    while len(p)>1 and not p[-1]:p.pop()
    return p


def remainder(a,b):
    a=list(a)
    while len(a)>=len(b) and any(a):
        shift=len(a)-len(b);factor=a[-1]/b[-1]
        for j,v in enumerate(b):a[shift+j]-=factor*v
        trim(a)
    return a


def main(a):
    assert not a.out.exists()
    d=json.loads(a.analysis.read_text('utf-8'))
    # u=omega_c*tau; all frequency ratios are exact rational numbers.
    aa=Q(20,95)**2;bb=Q(170,95)**2;A=aa+bb;B=aa*bb
    p=[Q(-1),A-3,3*B-A,B]
    sequence=[p,[p[1],2*p[2],3*p[3]]]
    while len(sequence[-1])>1:
        r=remainder(sequence[-2],sequence[-1]);assert any(r)
        sequence.append([-v for v in r])
    def variations(signs):return sum(x!=y for x,y in zip(signs,signs[1:]))
    zero=[1 if p[0]>0 else -1 for p in sequence if p[0]]
    infinity=[1 if p[-1]>0 else -1 for p in sequence]
    count=variations(zero)-variations(infinity);assert count==1
    def loss(log_tau):
        tau=math.exp(log_tau)*1e-9
        xl=2*math.pi*20e6*tau;xc=2*math.pi*95e6*tau;xh=2*math.pi*170e6*tau
        drop=1/(1+xl*xl)-1/(1+xh*xh)
        delta=3/drop
        return delta*xc/(1+xc*xc)
    lo,hi=math.log(.01),math.log(1e4);g=(math.sqrt(5)-1)/2
    x1=hi-g*(hi-lo);x2=lo+g*(hi-lo)
    for _ in range(150):
        if loss(x1)<loss(x2):hi=x2;x2=x1;x1=hi-g*(hi-lo)
        else:lo=x1;x1=x2;x2=lo+g*(hi-lo)
    minimum=loss((lo+hi)/2);tau=math.exp((lo+hi)/2)
    assert abs(minimum-d['minimum_polarization_epsilon_imag95'])<1e-12
    assert abs(tau-d['minimizing_tau_ns'])<1e-5
    errors=[]
    for case in d['cases']:
        pol=loss(math.log(case['tau_ns']))
        dc=case['sigma_DC_S_m']/(2*math.pi*95e6*8.8541878128e-12)
        total=pol+dc;eps=12
        alpha=2*math.pi*95e6/299792458*math.sqrt((math.sqrt(eps*eps+total*total)-eps)/2)
        db=-20*alpha*14/math.log(10)
        assert abs(db-case['two_way_7m_field_absorption_dB'])<1e-9
        errors.append(abs(db-case['two_way_7m_field_absorption_dB']))
    result=dict(status='PASS_INDEPENDENT_SCALAR_AND_EXACT_STURM_POSITIVE_DEBYE_BOUND',
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                analysis_sha256=hashlib.sha256(a.analysis.read_bytes()).hexdigest(),
                exact_positive_stationary_root_count=count,sturm_zero_signs=zero,sturm_infinity_signs=infinity,
                scalar_minimizing_tau_ns=tau,scalar_minimum_epsilon_imag95=minimum,
                maximum_absorption_recompute_error_dB=max(errors),
                scope='Mathematical positive-Debye model class only; no physical sample or FDTD certification.')
    a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
