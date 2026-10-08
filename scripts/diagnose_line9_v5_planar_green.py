"""Planar2D line-current angular-spectrum diagnostic; not a nonflat-field solver."""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.special import hankel2
from audit_line9_postprocessing import FREQ,inverse,weights
from review_line9_result_packages import indices,EPS0,C0
from hs_capsule_identity import sha256 as sha
from analyze_line9_basal_pair import plotting,save

MU0=1/(EPS0*C0*C0)


@lru_cache(None)
def nodes(order):
    return np.polynomial.legendre.leggauss(order)


def root_outgoing(value):
    q=np.sqrt(np.asarray(value,dtype=np.complex128))
    return np.where(q.imag>0,-q,q)


def integral(f,height_sum,separation,kernel,order=256,tail_exponent=36.):
    """exp(+iwt); ky outgoing Im<=0. Split light-line by sin/cosh substitutions.

    Reflected E/I_length=-omega*mu/(2pi*length)*integral_0^inf
    cos(q dx)*exp(-i ky0(hs+hr))*R(q)/ky0 dq.
    Returns that integral for all kernel columns, before physical prefactor.
    """
    assert height_sum>0 and separation>=0
    x,w=nodes(order);k0=2*np.pi*f/C0
    theta=(x+1)*np.pi/4;wt=w*np.pi/4
    umax=np.arcsinh(tail_exponent/(k0*height_sum));u=(x[None,:]+1)*umax[:,None]/2;wu=w[None,:]*umax[:,None]/2
    qp=k0[:,None]*np.sin(theta);kyp=k0[:,None]*np.cos(theta)
    qe=k0[:,None]*np.cosh(u);kye=-1j*k0[:,None]*np.sinh(u)
    fp=np.cos(qp*separation)*np.exp(-1j*kyp*height_sum)*wt
    fe=1j*np.cos(qe*separation)*np.exp(-1j*kye*height_sum)*wu
    return np.sum(fp[:,:,None]*kernel(qp,kyp),axis=1)+np.sum(fe[:,:,None]*kernel(qe,kye),axis=1)


def reflect_parts(q,ky0,k1,k2,k3,d1,d2):
    ky1=root_outgoing(k1[:,None]**2-q*q);ky2=root_outgoing(k2[:,None]**2-q*q);ky3=root_outgoing(k3[:,None]**2-q*q)
    r01=(ky0-ky1)/(ky0+ky1);r12=(ky1-ky2)/(ky1+ky2);r23=(ky2-ky3)/(ky2+ky3)
    p1=np.exp(-2j*ky1*d1);p2=np.exp(-2j*ky2*d2)
    full0=(r01+r12*p1)/(1+r01*r12*p1)
    sub=(r12+r23*p2)/(1+r12*r23*p2);full1=(r01+sub*p1)/(1+r01*sub*p1)
    once=(1-r01*r01)*r12*p1;twice=-(1-r01*r01)*r01*r12*r12*p1*p1
    return np.stack([r01,once,twice,full0-r01-once-twice,full1-full0],axis=-1)


def checks():
    # Constant reflection coefficients must reproduce the exact image line source.
    f=np.array([20e6,95e6,170e6]);height=16.;dx=1.3;k=2*np.pi*f/C0
    value=integral(f,height,dx,lambda q,ky:np.ones((*q.shape,1)),order=256)[:,0]
    expected=np.pi/2*hankel2(0,k*np.hypot(height,dx))
    err=float(np.linalg.norm(value-expected)/np.linalg.norm(expected));assert err<1e-10
    # Equal upper/lower media gives exactly zero reflection without fitting.
    q=k[:,None]*np.array([.1,.8,1.2])[None,:];ky=root_outgoing(k[:,None]**2-q*q)
    parts=reflect_parts(q,ky,k,k,k,7.3,7.)
    assert np.max(abs(parts))==0
    # At normal incidence, independent slab closed form fixes recurrence sign.
    k1=3*k;k2=2*k;k3=2*k;parts=reflect_parts(np.zeros((3,1)),k[:,None],k1,k2,k3,7.3,7.)[:,0,:]
    p=np.exp(-2j*k1*7.3);full=(-.5+.2*p)/(1-.1*p)
    err2=float(np.max(abs(parts[:,:4].sum(axis=1)-full)));assert err2<1e-12
    return dict(image_line_source_relative_L2=err,equal_media_zero_reflection=True,slab_closed_form_max_error=err2)


def main(a):
    assert not a.out.exists() and not a.numerical.exists(),'Fresh diagnostic required'
    m=json.loads((a.prepared/'manifest.json').read_text('utf-8'));db=a.prepared/'base_H1/geometries/line9_research_materials_v1_smoothed.json'
    assert sha(db)==m['material_sha256'];materials=json.loads(db.read_text('utf-8'))['materials'];g=m['geometry_diagnostic']
    assert len(g['boundaries'])==3 and [b['above'] for b in g['boundaries']]==[0,1,2]
    n=np.array([indices(materials,f) for f in FREQ]);k=2*np.pi*FREQ[:,None]/C0*n
    tx=m['groups'][0]['tx_m'];rx=m['groups'][0]['rx_m'];ys=g['surface_y_m'];hs=tx[1]-ys;hr=rx[1]-ys;dx=abs(rx[0]-tx[0]);d1=g['cover_base']['depth_m'];d2=g['basal_sand']['depth_m']-d1
    kernel=lambda q,ky:reflect_parts(q,ky,k[:,1],k[:,2],k[:,3],d1,d2)
    prefactor=-(2*np.pi*FREQ)*MU0/(2*np.pi*.025)
    coarse=integral(FREQ,hs+hr,dx,kernel,order=256)*prefactor[:,None]
    refined=integral(FREQ,hs+hr,dx,kernel,order=512)*prefactor[:,None]
    convergence=np.linalg.norm(refined-coarse,axis=0)/np.linalg.norm(refined,axis=0)
    assert np.max(convergence)<1e-7
    direct=-(2*np.pi*FREQ)*MU0/(4*.025)*hankel2(0,k[:,0]*np.hypot(dx,hs-hr))
    names=['direct_air','surface','cover_once','cover_twice','cover_higher','bottom_contrast'];response=np.column_stack([direct,refined]);metrics={}
    plt=plotting();fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for ax,window in zip(axes,['hann','blackman']):
        z,t=inverse(response,FREQ,weights(window,501));show=(t*1e9>=300)&(t*1e9<=450);gate=(t*1e9>=m['basal_gate_ns'][0])&(t*1e9<=m['basal_gate_ns'][1]);total0=z[:,:5].sum(axis=1)
        rows={}
        for j,label in enumerate(['空气直达有限带尾部','地表单次','覆盖层底单次','覆盖层二次往返','覆盖层更高次','底砂材料对比']):
            ax.semilogy(t[show]*1e9,abs(z[show,j]),lw=1,label=label)
            rows[names[j]]=dict(gate_L2_over_bottom=float(np.linalg.norm(z[gate,j])/np.linalg.norm(z[gate,5])),gate_peak_ns=float(t[np.flatnonzero(gate)[np.argmax(abs(z[gate,j]))]]*1e9))
        ax.axvspan(*m['basal_gate_ns'],color='red',alpha=.08);ax.set(title=window+'：平面二维解析诊断，尚非原非平场验证',xlabel='SFCW时间 / ns',ylabel='解析包络 / (V/m)/(A·m)');ax.legend(fontsize=8)
        rows['total_H0_over_bottom']=float(np.linalg.norm(total0[gate])/np.linalg.norm(z[gate,5]));metrics[window]=rows
    a.out.mkdir(parents=True);fig.savefig(a.out/'planar_v5_components.png',dpi=140);plt.close(fig)
    with h5py.File(a.numerical,'x') as h:
        h.create_dataset('frequency_Hz',data=FREQ);h.create_dataset('response',data=response);h.attrs['columns']=','.join(names)
    save(a.out/'analysis.json',dict(status='ANALYTIC_PLANAR_2D_COMPONENT_DIAGNOSTIC_NOT_NONFLAT_FDTD_CAUSAL_CERTIFICATION',script_sha256=sha(__file__),material_sha256=sha(db),prepared_manifest_sha256=sha(a.prepared/'manifest.json'),numerical_sha256=sha(a.numerical),
        tests=checks(),quadrature_256_512_relative_L2_by_reflected_component=convergence.tolist(),metrics=metrics,source_heights_m=[hs,hr],horizontal_offset_m=dx,cover_thickness_m=d1,mud_thickness_m=d2,
        convention='exp(+iwt),outgoingImky<=0,scalarEz invariantz; TE relative to y-normal interfaces; ideal line current, not3Dpoint dipole/portS21. NativeSFCWunits E/(I*.025m).',
        method='2D Fourier/angular spectrum derived from Helmholtz source equation and Ez/dEzdy continuity. Propagating q=k0sin(theta),evanescent q=k0cosh(u). No ray-only or single-angle approximation; integration tail exponent36.',
        source_reading='Chew Purdue Lecture35 sections35.1-35.2 angular-spectrum method; Lambot2007 sections1-2 global reflection recursion. Both sources discuss3D; scalar2D specialization derived here, not claimed copied from them.',
        sources=['https://engineering.purdue.edu/wcchew/ece604f20/Lecture%20Notes/Lect35.pdf','https://doi.org/10.1029/2007GL031459'],
        limits='Planar local-column medium, infinite horizontal extent and exact constitutive continuum, not the actual nonflat voxel model. Quadrature checks do not prove FDTD agreement, material/site accuracy or weak-target error floor. Cover_twice is an analytic path term, not an observed event label.'))
    print(json.dumps(dict(tests=checks(),quadrature_relative_max=float(max(convergence)),metrics=metrics)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['prepared','out','numerical']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
