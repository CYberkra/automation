"""Independent continuum 2D line-current/half-space Sommerfeld comparison.

Uses exact TE Fresnel coefficients, propagating and evanescent parts; no fit.
This is not a discrete-Yee/PML truth, nor an official gprMax 2D solver tool.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def reflection(f,c,nodes):
    mu=4*np.pi*1e-7;v=299792458.;eps0=1/(mu*v*v)
    omega=2*np.pi*f;k=omega/v
    mat=c['cover']
    er=mat['er']+mat['delta_er']/(1+1j*omega*mat['tau_s'])+mat['se']/(1j*omega*eps0)
    ks2=k*k*er
    z=np.polynomial.legendre.leggauss(nodes)
    theta=z[0]*np.pi/2;weights=z[1]*np.pi/2
    q=k[:,None]*np.sin(theta);ky=k[:,None]*np.cos(theta)
    ksoil=np.sqrt(ks2[:,None]-q*q+0j)
    ratio=(ky-ksoil)/(ky+ksoil)
    height=c['physical_source_logical_m'][1]+c['original_rx_logical_m'][1]-2*c['ground_y_m']
    offset=c['original_rx_logical_m'][0]-c['physical_source_logical_m'][0]
    prop=np.sum(ratio*np.exp(-1j*ky*height+1j*q*offset)*weights,axis=1)
    u=(z[0]+1)*1.5;uw=z[1]*1.5
    q=k[:,None]*np.cosh(u);ky=-1j*k[:,None]*np.sinh(u)
    ksoil=np.sqrt(ks2[:,None]-q*q+0j)
    ratio=(ky-ksoil)/(ky+ksoil)
    evan=np.sum(ratio*np.exp(-k[:,None]*height*np.sinh(u))*2*np.cos(q*offset)*uw,axis=1)
    return -mu*omega/(4*np.pi)*(prop+1j*evan)


def analyse(root):
    c=json.loads((root/'contract.json').read_text(encoding='utf-8'))
    d=np.load(root/'spectra_selected.npz');f=d['frequency_hz']
    reference=reflection(f,c,400)
    fine=reflection(f,c,800)
    quadrature=np.linalg.norm(fine-reference)/np.linalg.norm(fine)
    assert quadrature<1e-8
    native=d['full_cover_original']-d['full_air_original']
    tests={'full_native_500ns_scattered':native}
    eh=np.load(root/'eh_return_checks.npz')
    tests['full_plane_EH_return']=eh['full_space16_direct']
    tests['compact_plane_EH_return']=eh['compact_space16_direct']
    rows=[]
    for name,value in tests.items():
        low=f<=40e6
        rows.append(dict(path=name,relative_l2_to_continuum=float(np.linalg.norm(value-fine)/np.linalg.norm(fine)),
                         low20_40_relative_l2_to_continuum=float(np.linalg.norm((value-fine)[low])/np.linalg.norm(fine[low]))))
    (root/'continuum_flat_checks.json').write_text(json.dumps(dict(rows=rows,quadrature400_vs800_relative_l2=float(quadrature),
        assumptions=['Continuum infinite flat half-space','Ideal 2D z line current','No Yee dispersion/material-voxel placement/PML/finite-record truncation'],
        meaning='An independent asymptotic reference, not proof of a unique numerical error cause'),indent=2)+'\n',encoding='utf-8')
    np.savez_compressed(root/'continuum_flat_checks.npz',frequency_hz=f,continuum_scattered=fine,**tests)
    print(json.dumps(rows,indent=2));print('quadrature error',quadrature)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    analyse(p.parse_args().root)
