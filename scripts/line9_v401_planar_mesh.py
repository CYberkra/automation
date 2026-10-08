"""Prepare a narrow planar grid control; never alter nonflat production models.

The 80 m coarse pair must pass a declared comparison with the completed 210 m
pair before a separate 1.25 cm pair is frozen. Fixed physical PML = 2 m.
"""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import save


def main(a):
    assert not a.out.exists()
    old=json.loads((a.package/'manifest.json').read_text('utf-8'))
    records={g['id']:g for g in old['groups']}
    dx=a.spacing
    assert dx in [.025,.0125] and a.width==80
    if dx==.0125:
        assert a.gate is not None
        gate=json.loads(a.gate.read_text('utf-8'))
        assert gate['status']=='PASS_NARROW_PLANAR_DOMAIN_GATE'
        assert gate['coarse_spacing_m']==.025 and gate['width_m']==80
    a.out.mkdir(parents=True);groups=[]
    for tag in ['H0','H1']:
        ref=records['flat_ricker_'+tag]
        for key in ['input','geometry','material']:assert sha(a.package/ref[key])==ref[key+'_sha256']
        gp=a.out/('planar_'+tag)/'geometries';gp.mkdir(parents=True)
        with h5py.File(a.package/ref['geometry']) as src:
            col=src['data'][0,:,0]
            # A flat input is checked independently at both ends and midpoint.
            for k in [src['data'].shape[0]//2,src['data'].shape[0]-1]:np.testing.assert_array_equal(src['data'][k,:,0],col)
            factor=round(.025/dx);v=np.repeat(np.repeat(col,factor)[None,:,None],round(a.width/dx),axis=0)
            with h5py.File(gp/'full2d_compact.h5','x') as dst:
                for key,value in src.attrs.items():dst.attrs[key]=value
                dst.attrs['dx_dy_dz']=[dx]*3;dst.attrs['shape_nxyz']=v.shape
                dst.attrs['MeshControl']='Planar screening only, repeated physical interfaces, not a nonflat production domain replacement'
                src.copy('material_keys',dst)
                dst.create_dataset('data',data=v,compression='gzip',compression_opts=4)
        shutil.copyfile(a.package/ref['material'],gp/Path(ref['material']).name)
        inp=a.out/('planar_'+tag)/'cases/full2d_r0098/profile.in';inp.parent.mkdir(parents=True)
        lines=[]
        for line in (a.package/ref['input']).read_text('utf-8').splitlines():
            if line.startswith('#title:'):line=f'#title: planar {tag}; width80; mesh{dx}; fixed line current and 2m PML'
            elif line.startswith('#domain:'):line=f'#domain: 80 42.5 inf'
            elif line.startswith('#dx_dy_dz:'):line=f'#dx_dy_dz: {dx} {dx} {dx}'
            elif line.startswith('#pml_cells:'):line=f'#pml_cells: {round(2/dx)} {round(2/dx)} 0 {round(2/dx)} {round(2/dx)} 0'
            elif line.startswith('#hertzian_dipole:'):line=f'#hertzian_dipole: z 39.35 39.375 {dx/2} pulse'
            elif line.startswith('#rx:'):line=f'#rx: 40.65 39.425 {dx/2} full2d_r0098_rx1 Ez'
            lines.append(line)
        inp.write_text('\n'.join(lines)+'\n',encoding='utf-8')
        row=dict(ref,id='planar_'+tag,input=inp.relative_to(a.out).as_posix(),input_sha256=sha(inp),
                 geometry=(gp/'full2d_compact.h5').relative_to(a.out).as_posix(),geometry_sha256=sha(gp/'full2d_compact.h5'),
                 material=(gp/Path(ref['material']).name).relative_to(a.out).as_posix(),
                 native_shape=list(v.shape),spacing_m=[dx]*3,tx_m=[39.35,39.375,0],rx_m=[40.65,39.425,0],
                 expected_samples=int(np.ceil(1.2e-6/(dx/(299792458*np.sqrt(2)))))+1,
                 dt_s=float(dx/(299792458*np.sqrt(2))))
        groups.append(row)
    manifest=dict(old,groups=groups,dt_s=groups[0]['dt_s'],expected_samples=groups[0]['expected_samples'],
                  generator_sha256=sha(__file__),
                  width_m=80,spacing_m=dx,physical_pml_m=2,parent_manifest_sha256=sha(a.package/'manifest.json'),
                  preceding_domain_gate_sha256=sha(a.gate) if a.gate else None,
                  approval_basis='User authorizes 4.0.1 setup and continued research; bounded staged planar mesh test.',
                  processing='Exact501tones, preserved phase, Hann/Blackman, no AGC/taper/fit. Compare E/I_line; for previous display units divide ALL grids by common .025m reference, not actual dz.',
                  source_scaling_basis='z-Hertzian injection I*dz/(dx*dy*dz)=I/(dx*dy); integrated 2D line current remains I. Actual source SpatialScale changes with dz, so E/(I*dz) alone is not grid-invariant.',
                  domain_gate=dict(basal_delta_relative_L2_max=.005,basal_H0_relative_L2_max=.005,full_band_delta_relative_L2_max=.01,
                                   windows=['hann','blackman'],no_fit=True),
                  fine_prediction='Prior parameter-only Yee prediction: phase at95MHz about-3.156deg,170MHz about-17.938deg relative to continuum. No output fitting.',
                  limits='Infinite-planar mechanism/grid study only; not shrinking nonflat production geometry, field calibration or whole-line validation.')
    save(a.out/'manifest.json',manifest)
    print(json.dumps(dict(groups=2,shape=groups[0]['native_shape'],spacing=dx,samples=groups[0]['expected_samples'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','out','gate']:p.add_argument('--'+key,type=Path,required=key!='gate')
    p.add_argument('--spacing',type=float,required=True);p.add_argument('--width',type=float,default=80)
    main(p.parse_args())
