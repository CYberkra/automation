"""Independent input raster/actual material and raw acquisition audit; no solver."""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.out.exists(): raise ValueError('new output required')
    evidence=a.root/'artifacts/research_checks'
    capsules=('hs4_height_wavefield_pilot','hs4_height_wavefield_continuation','hs4_local_patch_controls','hs4_patch_refinement','hs4_patch_finest')
    rows=[]; unique_maps={}
    for name in capsules:
        folder=evidence/('2026-10-04_'+name)
        c=json.loads((folder/'execution_contract.json').read_text('utf-8'))
        for g in c['groups']:
            directory=Path(g['input']).parent
            raws=sorted(directory.glob('profile*.h5'))
            if not raws: continue  # Failed parser group/unattempted inputs are preserved.
            text=Path(g['input']).read_text('utf-8').splitlines()
            commands={}
            for line in text:
                if line.startswith('#'):
                    key,*values=line.split(); commands.setdefault(key,[]).append(values)
            spacing=np.array(list(map(float,commands['#dx_dy_dz:'][0])))
            domain=np.array(list(map(float,commands['#domain:'][0])))
            shape=np.rint(domain/spacing).astype(int)
            material=np.full(tuple(shape[::-1]),2,dtype=np.uint32)
            for words in commands['#box:']:
                limits=np.array(list(map(float,words[:6]))).reshape(2,3)/spacing
                if not np.allclose(limits,np.rint(limits),rtol=0,atol=1e-9): raise ValueError('off-grid box')
                lo,hi=np.rint(limits).astype(int)
                material[lo[2]:hi[2],lo[1]:hi[1],lo[0]:hi[0]]={'rock':3,'cover':4}[words[-1]]
            source=np.array(list(map(float,commands['#hertzian_dipole:'][0][1:4])))
            receiver=np.array(list(map(float,commands['#rx:'][0][:3])))
            srcstep=np.array(list(map(float,commands.get('#src_steps:',[['0','0','0']])[0])))
            rxstep=np.array(list(map(float,commands.get('#rx_steps:',[['0','0','0']])[0])))
            cfl=1/(299792458*np.sqrt(np.sum(1/spacing[[0,2]]**2)))
            pml=np.array(list(map(int,commands['#pml_cells:'][0])))
            if not np.allclose(pml[[0,2,3,5]]*spacing[[0,2,0,2]],[2,1,2,1],rtol=0,atol=1e-12): raise ValueError('physical PML differs')
            geometry_digest=hashlib.sha256(material.tobytes()).hexdigest()
            for raw in raws:
                stem=raw.stem; suffix=stem.removeprefix('profile'); index=int(suffix or '1')-1
                geom=directory/('hs4t2d_geom'+suffix+'.vtkhdf')
                with h5py.File(geom) as h:
                    actual=h['VTKHDF/CellData/Material'][:]
                    if not np.array_equal(actual,material): raise ValueError('actual material differs from independently rasterized input')
                with h5py.File(raw) as h:
                    if str(h.attrs['gprMax'])!='4.0.0' or not np.array_equal(h.attrs['dx_dy_dz'],spacing) or not np.isclose(float(h.attrs['dt']),cfl,rtol=1e-12,atol=0): raise ValueError('version/grid/CFL differs')
                    for key,position in [('srcs/src1',source+index*srcstep),('rxs/rx1',receiver+index*rxstep)]:
                        if not np.array_equal(h[key].attrs['GridPosition'],np.rint(position/spacing).astype(int)): raise ValueError('actual acquisition differs')
                    v=h['rxs/rx1/Ey'][:]
                    if v.dtype!=np.float64 or not np.isfinite(v).all() or len(v)!=int(h.attrs['Iterations']): raise ValueError('raw validity differs')
                rows.append({'capsule':name,'group':g['id'],'raw':raw.relative_to(a.root).as_posix(),'raw_sha256':sha256(raw),
                    'input_sha256':sha256(g['input']),'geometry':geom.relative_to(a.root).as_posix(),'geometry_sha256':sha256(geom),
                    'independent_material_array_sha256':geometry_digest,'native_spacing_m':spacing.tolist()})
            unique_maps[name+'/'+g['id']]=geometry_digest
    if len(rows)!=47: raise ValueError('unexpected trace count')
    a.out.parent.mkdir(parents=True,exist_ok=True)
    result={'status':'PASS','code_sha256':sha256(__file__),'successful_raw_traces':47,'trace_rows':rows,'material_maps':unique_maps,
        'scope':'Independently rasterized all actual geometries and checked main source/receiver/grid/time/float64. Parser failure/unattempted groups are excluded from successful count. Snapshots have a separate full audit. Duplicate maps may remain local with per-file hashes.'}
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('PASS:47 native traces and actual material maps independently audited')


if __name__=='__main__': main()
