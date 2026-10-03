"""Freeze minimal paired controls for HS4T2D geometry/B-scan attribution."""
import argparse
import json
from pathlib import Path

from hs_capsule_identity import sha256


ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_design'


def create_input(source, width, height, dz, halfspace, traces):
    shift=(width-12)/2
    lines=source.read_text('utf-8').splitlines()
    cover=[line.split() for line in lines if line.startswith('#box:') and line.split()[-1]=='cover']
    result=[]
    for line in lines:
        words=line.split()
        if line.startswith('#domain:'):
            line=f'#domain: {width:g} 0.05 33'
        elif line.startswith('#dx_dy_dz:'):
            line=f'#dx_dy_dz: 0.05 0.05 {dz:g}'
        elif line.startswith('#pml_cells:'):
            line=f'#pml_cells: 20 0 {round(1/dz)} 20 0 {round(1/dz)}'
        elif line.startswith(('#hertzian_dipole:', '#rx:')):
            start=2 if words[0]=='#hertzian_dipole:' else 1
            words[start]=f'{float(words[start])+shift:.12g}'
            words[start+2]=f'{12+height:g}'
            line=' '.join(words)
        elif line.startswith('#box:'):
            if words[-1]=='rock':
                line=f'#box: 0 0 0 {width:g} 0.05 12 rock'
            elif halfspace:
                continue
            else:
                words[1]=f'{float(words[1])+shift:.12g}'
                words[4]=f'{float(words[4])+shift:.12g}'
                line=' '.join(words)
        elif line.startswith('#geometry_view:'):
            if halfspace:
                result.append(f'#box: 0 0 0 {width:g} 0.05 12 cover')
            elif shift:
                result.extend([f'#box: 0 0 {cover[0][3]} {shift:g} 0.05 12 cover',
                               f'#box: {shift+12:g} 0 {cover[-1][3]} {width:g} 0.05 12 cover'])
            line=f'#geometry_view: 0 0 0 {width:g} 0.05 33 0.05 0.05 {dz:g} hs4t2d_geom n'
        result.append(line)
    if traces>1:
        result.extend(['#src_steps: 0.5 0 0','#rx_steps: 0.5 0 0'])
    return '\n'.join(result)+'\n'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new study directory required')
    base=json.loads((DESIGN/'execution_contract_cuda_cached_att4.json').read_text('utf-8'))
    specs=[('halfspace12',12,15,.05,True,13),
           ('wide36_rough',36,15,.05,False,13),('wide36_halfspace',36,15,.05,True,13),
           ('low2_rough',12,2,.05,False,13),('low2_halfspace',12,2,.05,True,13),
           ('zfine_rough',12,15,.025,False,1),('zfine_halfspace',12,15,.025,True,1)]
    args.out.mkdir(parents=True)
    groups=[]
    for name,width,height,dz,halfspace,traces in specs:
        source=DESIGN/'relief2p0'/('hs4t2d_t61.in' if traces==1 else 'hs4t2d_t01.in')
        folder=args.out/name
        folder.mkdir()
        target=folder/'profile.in'
        target.write_text(create_input(source,width,height,dz,halfspace,traces),encoding='utf-8')
        groups.append({'id':name,'input':str(target.resolve()),'sha256':sha256(target),
                       'source_input':str(source.relative_to(ROOT)),'source_sha256':sha256(source),
                       'width_m':width,'x_translation_m':(width-12)/2,'height_m':height,
                       'dz_m':dz,'halfspace':halfspace,'traces':traces,
                       'old_station_indices':[61] if traces==1 else list(range(1,122,10)),
                       'minimum_free_VRAM_GiB':3.5 if width==36 else 2.4 if dz==.025 else 1.4})
    c={'status':'FROZEN_APPROVED','approval_basis':'User 2026-10-03 requested to attempt to identify the cause; existing autonomous GPU research authorization. Minimal controlled development experiment, frozen before execution.',
       'purpose':'Separate reference/processing, lateral truncation, standoff and vertical-grid sensitivity; no training or physical acceptance',
       'backend':'CUDA','gpu_device':0,'precision':'double','python':base['python'],
       'max_runs':67,'max_wall_s':1800,'min_available_RAM_GiB':1,'groups':groups,
       'source_identities':base['source_identities'],
       'original_121_station_capsule':'artifacts/research_checks/2026-10-03_hs4t2d_relief08_raw',
       'padding':'original central 12m unchanged; extend endpoint cover depth into added sides. This changes exterior geometry and boundary placement, not a pure PML-coefficient test.',
       'halfspace_definition':'cover fills Z=0..12, retains flat ground and all source/material parameters; removes the rock/cover interface. Diagnostic contrast only, not clean training truth.',
       'refinement':'Z only 5cm -> 2.5cm; X/Y unchanged; PML physical thickness held at 1m. Sensitivity check, not a complete convergence certificate.',
       'processing':'same official 20-170MHz / 501 tones / Hann / 8x padding / 200ns tail taper; complex subtraction before reconstruction',
       'fixed_windows_ns':{'ground':[85,120],'underground':[160,220]},
       'low_height_comparison':'air delay align only for diagnostics, using 2*(15-2)/c; no per-trace tuning',
       'estimators':'window energy, raw complex spectral L2, unscaled complex differences, normalized phase correlation; no physical PASS threshold or label',
       'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in
                          ['scripts/prepare_hs4t2d_cause_controls.py','scripts/run_hs4t2d_cause_controls.py',
                           'scripts/run_hs4t2d_cause_controls.cmd','scripts/gprmax_cached_cuda_entry.py']}}
    (args.out/'execution_contract.json').write_text(json.dumps(c,indent=2)+'\n',encoding='utf-8')
    print('Frozen groups:',len(groups),'runs:',sum(g['traces'] for g in groups))


if __name__=='__main__':
    main()
