"""Freeze source-conditioning and full-transect Z refinement confirmation."""
import argparse
import json
from pathlib import Path

from hs_capsule_identity import sha256
from prepare_hs4t2d_cause_controls import ROOT, DESIGN, create_input


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new study required')
    args.out.mkdir(parents=True)
    base=json.loads((DESIGN/'execution_contract_cuda_cached_att4.json').read_text('utf-8'))
    groups=[]
    specs=[(f'ricker{width}_{role}',width,.05,role=='halfspace',1,'ricker95')
           for width in (12,24,36) for role in ('rough','halfspace')]
    specs.extend([(f'zfine13_{role}',12,.025,role=='halfspace',13,'impulse') for role in ('rough','halfspace')])
    for name,width,dz,halfspace,traces,condition in specs:
        source=DESIGN/'relief2p0'/('hs4t2d_t61.in' if traces==1 else 'hs4t2d_t01.in')
        text=create_input(source,width,15,dz,halfspace,traces)
        if condition=='ricker95':
            if text.count('#waveform: impulse 1 1 impulse')!=1:
                raise ValueError('unexpected source definition')
            text=text.replace('#waveform: impulse 1 1 impulse','#waveform: ricker 1 95e6 impulse')
        folder=args.out/name; folder.mkdir()
        target=folder/'profile.in'; target.write_text(text,encoding='utf-8')
        groups.append({'id':name,'input':str(target.resolve()),'sha256':sha256(target),
                       'source_input':str(source.relative_to(ROOT)),'source_sha256':sha256(source),
                       'width_m':width,'x_translation_m':(width-12)/2,'height_m':15,'dz_m':dz,
                       'halfspace':halfspace,'source_condition':condition,'traces':traces,
                       'old_station_indices':[61] if traces==1 else list(range(1,122,10)),
                       'minimum_free_VRAM_GiB':3.5 if width==36 else 2.4 if width==24 or dz==.025 else 1.4})
    c={'status':'FROZEN_APPROVED','approval_basis':'User requested to attempt to identify cause; confirm source conditioning and spatial Z-grid effects identified in first controlled study. Frozen before execution.',
       'purpose':'Separate instantaneous impulse grid texture from in-band transfer function; verify full B-scan under Z refinement',
       'backend':'CUDA','gpu_device':0,'precision':'double','python':base['python'],
       'max_runs':32,'max_wall_s':1200,'min_available_RAM_GiB':1,'groups':groups,
       'source_identities':base['source_identities'],
       'invariants':'geometry/material/source positions/time/PML same as declared counterparts; Ricker modifies waveform only; width extends endpoint materials and translates central domain; Z refinement keeps 1m physical PML',
       'processing':'load exact saved excitation, divide complex spectrum by that source; original 20-170MHz/501 tones/Hann/8x/200ns tail; no raw-amplitude comparison between impulse and Ricker',
       'scope':'source and full-transect vertical-grid sensitivity, not complete spatial convergence or physical validation',
       'code_identities':{str(ROOT/name):sha256(ROOT/name) for name in
                           ['scripts/prepare_hs4t2d_source_controls.py','scripts/prepare_hs4t2d_cause_controls.py',
                            'scripts/run_hs4t2d_cause_controls.py','scripts/run_hs4t2d_cause_controls.cmd',
                            'scripts/gprmax_cached_cuda_entry.py']}}
    (args.out/'execution_contract.json').write_text(json.dumps(c,indent=2)+'\n',encoding='utf-8')
    print('Frozen',len(groups),'groups /',sum(g['traces'] for g in groups),'runs')


if __name__=='__main__':
    main()
