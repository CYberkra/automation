"""Prepare source-only replacements of existing Line9 inputs; never solve."""
import argparse
import copy
import json
from pathlib import Path
import shutil

from gprMax.hash_cmds_file import get_user_objects
from hs_capsule_identity import sha256 as sha

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/'configs/research/line9_source_v0_2.json'


def replace_source(text):
    lines=text.splitlines(keepends=True)
    wave=[i for i,line in enumerate(lines) if line.startswith('#waveform:')]
    dipoles=[i for i,line in enumerate(lines) if line.startswith('#hertzian_dipole:')]
    if len(wave)!=1 or len(dipoles)!=1:
        raise ValueError('Exactly one waveform and one Hertzian dipole required')
    wi,di=wave[0],dipoles[0];w=lines[wi].split();d=lines[di].split()
    if len(w)!=5 or w[1]!='ricker' or len(d)!=6 or d[1]!='z' or d[-1]!=w[-1]:
        raise ValueError('Expected the existing z Ricker source without delayed start')
    for i in [wi,di]:
        end='\r\n' if lines[i].endswith('\r\n') else '\n' if lines[i].endswith('\n') else ''
        lines[i]=('#waveform: impulse 1 1 impulse' if i==wi else ' '.join(d[:-1]+['impulse']))+end
    return ''.join(lines)


def prepare(parent,out,group_ids):
    if out.exists():raise ValueError('Fresh package required; historical inputs stay read-only')
    design=json.loads(DESIGN.read_text('utf-8'));old=json.loads((parent/'manifest.json').read_text('utf-8'))
    by_id={g['id']:g for g in old['groups']}
    if not group_ids or len(set(group_ids))!=len(group_ids) or any(k not in by_id for k in group_ids):
        raise ValueError('Select unique existing group IDs')
    # Verify all selected inputs before creating output.
    for name in group_ids:
        for key in ['input','geometry','material']:
            if sha(parent/by_id[name][key])!=by_id[name][key+'_sha256']:
                raise ValueError('Parent identity changed: '+name+'/'+key)
    out.mkdir(parents=True);groups=[];checks=[]
    for name in group_ids:
        g=copy.deepcopy(by_id[name])
        for key in ['input','geometry','material']:
            dst=out/g[key];dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(parent/g[key],dst)
        card=out/g['input'];original=card.read_bytes().decode('utf-8');changed=replace_source(original)
        card.write_bytes(changed.encode('utf-8'))
        objects=get_user_objects(changed.splitlines(),input_dir=card.parent)
        before=original.splitlines();after=changed.splitlines()
        assert [v for v in before if not v.startswith(('#waveform:','#hertzian_dipole:'))]==[
            v for v in after if not v.startswith(('#waveform:','#hertzian_dipole:'))]
        g.update(input_sha256=sha(card),source_type='impulse',source_frequency_Hz=1.,source_amplitude_A=1.,
                 source_start_s=0.,source_frequency_slot_unused=True)
        g.pop('source_band_min_over_max',None)
        groups.append(g);checks.append(dict(id=name,only_waveform_and_source_ID_changed=True,
            geometry_sha256_unchanged=sha(out/g['geometry'])==g['geometry_sha256'],
            material_sha256_unchanged=sha(out/g['material'])==g['material_sha256'],
            official_parser_object_count=len(objects)))
    manifest=dict(status='PREPARED_OFFICIAL_IMPULSE_SOURCE_NOT_RUN',approval_basis=design['authorization'],
        execution_authorized=False,source_design_sha256=sha(DESIGN),generator_sha256=sha(__file__),
        parent_manifest_sha256=sha(parent/'manifest.json'),groups=groups,no_retry=True,
        dt_s=old['dt_s'],SFCW=design['sfcw'],source_only_checks=checks,
        limits=design['limits'],reused_solver_outputs=[])
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=manifest['status'],inputs=len(groups),new_solver_runs=0)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--groups',nargs='+')
    a=p.parse_args();design=json.loads(DESIGN.read_text('utf-8'))
    prepare(a.parent,a.out,a.groups or design['preparation_example_group_ids'])
