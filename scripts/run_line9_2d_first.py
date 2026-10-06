"""Freeze/run/audit the bounded Line9 2D pilot or coarse preview on its target machine."""
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import psutil

from build_pdf_profile_geometry import digest, save_json

ROOT=Path(__file__).resolve().parents[1]
TASK=ROOT/'configs/research/line9_2d_first_run_task_v0_1.json'


def resources(task):
    lines=subprocess.check_output(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total,memory.free',
        '--format=csv,noheader,nounits'],text=True).strip().splitlines()
    # Keep existing supervisor GPU0 mapping unambiguous.
    if len(lines)!=1:
        raise ValueError('This runner requires one visible physical GPU; freeze another mapping for multi-GPU hosts')
    name,uuid,driver,total,free=[v.strip() for v in lines[0].split(',')]
    r=dict(gpu=name,gpu_uuid=uuid,driver_version=driver,total_VRAM_bytes=int(total)*2**20,
        free_VRAM_bytes=int(free)*2**20,available_RAM_bytes=psutil.virtual_memory().available,
        total_RAM_bytes=psutil.virtual_memory().total)
    if r['available_RAM_bytes']<task['limits']['min_available_RAM_bytes'] or r['free_VRAM_bytes']<task['limits']['min_free_VRAM_bytes']:
        raise RuntimeError('Capacity rejected before capsule/solver creation: '+json.dumps(r))
    return r


def check_files(contract):
    if digest(TASK)!=contract['task_sha256']:
        raise ValueError('Task changed')
    for name,value in contract['code_identities'].items():
        if digest(Path(name))!=value: raise ValueError('Code changed')
    package=Path(importlib.util.find_spec('gprMax').origin).parent
    for name,value in contract['source_identities'].items():
        if digest(package/name)!=value: raise ValueError('Runtime changed')
    for name,value in contract['file_identities'].items():
        if digest(Path(name))!=value: raise ValueError('Geometry/material/input changed')
    for name,value in contract['toolchain']['binary_sha256'].items():
        if digest(Path(name))!=value: raise ValueError('Compiler/Python binary changed')


def audit(path,completed=False,group_ids=None):
    c=json.loads(path.read_text('utf-8'));check_files(c)
    rows=[]
    for g in c['groups']:
        if group_ids is not None and g['id'] not in group_ids: continue
        p=Path(g['input'])
        if digest(p)!=g['input_sha256']: raise ValueError('Frozen input changed')
        record=dict(id=g['id'],profile_x_m=g['profile_x_m'],acquisition_s_m=g['acquisition_s_m'],
                    input_sha256=digest(p))
        if completed:
            raw=p.with_suffix('.h5')
            with h5py.File(raw,'r') as h:
                data=h['rxs/rx1/Ez'][:]
                if str(h.attrs['gprMax'])!='4.0.0' or data.dtype!=np.float64 or not np.isfinite(data).all():
                    raise ValueError('Native version/dtype/finiteness mismatch')
                spacing=np.asarray(h.attrs['dx_dy_dz']);shape=np.asarray(h.attrs['nx_ny_nz'])
                np.testing.assert_array_equal(spacing,[.025]*3)
                np.testing.assert_array_equal(shape,[16000,3000,1])
                dt=float(h.attrs['dt']);iterations=int(h.attrs['Iterations'])
                if data.shape!=(iterations,) or not np.isclose(dt,.025/(299792458*np.sqrt(2)),rtol=1e-12,atol=0):
                    raise ValueError('Native samples/time grid mismatch')
                if not np.isclose((iterations-1)*dt,1.2e-6,rtol=0,atol=dt):
                    raise ValueError('Native time window mismatch')
                for key,position in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                    grid=np.asarray(h[key].attrs['GridPosition'])
                    np.testing.assert_array_equal(grid[:2],np.rint(np.asarray(position)[:2]/spacing[:2]).astype(int))
                    if grid[2] not in [0,1]: raise ValueError('Invariant-axis location mismatch')
                source=h['srcs/src1/excitation/samples'][:]
                if source.dtype!=np.float64 or not np.isfinite(source).all():
                    raise ValueError('Source samples not native float64 finite')
                attrs=h['srcs/src1/excitation'].attrs
                if str(attrs['Polarisation'])!='z' or float(attrs['SpatialScale'])!=.025 or float(attrs['WaveformAmplitude'])!=40 or float(attrs['WaveformFrequency'])!=100e6:
                    raise ValueError('Native source polarization/length/amplitude/frequency mismatch')
            record.update(raw_sha256=digest(raw),raw_path=str(raw.resolve()),dt_s=dt,
                iterations=iterations,dtype='float64')
        rows.append(record)
    if group_ids is not None and {r['id'] for r in rows}!=set(group_ids):
        raise ValueError('Unknown requested audit group')
    return dict(status='PASS',completed=completed,contract_sha256=digest(path),
        task_sha256=c['task_sha256'],stage=c['stage'],groups=rows,
        scope='Native execution/input identity only; not physics or record-tail acceptance')


def freeze(package,out,stage,pilots=None,station_ids=None):
    import gprMax
    if gprMax.__version__!='4.0.0' or out.exists(): raise ValueError('V4.0.0 and fresh capsule required')
    task=json.loads(TASK.read_text('utf-8'))
    hardware=resources(task)  # Reject before writing anything; no CUDA allocation.
    if sys.platform!='win32': raise ValueError('Use a separately frozen Linux supervisor; this runner is Windows')
    nvcc=shutil.which('nvcc');cl=shutil.which('cl')
    if nvcc is None or cl is None: raise ValueError('Verified CUDA/MSVC compiler environment required')
    nvcc_version=subprocess.check_output([nvcc,'--version'],text=True,encoding='utf-8',errors='replace')
    cl_version=subprocess.run([cl],capture_output=True,text=True,encoding='utf-8',errors='replace')
    toolchain=dict(nvcc_version=nvcc_version,MSVC_version=cl_version.stdout+cl_version.stderr,
        binary_sha256={str(Path(p).resolve()):digest(Path(p)) for p in [nvcc,cl,sys.executable]})
    manifest_path=package/'manifest.json'
    if digest(manifest_path)!=task['package_manifest_sha256']:
        raise ValueError('Wrong prepared package; copy the complete private portable archive')
    m=json.loads(manifest_path.read_text('utf-8'));byid={c['id']:c for c in m['cases']}
    native=Path(gprMax.__file__).parent
    for name,expected in m['installed_native_sources_sha256'].items():
        relative='grid/fdtd_grid.py' if name.replace('\\','/').endswith('/grid/fdtd_grid.py') else 'config.py'
        if digest(native/relative)!=expected:
            raise ValueError('Array allocator/config differs from resource model; re-audit before solving')
    ids=task['acquisition']['pilot_ids'] if stage=='pilots' else task['acquisition']['preview_pending_ids']
    if stage=='continuation':
        if not station_ids or len(set(station_ids))!=len(station_ids) or not set(station_ids)<=set(task['acquisition']['preview_ids']):
            raise ValueError('Continuation requires unique approved coarse-scan station IDs')
        ids=station_ids
    reuse={}
    if stage=='preview':
        if pilots is None: raise ValueError('Pilot capsule required')
        pc=pilots/'execution_contract.json';pv=pilots/'completed_verification.json'
        current=audit(pc,True)
        previous=json.loads(pv.read_text('utf-8'))
        if current!=previous: raise ValueError('Pilot verification/identity changed')
        review=pilots/'pilot_review.json';verdict=json.loads(review.read_text('utf-8'))
        if verdict['status']!='APPROVE_PREVIEW' or verdict['verification_sha256']!=digest(pv):
            raise ValueError('Pilot physical/processing review required before preview')
        for key in ['native_and_processing','record_tail','boundary_and_mesh_limits','resource_and_ETA']:
            if not isinstance(verdict.get(key),str) or not verdict[key].strip():
                raise ValueError('Pilot review must explain '+key)
        report=Path(verdict['analysis_report_path'])
        if digest(report)!=verdict['analysis_report_sha256']:
            raise ValueError('Pilot analysis changed')
        for preview_id,pilot_id in task['acquisition']['reuse_pilot_output'].items():
            row=next(r for r in current['groups'] if r['id']==pilot_id)
            reuse[preview_id]=dict(raw_path=row['raw_path'],raw_sha256=row['raw_sha256'],pilot_id=pilot_id)
    geo=package/'geometries/full2d.h5';db=package/'geometries/line9_research_materials_v1.json'
    if digest(geo)!=task['invariants']['geometry_sha256'] or digest(db)!=task['invariants']['material_database_sha256']:
        raise ValueError('Native model changed')
    for name in ids:
        if digest(package/byid[name]['input'])!=byid[name]['input_sha256']:
            raise ValueError('Prepared station input changed')
    out.mkdir(parents=True);(out/'geometries').mkdir()
    for p in [geo,db]: shutil.copyfile(p,out/'geometries'/p.name)
    groups=[]
    for name in ids:
        case=byid[name];p=out/case['input'];p.parent.mkdir(parents=True)
        shutil.copyfile(package/case['input'],p)
        groups.append(dict(id=name,input=str(p.resolve()),input_sha256=digest(p),
            profile_x_m=case['profile_x_m'],acquisition_s_m=case['acquisition_s_m'],tx_m=case['tx_m'],rx_m=case['rx_m']))
    runtime={p.relative_to(native).as_posix():digest(p) for p in native.rglob('*')
             if p.is_file() and p.suffix in ['.py','.pyd','.so']}
    sources=[ROOT/'scripts'/n for n in ['run_line9_2d_first.py','hs4_station_grid_controls.py',
        'gprmax_cached_cuda_entry.py','hs_capsule_identity.py','build_pdf_profile_geometry.py',
        'analyze_line9_2d_sfcw.py','run_line9_2d_v4.cmd']]
    c=dict(status='FROZEN_APPROVED',approval_basis=task['authorization'],stage=stage,
        task_sha256=digest(TASK),python=str(Path(sys.executable).resolve()),source_identities=runtime,
        code_identities={str(p):digest(p) for p in sources},
        file_identities={str(p.resolve()):digest(p) for p in out.rglob('*') if p.is_file()},
        groups=groups,max_runs=len(groups),no_retry=True,
        min_available_RAM_GiB=task['limits']['min_available_RAM_bytes']/2**30,
        min_free_VRAM_GiB=task['limits']['min_free_VRAM_bytes']/2**30,
        max_owned_RSS_GiB=24,min_system_available_during_run_GiB=1.5,
        max_group_wall_s=3600,max_batch_wall_s=18000 if stage=='pilots' else 172800,
        gpu_lock=str(Path(os.environ.get('GPRMAX_GPU_LOCK',ROOT/'artifacts/local_checks/hs4_gpu_exclusive.lock')).resolve()),hardware_at_freeze=hardware,
        toolchain=toolchain,
        reuse=reuse,private_scope='Input capsule contains site-derived geometry; do not upload to Git',
        no_3d=True,expected_preview_total=99 if stage=='preview' else None)
    if stage=='preview':
        c['pilot_review_identity']=dict(path=str(review.resolve()),sha256=digest(review))
        c['pilot_verification_identity']=dict(path=str(pv.resolve()),sha256=digest(pv))
        c['pilot_analysis_identity']=dict(path=str(report.resolve()),sha256=digest(report))
    save_json(out/'execution_contract.json',c)
    save_json(out/'preflight_verification.json',audit(out/'execution_contract.json',False))
    print(f'Frozen {len(groups)} new solves; {len(reuse)} pilot traces reused; no solver called')


def run(out):
    import hs4_station_grid_controls as supervisor
    c=json.loads((out/'execution_contract.json').read_text('utf-8'))
    task=json.loads(TASK.read_text('utf-8'));hardware=resources(task)
    if hardware['gpu_uuid']!=c['hardware_at_freeze']['gpu_uuid']: raise ValueError('Frozen GPU changed')
    if hardware['driver_version']!=c['hardware_at_freeze']['driver_version']: raise ValueError('Frozen GPU driver changed')
    if c['stage']=='preview':
        for name in ['pilot_review_identity','pilot_verification_identity','pilot_analysis_identity']:
            item=c[name]
            if digest(Path(item['path']))!=item['sha256']: raise ValueError('Pilot approval/evidence changed')
        for item in c['reuse'].values():
            if digest(Path(item['raw_path']))!=item['raw_sha256']: raise ValueError('Reused pilot output changed')
    supervisor.audit=audit
    supervisor.run(out/'execution_contract.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['freeze','run','verify','capacity'])
    p.add_argument('--package',type=Path);p.add_argument('--out',type=Path)
    p.add_argument('--stage',choices=['pilots','preview']);p.add_argument('--pilots',type=Path)
    a=p.parse_args()
    if a.action=='capacity': print(json.dumps(resources(json.loads(TASK.read_text('utf-8'))),indent=2))
    elif a.action=='freeze':
        if a.package is None or a.out is None or a.stage is None: p.error('freeze needs package,out,stage')
        freeze(a.package.resolve(),a.out.resolve(),a.stage,a.pilots.resolve() if a.pilots else None)
    elif a.action=='run': run(a.out.resolve())
    else: print(json.dumps(audit(a.out.resolve()/'execution_contract.json',True),indent=2))
