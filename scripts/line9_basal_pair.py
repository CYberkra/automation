"""Prepare/freeze/run two bounded Line9 basal-sand controls with passive snapshots."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import h5py
import numpy as np
import psutil

from hs_capsule_identity import sha256 as sha

ROOT=Path(__file__).resolve().parents[1]
FIELDS=['Ex','Ey','Ez','Hx','Hy','Hz']
CODE=['line9_basal_pair.py','hs4_station_grid_controls.py','hs_capsule_identity.py','gprmax_cached_cuda_entry.py','run_line9_2d_v4.cmd']


def save(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def basal_mask(data):
    """Only sandstone continuously connected to the model's lower edge."""
    if data.ndim!=3 or data.shape[2]!=1 or not np.all(data[:,0,0]==3):
        raise ValueError('Require 2D sandstone along the entire lower edge')
    non=data[:,:,0]!=3
    if not np.all(non.any(axis=1)):
        raise ValueError('Every column must have a non-sandstone contact')
    stop=np.argmax(non,axis=1)
    return (np.arange(data.shape[1])[None,:]<stop[:,None])[:,:,None]


def self_checks():
    data=np.array([[[3],[3],[2],[3],[2],[1],[0]],[[3],[2],[3],[2],[1],[1],[0]]],dtype=np.int16)
    mask=basal_mask(data); expected=np.zeros_like(mask);expected[0,:2]=True;expected[1,:1]=True
    np.testing.assert_array_equal(mask,expected)
    changed=data.copy();changed[mask]=2
    assert changed[0,3,0]==3 and changed[1,2,0]==3
    bad=data.copy();bad[0,0,0]=2
    try:basal_mask(bad)
    except ValueError:pass
    else:raise AssertionError('Missing basal contact must reject')
    return dict(preserve_sandstone_interbeds=True,reject_missing_basal_contact=True)


def prepare(package,review,out,evidence):
    if out.exists() or evidence.exists():raise ValueError('Fresh package/evidence required')
    checks=self_checks()
    audit=json.loads(review.read_text('utf-8'))
    if audit['package']!='line9_pkg2_x160-110_agl35_142st':raise ValueError('Wrong source package')
    row=next(r for r in audit['records'] if r['id']=='full2d_c0008')
    source=package/'cases'/row['id']/'profile.in';raw=source.with_suffix('.h5')
    geo=package/'geometries/full2d_compact.h5';db=package/'geometries/line9_research_materials_v1_smoothed.json'
    for path,expected in [(source,row['input_sha256']),(raw,row['native_sha256']),
                          (geo,audit['geometry_sha256']),(db,audit['materials_sha256'])]:
        if sha(path)!=expected:raise ValueError('Source identity changed: '+str(path))
    with h5py.File(raw) as h:
        dt=float(h.attrs['dt']);iterations=int(h.attrs['Iterations'])
        tx=h['srcs/src1'].attrs['Position'].tolist();rx=h['rxs/rx1'].attrs['Position'].tolist()
        assert h['rxs/rx1/Ez'].dtype==np.float32
        np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[4400,1600,1])
    # ~2ns through 480ns, then ~20ns. Passive observers, no solver coarsening.
    timeline=sorted(set(range(0,round(480e-9/dt)+1,34))|set(range(round(480e-9/dt)+34,iterations,340)))
    suffix=''.join(f'#snapshot: 0 0 0 110 40 0.025 0.1 0.1 0.025 {j} snap{j:05d}.h5\n' for j in timeline)
    original=source.read_bytes()
    if not original.endswith(b'\n'):raise ValueError('Expected terminal newline')
    with h5py.File(geo) as h:data=h['data'][:]
    mask=basal_mask(data)
    out.mkdir(parents=True);(out/'reference').mkdir()
    shutil.copyfile(raw,out/'reference/profile_fp32.h5')
    shutil.copyfile(source,out/'reference/profile_original.in')
    save(out/'reference/station_audit.json',row)
    groups=[]
    for name in ['H1','H0']:
        p=out/name/'cases'/row['id']/'profile.in';p.parent.mkdir(parents=True)
        gp=out/name/'geometries';gp.mkdir(parents=True)
        shutil.copyfile(geo,gp/geo.name);shutil.copyfile(db,gp/db.name)
        if name=='H0':
            with h5py.File(gp/geo.name,'r+') as h:
                changed=data.copy();changed[mask]=2;h['data'][:]=changed
        p.write_bytes(original+suffix.encode('ascii'))
        groups.append(dict(id=name,input=p.relative_to(out).as_posix(),input_sha256=sha(p),
                           geometry=(gp/geo.name).relative_to(out).as_posix(),geometry_sha256=sha(gp/geo.name),
                           material=(gp/db.name).relative_to(out).as_posix(),material_sha256=sha(db),
                           tx_m=tx,rx_m=rx))
    assert sha(out/groups[0]['input'])==sha(out/groups[1]['input'])
    with h5py.File(out/groups[1]['geometry']) as h:
        changed=h['data'][:];np.testing.assert_array_equal(changed[~mask],data[~mask]);assert np.all(changed[mask]==2)
    from gprMax.hash_cmds_file import get_user_objects
    counts=[len(get_user_objects((out/g['input']).read_text('utf-8').splitlines(),input_dir=(out/g['input']).parent)) for g in groups]
    history=int(1100*400*6*8*len(timeline))
    report=dict(status='PREPARED_NOT_RUN',calls_solver=False,self_checks=checks,source_git_base='4e34d180e7027bb22a862d6f7df9e34a7a0c33ed',
        approval_basis='User: 通过ssh连接我们那个机器然后开跑吧; two single-station FP64 controls with equal sparse snapshots.',
        source_audit_sha256=sha(review),source_input_sha256=sha(source),source_native_sha256=sha(raw),
        geometry_original_sha256=sha(geo),material_sha256=sha(db),script_sha256=sha(__file__),
        source_package=audit['package'],station_id=row['id'],chainage_m=row['chainage_m'],midpoint_agl_m=row['geometry']['midpoint_agl_m'],
        geometry_diagnostic=row['geometry'],groups=groups,dt_s=dt,iterations=iterations,
        domain_m=[110.,40.,.025],native_shape=[4400,1600,1],spacing_m=[.025]*3,
        pml_cells=[80,80,0,80,80,0],pml='HORIPML',source='z Hertzian first-half-step 40A impulse, frequency label 1Hz, spatial scale .025m',
        receiver='Ez',precision='native float64 required',time_window_s=8e-7,
        snapshot_iterations=timeline,snapshot_shape=[1100,400,1],snapshot_spacing_m=[.1,.1,.025],
        snapshot_six_field_history_bytes_per_solve=history,snapshot_fields=FIELDS,native_parser_counts=counts,
        changed_voxel_count=int(mask.sum()),changed_mask_sha256=__import__('hashlib').sha256(mask.tobytes()).hexdigest(),
        change='Only bottom-edge-connected material ID3 sandstone -> ID2 mudstone; all overlying interbeds retained, including side/bottom PML continuation of that region.',
        invariant_inputs_identical_except_external_geometry=True,
        processing=dict(frequency_start_Hz=20e6,frequency_step_Hz=300000.,frequency_count=501,frequency_stop_Hz=170e6,
            window='hann',zero_pad_factor=8,tail_taper=False,background='H1-H0 in complex domain only',AGC=False,per_trace_normalisation=False),
        limitations='This ablation includes propagation/interactions and changed lower PML material continuation; not a pure reflection separation or PML-amplitude control. Sparse pulse snapshots do not certify complete source-spectrum phase.')
    save(out/'manifest.json',report)
    evidence.mkdir(parents=True);save(evidence/'preparation.json',{k:v for k,v in report.items() if k not in ['geometry_diagnostic','snapshot_iterations']})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axs=plt.subplots(1,2,figsize=(14,5),layout='constrained')
    for ax,v,title in zip(axs,[data,changed],['H1：原始地质，保留底部砂岩','H0：仅底部连续砂岩替换为泥岩']):
        pic=ax.imshow(v[::4,::4,0].T,origin='lower',extent=[0,110,0,40],aspect='equal',cmap=ListedColormap(['#deeff8','#d5b577','#a599a6','#ead67b']),vmin=-.5,vmax=3.5)
        ax.plot(tx[0],tx[1],'^r',label='发射点');ax.plot(rx[0],rx[1],'vb',label='接收点')
        ax.set(title=title,xlabel='裁剪域局部x / m',ylabel='局部y / m');ax.legend(fontsize=8)
    cb=fig.colorbar(pic,ax=axs,ticks=[0,1,2,3]);cb.ax.set_yticklabels(['空气','粉质黏土','泥岩','砂岩'])
    fig.suptitle('X179.25m / 第8道；高航高；110×40m / 2.5cm；两次输入相同，仅导入体素不同；尚未求解')
    fig.savefig(evidence/'geometry_pair.png',dpi=140);plt.close(fig)
    print(json.dumps(dict(snapshot_count=len(timeline),six_field_history_GiB=history/2**30,changed_voxels=int(mask.sum())),ensure_ascii=False))


def resources():
    line=subprocess.check_output(['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total,memory.free','--format=csv,noheader,nounits'],text=True).strip().splitlines()
    if len(line)!=1:raise ValueError('One visible GPU required')
    name,uuid,driver,total,free=[s.strip() for s in line[0].split(',')]
    result=dict(gpu=name,gpu_uuid=uuid,driver_version=driver,total_VRAM_bytes=int(total)*2**20,
                free_VRAM_bytes=int(free)*2**20,available_RAM_bytes=psutil.virtual_memory().available)
    if result['free_VRAM_bytes']<12*2**30 or result['available_RAM_bytes']<20*2**30:
        raise RuntimeError('RAM/VRAM budget rejected: '+json.dumps(result))
    return result


def check_files(c):
    for category in ['code_identities','file_identities','binary_identities']:
        for p,expected in c[category].items():
            if sha(p)!=expected:raise ValueError(category+' changed: '+p)
    native=Path(importlib.util.find_spec('gprMax').origin).parent
    for p,expected in c['source_identities'].items():
        if sha(native/p)!=expected:raise ValueError('Installed runtime changed: '+p)


def audit(path,completed=False):
    c=json.loads(path.read_text('utf-8'));check_files(c);m=c['prepared_manifest']
    rows=[];sources=[]
    with h5py.File(Path(c['package'])/'reference/profile_fp32.h5') as ref:
        for g in c['groups']:
            p=Path(g['input']);row=dict(id=g['id'],input_sha256=sha(p))
            if completed:
                raw=p.with_suffix('.h5')
                with h5py.File(raw) as h:
                    if str(h.attrs['gprMax'])!='4.0.0':raise ValueError('Require V4.0.0')
                    for k in ['Iterations','nx_ny_nz','dx_dy_dz','dt']:
                        np.testing.assert_array_equal(h.attrs[k],ref.attrs[k])
                    v=h['rxs/rx1/Ez'][:];s=h['srcs/src1/excitation/samples'][:]
                    if v.dtype!=np.float64 or s.dtype!=np.float64 or not np.isfinite(v).all() or v.shape!=(m['iterations'],):
                        raise ValueError('Native float64, shape and finite fields required')
                    np.testing.assert_array_equal(s,ref['srcs/src1/excitation/samples'][:].astype(np.float64));sources.append(s)
                    for k in ref['srcs/src1/excitation'].attrs:
                        np.testing.assert_array_equal(h['srcs/src1/excitation'].attrs[k],ref['srcs/src1/excitation'].attrs[k])
                    for k in ['srcs/src1','rxs/rx1']:
                        np.testing.assert_array_equal(h[k].attrs['Position'],ref[k].attrs['Position'])
                    row.update(raw_sha256=sha(raw),raw_path=str(raw),dtype=str(v.dtype),iterations=len(v),dt_s=float(h.attrs['dt']))
                files=sorted(p.parent.glob('profile_snaps/snap*.h5'));snapshot_rows=[]
                if len(files)!=len(m['snapshot_iterations']):raise ValueError('Missing or extra snapshots')
                for file,j in zip(files,m['snapshot_iterations']):
                    with h5py.File(file) as h:
                        if int(h.attrs['iteration'])!=j or str(h.attrs['gprMax'])!='4.0.0':raise ValueError('Snapshot identity/time mismatch')
                        np.testing.assert_allclose(h.attrs['time'],j*m['dt_s'],atol=0,rtol=1e-14)
                        np.testing.assert_allclose(h.attrs['magnetic_time'],(j-.5)*m['dt_s'],atol=0,rtol=1e-14)
                        np.testing.assert_array_equal(h.attrs['origin'],[0.,0.,0.]);np.testing.assert_array_equal(h.attrs['dx_dy_dz'],m['snapshot_spacing_m'])
                        for key in FIELDS:
                            ds=h[key]
                            if ds.dtype!=np.float64 or list(ds.shape)!=m['snapshot_shape'] or not np.isfinite(ds[:]).all():raise ValueError('Snapshot validity failed')
                    snapshot_rows.append(dict(file=str(file),iteration=j,sha256=sha(file)))
                row.update(snapshot_count=len(files),snapshots=snapshot_rows)
            rows.append(row)
    if completed:np.testing.assert_array_equal(sources[0],sources[1])
    return dict(status='PASS_NATIVE_AND_OBSERVER_IDENTITY_NOT_PHYSICAL_ATTRIBUTION',completed=completed,
                contract_sha256=sha(path),groups=rows)


def freeze(package,out):
    import gprMax
    if out.exists() or gprMax.__version__!='4.0.0' or sys.platform!='win32':raise ValueError('Fresh Windows V4.0.0 capsule required')
    m=json.loads((package/'manifest.json').read_text('utf-8'));hardware=resources()
    if m['status']!='PREPARED_NOT_RUN' or len(m['groups'])!=2:raise ValueError('Require prepared two-case pair')
    nvcc=shutil.which('nvcc');cl=shutil.which('cl')
    if not nvcc or not cl:raise ValueError('CUDA/MSVC environment required')
    native=Path(gprMax.__file__).parent
    for g in m['groups']:
        for key in ['input','geometry','material']:
            if sha(package/g[key])!=g[key+'_sha256']:raise ValueError('Prepared package changed')
    if shutil.disk_usage(out.parent).free<20*2**30:raise RuntimeError('Require 20GiB free disk')
    with h5py.File(package/m['groups'][0]['geometry']) as h:data=h['data'][:]
    with h5py.File(package/m['groups'][1]['geometry']) as h:changed=h['data'][:]
    mask=basal_mask(data);np.testing.assert_array_equal(changed[~mask],data[~mask]);assert np.all(changed[mask]==2)
    history=m['snapshot_six_field_history_bytes_per_solve']
    # 6GiB reserved above six-field snapshot history for native/dispersive arrays and buffers.
    if hardware['free_VRAM_bytes']<history+6*2**30:raise RuntimeError('Explicit observer+native device budget rejected')
    groups=[dict(g,input=str((package/g['input']).resolve())) for g in m['groups']]
    c=dict(status='FROZEN_APPROVED',approval_basis=m['approval_basis'],python=str(Path(sys.executable).resolve()),
        package=str(package),prepared_manifest=m,groups=groups,max_runs=2,no_retry=True,
        source_identities={p.relative_to(native).as_posix():sha(p) for p in native.rglob('*') if p.is_file() and p.suffix in ['.py','.pyd','.so']},
        code_identities={str(ROOT/'scripts'/name):sha(ROOT/'scripts'/name) for name in CODE},
        file_identities={str(p.resolve()):sha(p) for p in package.rglob('*') if p.is_file()},
        binary_identities={str(Path(p).resolve()):sha(p) for p in [sys.executable,nvcc,cl]},
        hardware_at_freeze=hardware,toolchain=dict(nvcc=subprocess.check_output([nvcc,'--version'],text=True),
            cl=subprocess.run([cl],capture_output=True,text=True,errors='replace').stderr),
        solver_entrypoint='gprmax_cached_cuda_entry.py',additional_solver_args=[],
        min_available_RAM_GiB=20.,min_free_VRAM_GiB=12.,max_owned_RSS_GiB=24.,min_system_available_during_run_GiB=2.,
        max_group_wall_s=1800,max_batch_wall_s=3600,
        gpu_lock='E:/automation_djh/artifacts/local_checks/hs4_gpu_exclusive.lock',
        cancel_file=str(out/'USER_STOP'),lease_file=str(out/'session_heartbeat'),max_lease_age_s=120,
        snapshot_budget_basis='Native observers fit free device budget including six fields and6GiB native/dispersive reserve; no forced streaming or field-kernel edits.',
        limits='Matched lower-material contrast including interaction/PML continuation; not pure reflection, detection/SNR or physical validation.')
    out.mkdir(parents=True);save(out/'execution_contract.json',c)
    (out/'session_heartbeat').write_text('Active authorized SSH session\n',encoding='utf-8')
    save(out/'preflight_verification.json',audit(out/'execution_contract.json'))
    print(json.dumps(dict(status=c['status'],hardware=hardware,snapshot_count=len(m['snapshot_iterations']),max_runs=2),ensure_ascii=False))


def run(out):
    import hs4_station_grid_controls as supervisor
    c=json.loads((out/'execution_contract.json').read_text('utf-8'));hw=resources()
    if any(hw[k]!=c['hardware_at_freeze'][k] for k in ['gpu_uuid','driver_version']):raise ValueError('Frozen GPU/driver changed')
    check_files(c);supervisor.audit=audit;supervisor.run(out/'execution_contract.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['prepare','freeze','run','verify','checks'])
    for name in ['package','review','out','evidence']:p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.package.resolve(),a.review.resolve(),a.out.resolve(),a.evidence.resolve())
    elif a.action=='freeze':freeze(a.package.resolve(),a.out.resolve())
    elif a.action=='run':run(a.out.resolve())
    elif a.action=='verify':print(json.dumps(audit(a.out.resolve()/'execution_contract.json',True),ensure_ascii=False))
    else:print(json.dumps(self_checks()))
