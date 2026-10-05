"""Archive local benchmark identities, failures, resources and unfinished scope."""
import argparse
import json
from pathlib import Path
import datetime
import h5py
import numpy as np
import psutil
from hs_capsule_identity import sha256

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'artifacts/research_checks'
PREFIX='2026-10-05_local_uav_'


def save(p,v):
    if p.exists():
        raise ValueError('refuse overwrite '+str(p))
    p.write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def main(out):
    if out.exists():
        raise ValueError('fresh summary required')
    rows=[];native_count=0;gpu_times=[];cpu_times=[];ram_peaks=[]
    names=['free_space_r1','free_space_coarse_r2','free_space_fine_r3',
           'plane_r1','plane_r2','plane_r3','plane_r4','plane_r5','plane_fine_r6','plane_domain_r7']
    for name in names:
        p=BASE/(PREFIX+name);cp=p/'execution_contract.json'
        c=json.loads(cp.read_text('utf-8'))
        row=dict(capsule=str(p.relative_to(ROOT)).replace('\\','/'),contract_sha256=sha256(cp),native=[])
        # Frozen source copies preserve old attempts after implementation changes.
        for code,digest in c['code_identities'].items():
            if sha256(code)!=digest:
                candidate=p/('frozen_launcher.cmd' if code.endswith('.cmd') else 'frozen_runner.py')
                if not candidate.exists() or sha256(candidate)!=digest:
                    raise ValueError('missing matching frozen source: '+code+' in '+str(p))
        for group in c['groups']:
            inp=Path(group['input']);raw=inp.with_suffix('.h5')
            if sha256(inp)!=group['input_sha256']:
                raise ValueError('input changed')
            if raw.exists():
                native_count+=1
                with h5py.File(raw) as h:
                    fields={key:float(np.max(abs(h['rxs/rx1/'+key][:]))) for key in ('Ex','Ey','Ez')}
                    dtype=str(h['rxs/rx1/Ey'].dtype)
                    if dtype!='float64' or str(h.attrs['gprMax'])!='4.0.0':
                        raise ValueError('native V4 double required')
                row['native'].append(dict(id=group['id'],raw_sha256=sha256(raw),peak_fields_V_per_m=fields,dtype=dtype))
        if name in ('free_space_r1','free_space_fine_r3'):
            if row['native'] or (p/'execution.jsonl').exists():
                raise ValueError('preflight rejection must not consume attempt')
            row['status']='PREPARED_NOT_RUN_RESOURCE_PREFLIGHT_REJECTED'
            rejection=p/'capacity_rejection.json'
            if not rejection.exists():
                save(rejection,dict(status=row['status'],execution_attempt_consumed=False,calls_solver=False,
                     basis='Recorded runner rejected live RAM/VRAM before STARTED; exact rejected snapshot was not emitted',
                     min_available_RAM_GiB=c['min_available_RAM_GiB'],min_free_VRAM_GiB=c['min_free_VRAM_GiB']))
        elif name in ('plane_r1','plane_r2','plane_r4'):
            if row['native']:
                raise ValueError('initialisation failure should have no native output')
            row['status']='FAILED_INITIALISATION'
            row['reason']={'plane_r1':'OpenBLAS allocation under job commit cap; fixed by explicit thread environment',
                           'plane_r2':'Installed V4 rejects symmetry boundary in 2D',
                           'plane_r4':'psi90 projects onto dead Ex/Hy in TMy; psi0 required for Ey/Hx'}[name]
        elif name=='plane_r3':
            row['status']='NATIVE_COMPLETED_BUT_PHYSICALLY_INVALID_CONFIGURATION'
            row['reason']='psi90 excites Ex with incompatible symmetry faces and TFSF on domain faces; fields grow to >1e70 V/m from unit source'
            save(p/'diagnostic_rejection.json',dict(status=row['status'],reason=row['reason'],native=row['native'],
                 scope='Failed auxiliary input only; not evidence that historical dipole UAV model is unstable'))
        else:
            verified=json.loads((p/'completed_verification.json').read_text('utf-8'))
            if verified['status']!='PASS' or verified['contract_sha256']!=sha256(cp):
                raise ValueError('completed identity mismatch')
            for group,item in zip(c['groups'],verified['groups']):
                if sha256(Path(group['input']).with_suffix('.h5'))!=item['raw_sha256']:
                    raise ValueError('verified raw changed')
            row['status']='NATIVE_EXECUTION_IDENTITY_PASS_SEE_SEPARATE_ANALYTIC_DIAGNOSTIC'
        if 'free_space_coarse' in name:
            events=[json.loads(s) for s in (p/'execution.jsonl').read_text('utf-8').splitlines()]
            for e in events:
                if 'elapsed_s' in e:
                    gpu_times.append(e['elapsed_s']);ram_peaks.append(e['peak_owned_RSS_bytes'])
        elif name not in ('plane_r1','plane_r2','plane_r3','plane_r4') and name.startswith('plane_'):
            times=[];peaks=[]
            for group in c['groups']:
                s=json.loads((p/(group['id']+'_supervision.json')).read_text('utf-8'))
                times.append(s['wall_s']);peaks.append(s['peak_job_commit_bytes'])
            cpu_times+=times
            row['completed_wall_s_sum']=sum(times);row['peak_job_commit_bytes']=max(peaks)
        rows.append(row)
    analyses={}
    for name in ('free_space_coarse_results_r2','plane_results_r5','plane_results_fine_r6','plane_results_domain_r7','observation_r1','plane_reference_checks_r1','grid_comparison_r1'):
        p=BASE/(PREFIX+name)/'summary.json'
        analyses[name]=dict(path=str(p.relative_to(ROOT)).replace('\\','/'),sha256=sha256(p),status=json.loads(p.read_text('utf-8'))['status'])
    cells=round(16/.05)*round(16/.05)*round(24/.05)
    current=dict(available_RAM_bytes=psutil.virtual_memory().available,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    resource_estimate=dict(status='RESOURCE_ESTIMATE_ONLY_NOT_FROZEN_SOLVER_INPUTS',
         illustrative_domain_m=[16,16,24],spacing_m=[.05,.05,.05],cells=cells,
         host_basic_array_lower_bound_bytes=cells*94,
         device_fields_and_ID_lower_bound_bytes=cells*72,
         excluded='PML histories, dispersive states, boundary-grid enlargement, Python/CUDA context, construction intermediates',
         domain_converged=False,current_followup_snapshot=current,
         snapshot_scope='End-of-unit observation; not original preflight-rejection snapshot',
         conclusion='Full finite-height layered Green-function and 12-case layout/depth controls not executed; native air refinement resource-rejected')
    out.mkdir(parents=True)
    save(out/'summary.json',dict(status='COMPLETED_LOCAL_LIGHTWEIGHT_CHECKS_FULL_3D_UNFINISHED',capsules=rows,
         native_H5_count_including_invalid=native_count,valid_native_H5_count=native_count-5,
         successful_GPU_group_wall_s=gpu_times,GPU_peak_owned_RSS_bytes=max(ram_peaks),
         successful_CPU_group_wall_s_sum=sum(cpu_times),analyses=analyses,resource_estimate=resource_estimate,
         calls_training=False,field_fitted=False,code_sha256=sha256(__file__)))
    # Failed five-case capsule: actual conditional trace ensemble, not a spatial B-scan.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei'];plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    names_zh=dict(air='空气',halfspace='无损半空间',slab='无损三米平层',conductive='平层加电导',debye='平层加电导与Debye')
    for case in names_zh:
        with h5py.File(BASE/(PREFIX+'plane_r3')/(case+'.h5')) as h:
            magnitude=np.maximum.reduce([abs(h['rxs/rx1/'+k][:]) for k in ('Ex','Ey','Ez')])
            ax.semilogy(np.arange(len(magnitude))*float(h.attrs['dt'])*1e9,np.maximum(magnitude,1e-30),label=names_zh[case])
    ax.set(title='失败辅助配置：五场景同一接收点，非空间B-scan\n单位平面波激励下场值发散；不可用于物理分析',
           xlabel='时间(ns)',ylabel='三分量最大绝对场(V/m)')
    ax.legend();ax.grid(alpha=.2)
    fig.savefig(out/'invalid_plane_configuration.png',dpi=140);plt.close(fig)
    print(json.dumps(dict(status='ARCHIVED',native_count=native_count,resource_estimate=resource_estimate),indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    main(ap.parse_args().out)
