"""Portable HS4 task: CPU design, target-machine freeze, bounded GPU run, audit.

Never resumes a consumed historical execution contract. Paths in the design
are relative to its folder; execution contracts are freshly frozen on target.
"""
import argparse
from decimal import Decimal, ROUND_HALF_UP
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np
import psutil

from hs_capsule_identity import sha256
from hs4_station_grid_controls import commands, terminate_owned_tree

ROOT = Path(__file__).resolve().parents[1]
CHECKS = ROOT / 'artifacts/research_checks'
FINE = CHECKS / '2026-10-04_hs4_patch_finest'
OLD3D = CHECKS / '2026-10-02_halfspace_standard_hs'
STAGES = ('replay', 'centre1cm', '3d-centre', '3d-ends', '3d-xwide-centre')


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                               allow_nan=False) + '\n', encoding='utf-8')


def read(path):
    return json.loads(Path(path).read_text('utf-8'))


def version_runtime():
    if not (3, 11) <= sys.version_info[:2] <= (3, 13):
        raise ValueError('gprMax V4 requires Python 3.11-3.13')
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    import gprMax
    if str(gprMax.__version__) != '4.0.0':
        raise ValueError('gprMax V4.0.0 required')
    reference = read(FINE / 'execution_contract.json')['source_identities']
    identities = {p.relative_to(package).as_posix(): sha256(p)
                  for p in sorted(package.rglob('*')) if p.is_file()
                  and p.suffix in ('.py', '.pyd', '.so')}
    binary_changes, source_changes = [], []
    for name, digest in reference.items():
        actual = sha256(package / name) if (package / name).is_file() else None
        if actual != digest:
            (binary_changes if name.endswith(('.pyd', '.so')) else source_changes).append(name)
    if source_changes:
        raise ValueError('Reviewed V4 source differs/missing: ' + ', '.join(source_changes))
    return {'python': sys.executable, 'python_version': sys.version,
            'gprMax': '4.0.0', 'package': str(package), 'source_identities': identities,
            'historical_binary_changes': binary_changes,
            'comparison_scope': 'CROSS_BUILD_REPLAY_REQUIRED' if binary_changes else 'SAME_ARCHIVED_PACKAGE_FILES'}


def resources():
    free = int(subprocess.check_output(
        ['nvidia-smi', '--query-gpu=memory.free', '--format=csv,noheader,nounits'],
        text=True).splitlines()[0]) * 2**20
    return {'available_RAM_bytes': psutil.virtual_memory().available,
            'free_VRAM_bytes': free}


def case_record(folder, stage, name, text, source):
    target = folder / stage / name / 'profile.in'
    target.parent.mkdir(parents=True)
    target.write_text(text + '\n', encoding='utf-8')
    c = commands(target)
    dx = np.array(list(map(float, c['#dx_dy_dz:'][0])))
    cells = np.rint(np.array(list(map(float, c['#domain:'][0]))) / dx).astype(int)
    requested_tx = list(map(float, c['#hertzian_dipole:'][0][1:4]))
    requested_rx = list(map(float, c['#rx:'][0][:3]))
    tx = (np.rint(np.asarray(requested_tx)/dx)*dx).tolist()
    rx = (np.rint(np.asarray(requested_rx)/dx)*dx).tolist()
    active = cells > 1
    if any(not np.allclose(np.array(p)[active]/dx[active], np.rint(np.array(p)[active]/dx[active]), rtol=0, atol=1e-9)
           for p in (requested_tx, requested_rx)):
        raise ValueError('source/receiver must use actual grid coordinates')
    for words in c['#box:']:
        corners = np.array(list(map(float, words[:6]))).reshape(2, 3)
        if not np.allclose(corners/dx, np.rint(corners/dx), rtol=0, atol=1e-9):
            raise ValueError('off-grid geometry')
        if np.any(corners[0] < 0) or np.any(corners[1] > cells*dx) or np.any(corners[1] <= corners[0]):
            raise ValueError('invalid box bounds')
    return {'id': name, 'input': target.relative_to(folder).as_posix(),
            'input_sha256': sha256(target), 'source': source.relative_to(ROOT).as_posix(),
            'source_sha256': sha256(source), 'tx_m': tx, 'rx_m': rx,
            'requested_tx_m':requested_tx, 'requested_rx_m':requested_rx,
            'component': 'Ex' if stage.startswith('3d') else 'Ey',
            'spacing_m': dx.tolist(), 'cells': cells.tolist(),
            'core_device_array_lower_bound_bytes': int(np.prod(cells+1))*96}


def prepare(folder):
    if folder.exists():
        raise ValueError('new design directory required')
    folder.mkdir(parents=True)
    groups = {stage: [] for stage in STAGES}
    for stage in ('replay', 'centre1cm'):
        for variant in ('base', 'halfspace', 'crest', 'slope'):
            source = FINE / variant / 'profile.in'
            lines = source.read_text('utf-8').splitlines()
            if stage == 'centre1cm':
                lines = [('#dx_dy_dz: 0.01 0.05 0.01' if line.startswith('#dx_dy_dz:') else
                          '#pml_cells: 200 0 100 200 0 100' if line.startswith('#pml_cells:') else
                          '#geometry_view: 0 0 0 36 0.05 33 0.01 0.05 0.01 hs4t2d_geom n'
                          if line.startswith('#geometry_view:') else line) for line in lines]
            groups[stage].append(case_record(folder, stage, variant, '\n'.join(lines), source))
    source = OLD3D / 'hs4_col6_t13.in'
    originals = source.read_text('utf-8').splitlines()
    boxes = [line.split()[1:] for line in originals if line.startswith('#box:') and line.endswith(' cover')]
    if len(boxes) != 2304:
        raise ValueError('expected 48x48 archived surface')
    znew = [float((Decimal('9.05') + 2*(Decimal(w[2])-Decimal('9.05')))
                  .quantize(Decimal('.01'), rounding=ROUND_HALF_UP)) for w in boxes]
    profile = np.genfromtxt(CHECKS/'2026-10-03_hs4t2d_relief08_design/relief2p0/profile.csv',
                            delimiter=',', names=True)['grid_z_m']
    ztable = np.asarray(znew).reshape(48, 48)
    if not np.array_equal(ztable[6], profile):
        raise ValueError('3D bin6 is not the executed 0.8m 2D interface')
    for stage in ('3d-centre', '3d-ends', '3d-xwide-centre'):
        shift = 2 if stage == '3d-xwide-centre' else 0
        stations = [('left', 1), ('right', 25)] if stage == '3d-ends' else [('centre', 13)]
        for station, index in stations:
            for variant in ('rough08', 'flat', 'halfspace'):
                lines = []
                for line in originals:
                    if line.startswith(('#box:', '#geometry_view:', '#hertzian_dipole:', '#rx:')):
                        continue
                    if line.startswith('#domain:') and shift:
                        line = '#domain: 16 12 33'
                    lines.append(line)
                tx = 2.6 + .25*(index-1)
                lines += [f'#hertzian_dipole: x {1.6+shift:g} {tx:g} 27 impulse',
                          f'#rx: {1.6+shift:g} {tx+1.3:g} 27 {station} Ex',
                          f'#box: 0 0 0 {12+2*shift} 12 12 rock']
                if variant == 'rough08':
                    for w, z in zip(boxes, znew):
                        values = list(map(float, w[:6])); values[0] += shift; values[3] += shift; values[2] = z
                        lines.append('#box: ' + ' '.join(f'{v:g}' for v in values) + ' cover')
                    if shift:
                        for i in range(48):
                            y0, y1 = i*.25, (i+1)*.25
                            lines += [f'#box: 0 {y0:g} {ztable[0,i]:g} 2 {y1:g} 12 cover',
                                      f'#box: 14 {y0:g} {ztable[-1,i]:g} 16 {y1:g} 12 cover']
                else:
                    lines.append(f'#box: 0 0 {9.05 if variant == "flat" else 0:g} {12+2*shift} 12 12 cover')
                lines.append(f'#geometry_view: 0 0 0 {12+2*shift} 12 33 0.05 0.05 0.05 hs4_geom n')
                name = station + '_' + variant
                groups[stage].append(case_record(folder, stage, name, '\n'.join(lines), source))
    result = {'status': 'PREPARED_NOT_RUN', 'solver_executed': False,
              'generator_sha256': sha256(__file__), 'stages': groups,
              '3d_bin6_relief_m': float(np.ptp(ztable[6])), '3d_global_z_range_m': [float(ztable.min()), float(ztable.max())],
              '3d_global_relief_m': float(np.ptp(ztable)),
              '3d_polarization': 'ideal x electric point dipole; Ex receiver; inline baseline along scan Y',
              'equipment_mismatch': 'Project instrument baseline is cross-track; this inherited inline setup is a mechanism control, not verified equipment model.',
              'geometry': '2*(archived Z-9.05)+9.05; full 48x48 surface; bin6 equals executed2D0.8m transect. Xwide translates entire original interior+2m and continues edge profiles into both new margins.',
              'resource_floors_GiB': {'replay': [2.25, 3.5], 'centre1cm': [3.75, 4.3],
                                     '3d-centre': [8, 5], '3d-ends': [8, 5], '3d-xwide-centre': [10, 6.5]},
              'processing': 'native source normalization,501tones20-170MHz,200ns tail,Hann/8x; complex difference before magnitude; fixed160-180/180-220/160-220ns. No imaging or training.'}
    save(folder/'task.json', result)
    print(json.dumps({'stages': {k: len(v) for k,v in groups.items()}, '3d_bin6_relief_m':result['3d_bin6_relief_m']}))


def verify_design(folder):
    task = read(folder/'task.json')
    for groups in task['stages'].values():
        for g in groups:
            if sha256(folder/g['input']) != g['input_sha256'] or sha256(ROOT/g['source']) != g['source_sha256']:
                raise ValueError('design/source changed')
    return task


def freeze(design, stage, out, replay):
    task = verify_design(design)
    runtime = version_runtime()
    if out.exists():
        raise ValueError('new execution capsule required')
    if stage != 'replay':
        if replay is None:
            raise ValueError('run and analyze target-machine replay first; --replay-report required')
        report = read(replay)
        if report['stage'] != 'replay' or report['runtime_source_identities'] != runtime['source_identities']:
            raise ValueError('replay belongs to another runtime')
    out.mkdir(parents=True)
    groups = []
    for g in task['stages'][stage]:
        g = dict(g); target = out/g['id']/'profile.in'; target.parent.mkdir()
        target.write_bytes((design/g['input']).read_bytes()); g['input'] = target.relative_to(out).as_posix()
        groups.append(g)
    ram, vram = task['resource_floors_GiB'][stage]
    save(out/'execution_contract.json', {'status':'FROZEN_APPROVED', 'stage':stage, 'groups':groups,
         'approval_basis':'User2026-10-04 autonomous investigation and explicit complete cross-PC Git handoff.',
         'runtime':runtime, 'design_sha256':sha256(design/'task.json'), 'max_runs':len(groups),
         'replay_report_sha256':sha256(replay) if replay else None,
         'code_identities':{name:sha256(ROOT/name) for name in
                            ('scripts/hs4_cross_pc.py','scripts/hs4_station_grid_controls.py','scripts/gprmax_cached_cuda_entry.py')},
         'min_available_RAM_GiB':ram, 'min_free_VRAM_GiB':vram,
         'max_owned_RSS_GiB':ram, 'min_system_available_during_run_GiB':.25,
         'max_group_wall_s':7200, 'max_batch_wall_s':7200*len(groups), 'no_retry':True})
    print('Frozen on target:', stage, len(groups), 'traces; resource floors', ram, vram)


def audit(out, require_geometry=True):
    c = read(out/'execution_contract.json'); records = []
    for g in c['groups']:
        p = out/g['input']
        if sha256(p) != g['input_sha256']:
            raise ValueError('input changed')
        dx = np.asarray(g['spacing_m']); shape=np.asarray(g['cells']); raw=p.with_suffix('.h5')
        with h5py.File(raw) as h:
            v = h['rxs/rx1/'+g['component']][:]
            dt = float(h.attrs['dt']); expected_dt=1/(299792458*np.sqrt(np.sum(1/dx[shape>1]**2)))
            if str(h.attrs['gprMax']) != '4.0.0' or not np.array_equal(h.attrs['dx_dy_dz'],dx) or not np.array_equal(h.attrs['nx_ny_nz'],shape):
                raise ValueError('native version/grid differs')
            if not np.isclose(dt,expected_dt,rtol=1e-12,atol=0) or v.dtype != np.float64 or not np.isfinite(v).all() or len(v)!=int(h.attrs['Iterations']):
                raise ValueError('native dtype/time/finite check failed')
            for key, position in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                if not np.array_equal(h[key].attrs['GridPosition'],np.rint(np.asarray(position)/dx).astype(int)):
                    raise ValueError('actual pose differs')
        geom=p.parent/('hs4_geom.vtkhdf' if g['component']=='Ex' else 'hs4t2d_geom.vtkhdf')
        record={'id':g['id'],'raw_sha256':sha256(raw),'dt_s':dt,'dtype':str(v.dtype),'actual_pose_check':'PASS'}
        if require_geometry:
            boxes=[]
            for words in commands(p)['#box:']:
                limits=np.rint(np.asarray(list(map(float,words[:6]))).reshape(2,3)/dx).astype(int)
                boxes.append((limits, {'rock':3,'cover':4}[words[-1]]))
            with h5py.File(geom) as h:
                actual=h['VTKHDF/CellData/Material']
                if actual.shape != tuple(shape[::-1]):
                    raise ValueError('material output shape differs')
                for z0 in range(0,int(shape[2]),16):
                    z1=min(z0+16,int(shape[2])); expected=np.full((z1-z0,shape[1],shape[0]),2,dtype=np.uint32)
                    for (lo,hi), material in boxes:
                        a,b=max(z0,int(lo[2])),min(z1,int(hi[2]))
                        if a<b: expected[a-z0:b-z0,lo[1]:hi[1],lo[0]:hi[0]]=material
                    if not np.array_equal(actual[z0:z1],expected):
                        raise ValueError('actual full material map differs')
            record['geometry_sha256']=sha256(geom)
        records.append(record)
    return {'status':'PASS','stage':c['stage'],'contract_sha256':sha256(out/'execution_contract.json'),
            'runtime_source_identities':c['runtime']['source_identities'],'groups':records,
            'scope':'Execution/raw/geometry evidence; no absolute convergence, finite-antenna or field acceptance.'}


def run(out):
    c=read(out/'execution_contract.json'); current=version_runtime()
    if c['status']!='FROZEN_APPROVED' or c['max_runs']!=len(c['groups']) or c['runtime']['source_identities']!=current['source_identities'] or Path(c['runtime']['python']).resolve()!=Path(sys.executable).resolve():
        raise ValueError('frozen target runtime required')
    for name,digest in c['code_identities'].items():
        if sha256(ROOT/name)!=digest: raise ValueError('frozen code changed')
    for g in c['groups']:
        if sha256(out/g['input'])!=g['input_sha256']: raise ValueError('frozen input changed')
    log=out/'execution.jsonl'
    if log.exists() or any((out/g['input']).with_suffix('.h5').exists() for g in c['groups']):
        raise ValueError('consumed attempt; preserve it and freeze a new justified task')
    def capacity():
        r=resources()
        if r['available_RAM_bytes']<c['min_available_RAM_GiB']*2**30 or r['free_VRAM_bytes']<c['min_free_VRAM_GiB']*2**30:
            raise RuntimeError('insufficient live RAM/VRAM; no solver started')
        return r
    capacity()
    lockpath=ROOT/'artifacts/local_checks/hs4_gpu_exclusive.lock'; lockpath.parent.mkdir(parents=True,exist_ok=True)
    with lockpath.open('a+b') as lock:
        if os.name=='nt':
            import msvcrt
            lock.seek(0); msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        def event(value):
            with log.open('a',encoding='utf-8') as f: f.write(json.dumps(value,allow_nan=False)+'\n')
        batch=time.monotonic(); process=None
        event({'status':'STARTED','contract_sha256':sha256(out/'execution_contract.json'),'unix_s':time.time()})
        try:
            for g in c['groups']:
                p=(out/g['input']).resolve(); r=capacity()
                command=[sys.executable,str(ROOT/'scripts/gprmax_cached_cuda_entry.py'),str(p),'-gpu','0','-gpu_precision','double','--hide-progress-bars']
                env=os.environ.copy(); env['HS4_CUDA_CACHE_LOG']=str(p.parent/'cuda_cache.jsonl')
                event({'status':'STARTED','group':g['id'],**r})
                with (p.parent/'stdout.log').open('xb') as stdout,(p.parent/'stderr.log').open('xb') as stderr:
                    process=subprocess.Popen(command,cwd=p.parent,env=env,stdout=stdout,stderr=stderr)
                    start=time.monotonic(); peak=0
                    while process.poll() is None:
                        try:
                            parent=psutil.Process(process.pid); owned=[parent]+parent.children(recursive=True); rss=0
                            for child in owned:
                                try: rss+=child.memory_info().rss
                                except psutil.NoSuchProcess: pass
                            peak=max(peak,rss)
                        except psutil.NoSuchProcess: pass
                        if peak>c['max_owned_RSS_GiB']*2**30 or psutil.virtual_memory().available<.25*2**30 or time.monotonic()-start>c['max_group_wall_s'] or time.monotonic()-batch>c['max_batch_wall_s']:
                            terminate_owned_tree(process); raise RuntimeError('owned solver resource/wall guard')
                        time.sleep(.5)
                    if process.returncode: raise RuntimeError('solver failed; no automatic retry')
                event({'status':'COMPLETED','group':g['id'],'raw_sha256':sha256(p.with_suffix('.h5')),'elapsed_s':time.monotonic()-start,'peak_owned_RSS_bytes':peak})
                print('Completed',g['id'],flush=True); process=None
            verification=audit(out); save(out/'completed_verification.json',verification)
            event({'status':'COMPLETED','traces':len(c['groups']),'verification_sha256':sha256(out/'completed_verification.json')})
        except BaseException as exc:
            if process is not None and process.poll() is None: terminate_owned_tree(process)
            event({'status':'FAILED','error':str(exc)}); raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('prepare','preflight','freeze','run','verify'))
    p.add_argument('--design',type=Path); p.add_argument('--stage',choices=STAGES)
    p.add_argument('--out',type=Path,required=True); p.add_argument('--capsule',type=Path)
    p.add_argument('--replay-report',type=Path); p.add_argument('--execute',action='store_true')
    a=p.parse_args()
    if a.action=='prepare': prepare(a.out)
    elif a.action=='freeze': freeze(a.design,a.stage,a.out,a.replay_report)
    elif a.action=='run':
        if not a.execute: raise ValueError('--execute required')
        run(a.out)
    elif a.action=='verify':
        if a.out.exists(): raise ValueError('new audit output required')
        save(a.out,audit(a.capsule))
    else:
        if a.out.exists(): raise ValueError('new preflight output required')
        runtime=version_runtime(); save(a.out,{'runtime':runtime,**resources(),
            'compilers':{name:__import__('shutil').which(name) for name in ('nvcc','cl','gcc')},
            'gpu':subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader'],text=True).strip()})


if __name__=='__main__':
    main()
