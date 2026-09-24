"""Static double-precision array budget for the bounded dz4cm pair; no solve."""
import argparse, hashlib, json, math
from pathlib import Path
import psutil

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--source-root',type=Path,required=True)
a=p.parse_args()
rows=[]
for nz,tz in ((480,20),(600,25)):
    nx=ny=480;tx=ty=20
    nodes=(nx+1)*(ny+1)*(nz+1);cells=nx*ny*nz
    # Official PML.initialise_field_arrays: four double histories per face,
    # including Yee padding; one CFS term, two faces per axis.
    px=(2*tx+1)*(ny*(nz+1)+(ny+1)*nz)
    py=(2*ty+1)*(nx*(nz+1)+(nx+1)*nz)
    pz=(2*tz+1)*(nx*(ny+1)+(nx+1)*ny)
    phi=2*(px+py+pz)*8
    gpu=72*nodes+phi
    host=22*cells+72*nodes+phi
    dz=24/nz;dt=1/(299792458*math.sqrt(2/.05**2+1/dz**2))
    rows.append(dict(dz_m=dz,cells=cells,dt_s=dt,steps_estimate=math.ceil(800e-9/dt)+1,
                     PML_history_bytes=phi,GPU_main_arrays_bytes=gpu,host_main_arrays_bytes=host,
                     GPU_main_arrays_GiB=gpu/2**30,host_main_arrays_GiB=host/2**30))
assert rows[1]['cells']==138240000 and rows[1]['GPU_main_arrays_GiB']<12
assert all(abs(z/.04-round(z/.04))<1e-9 for z in (4,19,24))
r=dict(solver_invoked=False,rows=rows,available_RAM_bytes=psutil.virtual_memory().available,
       source_sha256={s:hashlib.sha256((a.source_root/'gprMax'/s).read_bytes()).hexdigest()
                      for s in ('grid/fdtd_grid.py','grid/cuda_grid.py','pml.py')},
       exclusions=['small source/receiver/material/PML coefficient arrays','CUDA context/compiler overhead','host build temporaries and driver allocations'],
       empirical_host_peak_GiB=22430621696/2**30,
       scaled_host_peak_GiB=22430621696/2**30*1.25,
       limits=dict(job_commit_GiB=36,minimum_available_RAM_GiB=36,minimum_available_VRAM_GiB=13,wall_minutes_each=30),
       PML='Fixed sigma_max=0.21235349838321013 on every face in both models; same continuous quartic profile and 1m thickness as prior fixed-PML. Vertical samples change20->25; dt changes; not bitwise cross-grid coefficient equivalence.')
with a.output.open('x') as f:json.dump(r,f,indent=2)
print(json.dumps(r))
