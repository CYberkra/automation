"""Static axial dispersion and reduced-grid resource estimates; no solver launch."""
import argparse,json,math,hashlib
from pathlib import Path
from scipy.constants import c

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
f=170e6;omega=2*math.pi*f
def axial_phase_error(dy,dz):
 dt=1/(c*math.sqrt(dy**-2+dz**-2));phase=0
 for er,length in [(1,15),(16,3),(9,16.75)]:
  kh=2/dz*math.asin(math.sqrt(er)*dz/(c*dt)*math.sin(omega*dt/2));phase-=2*(kh-omega*math.sqrt(er)/c)*length
 return math.degrees(phase)
base_dt=1/(c*math.sqrt(2/.0125**2));rows=[]
for dy,dz in [(.025,.025),(.0125,.0125),(.00625,.00625),(.005,.005),(.003125,.003125),(.0125,.00625),(.0125,.005),(.0125,.003125)]:
 ny,nz=round(32/dy),round(50/dz);dt=1/(c*math.sqrt(dy**-2+dz**-2));cells=ny*nz;gpu=72*2*(ny+1)*(nz+1)/2**30
 aligned=all(math.isclose(v/dy,round(v/dy),abs_tol=1e-8) for v in [15.35,16.65,14,18]) and all(math.isclose(v/dz,round(v/dz),abs_tol=1e-8) for v in [45,30,27,9.75,10.25])
 rows.append({'dy_m':dy,'dz_m':dz,'cells':cells,'dt_s':dt,'main_field_ID_GiB':gpu,'geometry_aligned':aligned,'axial_lossless_170MHz_phase_error_deg':axial_phase_error(dy,dz),'wall_s_work_scaled_from_FINE_target':122.187*cells/10240000*base_dt/dt,'main_arrays_below_13GiB_screen':gpu<13})
assert all(row['geometry_aligned'] for row in rows)
prediction=rows[1]['axial_lossless_170MHz_phase_error_deg']-rows[0]['axial_lossless_170MHz_phase_error_deg'];assert abs(prediction-63.604309166983626)<1e-8
result={'rows':rows,'checks':{'all_geometries_exact':True,'previous_phase_prediction_reproduced':True},'solver_invoked':False,'execution_authorized_by_this_file':False,
 'limitations':['Memory is main six double field/six uint32 ID arrays with two invariant-axis Yee nodes only; excludes PML, coefficients, CUDA context and host construction.','Runtime is cell-step scaling from one measured run, not benchmark or guaranteed ETA.','Axial lossless phase excludes target scattering, oblique propagation and conduction; anisotropic mesh still has lateral dispersion.','13GiB is only a preliminary main-array screen, not a physical or runtime acceptance threshold.','No user phase tolerance invented. Full20-170MHz remains in scope.'],
 'decision':'Assess directional depth refinement as an affordable diagnostic; do not assert full convergence from axial prediction.', 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
with a.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
print(json.dumps(result))
