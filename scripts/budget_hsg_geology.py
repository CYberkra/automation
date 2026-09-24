"""Reject oversized HSG candidates using rigorous field/ID storage lower bounds."""
import argparse,hashlib,json,math
from pathlib import Path
from budget_lossy_geology import memory
from study_layer_loss_budget import propagation

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-root',type=Path,required=True)
p.add_argument('--installed-root',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
files=['subgrids/grid.py','subgrids/user_objects.py','subgrids/cuda_subgrid_hsg.py','subgrids/cuda_subgrid_updates.py','subgrids/cuda_precursor_nodes.py']
hashes={n:hashlib.sha256((a.source_root/n).read_bytes()).hexdigest() for n in files}
assert all(hashlib.sha256((a.installed_root/n).read_bytes()).hexdigest()==h for n,h in hashes.items())
rows=[]
for label,coarse,extent in [('deep20m',.1,(20,20,22)),('wide_deep',.1,(22,22,22)),('cover_only_coarse_rock',.1,(20,20,3)),('cover_only_resolved_rock',.05,(20,20,3)),('local_mechanism_only',.05,(4,4,4))]:
 ratio=3;sep=3;pml=6;gap=ratio//2+2;pad=sep*ratio+gap+pml
 fine=coarse/ratio
 working=[round(d/fine) for d in extent]
 assert all(math.isclose(d/fine,n,abs_tol=1e-8) for d,n in zip(extent,working))
 total=[n+2*pad for n in working]
 # CUDA subgrid inherits full six fields and six uint32 ID arrays.
 sub_lower=72*math.prod(n+1 for n in total)
 domain=(24,24,28) if label=='local_mechanism_only' else (24,24,44)
 main=memory(domain,(coarse,coarse,coarse),8)
 rock_wavelength=2*math.pi/propagation(170e6,9,.001).imag
 wet_wavelength=2*math.pi/propagation(170e6,20,.02).imag
 rows.append(dict(label=label,main_domain_m=domain,full_depth20_candidate=label!='local_mechanism_only',main_spacing_m=coarse,ratio=ratio,fine_spacing_m=fine,working_extent_m=extent,
   pad_cells_per_face=pad,total_subgrid_cells_xyz=total,main_arrays_GiB=main['GPU_main_arrays_GiB'],
   subgrid_fields_ID_only_GiB=sub_lower/2**30,combined_lower_bound_GiB=main['GPU_main_arrays_GiB']+sub_lower/2**30,
   unrefined_rock_cells_per_wavelength=rock_wavelength/coarse,fine_wet_cells_per_wavelength=wet_wavelength/fine,
   excluded=['subgrid PML histories','precursor device histories and interpolation/filter buffers','coefficient tables','CUDA context/compiler overhead'],
   placement_validated=False))
assert rows[0]['pad_cells_per_face']==18
assert rows[0]['combined_lower_bound_GiB']>16
assert rows[2]['unrefined_rock_cells_per_wavelength']<10
# Local cubic candidate: exact PML history count, conservative linear precursor
# buffer bound (24 faces/components; use padded subgrid width also for scratch).
n=276;t=6
local_pml_bytes=2*3*(2*t+1)*(2*n*(n+1))*8
local_precursor_bound=24*(3*(240+1)**2+(n+1)**2)*8+24*4*(n+1)*8
rows[-1]['subgrid_PML_history_GiB']=local_pml_bytes/2**30
rows[-1]['linear_precursor_buffer_bound_GiB']=local_precursor_bound/2**30
rows[-1]['accounted_arrays_with_precursor_bound_GiB']=rows[-1]['combined_lower_bound_GiB']+(local_pml_bytes+local_precursor_bound)/2**30
rows[-1]['still_excluded']=['small coefficient tables','CUDA runtime/compiler/context','host build temporaries','unverified allocation branches']
r=dict(solver_invoked=False,large_grid_allocated=False,installed_sources_match=True,source_sha256=hashes,
 rows=rows,checks=dict(default_padding18=True,deep_candidate_exceeds16GiB_before_extra_arrays=True,cover_only_coarse_rock_fails_sampling_screen=True),
 scope='Necessary storage bounds, not a complete HSG memory estimate or proof against all configurations. No placement/material-coupling accuracy certification.',
 next='Do not run oversized full-depth candidates. Need a smaller explicitly scoped mechanism model or verified symmetry/equivalent-source reduction or larger GPU before full-depth production.')
with a.output.open('x',encoding='utf-8') as f:json.dump(r,f,indent=2)
print(json.dumps(r))
