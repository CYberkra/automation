"""Allocation lower bound for full-profile extruded3D; deliberately no solver."""
import argparse,json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    import gprMax
    source=Path(gprMax.__file__).parent
    files=['grid/cuda_grid.py','grid/fdtd_grid.py','config.py']
    rows=[]
    for dl in [.025,.05]:
        for width in [10.,20.]:
            n=np.rint(np.array([210.,42.5,width])/dl).astype(int);nodes=int(np.prod(n+1))
            # Installed CUDA copies six float64 E/H arrays in full; no PML,
            # material ID or Debye-history estimate needed for this lower bound.
            lower=6*8*nodes
            rows.append({'spacing_m':dl,'domain_m':[210.,42.5,width],'cells_xyz':n.tolist(),
                         'field_nodes_per_component':nodes,'six_fields_only_bytes':lower,'six_fields_only_GiB':lower/2**30,
                         'fits_16GiB_even_before_ID_Debye_PML':lower<16*2**30,
                         'mesh_source_receiver_geometry_equivalent_to_current':dl==.025})
    result={'status':'FULL_PROFILE3D_FIELD_ARRAY_LOWER_BOUND_ONLY_NO_EXECUTION',
            'script_sha256':sha(__file__),'installed_gprMax':gprMax.__version__,
            'allocation_source_identities':{p:sha(source/p) for p in files},'rows':rows,
            'precision':'float64','calls_solver':False,
            'limits':'Lower bound from installed six component field allocations. IDs, dispersion, PML, histories and compilation increase capacity further. 5cm changes grid and Rx y alignment; not approved equivalent replacement. No cheap cropped3D claimed equivalent.'}
    assert not a.out.exists();a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
