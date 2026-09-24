"""Hypothetical geology propagation and dense-grid budgets; no solver or field data."""
import argparse, hashlib, json, math
from pathlib import Path
from study_layer_loss_budget import propagation, C0


def memory(domain, grid, real_bytes):
    n=[round(d/h) for d,h in zip(domain,grid)]
    assert all(math.isclose(d/h,k,abs_tol=1e-8) for d,h,k in zip(domain,grid,n))
    nx,ny,nz=n;tx,ty,tz=[math.ceil(1/h) for h in grid]
    cells=math.prod(n);nodes=math.prod(k+1 for k in n)
    histories=2*((2*tx+1)*(ny*(nz+1)+(ny+1)*nz)
                 +(2*ty+1)*(nx*(nz+1)+(nx+1)*nz)
                 +(2*tz+1)*(nx*(ny+1)+(nx+1)*ny))
    gpu=(6*real_bytes+24)*nodes+histories*real_bytes
    return dict(domain_m=domain,grid_m=grid,cells=cells,PML_cells=[tx,ty,tz],
                precision='double' if real_bytes==8 else 'single',
                GPU_main_arrays_GiB=gpu/2**30,host_main_arrays_GiB=(gpu+22*cells)/2**30,
                baseline_1p3m_exact_on_y_grid=math.isclose(1.3/grid[1],round(1.3/grid[1]),abs_tol=1e-8))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    materials=[('cover',16,.01),('sandstone',9,.001),('wet_patch',20,.02)]
    spectra=[]
    for name,er,sigma in materials:
        for f in (20e6,80e6,170e6):
            g=propagation(f,er,sigma);lam=2*math.pi/g.imag
            spectra.append(dict(material=name,er=er,sigma_S_m=sigma,frequency_Hz=f,
                alpha_Np_m=g.real,phase_wavelength_m=lam,
                ten_cells_max_spacing_m=lam/10,
                absorption_dB_per_m_one_way=20/math.log(10)*g.real))
    loss=[]
    for sigma in (.001,.01,.05):
        for f in (20e6,80e6,170e6):
            cover=20/math.log(10)*6*propagation(f,16,sigma).real
            rock=20/math.log(10)*34*propagation(f,9,.001).real
            loss.append(dict(cover_sigma_S_m=sigma,frequency_Hz=f,cover_3m_two_way_dB=cover,
                             rock_17m_two_way_dB=rock,total_absorption_dB=cover+rock))
    rows=[]
    for domain in ((24,24,44),(24,24,48)):
        for grid in ((.05,.05,.04),(.04,.04,.04),(.04,.025,.04),(.025,.025,.025)):
            for real in (8,4):
                r=memory(domain,grid,real)
                r['minimum_cells_per_phase_wavelength_170MHz']={name:min(2*math.pi/propagation(170e6,er,sigma).imag/h for h in grid) for name,er,sigma in materials}
                r['below_13GiB_main_array_budget']=r['GPU_main_arrays_GiB']<13
                r['all_materials_at_least10_cells']=min(r['minimum_cells_per_phase_wavelength_170MHz'].values())>=10
                r['passes_initial_dense_screen']=r['below_13GiB_main_array_budget'] and r['all_materials_at_least10_cells'] and r['baseline_1p3m_exact_on_y_grid']
                rows.append(r)
    # Reproduce completed dz4 main-array budget before trusting new sizes.
    check=memory((24,24,24),(.05,.05,.04),8)
    assert math.isclose(check['GPU_main_arrays_GiB'],10.379942424595356,rel_tol=1e-14)
    assert propagation(80e6,9,0).real==0 and all(x['total_absorption_dB']>0 for x in loss)
    assert not any(r['passes_initial_dense_screen'] for r in rows)
    result=dict(status='hypothetical_budget_not_execution_contract',solver_invoked=False,field_data_used=False,
        materials=spectra,loss_20m=loss,dense_grids=rows,
        geometry='Flat domain24x24x44: groundz24, sourcez39, target20m depth z4; 1m PML leaves3m bottom clearance. 48m vertical extent is a slope-envelope resource scenario only, not a frozen slope mesh.',
        slope_relief_24m_at10deg_m=24*math.tan(math.radians(10)),
        lossless_vertical_arrival_ns=dict(ground=30/C0*1e9,cover_base=(30+6*4)/C0*1e9,depth20=(30+6*4+34*3)/C0*1e9),
        absorption_exclusions=['interface transmission/reflection','spreading','antenna coupling','scattering','noise/dynamic range','dielectric relaxation'],
        memory_exclusions=['runtime/context','small coefficient/source/receiver arrays','construction temporaries','dispersive auxiliary fields','snapshots'],
        notes=['10 cells per phase wavelength is a screening heuristic, not accuracy certification.','Single precision keeps uint32 material IDs; GPU bytes do not halve.','Material values and10deg slope are hypotheses, not site estimates.','Arrival estimates ignore loss/dispersion/obliquity; old400ns taper would attenuate a20m return and must be redesigned.'],
        checks=dict(prior_dz4_memory_reproduced=True,zero_loss_limit=True,positive_absorption=True,no_dense_candidate_passes=True),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        propagation_script_sha256=hashlib.sha256(Path('scripts/study_layer_loss_budget.py').read_bytes()).hexdigest())
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result))


if __name__=='__main__':main()
