"""Freeze the coarse station schedule and evidence-based ETA; no solver."""
import argparse
import json
import math
import statistics
from pathlib import Path

from build_pdf_profile_geometry import digest, save_json

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT/'configs/research/line9_2d_first_run_task_v0_1.json'
REFERENCE = ROOT/'artifacts/research_checks/2026-10-05_slope_bscan_2d_fine2m_r1'


def main(package, evidence):
    if TASK.exists() or evidence.exists():
        raise ValueError('Fresh task and evidence paths required')
    manifest = json.loads((package/'manifest.json').read_text('utf-8'))
    geometry = manifest['geometries']['full2d']
    indices = sorted(set(range(0,391,4)) | {390})
    selected = [f'full2d_s{i:04d}' for i in indices]
    cases = {c['id']:c for c in manifest['cases']}
    pilots = [cases[name] for name in manifest['full2d_pilot_ids']]
    reused = {c['id']:next((p['id'] for p in pilots if p['tx_m']==c['tx_m'] and p['rx_m']==c['rx_m']),None)
              for c in [cases[n] for n in selected]}
    reused = {k:v for k,v in reused.items() if v is not None}
    pending = [n for n in selected if n not in reused]
    assert len(selected)==99 and len(pending)==97 and len(reused)==2
    distances = [cases[n]['acquisition_s_m'] for n in selected]
    assert distances[0]==0 and distances[-1]==195
    assert all(b-a==2 for a,b in zip(distances[:-2],distances[1:-1])) and distances[-1]-distances[-2]==1
    old = json.loads((REFERENCE/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (REFERENCE/'execution.jsonl').read_text('utf-8').splitlines()]
    measured = [r['elapsed_s'] for r in events if r['status']=='COMPLETED' and 'group' in r]
    assert len(measured)==26 and old['status']=='PASS'
    stdout = REFERENCE/'slope_rough_s01/stdout.log'
    assert 'NVIDIA GeForce RTX 4090 Laptop GPU' in stdout.read_text('utf-8')
    old_cells = 2880*2520
    old_steps = old['groups'][0]['iterations']
    dt = .025/(299792458*math.sqrt(2))
    # Same V4 nearest iteration convention; future native output is authoritative.
    new_steps = round(1.2e-6/dt)+1
    factor = geometry['cells']*new_steps/(old_cells*old_steps)
    seconds = statistics.median(measured)*factor
    eta = dict(status='COMPUTE_WORK_EXTRAPOLATION_NOT_MEASURED_NEW_MODEL',
        reference_gpu='NVIDIA GeForce RTX 4090 Laptop GPU', reference_completed_traces=26,
        reference_median_wall_s=statistics.median(measured),reference_cells=old_cells,
        reference_iterations=old_steps,new_cells=geometry['cells'],new_CFL_dt_s=dt,
        new_initial_iterations=new_steps,cell_iteration_ratio=factor,
        nominal_seconds_per_new_trace=seconds,nominal_full391_hours=391*seconds/3600,
        nominal_preview99_hours=99*seconds/3600,nominal_preview_with_pilots102_hours=102*seconds/3600,
        uncertainty='Cell-iteration scaling only; different material count/PML/import costs and GPU clocks. Desktop4090/4090D may be faster but no measured conversion factor. Replace with this model pilot timings.',
        source_sha256={str(p.relative_to(ROOT)):digest(p) for p in [REFERENCE/'execution.jsonl',REFERENCE/'completed_verification.json',stdout]})
    task = dict(schema='line9-first-2d-task/1',status='APPROVED_TASK_NOT_MACHINE_EXECUTION_CONTRACT',
        date='2026-10-06',authorization='User: detailed Git task including SFCW, first run2D; continue and moderately downsample first scan',
        role='Line9 site adaptation, not independent blind validation',
        package_manifest_sha256=digest(package/'manifest.json'),
        portable_archive_sha256='c0ca87b431b48538cf8baa728483a8d6cc69de1a5165c8f5c70cf950aadd925b',
        acquisition=dict(direction='X220 to25; s=220-X',full_station_count=391,preview_station_count=99,
            nominal_preview_spacing_m=2,final_interval_m=1,pilot_ids=[p['id'] for p in pilots],
            preview_ids=selected,reuse_pilot_output=reused,preview_pending_ids=pending,
            maximum_unique_FDTD_runs_initial=102,full391_is_not_initial_batch=True),
        invariants=dict(solver_version='4.0.0',domain_m=[400,75,.025],grid_xyz_m=[.025]*3,
            raw_field_dtype='float64',backend='CUDA',CPU_fallback=False,time_window_s=1.2e-6,
            waveform='Ricker100MHz,40A on z current element0.025m, nominal1Am',
            polarization='z',receiver='rxs/rx1/Ez',midpoint_AGL_m=15,baseline_m=1.3,
            baseline_role='ALONG_TRACK_2D_SURROGATE',geometry_sha256=geometry['sha256'],
            material_database_sha256=manifest['material_database_sha256'],
            PML='HORIPML;80/80/0/80/80/0 cells;2m on active axes; averaging n'),
        SFCW=dict(min_hz=20e6,max_hz=170e6,step_hz=.3e6,count=501,
            estimator='V4 official direct_frequency_response; complex Y/X/source.spatial_scale',
            source_floor_db=-100,all_tones_valid_required=True,source_length_m=.025,
            phase_and_native_time_offsets_preserved=True,extra_time_shift_s=0,
            windows=['rectangular','hann'],zero_pad_factor=8,tail_taper_fraction_primary=0,
            tail_sensitivity_ns=20,frequency_downsampling=False,gain=False,background_removal=False,
            quantity='Field response per current moment; not antenna-port S21',
            signed_display='Real part of remodulated complex profile; complex array retained; not envelope used as signed data'),
        pilot_review=dict(required='Native audits +501valid tones +independent DFT/inverse +Chinese total-field panels +record-tail diagnostics +timing and memory. Freeze preview only after pilot_review.json APPROVE_PREVIEW.',
            automatic_physical_acceptance=False,tail_is_diagnostic_not_complete_record_proof=True),
        limits=dict(pilot_max_seconds_per_trace=3600,pilot_max_batch_seconds=18000,
            preview_max_seconds_per_trace=3600,preview_max_batch_seconds=172800,no_retry=True,
            min_available_RAM_bytes=14968955386,min_free_VRAM_bytes=12522413562,
            max_owned_RSS_GiB=24,min_system_available_during_run_GiB=1.5,
            concurrent_solver_count=1,no_3d_in_this_task=True),
        sampling_limit='2m is coarse morphology preview, not spatial-Nyquist certification or imaging-ready sampling; no interpolated columns counted as computed.',
        estimated_cost=eta)
    TASK.parent.mkdir(parents=True,exist_ok=True)
    save_json(TASK,task)
    evidence.mkdir(parents=True)
    save_json(evidence/'planning_evidence.json',dict(task_sha256=digest(TASK),calls_solver=False,
        script_sha256=digest(Path(__file__)),schedule_checks='PASS99 preview +5pilots -2reused=102 unique',estimated_cost=eta))
    print(json.dumps(dict(stations=99,unique_total=102,full391_hours=eta['nominal_full391_hours'],
        preview99_hours=eta['nominal_preview99_hours'],reference_GPU=eta['reference_gpu']),ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,required=True)
    p.add_argument('--evidence',type=Path,required=True)
    a=p.parse_args();main(a.package,a.evidence)
