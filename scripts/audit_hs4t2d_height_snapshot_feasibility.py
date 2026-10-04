"""Read-only V4 identity and dense-snapshot storage estimate; never run FDTD."""
import argparse
import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('new output path required')
    old_contract = ROOT/'artifacts/research_checks/2026-10-03_hs4t2d_relief08_snapshot/execution_contract.json'
    old = json.loads(old_contract.read_text('utf-8'))
    package = Path(importlib.util.find_spec('gprMax').origin).parent
    checked = {}
    for name in ('snapshots.py', 'updates/cuda_updates.py', 'user_objects/cmds_output.py',
                 'cuda_opencl/knl_snapshots.py', 'utilities/host_info.py'):
        digest = sha256(package/name)
        if digest != old['source_identities'][name]:
            raise ValueError('snapshot implementation differs from old frozen runtime: '+name)
        checked[name] = digest
    raw = ROOT/'artifacts/research_checks/2026-10-04_hs4t2d_joint_grid_centre/centre_rough/profile.h5'
    with h5py.File(raw, 'r') as h:
        dt = float(h.attrs['dt'])
        n = len(h['rxs/rx1/Ey'])
        if str(h.attrs['gprMax']) != '4.0.0' or h['rxs/rx1/Ey'].dtype != np.float64:
            raise ValueError('V4 double reference required')
    hop = 10
    frames = len(range(0, n, hop))
    extent = [12, 0, 8, 24, .05, 28]
    spacing = [.1, .05, .1]
    shape = [round((extent[k+3]-extent[k])/spacing[k]) for k in range(3)]
    cells = int(np.prod(shape))
    frame_bytes = cells*6*8
    history_bytes = frame_bytes*frames
    report = {
        'status': 'DESIGN_FEASIBILITY_ONLY_NOT_EXECUTION_CONTRACT',
        'calls_solver': False,
        'script_sha256': sha256(Path(__file__)),
        'reference_raw': raw.relative_to(ROOT).as_posix(),
        'reference_raw_sha256': sha256(raw),
        'old_snapshot_contract_sha256': sha256(old_contract),
        'current_runtime_matches_old_frozen_snapshot_sources': checked,
        'fdtd_spacing_m': [.025, .05, .025],
        'dt_s': dt, 'native_iterations': n,
        'snapshot_native_hop': hop, 'snapshot_period_ns': hop*dt*1e9,
        'snapshot_temporal_nyquist_MHz': 1/(2*hop*dt)/1e6,
        'snapshot_frames': frames,
        'first_E_time_ns': 0, 'last_E_time_ns': (frames-1)*hop*dt*1e9,
        'extent_m': extent, 'output_spacing_m': spacing, 'output_shape_xyz': shape,
        'components_stored': ['Ex', 'Ey', 'Ez', 'Hx', 'Hy', 'Hz'],
        'six_field_frame_bytes': frame_bytes,
        'six_field_history_bytes_per_trace': history_bytes,
        'six_field_history_GiB_per_trace': history_bytes/2**30,
        'four_trace_payload_GiB_excluding_HDF5_headers': 4*history_bytes/2**30,
        'initial_plus_replacement_host_buffers_GiB': 2*history_bytes/2**30,
        'budget_scope': 'Array payload estimate only, not measured process peak. Add native model, Python, receiver probes and file overhead. Serial traces; release one solver before next.',
        'gpu_storage_scope': 'Actual V4 CUDA buffers allocate all six fields. Auto gpu-to-cpu policy can use one GPU frame but retains host history; do not assume streaming to disk or a CLI flag.',
        'required_pre_execution_checks': [
            'Freeze new inputs, source/runtime/analysis identities, run count, resource limits and no-retry contract.',
            'Use Ricker95 only after fine-grid full-complex transfer equivalence to archived impulse is checked; earlier source control was coarse.',
            'Native receiver and co-located probe histories must test hop10 versus hop5/full-native spectra; Nyquist alone does not certify no aliasing.',
            'Match snapshot cell-centre interpolation with native neighbouring E/H probes; do not compare an interpolated cell centre to an unrelated Yee receiver node.',
            'Use separate native E/H time levels in DTFT; collocate time before Poynting direction diagnostics.',
            'Reproduce current 501-tone source normalization, physical 200ns tail taper, Hann and carrier before receiver/B-scan attribution.',
            'Minimum planned available host RAM 3.5GiB and free VRAM 3.5GiB; recheck live and verify actual peak in pilot.',
            'ROI excludes outer PML and outside-profile propagation; absence of a visible return cannot exclude boundary reflections.'
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('status', 'snapshot_frames',
        'snapshot_period_ns', 'six_field_history_GiB_per_trace',
        'initial_plus_replacement_host_buffers_GiB')}, ensure_ascii=False))


if __name__ == '__main__':
    main()
