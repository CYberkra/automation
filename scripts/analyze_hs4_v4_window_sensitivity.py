"""CPU-only tail/Hann/rectangular sensitivity of completed material controls.

Diagnostic processing variants; the frozen production Hann chain stays intact.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from hs_capsule_identity import sha256
from hs4_v4_factor_controls import ROOT, CENTRE
from analyze_hs4_v4_factor_controls import compare
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response, reconstruct_time_response


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise ValueError('new diagnostic directory required')
    study = ROOT/'artifacts/research_checks/2026-10-04_hs4_v4_material_controls_r1'
    v = json.loads((study/'completed_verification.json').read_text('utf-8'))
    if not v['completed'] or v['status'] != 'PASS' or v['contract_sha256'] != sha256(study/'execution_contract.json'):
        raise ValueError('completed material study required')
    verify_official_runtime()
    native = {}
    for tag in ('debye_n', 'matched95', 'low_sigma'):
        for role in ('rough', 'halfspace'):
            p = CENTRE/f'centre_{role}'/'profile.h5' if tag == 'debye_n' else study/f'{tag}_centre_{role}'/'profile.h5'
            if tag == 'debye_n':
                digest = json.loads((p.parent/'audit.json').read_text('utf-8'))[0]['sha256']
            else:
                digest = next(g['raw_sha256'] for g in v['groups'] if g['id'] == f'{tag}_centre_{role}')
            if sha256(p) != digest:
                raise ValueError('native identity differs')
            native[tag, role] = (load_source(p), load_receiver(p, receiver_path='/rxs/rx1', component='Ey'))
    variants = [(200, 'hann'), (100, 'hann'), (0, 'hann'), (200, 'rectangular')]
    rows, arrays = {}, {}
    for tail_ns, window in variants:
        key = f'tail{tail_ns}_{window}'
        products = {}
        for tag in ('debye_n', 'matched95', 'low_sigma'):
            responses = []
            for role in ('rough', 'halfspace'):
                s, r = native[tag, role]
                fraction = 0 if tail_ns == 0 else (round(tail_ns*1e-9/r.dt)-.25)/len(r.samples)
                f = direct_frequency_response(s, r, FREQ, tail_taper_fraction=fraction)
                if not f.source_valid.all():
                    raise ValueError('all501 source bins required')
                responses.append(f)
            product = reconstruct_time_response(replace(responses[0], response=responses[0].response-responses[1].response), window=window, zero_pad_factor=8)
            products[tag] = product
            arrays[key+'_'+tag+'_complex_envelope'] = product.complex_envelope
        time = products['debye_n'].time*1e9
        rows[key] = {'conductivity_comparison': compare(products['matched95'], products['low_sigma'], time),
                     'spectrum_comparison': compare(products['debye_n'], products['matched95'], time)}
    report = {'status': 'COMPLETED_CPU_WINDOW_SENSITIVITY', 'code_sha256': sha256(__file__),
              'contract_sha256': v['contract_sha256'], 'solver_called': False, 'variants': rows,
              'production_parameters_changed': False, 'physical_attribution_certified': False,
              'scope': 'Checks whether material-effect conclusion depends on tail/window; no adaptive tuning or instrument-window certification.'}
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out/'comparison_arrays.npz', time_ns=time, **arrays)
    (args.out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({key: {name: values['underground'] for name, values in row.items()} for key, row in rows.items()}, indent=2))


if __name__ == '__main__':
    main()
