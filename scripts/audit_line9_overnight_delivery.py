"""Independent delivery summary: input bytes, trace counts, tail advisory and peak swaps."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from hs4_station_grid_controls import commands
from line9_v401_version_controls import save


def main(a):
    result_path = a.public/'delivery_verification.json'
    assert not result_path.exists()
    analysis = json.loads((a.public/'analysis.json').read_text('utf-8'))
    c = json.loads((a.raw/'execution_contract.json').read_text('utf-8'))
    v = json.loads((a.raw/'completed_verification.json').read_text('utf-8'))
    events = [json.loads(s) for s in (a.raw/'execution.jsonl').read_text('utf-8-sig').splitlines()]
    done = [r for r in events if r['status'] == 'COMPLETED' and 'group' in r]
    assert len(done) == len(v['groups']) == len(c['groups']) == 42
    assert not any(r['status'] in ['FAILED','GUARD_STOP'] for r in events)
    assert v['contract_sha256'] == analysis['contract_sha256'] == sha(a.raw/'execution_contract.json')
    assert events[-1]['verification_sha256'] == sha(a.raw/'completed_verification.json')
    m = c['study_manifest']; templates = json.loads((a.parent/'manifest.json').read_text('utf-8'))
    templates = {g['role']: g for g in templates['groups']}
    invariants = ['#domain_mode:','#domain:','#dx_dy_dz:','#time_window:','#omp_threads:',
                  '#pml_cells:','#pml_formulation:','#dispersive_averaging:','#waveform:']
    for g in m['groups']:
        old = templates[g['role']]
        x, y = commands(a.package/g['input']), commands(a.parent/old['input'])
        assert all(x[k] == y[k] for k in invariants)
        assert g['geometry_sha256'] == old['geometry_sha256']
        assert g['material_sha256'] == old['material_sha256']
        for key in ['input','geometry','material']:
            assert sha(a.package/g[key]) == g[key+'_sha256']
        assert sha(a.raw/'inputs'/g['id']/'profile.in') == g['input_sha256']
    with h5py.File(a.numerical) as h:
        assert sha(a.numerical) == analysis['numerical_sha256']
        np.testing.assert_array_equal(h['profile_x_m'][:], np.arange(190.,169.99,-.5))
        assert h['H1_computed_mask'][:].all() and h['H1_response'].shape == (501,41)
        assert h['anchor_delta_response'].shape == (501,3)
        summaries = {}
        t = h['time_s'][:]*1e9; deep = (t>=250)&(t<=550)
        for window in ['hann','blackman']:
            z = h[f'{window}/H1_complex_bandpass'][:]
            peak_ns = t[deep][np.argmax(abs(z[deep]),axis=0)]
            np.testing.assert_array_equal(peak_ns, [r['deep_peak_ns'] for r in analysis['metrics'][window]['H1_peaks']])
            increases = np.flatnonzero(np.diff(peak_ns)>0)
            summaries[window] = dict(anchor_delta_peak_ns=[r['delta_peak_ns'] for r in analysis['metrics'][window]['matched_anchors']],
                strongest_peak_time_increases=[dict(from_profile_x_m=m['stations'][i]['profile_x_m'],
                    to_profile_x_m=m['stations'][i+1]['profile_x_m'],from_peak_ns=float(peak_ns[i]),
                    to_peak_ns=float(peak_ns[i+1])) for i in increases],
                strongest_peak_not_a_certified_basal_picker=True)
    rows = analysis['native_audits']; tail = [r['receiver_final_5pct_peak_relative_db'] for r in rows]
    result = dict(status='PASS_COMPLETED_BATCH_AND_NUMERICAL_DELIVERY_WITH_TAIL_ADVISORY',
        script_sha256=sha(__file__),completed_new_runs=42,reused_runs=2,H1_stations=41,matched_anchors=3,
        sum_case_wall_s=sum(r['elapsed_s'] for r in done),batch_start_unix_s=events[0]['unix_s'],
        peak_owned_RSS_GiB=max(r['peak_owned_RSS_bytes'] for r in done)/2**30,
        native_identity_precision='float64',all_42_input_invariants_match_parent=True,
        archive_sha256=sha(Path(str(a.raw)+'.zip')),retrieved_files_checked=analysis['retrieved_files_checked'],
        official_CLI_complex_profile_relative_L2=analysis['official_CLI_complex_profile_relative_L2'],
        max_independent_DFT_relative_L2=max(r['relative_L2'] for r in rows),
        tail_relative_db_range=[min(tail),max(tail)],native_records_above_minus60_db=sum(v>-60 for v in tail),
        H1_records_above_minus60_db=sum(r['role']=='H1' and r['receiver_final_5pct_peak_relative_db']>-60 for r in rows),
        summaries=summaries,limits='Tail advisory not satisfied at every station. Strongest wide-window peak can swap between returns; do not label all maxima basal. Three anchor contrasts are not a full-line clean reference. Spatial/mesh/time convergence and field calibration not certified.')
    save(result_path,result); print(json.dumps(result,ensure_ascii=False))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['public','raw','numerical','package','parent']:
        p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
