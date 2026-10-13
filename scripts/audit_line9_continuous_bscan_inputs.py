"""Read-only compatibility audit for the requested51-trace mesh comparison."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from study_line9_v31_smoothing_domain import load, save, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['old-bscan','smoothed','coarse','out']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    original=load(a.old_bscan/'resolved_contract.json');done=load(a.old_bscan/'completion.json')
    fine=load(a.smoothed/'contract.json');coarse=load(a.coarse/'contract.json')
    rows=[];fine_reusable=[];coarse_reusable=[]
    for r in done['receipts']:
        s=r['station_global_x_m']
        if s>180+1e-8:
            continue
        assert s>=150-1e-8
        raw=a.old_bscan/'results'/f"line9_v31_bscan{r['trace_index']}.h5"
        assert sha(raw)==r['output_sha256']
        y=np.array([r['tx_local_m'][1],r['rx_local_m'][1]])
        common=np.floor(y/.1+.5+1e-9)*.1
        rows.append(dict(station_m=s,old_raw_sha256=sha(raw),old_fine_y_m=y.tolist(),
            both_old_y_on_coarse_nodes=bool(np.all(np.isclose(y/(1/30),np.round(y/(1/30)),rtol=0,atol=1e-7))),
            proposed_common_y_m=common.tolist(),proposed_y_change_m=(common-y).tolist(),
            proposed_common_domain_x_origin_m=s-25.05,proposed_common_local_tx_rx_x_m=[24.4,25.7],
            proposed_actual_global_tx_rx_x_m=[s-.65,s+.65],proposal_only_not_authorized=True))
    assert len(rows)==51 and np.allclose(np.diff([r['station_m'] for r in rows]),.6,rtol=0,atol=1e-9)
    for case in fine['cases']:
        if not(case['id'].startswith('full_x') and case['averaging']=='y' and not case['geometry_only'] and case['station_m']<=180):
            continue
        rec=load(a.smoothed/'cases'/case['id']/'receipt.json')
        assert sha(a.smoothed/'cases'/case['id']/'impulse.h5')==rec['raw_sha256']
        fine_reusable.append(dict(station_m=case['station_m'],raw_sha256=rec['raw_sha256'],
            build_solver_seconds=rec['official_build_and_solver_seconds'],strict_proposed_coordinates_match=False))
    for case in coarse['cases']:
        rec=load(a.coarse/'cases'/case['id']/'receipt.json')
        assert sha(a.coarse/'cases'/case['id']/'impulse.h5')==rec['raw_sha256']
        coarse_reusable.append(dict(station_m=case['station_m'],raw_sha256=rec['raw_sha256'],
            build_solver_seconds=rec['official_build_and_solver_seconds'],strict_proposed_coordinates_match=False))
    assert len(fine_reusable)==4 and len(coarse_reusable)==3
    for r in rows:
        for h in [.025,1/30]:
            assert np.allclose(np.array(r['proposed_common_y_m'])/h,np.round(np.array(r['proposed_common_y_m'])/h),atol=1e-8,rtol=0)
            assert np.allclose(np.array([24.4,25.7])/h,np.round(np.array([24.4,25.7])/h),atol=1e-8,rtol=0)
    summary=dict(status='AWAITING_ACTUAL_COORDINATE_DECISION_NO_SOLVER_STARTED',calls_solver=False,
        new_fdtd=0,original_fine_records_verified=51,smoothed_fine_existing=4,coarse_existing=3,
        both_old_fine_y_on_coarse_nodes=sum(r['both_old_y_on_coarse_nodes'] for r in rows),
        nominal_stations_m=[150,180],station_step_m=.6,trace_count=51,
        proposed_strict_new_solves=102,proposed_max_vertical_change_m=max(abs(x) for r in rows for x in r['proposed_y_change_m']),
        median_fine_build_solver_seconds=float(np.median([r['build_solver_seconds'] for r in fine_reusable])),
        median_coarse_build_solver_seconds=float(np.median([r['build_solver_seconds'] for r in coarse_reusable])),
        proposal_execution_authorized=False,source_inputs_modified=False)
    summary['linear_solver_eta_minutes']=51*(summary['median_fine_build_solver_seconds']+summary['median_coarse_build_solver_seconds'])/60
    save(a.out/'compatibility_audit.json',dict(summary=summary,stations=rows,fine_existing=fine_reusable,coarse_existing=coarse_reusable,
        source_contract_sha256=sha(a.old_bscan/'resolved_contract.json'),audit_script_sha256=sha(__file__)))
    print(summary)


if __name__=='__main__':
    main()
