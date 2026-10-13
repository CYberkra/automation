"""Audit the bounded 45x45m follow-up against both frozen references."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

from analyze_line9_v31_smoothing_domain import (
    audit, coefficients, compare, independent_check, plot_sparse_compare, read_sf,
)
from study_line9_v31_smoothing_domain import load, save, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mid',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.root=a.root.resolve();a.mid=a.mid.resolve();a.out=a.out.resolve()
    a.out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    parent,results,receipts,parent_checks=audit(a.root,'B')
    contract=load(a.mid/'contract.json')
    assert contract['parent_contract_sha256']==sha(a.root/'contract.json')
    assert contract['domain_screen']==parent['domain_screen']
    assert len(contract['cases'])==7
    complete=load(a.mid/'stage_B_complete.json')
    assert complete['status']=='COMPLETE' and len(complete['receipts'])==7
    rows=[];wide_rows=[];checks={};refs=[];cands=[];stations=[]
    for r in complete['receipts']:
        assert r['status']=='PASS'
        for f in r['files']:
            path=a.mid/f['path'];assert path.stat().st_size==f['bytes'] and sha(path)==f['sha256']
        case=next(c for c in contract['cases'] if c['id']==r['case'])
        s=case['station_m'];refid=f'full_x{s:.1f}_y'
        candidate=read_sf(a.mid/'cases'/case['id']/'sfcw.h5')
        checks[case['id']]=independent_check(a.mid/'cases'/case['id']/'impulse.h5',candidate)
        ref=results[refid]
        coef=coefficients(a.mid,case['id']);orig=coefficients(a.root,refid)
        assert np.array_equal(coef['epsilon'],orig['epsilon']) and np.array_equal(coef['sigma'],orig['sigma'])
        measurement=compare(candidate,ref,case,parent['domain_screen'])
        rows.append(dict(case=case['id'],station_m=s,**measurement,
                         same_retained_surface_coefficients=True,
                         fdtd_speedup=receipts[refid]['official_fdtd_seconds']/r['official_fdtd_seconds'],
                         whole_run_speedup=receipts[refid]['official_build_and_solver_seconds']/r['official_build_and_solver_seconds']))
        if s in [150.,174.,198.]:
            wide=results[f'wide70_x{s:.1f}']
            wide_rows.append(dict(station_m=s,**compare(candidate,wide,case,parent['domain_screen'])))
        refs.append(ref['real']);cands.append(candidate['real']);stations.append(s)
    summary=dict(domain_m=[45.,45.],station_count=7,
                 same50_reference_all_pass=all(r['passes_frozen_domain_screen'] for r in rows),
                 direct70_reference_all_pass=all(r['passes_frozen_domain_screen'] for r in wide_rows),
                 passes_both_references=all(r['passes_frozen_domain_screen'] for r in rows+wide_rows),
                 max_same50_spectrum_error=max(r['complex_band_relative_l2'] for r in rows),
                 max_same50_window_error=max(w['complex_relative_l2'] for r in rows for w in r['windows']),
                 max_direct70_window_error=max(w['complex_relative_l2'] for r in wide_rows for w in r['windows']),
                 median_fdtd_speedup=float(np.median([r['fdtd_speedup'] for r in rows])),
                 median_whole_run_speedup=float(np.median([r['whole_run_speedup'] for r in rows])),
                 original50_guard_result_preserved='failed at 174m; not rewritten',
                 screening_tolerances_changed=False)
    save(a.out/'mid_comparisons.json',dict(summary=summary,against50=rows,against70=wide_rows,
                                         parent_contract_sha256=sha(a.root/'contract.json'),contract_sha256=sha(a.mid/'contract.json')))
    save(a.out/'audit.json',dict(status='HASHES_AND_INDEPENDENT_DFT_INVERSE_PASS',checks=checks,
                               parent_outputs_checked=len(parent_checks),mid_outputs_checked=len(checks),
                               total_new_solved_outputs=len(parent_checks)+len(checks),
                               script_sha256=sha(Path(__file__))))
    plot_sparse_compare(a.out,'mid45x45_sparse_bscan.png',stations,np.array(refs).T,np.array(cands).T,
                        ref['t'],'原生平滑 y｜50×50 m → 45×45 m｜只显示已计算7站')
    fig,ax=plt.subplots(figsize=(11,4.8),constrained_layout=True)
    ax.plot(stations,[100*max(w['complex_relative_l2'] for w in r['windows']) for r in rows],'-o',label='45×45 对 50×50：最差窗口')
    ax.plot([r['station_m'] for r in wide_rows],[100*max(w['complex_relative_l2'] for w in r['windows']) for r in wide_rows],'-s',label='45×45 对 70×50：最差窗口')
    ax.axhline(1,color='black',ls='--',label='保持原1%窗口容差')
    ax.set(xlabel='站位（m）',ylabel='复响应相对差（%）',title='中间尺寸直接核对两个有限域参照；不改原未通过结果')
    ax.grid(alpha=.2);ax.legend();fig.savefig(a.out/'mid_reference_errors.png',dpi=170);plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
