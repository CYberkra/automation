"""Audit the authorized 20–110MHz/5cm trial; reuse existing 2.5cm FDTD."""
from __future__ import annotations

import argparse
from pathlib import Path
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

from study_line9_v31_smoothing_domain import load, save, sha
from analyze_line9_v31_smoothing_domain import (
    read_sf, independent_check, compare, coefficients, plot_sparse_compare, relative)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'out']:
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DejaVu Sans'],
                         'axes.unicode_minus': False})
    c = load(a.root/'contract.json')
    complete = load(a.root/'stage_A_complete.json')
    start = load(a.root/'execution_start.json')
    assert complete['status'] == 'COMPLETE' and complete['cases'] == c['max_new_solves'] == 3
    assert sha(a.root/'contract.json') == start['contract_sha256']
    assert start['gprMax'] == c['gprMax'] == '4.0.1' and '4090' in start['gpu']
    assert len(complete['reused_references']) == 3
    assert c['record_ns'] == 1100 and c['official_sfcw']['steps'] == 301
    hashed = 0
    for r in complete['receipts'] + complete['reused_references']:
        assert r['status'] == 'PASS'
        for item in r['files']:
            f = a.root/item['path']
            assert f.stat().st_size == item['bytes'] and sha(f) == item['sha256']
            hashed += 1
    for item in load(a.root/'input_manifest.json')['files']:
        f = a.root/item['path']
        assert f.stat().st_size == item['bytes'] and sha(f) == item['sha256']
    rows, checks, refs, news, wides = [], {}, {}, {}, {}
    for case in c['cases']:
        cid = case['id']
        rdef = next(r for r in c['reuse_references'] if r['id'] == case['baseline_case'])
        rd = a.root/'references'/rdef['id']
        nd = a.root/'cases'/cid
        receipt = load(nd/'receipt.json')
        old_receipt = load(rd/'original_receipt.json')
        assert sha(rd/'impulse.h5') == rdef['raw_sha256'] == old_receipt['raw_sha256']
        assert sha(rd/'sfcw_20_170.h5') == rdef['wide_sfcw_sha256'] == old_receipt['sfcw_sha256']
        assert sha(rd/'original_solver.log') == old_receipt['log_sha256']
        assert sha(a.root/case['geometry']) == case['geometry_sha256']
        assert sha(a.root/rdef['case']['geometry']) == rdef['case']['geometry_sha256']
        for d, definition in [(nd, case), (rd, rdef['case'])]:
            with h5py.File(d/'impulse.h5') as h:
                y = h['rxs/rx1/Ez']; x = h['srcs/src1/excitation']
                assert y.dtype == np.dtype('float32') and y.shape == (definition['expected_samples'],)
                assert np.isfinite(y[:]).all() and np.max(np.abs(y[:])) > 0
                assert np.count_nonzero(x['samples'][:]) == 1 and np.max(x['samples'][:]) == 1
                assert np.isclose(y.attrs['SampleInterval'], definition['dt_s'], rtol=1e-12, atol=0)
                assert y.attrs['TimeSampleOffset'] == 0
                assert x.attrs['TimeSampleOffset'] == float(h.attrs['dt'])/2
                for group, key in [('srcs/src1', 'tx_m'), ('rxs/rx1', 'rx_m')]:
                    assert np.allclose(h[group].attrs['Position'], definition[key], rtol=0, atol=1e-8)
        log = (nd/'solver.log').read_text(encoding='utf-8')
        assert 'CUDA backend using Device 0' in log and 'Interface averaging is enabled' in log
        sf = read_sf(nd/'sfcw.h5', frequency_spec=(20e6, 110e6, 301))
        ref = read_sf(rd/'sfcw.h5', frequency_spec=(20e6, 110e6, 301))
        wide = read_sf(rd/'sfcw_20_170.h5')
        checks[cid] = independent_check(nd/'impulse.h5', sf)
        checks[rdef['id']+'_narrow'] = independent_check(rd/'impulse.h5', ref)
        checks[rdef['id']+'_wide'] = independent_check(rd/'impulse.h5', wide)
        shared_error = relative(ref['H'], wide['H'][:301])
        assert shared_error < 1e-10
        co = coefficients(a.root, cid)
        assert co['mixed_ez_edge_count'] > 0
        stats = compare(sf, ref, rdef['case'], c['domain_screen'])
        env_windows = []
        for lo, hi in c['domain_screen']['common_windows_ns']:
            mask = (ref['t'] >= lo) & (ref['t'] < hi)
            env_windows.append(dict(lo_ns=lo, hi_ns=hi,
                magnitude_relative_l2=relative(np.abs(sf['env'][mask]), np.abs(ref['env'][mask]))))
        rows.append(dict(case=cid, station_m=case['station_m'], **stats,
            envelope_magnitude_windows=env_windows, shared_narrow_wide_spectrum_l2=shared_error,
            tail_db=receipt['tail_db'], reference_tail_db=old_receipt['tail_db'],
            fdtd_seconds=receipt['official_fdtd_seconds'],
            build_solver_seconds=receipt['official_build_and_solver_seconds'],
            fdtd_speedup=old_receipt['official_fdtd_seconds']/receipt['official_fdtd_seconds'],
            build_solver_speedup=old_receipt['official_build_and_solver_seconds']/receipt['official_build_and_solver_seconds'],
            tx_coordinate_delta_m=case['tx_coordinate_delta_m'],
            rx_coordinate_delta_m=case['rx_coordinate_delta_m'],
            mixed_ez_edge_count=co['mixed_ez_edge_count'],
            min_center_sandstone_interval_m=min(v['thickness_m'] for v in case['center_column_intervals'] if v['material_id']==3)))
        refs[case['station_m']], news[case['station_m']], wides[case['station_m']] = ref, sf, wide
    summary = dict(new_fdtd=3, reused_fdtd=3,
        all_three_pass=all(r['passes_frozen_domain_screen'] for r in rows),
        max_spectrum_relative_l2=max(r['complex_band_relative_l2'] for r in rows),
        max_window_relative_l2=max(w['complex_relative_l2'] for r in rows for w in r['windows']),
        max_phase_rms_degrees=max(r['strong_bin_phase_rms_degrees'] for r in rows),
        median_build_solver_speedup=float(np.median([r['build_solver_speedup'] for r in rows])),
        total_fdtd_seconds=sum(r['fdtd_seconds'] for r in rows),
        total_build_solver_seconds=sum(r['build_solver_seconds'] for r in rows),
        batch_wall_seconds=complete['elapsed_since_first_stage_s'],
        narrow_tones=301, wide_tones=501, sample_spacing_ns=float(ref['t'][1]-ref['t'][0]),
        nominal_bandwidth_resolution_ratio=150/90,
        min_cells_per_wavelength_at110MHz=299792458/(110e6*np.sqrt(18)*.05),
        timing_scope='Single current run vs reused earlier GPU baselines; not repeated benchmark',
        current_model_modified=False, default_grid_modified=False,
        attribution='Practical replacement includes vertical TX/RX quantization up to2.5cm. Finite2.5cm baseline is not analytic truth. Tail remains unresolved.')
    save(a.out/'comparisons.json', dict(summary=summary, cases=rows, thresholds=c['domain_screen']))
    save(a.out/'audit.json', dict(status='PASS', receipt_files_verified=hashed, checks=checks,
        contract_sha256=sha(a.root/'contract.json'), analysis_script_sha256=sha(Path(__file__)),
        new_fdtd=3, reused_fdtd=3, input_manifest_verified=True))
    stations = sorted(refs)
    plot_sparse_compare(a.out, 'narrow_grid_sparse_bscan.png', stations,
        np.array([refs[s]['real'] for s in stations]).T,
        np.array([news[s]['real'] for s in stations]).T, ref['t'],
        '20–110 MHz／301点｜2.5cm旧基准 → 5cm候选｜平滑y／1100ns｜仅三站',
        reference_events={r['case']['station_m']:r['case']['vertical_reference_events'] for r in c['reuse_references']})
    fig, axes = plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for j,s in enumerate(stations):
        r,q=refs[s],news[s]
        for ax,key,label in [(axes[0,j],'real','带符号实响应'),(axes[1,j],'env','复包络模')]:
            for sf,name in [(r,'2.5cm旧正演'),(q,'5cm新正演')]:
                y=sf[key] if key=='real' else np.abs(sf[key])
                mask=(sf['t']>=100)&(sf['t']<=400)
                ax.plot(sf['t'][mask],y[mask],label=name,lw=1)
            ax.set(xlim=(100,400),xlabel='时间（ns）',ylabel=label,title=f'{s:g}m｜20–110MHz')
            ax.legend();ax.grid(alpha=.2)
    fig.suptitle('同频带对照｜没有增益、时移或逐道归一化；部分收发高度相差2.5cm')
    fig.savefig(a.out/'narrow_ascan_comparison.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for j,s in enumerate(stations):
        r,q=refs[s],news[s];ratio=q['H']/r['H']
        axes[0,j].plot(r['f']/1e6,20*np.log10(np.maximum(np.abs(ratio),1e-15)))
        axes[1,j].plot(r['f']/1e6,np.angle(ratio,deg=True))
        axes[0,j].set(title=f'{s:g}m｜5cm／2.5cm',ylabel='幅度比（dB）')
        axes[1,j].set(xlabel='频率（MHz）',ylabel='相位差（°，主值）')
    for ax in axes.flat:ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.2)
    fig.suptitle('官方direct复谱｜20–110MHz／301点；没有响应拟合或幅度校准')
    fig.savefig(a.out/'narrow_spectrum_comparison.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8),constrained_layout=True)
    for j,s in enumerate(stations):
        for sf,name in [(wides[s],'20–170MHz／501点'),(refs[s],'20–110MHz／301点')]:
            mask=(sf['t']>=100)&(sf['t']<=400)
            axes[0,j].plot(sf['t'][mask],sf['real'][mask],label=name,lw=1)
            axes[1,j].plot(sf['t'][mask],np.abs(sf['env'][mask]),label=name,lw=1)
        for ax in axes[:,j]:ax.set(xlim=(100,400),xlabel='时间（ns）',title=f'{s:g}m｜同一2.5cm原始数据');ax.legend();ax.grid(alpha=.2)
        axes[0,j].set_ylabel('官方带符号实响应');axes[1,j].set_ylabel('官方复包络模')
    fig.suptitle('仅改变处理频带｜各用官方Gaussian窗和1/N逆变换；不同频带峰值不能视为同一标定幅度')
    fig.savefig(a.out/'bandwidth_only_comparison.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4.3),constrained_layout=True)
    case=c['cases'][0];rd=c['reuse_references'][0]['case']
    cmap=ListedColormap(['#e5eef3','#d5ad69','#87786c','#c0c7ba']);norm=BoundaryNorm(np.arange(-.5,4.5),4)
    for ax,definition in zip(axes,[rd,case]):
        dl=definition['grid_m'];ox=definition['origin_global_x_elevation_m'][0]
        with h5py.File(a.root/definition['geometry']) as h:data=h['data'][:,:,0]
        ix0,ix1=np.floor((np.array([148.,152.])-ox)/dl).astype(int)
        iy0,iy1=np.floor((np.array([433.4,434.3])-410)/dl).astype(int)
        im=ax.imshow(data[ix0:ix1,iy0:iy1].T,origin='lower',aspect='auto',cmap=cmap,norm=norm,interpolation='nearest',
            extent=[ox+ix0*dl,ox+ix1*dl,410+iy0*dl,410+iy1*dl])
        ax.set(title=f'{dl*100:g}cm｜原材料体素',xlabel='剖面里程（m）',ylabel='高程（m）')
    cb=fig.colorbar(im,ax=axes,ticks=[0,1,2,3],shrink=.8);cb.ax.set_yticklabels(['空气','粉质粘土','泥岩','砂岩'])
    fig.suptitle('150m薄层接近区｜同一连续地层重离散；平滑不能恢复消失的几何细节')
    fig.savefig(a.out/'thin_segment_geometry.png',dpi=180);plt.close(fig)
    print(summary)


if __name__ == '__main__':
    main()
