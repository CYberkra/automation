"""CPU material/return balance using frozen factorial JSON, no solver or fitting."""
import argparse
import json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_planar_point_layer_green import layer_parts, point_integral
from diagnose_line9_v5_planar_green import integral, MU0
from review_line9_result_packages import indices, C0
from analyze_line9_v401_version_controls import FREQ

ROOT = Path(__file__).resolve().parents[1]
IDS = ['span3_dc003', 'span3_dc0003', 'span03_dc003', 'span03_dc0003']
LABELS = ['原高跨度/高DC', '只降DC', '只降极化跨度', '两项均降']


def read(path):
    return json.loads(path.read_text('utf-8'))


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def main(a):
    assert not a.out.exists() and not a.arrays.exists()
    oldpub = ROOT/'artifacts/research_checks/2026-10-08_v401_polarization_dc_r1'
    dcpublic = ROOT/'artifacts/research_checks/2026-10-08_v401_nonflat_material_r1'
    old = read(oldpub/'analysis.json')
    pm = read(a.factorial/'manifest.json'); dm = read(a.dc/'manifest.json')
    assert pm == read(oldpub/'preparation.json') and dm == read(dcpublic/'preparation.json')
    assert read(oldpub/'independent_audit.json')['status'].startswith('PASS')
    groups = {g['id']: g for g in pm['groups']}
    dcgroups = {g['id']: g for g in dm['groups']}
    records = []; spectra = []
    for label in IDS:
        g = dcgroups['baseline_H1' if label=='span3_dc003' else 'low_H1'] if label.startswith('span3_') else groups[label+'_H1']
        folder = a.dc if label.startswith('span3_') else a.factorial
        path = folder/g['material']; assert sha(path)==g['material_sha256']
        db = read(path)['materials']
        records.append({'id': label, 'material_sha256': sha(path), 'materials': db})
        spectra.append(np.array([indices(db, f)**2 for f in FREQ]))
    for rec in records[1:]:
        for key in ['material_000_air', 'material_001_cover', 'material_003_sandstone']:
            assert rec['materials'][key] == records[0]['materials'][key]
    geo = pm['geometry_diagnostic']; tx = pm['groups'][0]['tx_m']; rx = pm['groups'][0]['rx_m']
    height = tx[1]+rx[1]-2*geo['surface_y_m']; offset = abs(tx[0]-rx[0])
    thickness = geo['cover_base']['depth_m']
    contract = {'status': 'FROZEN_CPU_MATERIAL_RETURN_BALANCE_NOT_FDTD', 'script_sha256': sha(__file__),
                'source_analysis_sha256': sha(oldpub/'analysis.json'), 'materials': records,
                'source_preparation_sha256': [sha(oldpub/'preparation.json'), sha(dcpublic/'preparation.json')],
                'helper_sha256': {n: sha(Path(__file__).with_name(n)) for n in
                    ['line9_planar_point_layer_green.py', 'diagnose_line9_v5_planar_green.py', 'review_line9_result_packages.py']},
                'air_height_sum_m': height, 'offset_m': offset, 'cover_thickness_m': thickness,
                'frequency_Hz': FREQ.tolist(), 'quadrature_orders': [256, 512], 'tail_exponent': 36,
                'windows': ['hann', 'blackman'], 'comparison': 'Within each source dimension, weighted fullband component norms relative to its original material',
                'limits': 'Planar material sensitivity only, no nonflat ray/bounce identity or calibrated site loss; actual H0 is geological response, not pure clutter'}
    a.out.mkdir(parents=True); save(a.out/'contract.json', contract)
    responses = []; checks = []; interface = []
    k = 2*np.pi*FREQ/C0
    for label, e in zip(IDS, spectra):
        kernel = lambda q, ky: layer_parts(q, ky, k, e[:, 1], e[:, 2], thickness)
        vals = []
        for order in [256, 512]:
            line = -FREQ[:, None]*MU0*integral(FREQ, height, offset, lambda q, ky: kernel(q, ky)[0], order=order)
            point = point_integral(FREQ, height, offset, np.pi/2, kernel, order=order)
            vals.append(np.stack([line, point], axis=1))
        er = np.linalg.norm(vals[1]-vals[0], axis=0)/np.linalg.norm(vals[1], axis=0)
        assert np.max(er)<1e-6; checks.append(er.tolist()); responses.append(vals[1])
        n1 = np.sqrt(e[:, 1]); n2 = np.sqrt(e[:, 2]); r = (n1-n2)/(n1+n2)
        interface.append({'id': label, 'epsilon_mud_real_20_95_170': e[[0,250,500],2].real.tolist(),
                          'epsilon_mud_loss_20_95_170': (-e[[0,250,500],2].imag).tolist(),
                          'normal_cover_mud_reflection_abs_20_95_170': abs(r[[0,250,500]]).tolist()})
    response = np.stack(responses, axis=1)
    metrics = {}; balance = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        norms = np.linalg.norm(response*w[:,None,None,None], axis=0)
        relative = norms/norms[0]; metrics[window] = {'component_norm_by_case_dimension': norms.tolist(),
                                                     'relative_to_original_by_case_dimension': relative.tolist()}
        balance[window] = {}
        for gate in ['basal','deep']:
            cases = old['metrics'][window][gate]['cases']; base = cases[IDS[0]]; rows = []
            for label in IDS:
                cur = cases[label]; hg = cur['H0_norm']/base['H0_norm']; dg = cur['delta_norm']/base['delta_norm']
                improvement = base['H0_over_delta']/cur['H0_over_delta']
                assert abs(improvement-dg/hg)<1e-10*improvement
                rows.append({'id': label, 'H0_remaining_amplitude_norm_fraction': hg,
                             'bottom_delta_amplitude_norm_gain': dg, 'H0_over_delta_ratio_improvement': improvement,
                             'deep_gain_dB': float(20*np.log10(dg)), 'H0_reduction_dB': float(-20*np.log10(hg)),
                             'ratio_improvement_dB': float(20*np.log10(improvement))})
            balance[window][gate] = rows
    a.arrays.parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.arrays, response=response, frequency_Hz=FREQ, epsilon=np.stack(spectra,axis=1))
    result = {'status': 'CPU_MATERIAL_RETURN_BALANCE_NO_NEW_SOLVER', 'script_sha256': sha(__file__),
              'contract_sha256': sha(a.out/'contract.json'), 'arrays_sha256': sha(a.arrays),
              'quadrature256_512_relative_L2': checks, 'normal_interfaces': interface,
              'planar_component_relative_metrics': metrics, 'actual_factorial_norm_balance': balance,
              'new_gprmax_runs': 0, 'site_parameters_changed': False, 'limits': contract['limits']}
    save(a.out/'analysis.json', result)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']; plt.rcParams['axes.unicode_minus']=False
    fig, axs = plt.subplots(1,3, figsize=(17,5), layout='constrained'); x=np.arange(4)
    for dim, label in enumerate(['二维TE线源','三维横向点偶极']):
        v=np.array(metrics['hann']['relative_to_original_by_case_dimension'])[:,dim]
        axs[0].plot(x,v[:,1],'-o',label=label+'：一次项'); axs[0].plot(x,v[:,2],'--o',label=label+'：两次项')
    axs[0].set(ylabel='全带场幅范数 / 同维度原配方',title='无限平层：界面返回也会随电性变化');axs[0].legend(fontsize=8)
    b=balance['hann']['basal']
    axs[1].plot(x,[r['H0_remaining_amplitude_norm_fraction'] for r in b],'-o',label='原起伏H0地质响应')
    axs[1].plot(x,[r['bottom_delta_amplitude_norm_gain'] for r in b],'-o',label='原起伏底砂替换差场');axs[1].set(yscale='log',ylabel='场幅范数 / 原配方',title='已完成原生对照：底砂窗，未重新求解');axs[1].legend(fontsize=8)
    axs[2].bar(x-.15,[r['deep_gain_dB'] for r in b],.3,label='底砂差场增强');axs[2].bar(x+.15,[r['H0_reduction_dB'] for r in b],.3,label='H0地质响应减弱')
    axs[2].set(ylabel='20log10 场幅范数变化 / dB',title='相对改善的两个因素；不是能量占比');axs[2].legend(fontsize=8)
    for ax in axs: ax.set_xticks(x,LABELS,rotation=15);ax.grid(alpha=.2)
    fig.suptitle('低损耗图更清楚：浅层返回消失了，还是深部回波强了？\n精确20–170MHz / 501点 / Hann；CPU平层参考与既有原生对照分开，无AGC/拟合/正式材料替换')
    fig.savefig(a.out/'material_return_balance.png',dpi=140);plt.close(fig)
    print(json.dumps({'status':result['status'],'quadrature_max_error':float(np.max(checks)),
                      'combined_change_hann_basal':b[-1]},ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['factorial','dc','out','arrays']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
