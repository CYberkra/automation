"""CPU free-space line/dipole comparison; not a3D geological run or instrument fit."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.special import hankel2

from analyze_line9_v401_version_controls import FREQ, inverse
from hs_capsule_identity import sha256 as sha

C0 = 299792458.
EPS0 = 8.8541878128e-12
MU0 = 1 / (EPS0 * C0**2)
REFERENCE_M = .025
GATES = dict(early=[0., 20.], wide=[300., 450.],
             inherited_basal=[345.6180971390552, 374.6081170991351], late=[450., 1100.])


def point_dipole(frequency, radius_xy, z):
    """Ez/(I dl), exp(+iwt), z-oriented dipole at a general complex z offset.

    Radiation,induction and electrostatic terms are all retained. The complex
    extension is used only for a convergent line-source integral,not a geometry.
    """
    omega = 2 * np.pi * frequency
    r = np.sqrt(radius_xy**2 + z * z + 0j)
    nz2 = z * z / (r * r)
    return np.exp(-1j * omega * r / C0) / (4 * np.pi * EPS0) * (
        (3 * nz2 - 1) * (1 / (1j * omega * r**3) + 1 / (C0 * r**2))
        + (nz2 - 1) * 1j * omega / (C0**2 * r))


def integrated_line(frequency, radius_xy, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    rotation = np.exp(-1j * np.pi / 4)
    maximum = 40 / (2 * np.pi * frequency / C0 * np.sin(np.pi / 4))
    u = (nodes[None, :] + 1) * maximum[:, None] / 2
    z = rotation * u
    return 2 * rotation * np.sum(point_dipole(frequency[:, None], radius_xy, z)
                                 * weights[None, :] * maximum[:, None] / 2, axis=1)


def summarize(profiles, time_ns):
    early = (time_ns >= 0) & (time_ns <= 20)
    norm = np.linalg.norm
    norms = norm(profiles[early], axis=0)
    assert np.all(norms > 0)
    result = {}
    for name, bounds in GATES.items():
        keep = (time_ns >= bounds[0]) & (time_ns <= bounds[1])
        result[name] = dict(gate_ns=bounds, norm_over_own_early= (norm(profiles[keep], axis=0) / norms).tolist(),
            peak_ns=[float(time_ns[keep][np.argmax(abs(profiles[keep, j]))]) for j in range(profiles.shape[1])])
    return result, norms


def main(a):
    assert not a.out.exists() and not a.numerical.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    p = read(a.public / 'analysis.json')
    audit = read(a.public / 'independent_audit.json')
    assert audit['analysis_sha256'] == sha(a.public / 'analysis.json') and audit['status'].startswith('PASS')
    assert p['numerical_sha256'] == sha(a.response)
    native_records = {r['id']: r['native_sha256'] for r in p['native']}
    native_ids = ['high_x19000_H0', 'low_x19000_H0']
    positions = []
    for sid in native_ids:
        path = a.source / (sid + '.h5')
        assert sha(path) == native_records[sid]
        with h5py.File(path) as h:
            assert h.attrs['gprMax'] == '4.0.1' and h['rxs/rx1/Ez'].dtype == np.float64
            positions.append((np.array(h['srcs/src1'].attrs['Position']), np.array(h['rxs/rx1'].attrs['Position'])))
    for tx, rx in positions:
        np.testing.assert_array_equal(tx, positions[0][0])
        np.testing.assert_array_equal(rx, positions[0][1])
        assert tx[2] == rx[2] == 0
    radius = float(np.linalg.norm(positions[0][1] - positions[0][0]))
    k = 2 * np.pi * FREQ / C0
    line = -2 * np.pi * FREQ * MU0 / 4 * hankel2(0, k * radius)
    orders = [256, 512]
    integrals = [integrated_line(FREQ, radius, order) for order in orders]
    convergence = [float(np.linalg.norm(q - line) / np.linalg.norm(line)) for q in integrals]
    assert max(convergence) < 1e-8, convergence
    dipole = point_dipole(FREQ, radius, np.zeros(501))
    with h5py.File(a.response) as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:], FREQ)
        ids = h['ids'].asstr()[:].tolist()
        native = h['response'][:, [ids.index(sid) for sid in native_ids]]
    # Only the2D line is divided by the established .025m bookkeeping reference.
    # Physical source distributions differ; never treat their absolute ratio as
    # a3D/2D source-power,antenna gain or sim/field calibration.
    z = np.column_stack([native, line / REFERENCE_M, dipole])
    names = ['native_high_H0', 'native_low_H0', 'analytic_2D_line_reference', 'analytic_3D_point_dipole']
    labels = ['原高损耗H0（仍有真实覆盖层）', '低损耗诊断H0', '解析二维自由空间线源', '解析三维自由空间短偶极子']
    a.out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    metrics = {}
    early_errors = {}
    for window, w in [('hann', np.hanning(501)), ('blackman', np.blackman(501))]:
        w = w / w.mean()
        profiles, t = inverse(z, w)
        tn = t * 1e9
        metrics[window], norms = summarize(profiles, tn)
        early = (tn >= 0) & (tn <= 20)
        early_errors[window] = [float(np.linalg.norm(profiles[early, j] - profiles[early, 2]) / np.linalg.norm(profiles[early, j])) for j in [0, 1]]
        normalized = profiles / norms[None, :]
        fig, axes = plt.subplots(2, 2, figsize=(15, 10), layout='constrained')
        keep = (tn >= 0) & (tn <= 80)
        for j in [0, 1, 2]:
            axes[0, 0].plot(tn[keep], profiles[keep, j].real, lw=1, label=labels[j])
        axes[0, 0].set(title='原H0与二维自由空间解析：共同物理尺度，未拟合', ylabel='(V/m)/(A·m)；二维参考', xlabel='SFCW时间 / ns')
        for j in range(4):
            axes[0, 1].plot(tn[keep], normalized[keep, j].real, lw=1, label=labels[j])
        axes[0, 1].set(title='只比较形状：各自0–20ns复数L2归一', xlabel='SFCW时间 / ns', ylabel='各自早段归一复场的实部（不同尺度）')
        keep = (tn >= 300) & (tn <= 450)
        for j in range(4):
            axes[1, 0].semilogy(tn[keep], abs(normalized[keep, j]), lw=1, label=labels[j])
        axes[1, 0].axvspan(*GATES['inherited_basal'], color='green', alpha=.1)
        axes[1, 0].set(title='晚时窗：H0含地层，解析两项只有自由空间直达', xlabel='SFCW时间 / ns', ylabel='各自早段L2归一复包络幅值')
        for j in range(4):
            keys = ['wide', 'inherited_basal', 'late']
            axes[1, 1].plot(range(3), [metrics[window][key]['norm_over_own_early'][j] for key in keys], marker='o', label=labels[j])
        axes[1, 1].set_yscale('log')
        axes[1, 1].set_xticks(range(3), ['300–450ns', '原底砂到时窗', '450–1100ns'])
        axes[1, 1].set(title='各窗复数范数/各自早段范数；非能量占比', ylabel='L2范数比（各自尺度）')
        for ax in axes.flat:
            ax.legend(fontsize=8)
        fig.suptitle(f'190m收发间距{radius:.6f}m / {window} / 精确20–170MHz、501点 / CPU解析核查，无新FDTD\n三维只有理想自由空间偶极，未含地表、地层或真实天线；归一图仅比较形状，不比较功率、不改生产数据。')
        fig.savefig(a.out / f'direct_dimensionality_{window}.png', dpi=130)
        plt.close(fig)
    with h5py.File(a.numerical, 'x') as h:
        h['frequency_Hz'] = FREQ
        h['response'] = z
        h['integrated_3D_line_response_E_over_I'] = np.column_stack(integrals)
        h['ids'] = np.array(names, dtype=h5py.string_dtype('utf-8'))
    report = dict(status='CPU_FREE_SPACE_DIMENSIONALITY_DIAGNOSTIC_NOT_3D_GEOLOGY_OR_ANTENNA_VALIDATION',
        script_sha256=sha(__file__), analysis_sha256=sha(a.public / 'analysis.json'),
        independent_source_audit_sha256=sha(a.public / 'independent_audit.json'), response_sha256=sha(a.response),
        numerical_sha256=sha(a.numerical), native={sid: native_records[sid] for sid in native_ids},
        source_tx_m=positions[0][0].tolist(), receiver_m=positions[0][1].tolist(),
        air_distance_m=radius, air_travel_time_ns=radius / C0 * 1e9, response_columns=names,
        line_integral_orders=orders, line_integral_contour_angle_rad=-np.pi / 4,
        line_integral_tail_exponent=40., integrated_3D_vs_Hankel_relative_L2=convergence,
        own_early_gate_ns=GATES['early'], metrics=metrics,
        native_H0_vs_2D_direct_unfitted_early_relative_L2=early_errors,
        calls_solver=False, calls_CUDA=False, production_data_changed=False,
        primary_sources=['https://docs.gprmax.com/en/latest/input_hash_cmds.html#hertzian-dipole',
                         'https://docs.gprmax.com/en/latest/comparisons_analytical.html#hertzian-dipole-in-free-space',
                         'https://raw.githubusercontent.com/gprMax/gprMax/master/testing/analytical_solutions.py'],
        reading='Official source semantics and analytical-comparison Hertzian sections; official analytical_solutions.py Hertzian function,not other validation cases. General3D frequency formula derived from its complete time-domain radiation/induction/electrostatic form; line integral cross-check independent of native solver.',
        limits='Infinite2D line versus infinitesimal3D point source distributions differ. Absolute2D/3D ratio has no instrument meaning. Shape normalization is CPU diagnostic only,not amplitude fitting,production normalization or field calibration.3D has no ground,layers,antenna/cables/ports. Native H0 includes real cover reflections. Finite-band tails and complex interference remain; mathematical checks are not3D FDTD convergence or unique late-event attribution.')
    (a.out / 'analysis.json').write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=report['status'], line_integral_errors=convergence, metrics=metrics)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source', 'public', 'response', 'out', 'numerical']:
        p.add_argument('--' + key, type=Path, required=True)
    main(p.parse_args())
