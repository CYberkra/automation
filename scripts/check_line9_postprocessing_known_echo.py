"""Known inclined echo with strong early arrivals: array-only processing counterexample."""
import argparse
from pathlib import Path

import numpy as np

from audit_line9_postprocessing import FREQ, inverse, weights, rel, save, sha


def main(out):
    if out.exists():
        raise ValueError('Use fresh output directory')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    positions = np.linspace(0, 13.8, 70)
    delays = np.linspace(290e-9, 258e-9, 70)
    target = 1e-4 * np.exp(-2j * np.pi * FREQ[:, None] * delays[None, :])
    early = (np.exp(-2j * np.pi * FREQ * 4.3e-9) +
             .3 * np.exp(-2j * np.pi * FREQ * 8e-9))[:, None]
    combined = early + target
    findings = {}
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    for column, name in enumerate(['gaussian', 'hann']):
        w = weights(name, 501)
        isolated, time = inverse(target, FREQ, w)
        total, _ = inverse(combined, FREQ, w)
        background, _ = inverse(np.broadcast_to(early, target.shape), FREQ, w)
        assert rel(total - background, isolated) < 1e-10
        error_ns = np.max(abs(time[np.argmax(abs(isolated), axis=0)] - delays)) * 1e9
        assert error_ns <= (time[1] * 1e9) / 2
        near = abs(time[:, None] - delays[None, :]) <= 10e-9
        findings[name] = dict(isolated_target_max_delay_error_ns=float(error_ns),
                              exact_complex_difference_relative_L2=rel(total - background, isolated),
                              early_leakage_over_target_near_true_interface_complex_L2=float(np.linalg.norm(background[near]) / np.linalg.norm(isolated[near])))
        for row, value, title in [(0, isolated, '仅已知倾斜弱界面'), (1, total, '叠加直达及地表强回波')]:
            keep = time <= 600e-9
            # Same reference1 across all panels, no renormalisation of the weak echo.
            db = 20 * np.log10(np.maximum(abs(value[keep]), 1e-12))
            ax = axes[row, column]
            pic = ax.imshow(db, origin='upper', aspect='auto', extent=[0, 13.8, time[keep][-1] * 1e9, 0],
                            cmap='gray_r', vmin=-110, vmax=-50, interpolation='nearest')
            ax.plot(positions, delays * 1e9, '--', color='#1ec766', lw=1, label='已知真值258–290ns')
            ax.legend(fontsize=8, loc='lower right')
            ax.set(xlabel='合成横距 / m', ylabel='时间 / ns', title=name + '：' + title)
            fig.colorbar(pic, ax=ax, label='幅度dB / 固定参考1')
    fig.suptitle('纯数组反例：同一倾斜界面，幅度1e−4；早时幅度1及0.3\n20–170MHz/0.3MHz/501点、无AGC/叠道；非FDTD结果，也非实际三包的目标真值')
    fig.savefig(out / 'known_inclined_echo_counterexample.png', dpi=135)
    plt.close(fig)
    save(out / 'known_echo_check.json', dict(calls_solver=False, calls_training=False,
        script_sha256=sha(__file__), reconstruction_script_sha256=sha(Path(__file__).with_name('audit_line9_postprocessing.py')),
        declared_echo_amplitudes=[1., .3, 1e-4], early_delays_ns=[4.3, 8.],
        target_delay_limits_ns=[258., 290.], results=findings,
        limits='Array-only assumed echoes. Demonstrates masking by finite-band reconstruction, not actual physical target detection or a recovered gprMax waveform.'))
    print(findings)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    main(args.out)
