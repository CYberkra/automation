"""Visualization report for the 5 cm free-space dipole batch on the RTX4090 machine.

Reads the accepted analysis capsule and the native H5s, renders a Chinese-
annotated figure: left = FDTD vs analytic continuum dipole spectrum (magnitude
and phase) for one representative receiver; right = per-receiver complex
relative L2 and max phase error bars plus reciprocity/symmetry numbers.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

from hs_capsule_identity import sha256
from analyze_uav_local_free_space import dipole_field, transfer, FREQ

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

STUDY = ROOT / 'artifacts/research_checks/2026-10-05_local_uav_free_space_fine_4090_r3'
ANALYSIS = ROOT / 'artifacts/research_checks/2026-10-05_local_uav_free_space_fine_4090_r3_analysis'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    contract = json.loads((STUDY/'execution_contract.json').read_text('utf-8'))
    summary = json.loads((ANALYSIS/'summary.json').read_text('utf-8'))
    audit = json.loads((STUDY/'completed_verification.json').read_text('utf-8'))
    identities = {r['id']: r['raw_sha256'] for r in audit['groups']}

    groups = {g['id']: g for g in contract['groups']}
    g = groups['d0.05_x_forward']
    raw = Path(g['input']).with_suffix('.h5')
    if sha256(raw) != identities[g['id']]:
        raise ValueError('raw changed')
    value, _, _ = transfer(raw, 1, 'Ex', g['spacing_m'])
    exact = dipole_field(np.array(g['receivers_m'][0])-np.array(g['tx_m']), 'x')[:, 0]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), layout='constrained')
    ax = axes[0]
    ax.plot(FREQ/1e6, 20*np.log10(np.abs(value)), label='FDTD 求解（V4.0.0/CUDA float64）')
    ax.plot(FREQ/1e6, 20*np.log10(np.abs(exact)), '--', label='解析连续偶极（Hertzian 精确解）')
    ax.set_xlabel('频率 (MHz)'); ax.set_ylabel('|E/(I·dl)| (dB)')
    ax.set_title('空气域 5 cm 网格：x 极化源→+x 1.3 m 接收，幅值谱对比\n（源为实际时钟 Ricker 反卷积，无经验相位/幅度拟合）', fontsize=10)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    axp = ax.twinx()
    axp.plot(FREQ/1e6, np.angle(value/exact)*180/np.pi, 'r-', lw=0.8, alpha=0.7, label='相位差')
    axp.set_ylabel('相位差 (°)', color='r'); axp.tick_params(axis='y', colors='r')
    axp.set_ylim(-1, 1)

    ax2 = axes[1]
    labels = [f"{m['case'].replace('d0.05_','')}\nrx{m['receiver']}" for m in summary['metrics']]
    x = np.arange(len(labels))
    l2 = [100*m['complex_relative_L2_to_continuum'] for m in summary['metrics']]
    ph = [m['maximum_phase_error_deg'] for m in summary['metrics']]
    b = ax2.bar(x-0.2, l2, width=0.4, label='复谱相对 L2 误差 (%)')
    ax2.axhline(5, color='gray', ls=':', lw=1, label='声明诊断门 5%')
    for xi, v in zip(x-0.2, l2):
        ax2.text(xi, v+0.01, f'{v:.2f}', ha='center', fontsize=7)
    ax2.set_xticks(x); ax2.set_xticklabels(labels, fontsize=6, rotation=45, ha='right')
    ax2.set_ylabel('复谱相对 L2 误差 (%)'); ax2.set_ylim(0, 0.6)
    ax2b = ax2.twinx()
    ax2b.bar(x+0.2, ph, width=0.4, color='r', alpha=0.6, label='最大相位误差 (°)')
    for xi, v in zip(x+0.2, ph):
        ax2b.text(xi, v+0.005, f'{v:.3f}', ha='center', fontsize=7, color='r')
    ax2b.set_ylabel('最大相位误差 (°)', color='r'); ax2b.tick_params(axis='y', colors='r')
    ax2b.set_ylim(0, 5)
    rec = summary['symmetry_reciprocity_relative_L2']
    ax2.set_title('10 个接收点全部过诊断门（5%/5°）\n'
                  f"互易/对称差：轴向 {rec['d0.05_axial']:.1e}、侧向 {rec['d0.05_broadside']:.1e}、"
                  f"reverse_x {rec['reverse_x_reciprocity']:.1e}、reverse_y {rec['reverse_y_reciprocity']:.1e}", fontsize=9)
    ax2.legend(loc='upper left', fontsize=8); ax2b.legend(loc='upper right', fontsize=8)
    ax2.grid(alpha=0.3)
    fig.suptitle('UAV 空气基准·5 cm 细网格批（RTX4090 本机，4 道，对比 3060 机 10 cm：L2 0.64–1.11%→本批 0.16–0.27%，相位 0.564°→0.137°）',
                 fontsize=11)
    a.out.mkdir(parents=True)
    fig.savefig(a.out/'free_space_fine_4090_report.png', dpi=140)
    plt.close(fig)
    print(str(a.out/'free_space_fine_4090_report.png'))


if __name__ == '__main__':
    main()
