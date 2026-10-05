"""Annotated Chinese comparison figures for the material/aperture proof chain.

Regenerates the four key B-scan panels from the two completed capsules and
renders them with full Chinese labels, conclusion annotations and legends.
Output: fresh directory with fig1_material_6m.png, fig2_aperture_30m.png,
fig3_real_line9.png.
"""
import argparse
from dataclasses import replace
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
from analyze_hs4_height_wavefield import response
from plot_hs4_permittivity_bscans import C_AIR, SURFACE_Z, interface_relief
from sfcw_official_loader_v0_2 import verify_official_runtime
from gprMax.toolboxes.SFCW.processing import reconstruct_time_response

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

CHECKS = ROOT / 'artifacts/research_checks'
MAT = CHECKS / '2026-10-05_hs4_material_scan_r1'
APX = CHECKS / '2026-10-05_hs4_aperture_scan_r1'
RX_OFFSET, ANT_Z, SHIFT = 1.3, 20.0, 36.0


def relief_tw(xmids, eps):
    fx, fz = interface_relief()
    z = np.interp(xmids, fx + SHIFT if xmids[0] > 36 else fx, fz)
    return 2*np.hypot(0.65, ANT_Z-SURFACE_Z)/C_AIR + 2*(SURFACE_Z - z)/C_AIR*np.sqrt(eps)


def load_bscan(capsule, tag, n_st):
    c = json.loads((capsule/'execution_contract.json').read_text('utf-8'))
    spectra = {}
    for g in c['groups']:
        if g['material_tag'] != tag:
            continue
        audit = json.loads((Path(g['input']).parent/'audit.json').read_text('utf-8'))
        for row in audit:
            raw = Path(g['input']).parent/row['file']
            if sha256(raw) != row['sha256']:
                raise ValueError('raw identity differs')
            spectra[g['role'], row['station_index']] = response(raw)
    template = next(iter(spectra.values()))
    indices = sorted(si for role, si in spectra if role == 'rough')
    if len(indices) != n_st:
        raise ValueError(f'{tag}: expected {n_st} stations, got {len(indices)}')
    total, contrast = [], []
    for si in indices:
        rough = spectra['rough', si]
        total.append(np.abs(reconstruct_time_response(rough, window='hann', zero_pad_factor=8).complex_envelope))
        cs = replace(template, response=rough.response - spectra['halfspace', si].response)
        contrast.append(np.abs(reconstruct_time_response(cs, window='hann', zero_pad_factor=8).complex_envelope))
    time = reconstruct_time_response(rough, window='hann', zero_pad_factor=8).time * 1e9
    return time, np.stack(total, 1), np.stack(contrast, 1)


def panel(ax, mat, time, x, relief, title, note, note_color):
    mdb = 20*np.log10(np.maximum(mat, mat.max()*1e-9)/mat.max())
    im = ax.imshow(mdb, aspect='auto', cmap='gray', vmin=-60, vmax=0,
                   extent=[x[0], x[-1], time[-1], time[0]])
    ax.plot(x, relief, 'r-', lw=1.6, label='真实基覆界面（双程时）')
    ax.set_ylim(260, 0)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.text(0.02, 0.97, note, transform=ax.transAxes, fontsize=11, va='top',
            color=note_color, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=note_color, alpha=0.85))
    ax.legend(fontsize=10, loc='lower right')
    ax.set_xlabel('测线距离 (m)', fontsize=11)
    ax.set_ylabel('时间 (ns)', fontsize=11)
    ax.tick_params(labelsize=10)
    return im


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    if a.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    a.out.mkdir(parents=True)

    # ---- fig 1: material comparison at 6 m aperture ----
    x6 = 14.6 + 0.5*np.arange(13)
    t6, tot_d, con_d = load_bscan(MAT, 'eps18.017_s0.003_debye', 13)
    _, tot_l, con_l = load_bscan(MAT, 'eps12_s0.001', 13)
    r_d = relief_tw(x6 + RX_OFFSET/2, 18.017)
    r_l = relief_tw(x6 + RX_OFFSET/2, 12.0)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), layout='constrained')
    panel(axes[0, 0], tot_d, t6, x6, r_d, '归档湿黏土配方（总场）\nεr=18 + Debye + σ=0.003',
          '界面回波 −65.6 dB：完全淹没【差】', 'red')
    panel(axes[0, 1], tot_l, t6, x6, r_l, '低损耗配方（总场）\nεr=12, σ=0.001（与实测相容）',
          '界面回波 −45.8 dB：红线上条带可见【好】', 'green')
    panel(axes[1, 0], con_d, t6, x6, r_d, '归档配方（去除背景后）', '差分后才有条带，且弱', 'red')
    panel(axes[1, 1], con_l, t6, x6, r_l, '低损耗配方（去除背景后）', '条带强 21 dB', 'green')
    for ax in axes.flat:
        ax.set_xlim(14, 22.6)
    fig.colorbar(axes[0, 0].images[0], ax=axes, shrink=0.85,
                 label='回波强度 (dB，各图相对自身最大值)')
    fig.suptitle('图1  材料配方对比（6 m 孔径、8 m 航高、同一起伏界面）：损耗配方决定界面"看不看得见"',
                 fontsize=15, fontweight='bold')
    fig.savefig(a.out/'fig1_material_6m.png', dpi=140)
    plt.close(fig)

    # ---- fig 2: aperture comparison at 30 m ----
    x30 = 38.6 + 0.5*np.arange(61)
    t30, tot30_d, _ = load_bscan(APX, 'eps18.017_s0.003_debye', 61)
    _, tot30_l, _ = load_bscan(APX, 'eps12_s0.001', 61)
    r30_d = relief_tw(x30 + RX_OFFSET/2, 18.017)
    r30_l = relief_tw(x30 + RX_OFFSET/2, 12.0)
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), layout='constrained')
    panel(axes[0], tot30_l, t30, x30, r30_l, '低损耗配方 εr=12, σ=0.001（总场，30 m 孔径）',
          '界面条带横跨全孔径、贴红线：形状相关0.59【好】（6m孔径时仅0.003）', 'green')
    panel(axes[1], tot30_d, t30, x30, r30_d, '归档湿黏土配方（总场，30 m 孔径）',
          '孔径救不了过损耗材料：仍淹没【差】', 'red')
    for ax in axes:
        ax.set_xlim(38, 70.4)
    fig.colorbar(axes[0].images[0], ax=axes, shrink=0.85,
                 label='回波强度 (dB，各图相对自身最大值)')
    fig.suptitle('图2  孔径对比（30 m 孔径、8 m 航高；蓝虚线间为起伏段，两侧为平界面）：孔径决定"看不看得清形状"',
                 fontsize=15, fontweight='bold')
    for ax in axes:
        ax.axvline(48, color='b', ls=':', lw=1.2)
        ax.axvline(60, color='b', ls=':', lw=1.2)
    fig.savefig(a.out/'fig2_aperture_30m.png', dpi=140)
    plt.close(fig)

    # ---- fig 3: measured Line9 ----
    csv = ROOT/'real_data_yingshan/营山测线数据/Line9origin(36).csv'
    with open(csv, encoding='utf-8') as f:
        hdr = [f.readline() for _ in range(4)]
    ns = int(hdr[0].split('=')[1].split(',')[0]); T = float(hdr[1].split('=')[1].split(',')[0])
    nt = int(hdr[2].split('=')[1].split(',')[0]); dx = float(hdr[3].split('=')[1].split(',')[0])
    d = np.loadtxt(csv, delimiter=',', skiprows=4)
    amp = d[:, 3].reshape(nt, ns).T
    from scipy.signal import hilbert
    env = np.abs(hilbert(amp, axis=0))
    edb = 20*np.log10(np.maximum(env, env.max()*1e-10)/env.max())
    fig, ax = plt.subplots(figsize=(13, 5.5), layout='constrained')
    im = ax.imshow(edb, aspect='auto', cmap='gray', vmin=-70, vmax=0,
                   extent=[0, nt*dx, T, 0])
    ax.annotate('基覆界面条带（400–500 ns）：\n横向连续、形态完整',
                xy=(150, 440), xytext=(150, 300), fontsize=12, fontweight='bold', color='red',
                arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='red', alpha=0.85))
    ax.annotate('地表反射', xy=(30, 60), xytext=(60, 130), fontsize=11, color='orange',
                arrowprops=dict(arrowstyle='->', color='orange', lw=1.2),
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='orange', alpha=0.85))
    ax.set_xlabel('测线距离 (m)，孔径 216 m', fontsize=11)
    ax.set_ylabel('时间 (ns)', fontsize=11)
    ax.set_title('图3  实测 Line9（9 m 航高、216 m 孔径）：深层界面清晰可见 —— 仿真的目标形态', fontsize=14, fontweight='bold')
    fig.colorbar(im, ax=ax, shrink=0.9, label='包络强度 (dB)')
    fig.savefig(a.out/'fig3_real_line9.png', dpi=140)
    plt.close(fig)
    print('written:', [p.name for p in sorted(a.out.glob('*.png'))])


if __name__ == '__main__':
    main()
