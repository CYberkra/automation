"""Render audited original/flat native Ricker fields, retaining all timing conventions."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from PIL import Image
    assert not a.out.exists()
    read = lambda p: json.loads(p.read_text('utf-8'))
    c = read(a.execution / 'execution_contract.json')
    v = read(a.execution / 'completed_verification.json')
    m = c['study_manifest']
    assert v['completed'] and v['contract_sha256'] == sha(a.execution / 'execution_contract.json')
    names = ['base_snapshot_H0', 'flat_snapshot_H0']
    rows = {r['id']: r for r in v['groups']}
    assert rows[names[0]]['observer_receiver_and_source_bitwise_equal']
    dt = m['dt_s']
    selected = m['snapshot_iterations'][::2]
    observations = {}
    fields = {}
    native = {}
    geometries = {}
    geometries_by_group = {g['id']: g for g in c['groups']}
    receipts = []
    for name in names:
        row = rows[name]
        assert len(row['snapshots']) == 250
        path = Path(row['native_path'])
        assert sha(path) == row['native_sha256']
        with h5py.File(path) as h:
            native[name] = h['rxs/rx1/Ez'][:]
            assert h['rxs/rx1/Ez'].attrs['TimeSampleOffset'] == 0
        observations[name] = {r['iteration']: r for r in row['snapshots']}
        arrays = []
        for j in selected:
            record = observations[name][j]
            path = Path(record['file'])
            assert sha(path) == record['sha256']
            with h5py.File(path) as h:
                assert h.attrs['gprMax'] == '4.0.1' and h.attrs['iteration'] == j
                assert abs(h.attrs['time'] - j * dt) < 1e-20
                assert abs(h.attrs['magnetic_time'] - (j - .5) * dt) < 1e-20
                x = h['Ez'][:, :, 0]
                assert x.dtype == np.float64 and x.shape == (350, 315) and np.isfinite(x).all()
                arrays.append(x.T)
            receipts.append(dict(group=name, iteration=j, snapshot_sha256=record['sha256']))
        fields[name] = arrays
        g = geometries_by_group[name]
        path = Path(c['package']) / g['geometry']
        assert sha(path) == g['geometry_sha256']
        with h5py.File(path) as h:
            material = h['data'][5800:7200, :, 0]
        assert material.shape == (1400, 1700) and set(np.unique(material)) == {0, 1, 2}
        # Both observed H0 models have contiguous mud/cover/air columns.
        assert (np.diff(material.astype(np.int16), axis=1) <= 0).all()
        geometries[name] = (np.count_nonzero(material != 0, axis=1) * .025,
                            np.count_nonzero(material == 2, axis=1) * .025)
    with h5py.File(Path(c['package']) / m['baseline']['native']) as h:
        assert h['rxs/rx1/Ez'][:].tobytes() == native[names[0]].tobytes()
        source = h['srcs/src1/excitation/samples'][:]
        source_offset = float(h['srcs/src1/excitation'].attrs['TimeSampleOffset'])
        source_peak_ns = float((np.argmax(abs(source)) * dt + source_offset) * 1e9)
    differences = [x - y for x, y in zip(fields[names[0]], fields[names[1]])]
    common = max(float(abs(x).max()) for name in names for x in fields[name])
    separate = max(float(abs(x).max()) for x in differences)
    assert common > 0 and separate > 0
    a.out.mkdir(parents=True)
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(figsize=(13, 6), layout='constrained')
    grid = fig.add_gridspec(2, 3, height_ratios=[3, 1])
    axes = [fig.add_subplot(grid[0, k]) for k in range(3)]
    titles = ['原非平几何H0：保留覆盖层底', '全横向平层H0：绝对收发位置不变', '原场−平层：全部几何差场，非纯多次波']
    images = []
    xx = 145 + (np.arange(1400) + .5) * .025
    for k, ax in enumerate(axes):
        limit = common if k < 2 else separate
        im = ax.imshow(np.zeros((315, 350)), origin='lower', extent=[145, 180, 10, 41.5], cmap='gray',
                       aspect='equal', norm=SymLogNorm(limit * 1e-5, vmin=-limit, vmax=limit, base=10))
        images.append(im)
        geo = geometries[names[min(k, 1)]] if k < 2 else geometries[names[0]]
        ax.plot(xx, geo[0], color='lime', lw=.7)
        ax.plot(xx, geo[1], color='cyan', lw=.7)
        ax.axvspan(154, 164, facecolor='none', edgecolor='magenta', lw=.6, hatch='//', alpha=.35)
        ax.axvspan(164, 174, facecolor='none', edgecolor='orange', lw=.6, hatch='\\', alpha=.35)
        g = geometries_by_group[names[0]]
        ax.scatter([g['tx_m'][0], g['rx_m'][0]], [g['tx_m'][1], g['rx_m'][1]], c=['orange', 'red'], s=14)
        ax.set(xlabel='模型local x / m（剖面里程=x+20m）', ylabel='模型y / m', title=titles[k],
               xlim=(145, 180), ylim=(10, 41.5))
        fig.colorbar(im, ax=ax, shrink=.8, label='原生Ez / V/m；固定对称对数灰度')
    graph = fig.add_subplot(grid[1, :])
    tn = np.arange(len(native[names[0]])) * dt * 1e9
    for name, label in zip(names, titles):
        graph.plot(tn, native[name], lw=.8, label=label)
    keep = (tn >= 250) & (tn <= 450)
    limit = max(float(abs(native[name][keep]).max()) for name in names)
    graph.set(xlim=(250, 450), ylim=(-1.05 * limit, 1.05 * limit), xlabel='原生时间 / ns（含源延迟）', ylabel='接收Ez / V/m')
    graph.legend(fontsize=8)
    cursor = graph.axvline(0, color='red', lw=1)
    title = fig.suptitle('')
    frames = []
    static = []
    targets = {min(selected, key=lambda j: abs(j * dt * 1e9 - t)) for t in [100, 150, 200, 250, 300, 330, 360, 400, 450]}
    for index, j in enumerate(selected):
        for im, values in zip(images, [fields[names[0]][index], fields[names[1]][index], differences[index]]):
            im.set_data(values)
        cursor.set_xdata([j * dt * 1e9] * 2)
        title.set_text(f'190m原生Ricker波场 t={j * dt * 1e9:.2f}ns，非501点SFCW；源峰{source_peak_ns:.3f}ns\n绿线地表、青线覆盖层底；均去深部底砂。紫/橙框仅标另两控制的远/近ROI，本动画不替换这些ROI。')
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()).convert('P', palette=Image.Palette.ADAPTIVE, colors=128))
        if j in targets:
            path = a.out / f'wavefield_{j:05d}.png'
            fig.savefig(path, dpi=120)
            static.append(dict(file=path.name, iteration=j, sha256=sha(path)))
    gif = a.out / 'cover_wavefield.gif'
    frames[0].save(gif, save_all=True, append_images=frames[1:], duration=90, loop=0, optimize=True)
    plt.close(fig)
    result = dict(status='AUDITED_ORIGINAL_FLAT_NATIVE_RICKER_MOVIE_NOT_SFCW_OR_UNIQUE_PATH',
        renderer_sha256=sha(__file__), contract_sha256=sha(a.execution / 'execution_contract.json'),
        verification_sha256=sha(a.execution / 'completed_verification.json'), frames=len(frames),
        common_abs_limit_V_per_m=common, difference_abs_limit_V_per_m=separate,
        symlog_linear_threshold_fraction=1e-5, source_peak_native_ns=source_peak_ns,
        electric_time_step_ns=68 * dt * 1e9, magnetic_time_offset_s=-dt / 2,
        source_time_offset_s=source_offset, rendered_last_electric_time_ns=selected[-1] * dt * 1e9,
        selected_snapshot_hashes=receipts, static_frames=static, gif_bytes=gif.stat().st_size, gif_sha256=sha(gif),
        observer_receiver_bitwise_equal=True,
        limits='Native Ricker movie, not SFCW. Green/cyan indicate actual ground/cover interfaces for each geometry. No bottom sandstone remains. Regional far/near controls are not animated; their outlines are guides only. Difference includes entire geometry change, not pure scattering order.')
    (a.out / 'render_receipt.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], frames=len(frames), gif_bytes=result['gif_bytes'])))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
