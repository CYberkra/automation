"""Read-only categorical preview of a legacy ZYX geometry and companion table.

Output is a private derived audit: use artifacts/local_checks, not research_checks.
No axis conversion is written to the source; no solver is called.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def preview(geometry, materials, out):
    if out.exists():
        raise ValueError('refuse to overwrite a previous preview')
    rows = []
    for line in materials.read_text('utf-8-sig').splitlines():
        fields = line.split()
        if fields and fields[0] == '#material:':
            if len(fields) != 6:
                raise ValueError('expected legacy scalar material definitions')
            rows.append(dict(id=len(rows), name=fields[5],
                             epsilon_r=float(fields[1]), sigma_S_m=float(fields[2]),
                             mu_r=float(fields[3]), sigma_magnetic=float(fields[4])))
        elif fields and not fields[0].startswith('##'):
            raise ValueError('unreviewed companion material command')
    if [r['name'] for r in rows] != ['air', 'soil', 'mudstone', 'sandstone']:
        raise ValueError('this categorical palette requires the reviewed four-material order')

    geometry_hash = digest(geometry)
    materials_hash = digest(materials)
    with h5py.File(geometry, 'r') as handle:
        data = handle['data']
        native_dtype = str(data.dtype)
        size = tuple(int(np.asarray(handle.attrs[k]).item()) for k in ('nx', 'ny', 'nz'))
        spacing = np.asarray(handle.attrs['dx_dy_dz'], dtype=float)
        origin = np.array([np.asarray(handle.attrs[k]).item() for k in ('x0', 'y0', 'z0')])
        if data.shape != size[::-1] or data.ndim != 3:
            raise ValueError('ZYX storage inference requires shape matching nz,ny,nx')
        if any(v <= 0 for v in size) or not np.all(spacing > 0):
            raise ValueError('invalid geometry dimensions')
        if (size != (7000, 9000, 2) or
                not np.array_equal(spacing, [.005, .005, .002]) or np.any(origin != 0)):
            raise ValueError('preview layout is specific to the reviewed 35x45m Line9 geometry')
        image = np.empty((size[1]//5, size[0]//5), dtype=np.uint8)
        counts = np.zeros((size[2], len(rows)), dtype=np.int64)
        different = 0
        # Bounded-memory reading: at most two 500-row blocks in memory.
        for start in range(0, size[1], 500):
            end = min(start+500, size[1])
            first = data[0, start:end, :]
            second = data[1, start:end, :]
            for plane, block in enumerate((first, second)):
                if block.min() < 0 or block.max() >= len(rows):
                    raise ValueError('unknown material ID')
                counts[plane] += np.bincount(block.ravel(), minlength=len(rows))
            different += np.count_nonzero(first != second)
            image[start//5:end//5] = first[::5, ::5]
    if digest(geometry) != geometry_hash or digest(materials) != materials_hash:
        raise ValueError('source changed while reading')
    out.mkdir(parents=True)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch

    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei']
    plt.rcParams['axes.unicode_minus'] = False
    extent = [origin[0], origin[0]+size[0]*spacing[0],
              origin[1], origin[1]+size[1]*spacing[1]]
    colors = ['#dceff9', '#b97445', '#7d8494', '#e4c982']
    names = ['空气', '土层', '泥岩', '砂岩']
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-.5, len(rows)+.5), len(rows))
    fig = plt.figure(figsize=(12, 10), layout='constrained')
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 1.05], height_ratios=[1.2, 1])
    full = fig.add_subplot(grid[:, 0])
    detail = fig.add_subplot(grid[0, 1])
    for ax in (full, detail):
        ax.imshow(image, cmap=cmap, norm=norm, origin='lower', extent=extent,
                  interpolation='nearest', aspect='equal')
        ax.set_xlabel('x：模型坐标 (m)')
        ax.set_ylabel('y：模型坐标 (m)')
    full.set_title('全剖面：35 m × 45 m（比例一致）')
    detail.set(xlim=extent[:2], ylim=(18, 38), title='局部：层间界面与 ID0 区域')
    note = fig.add_subplot(grid[1, 1])
    note.axis('off')
    text = '\n'.join([
        '材料名称按随附 materials.txt 行序解释：',
        *[f'ID {r["id"]}  {names[r["id"]]}：εr={r["epsilon_r"]:g}，σ={r["sigma_S_m"]:g} S/m' for r in rows],
        '',
        '网格：dx=dy=5 mm，dz=2 mm',
        '单元数：nx=7000，ny=9000，nz=2',
        '薄方向厚度：4 mm；图为第一个 z 切片',
        f'两个 z 切片不同单元数：{different:,}',
        '',
        '蓝色包括外部空气及内部 ID0；原样展示，未填补。',
        '上边缘可见 ID1 薄带，不自动视为地表。',
        '模型坐标不是现场绝对高程，未标注推测航高。',
        '显示每5格取1格（25 mm）；不修改原始体素。',
    ])
    note.text(0, .98, text, va='top', fontsize=10, linespacing=1.6)
    fig.legend(handles=[Patch(facecolor=c, label=f'ID{i} {n}')
                        for i, (c, n) in enumerate(zip(colors, names))],
               loc='outside lower center', ncol=4, frameon=False)
    fig.suptitle('unified_geometry.h5：原始材料几何二维示意\n按属性推断 Z–Y–X 存储；材料表配对待核验；无天线或飞行轨迹', fontsize=14)
    figure = out/'unified_geometry_2d.png'
    fig.savefig(figure, dpi=160)
    plt.close(fig)
    summary = dict(status='READ_ONLY_PRIVATE_GEOMETRY_PREVIEW_NOT_SOLVER_VALIDATION',
                   geometry_sha256=geometry_hash, materials_sha256=materials_hash,
                   script_sha256=digest(Path(__file__)),
                   dimensions_xyz_m=(np.array(size)*spacing).tolist(),
                   cells_xyz=list(size), spacing_xyz_m=spacing.tolist(),
                   origin_xyz_m=origin.tolist(), native_array_shape_zyx=list(size[::-1]),
                   native_dtype=native_dtype, material_mapping='Companion line order, not independently certified',
                   materials=rows, counts_by_z_and_material=counts.tolist(),
                   z_slice_differing_cells=int(different), display_decimation_xy=5,
                   calls_solver=False, originals_unchanged=True,
                   figure_sha256=digest(figure), source_relative_height_m=None)
    (out/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=Path, required=True)
    parser.add_argument('--materials', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    preview(args.geometry, args.materials, args.out)
