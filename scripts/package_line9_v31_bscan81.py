"""Package the historical 81-trace Line9 run; no solver or processing calls."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

import h5py
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--public', type=Path, required=True)
    a = ap.parse_args()
    repo = Path(__file__).resolve().parents[1]
    base = a.source.resolve()
    retrieved, stage = base / 'retrieved_r2', base / 'gpu_stage_r2'
    c = load(retrieved / 'resolved_contract.json')
    done = load(retrieved / 'completion.json')
    receipts = load(retrieved / 'receipts.json')
    ledger = {r['path'].replace('\\', '/'): r for r in done['files']}
    assert done['trace_count'] == c['trace_count'] == len(receipts) == 81
    assert done['gprMax'] == '4.0.1' and done['native_dtype'] == 'float32'
    assert c['official_sfcw']['steps'] == 501 and c['official_sfcw']['f_stop_hz'] == 170_000_000
    assert done['receipts'] == receipts
    assert not a.out.exists() and not a.public.exists()
    a.out.mkdir(parents=True)
    a.public.mkdir(parents=True)
    package = a.out / 'package'
    package.mkdir()
    bindings = []

    def copy(src, target, historical_key=None):
        dst = package / target
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        digest = sha(dst)
        assert digest == sha(src)
        if historical_key:
            record = ledger[historical_key]
            assert digest == record['sha256'] and dst.stat().st_size == record['bytes'], historical_key
            bindings.append(dict(package_path=target, historical_path=historical_key, sha256=digest))
        return dst

    # Initial staged package is retained separately from the per-trace inputs.
    for src in sorted((stage / 'package').rglob('*')):
        if not src.is_file() or '__pycache__' in src.parts:
            continue
        rel = src.relative_to(stage / 'package').as_posix()
        mutable = rel in ['config/pilot_ascan.json', 'model/line9_pilot_geometry.h5']
        copy(src, 'source_package_staged/' + rel, None if mutable else 'package/' + rel)
    assert sha(stage / 'original.zip') == c['archive_sha256']
    copy(stage / 'original.zip', 'original_v31_package.zip', 'original.zip')
    copy(stage / 'run_line9_v31_bscan_gpu.py', 'historical_execution/run_line9_v31_bscan_gpu.py', 'run_line9_v31_bscan_gpu.py')
    copy(stage / 'worker_patch.json', 'historical_execution/worker_patch.json', 'worker_patch.json')
    copy(retrieved / 'candidate_contract.json', 'historical_execution/candidate_contract.json', 'candidate_contract.json')
    for name in ['candidate_plans.json', 'candidate_geometry_manifest.json']:
        copy(stage / name, 'historical_execution/' + name, name)
    for name in ['attempt.json', 'resolved_contract.json', 'receipts.json']:
        copy(retrieved / name, 'execution_evidence/' + name, name)
    copy(retrieved / 'completion.json', 'execution_evidence/completion.json')
    for src in sorted((retrieved / 'results').glob('*.h5')):
        copy(src, 'results/' + src.name, 'package/results/' + src.name)
    assert len(list((package / 'results').glob('*.h5'))) == 83

    geometry_rows, trace_rows = [], []
    merged_path = package / 'results/line9_v31_bscan_merged.h5'
    with h5py.File(merged_path) as merged:
        matrix = merged['rxs/rx1/Ez'][:]
        assert matrix.shape == (19638, 81) and matrix.dtype == np.dtype('float32')
        assert np.isfinite(matrix).all()
        for index, (trace, receipt) in enumerate(zip(c['trace_manifest'], receipts), 1):
            assert receipt['trace_index'] == index and receipt['status'] == 'PASS'
            assert receipt['station_global_x_m'] == trace['station_global_x_m']
            geometry = trace['geometry_file']
            dst = copy(base / 'gpu_candidate_geometries' / geometry,
                       'model_geometries_used/' + geometry, 'candidate_geometries/' + geometry)
            assert sha(dst) == trace['geometry_sha256'] == receipt['geometry_sha256']
            geometry_rows.append(dict(trace_index=index, **trace))
            native = package / 'results' / f'line9_v31_bscan{index}.h5'
            assert sha(native) == receipt['output_sha256']
            with h5py.File(native) as h:
                y = h['rxs/rx1/Ez'][:]
                assert np.array_equal(y, matrix[:, index - 1])
                assert y.dtype == np.dtype('float32') and y.shape == (19638,)
                assert h.attrs['gprMax'] == '4.0.1'
                assert np.isclose(h.attrs['dt'], receipt['dt_s'], rtol=1e-12, atol=0)
                assert np.allclose(h['srcs/src1'].attrs['Position'][:2], trace['tx_local_m'], atol=1e-8, rtol=0)
                assert np.allclose(h['rxs/rx1'].attrs['Position'][:2], trace['rx_local_m'], atol=1e-8, rtol=0)
                excitation = h['srcs/src1/excitation']
                source = excitation['samples'][:]
                assert np.count_nonzero(source) == 1 and np.max(source) == 1
                assert excitation.attrs['TimeSampleOffset'] == h.attrs['dt'] / 2
                peak, tail = np.max(np.abs(y)), np.max(np.abs(y[int(.95 * len(y)):]))
                tail_db = float(20 * np.log10(float(tail) / float(peak)))
            trace_rows.append(dict(index=index, station_m=trace['station_global_x_m'], tail_db=tail_db))

    sfcw = package / 'results/line9_v31_bscan_sfcw.h5'
    assert sha(merged_path) == done['official_merge']['sha256']
    assert sha(sfcw) == done['official_sfcw']['sha256']
    with h5py.File(sfcw) as h:
        assert np.array_equal(h['frequency'][:], np.linspace(20e6, 170e6, 501))
        assert h['response'].shape == (501, 81) and np.isfinite(h['response'][:]).all()
        assert h.attrs['Method'] == 'direct' and h.attrs['TailTaperFraction'] == 0
        assert h['time_response'].attrs['Window'] == 'gaussian'
        assert h['time_response'].attrs['ZeroPadFactor'] == 4 and h['time_response'].attrs['TimeShift'] == 0
    save(package / 'model_geometries_used/stations.json', geometry_rows)

    public_old = repo / 'artifacts/research_checks/2026-10-13_line9_v31_bscan150_198_gpu_r2'
    for src in sorted(public_old.iterdir()):
        if src.is_file():
            copy(src, ('figures/' if src.suffix == '.png' else 'historical_audit/') + src.name)
    copy(repo / 'docs/research/2026-10-13_line9_v31_bscan150_198_gpu.md', 'report_original.md')
    for name in ['plot_bscan.py', 'plot_bscan_windowed.py', 'audit_retrieved_r2.py']:
        copy(base / name, 'historical_analysis_code/' + name)
    copy(Path(__file__), 'package_code/package_line9_v31_bscan81.py')
    audit = dict(status='PACKAGING_AND_HISTORICAL_HASH_BINDING_PASS', new_solves=0,
        new_sfcw_processing=0, trace_count=81, used_geometry_count=81, result_h5_count=83,
        native_dtype='float32', native_shape=[19638, 81], official_response_shape=[501, 81],
        grid_m=.025, averaging='n', record_ns=1100, frequency_hz=[20_000_000, 170_000_000],
        source_and_receiver_positions_verified=True, native_columns_bit_identical_to_merge=True,
        completion_bound_files=len(bindings), file_bindings=bindings, traces=trace_rows,
        tail_screen_passed=False, missing_local_historical_logs='81 per-trace solver logs and official merge log; hashes retained in original completion ledger/receipts. Official SFCW log included.',
        historical_naming_note='Resolved task_id retains_r1; executed successful root is_r2. Original contracts unchanged.',
        staged_mutable_files_note='Initial pilot config and pilot geometry are staged snapshots, not runtime-end snapshots. Actual81 geometries/poses are in model_geometries_used/stations.json and resolved_contract.json.',
        prepared_count_note='Historical candidate prepared97 geometries; this archive includes the81 actually used. Candidate planning metadata preserved unchanged.',
        raw_pdf_included=False, geometry_source_vectors_and_calibration_included=True)
    save(package / 'packaging_audit.json', audit)
    save(a.public / 'packaging_audit.json', audit)
    (package / 'README.md').write_text('''# 九号测线150–198米：81道模型与结果

这是2026-10-13在ROG完成的历史批次。2.5厘米网格、50×50米局部域、5米名义离地、PML120、1100ns、原生float32；gprMax4.0.1 GPU0 single。官方impulse/direct SFCW为20–170MHz/501点/0.3MHz/Gaussian/补零4。介电平滑为n。

双击“查看Bscan.html”看四张灰度图。原始波形和数值响应在results目录，图像仅作展示。模型不是B-scan截图。

- source_package_staged/model/line9_350m_voxels.h5：350米母模型；contacts_m.json及材料JSON同目录。
- model_geometries_used：实际81个局部几何H5，stations.json逐道列出站位、收发坐标和哈希。
- source_package_staged：当时暂存的构建/正演/SFCW脚本、配置与重建矢量。pilot配置和pilot几何是初始暂存快照，运行时被逐道替换；实际输入以execution_evidence/resolved_contract.json和81个已核验几何为准。
- original_v31_package.zip：原V3.1提供包，保持字节不变；其中单道初始配置不是81道运行契约。
- results：line9_v31_bscan1.h5至81.h5、line9_v31_bscan_merged.h5、line9_v31_bscan_sfcw.h5，共83份结果H5。原生合并Ez形状19638×81；官方复频谱501×81。
- execution_evidence：原完成记录、81份回执和已解析运行契约。历史task_id带r1而实际成功目录为r2，差异说明在historical_audit中，未改写原文件。
- historical_execution/historical_analysis_code：当时运行与出图代码。运行脚本含按首道耗时选择道数的历史逻辑；再次执行不保证选到同一道数，也会启动仿真。此包只交付资料，不自动执行。
- historical_audit、packaging_audit.json：历史审计、来源绑定和本次装包核查。逐道求解日志及合并日志本地未留存，原哈希保留；官方SFCW日志已包含。
- manifest.json：除自身外所有包内文件的SHA-256和长度。

本包包含重建矢量、比例尺配置和模型，不含原始PDF，也不含求解器环境或凭据。历史分析脚本保留当时路径，不声称解压后可直接运行。物理结论见report_original.md：81道尾部均未满足−60dB建议初筛，不能称时窗收敛或实测验证。当前20–110MHz/3.333厘米配置属于之后的另一批，不混入本包。
''', encoding='utf-8')
    figures = [
        ('line9_v31_native_impulse_bscan.png', '原生冲激Ez：全记录'),
        ('line9_v31_native_impulse_bscan_windowed.png', '原生冲激Ez：分时窗显示'),
        ('line9_v31_official_direct_sfcw_bscan.png', '官方direct SFCW：全周期显示'),
        ('line9_v31_official_direct_sfcw_bscan_windowed.png', '官方direct SFCW：分时窗显示')]
    html = '<!doctype html><meta charset="utf-8"><title>九号测线150–198米81道B-scan</title><style>body{font:16px/1.7 sans-serif;max-width:1280px;margin:30px auto;padding:20px}img{width:100%;height:auto}</style><h1>九号测线150–198米／81道／2.5厘米</h1><p>历史批次：20–170MHz、501点，平滑n，1100ns，5米名义离地。各图的色标及窗口说明保留在图中；分时窗面板不可直接用亮度比较幅度。SFCW周期约3333ns不代表1100ns以后有额外传播记录。</p><p><a href="README.md">文件说明</a> · <a href="report_original.md">原报告</a></p>'
    for name, title in figures:
        assert (package / 'figures' / name).is_file()
        html += f'<h2>{title}</h2><a href="figures/{name}"><img src="figures/{name}" alt="{title}"></a>'
    (package / '查看Bscan.html').write_text(html, encoding='utf-8')
    manifest = dict(schema='line9-private-package/1', files=[
        dict(path=f.relative_to(package).as_posix(), sha256=sha(f), bytes=f.stat().st_size)
        for f in sorted(package.rglob('*')) if f.is_file()])
    save(package / 'manifest.json', manifest)
    archive = a.out / 'Line9_150-198m_81道_2.5cm_模型与仿真结果完整包.zip'
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for f in sorted(package.rglob('*')):
            if f.is_file():
                z.write(f, f.relative_to(package).as_posix())
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist()) == len(manifest['files']) + 1
        for row in manifest['files']:
            data = z.read(row['path'])
            assert len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256']
        assert z.read('manifest.json') == (package / 'manifest.json').read_bytes()
    index = dict(status='ZIP_ALL_MEMBERS_HASH_VERIFIED', date='2026-10-16',
        zip_path=str(archive.resolve()), sha256=sha(archive), bytes=archive.stat().st_size,
        archive_members=len(manifest['files']) + 1, used_geometries=81, result_h5_count=83,
        mother_geometry_included=True, models_and_results_hash_bound_to_historical_run=True,
        new_solves=0, native_data_or_model_modified=False,
        packaging_script_sha256=sha(__file__), private_data_uploaded=False)
    save(a.public / 'delivery_index.json', index)
    print(json.dumps(index, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
