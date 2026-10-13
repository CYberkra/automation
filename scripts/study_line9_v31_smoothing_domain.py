"""Prepare and run the authorized V3.1 smoothing and domain comparisons.

Private geometry stays in the supplied study folder. All FDTD and SFCW work
uses gprMax 4.0.1 without modifying its numerical kernels.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np

C = 299792458.0
DB = "line9_wang2026_constant_v2"
ANCHORS = [150.0, 158.4, 162.0, 174.0, 186.0, 190.2, 198.0]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def duration_seconds(text):
    units = {"hour": 3600.0, "minute": 60.0, "second": 1.0}
    parts = re.findall(r"([0-9.]+)\s+(hour|minute|second)s?", text)
    if not parts:
        raise ValueError(f"No duration in {text!r}")
    return sum(float(value) * units[unit] for value, unit in parts)


def surface_polyline(package):
    sel = load(package / "config/selection_pdf_reconstructed.json")
    vec = load(package / "reference/line9_vectors_extracted.json")
    assert vec["source_pdf_sha256"] == sel["source_pdf_sha256"]
    spec = importlib.util.spec_from_file_location("original_v31_builder", package / "01_build_geometry.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    p = mod.join_contact(vec["paths"], sel["contacts"]["surface"], sel["contact_join_tolerance_pt"])
    cal = sel["calibration"]
    sx = (p[-1, 0] - p[0, 0]) / cal["profile_width_m"]
    sy = cal["vertical_reference"]["vertical_units_per_m"]
    collar_y = np.interp(cal["borehole_pdf_x_pt"], p[:, 0], p[:, 1])
    return np.column_stack(((p[:, 0] - p[0, 0]) / sx,
                            (p[:, 1] - collar_y) / sy + cal["borehole_collar_elevation_m"]))


def write_geometry(path, data, attrs, dl, origin_x, bottom=410.0):
    with h5py.File(path, "w") as h:
        h.create_dataset("data", data=data, compression="gzip", compression_opts=4)
        h.create_dataset("material_keys", data=np.asarray(["air", "cover", "mudstone", "sandstone"], dtype="S"))
        for k, v in attrs.items():
            h.attrs[k] = v
        h.attrs["dx_dy_dz"] = [dl, dl, dl]
        h.attrs["shape_nxyz"] = data.shape
        h.attrs["origin_xyz"] = [0.0, 0.0, 0.0]
        h.attrs["global_profile_x_origin_m"] = origin_x
        h.attrs["global_profile_x_end_m"] = origin_x + data.shape[0] * dl
        h.attrs["elevation_offset_m"] = bottom
        h.attrs["MaterialDatabase"] = DB


def vertical_events(data, tx, rx, dl):
    ix = data.shape[0] // 2
    col = data[ix, :, 0]
    k = int(np.flatnonzero(col != 0)[-1])
    surface = (k + 1) * dl
    t = ((tx[1] - surface) + (rx[1] - surface)) / C * 1e9
    events = [{"label": "surface_vertical_reference", "time_ns": float(t)}]
    eps = {1: 12.0, 2: 18.0, 3: 9.0}
    previous = int(col[k])
    for j in range(k, -1, -1):
        t += 2 * dl * math.sqrt(eps[int(col[j])]) / C * 1e9
        if j and col[j - 1] != col[j]:
            events.append({"label": f"{previous}_to_{int(col[j-1])}_vertical_reference", "time_ns": float(t)})
            previous = int(col[j - 1])
    return events


def prepare(a):
    root = a.root.resolve()
    if root.exists():
        raise FileExistsError(root)
    root.mkdir(parents=True)
    geomdir = root / "geometries"
    geomdir.mkdir()
    package = a.source_package.resolve()
    bscan = a.bscan_root.resolve()
    orig = load(bscan / "retrieved_r2/resolved_contract.json")
    completion = load(bscan / "retrieved_r2/completion.json")
    assert orig["trace_count"] == 81 and orig["record_ns"] == 1100.0
    assert orig["gpu_precision"] == "single"
    assert sha(a.mother) == orig["mother_geometry_sha256"]
    shutil.copy2(package / f"model/{DB}.json", geomdir / f"{DB}.json")
    poly = surface_polyline(package)
    trace_map = {round(float(t["station_global_x_m"]), 4): t for t in orig["trace_manifest"]}
    receipt_map = {round(float(t["station_global_x_m"]), 4): t for t in completion["receipts"]}
    cases = []
    base_info = {}

    def case(cid, station, geo, data, attrs, dl, averaging, stage, family, geometry_only=False):
        rec = receipt_map[round(station, 4)]
        old = trace_map[round(station, 4)]
        xorigin = float(attrs["global_profile_x_origin_m"])
        bottom = float(attrs["elevation_offset_m"])
        tx = [station - 0.65 - xorigin, rec["tx_local_m"][1] + 410.0 - bottom, 0.0]
        rx = [station + 0.65 - xorigin, rec["rx_local_m"][1] + 410.0 - bottom, 0.0]
        domain = [data.shape[0] * dl, data.shape[1] * dl]
        for p in (tx, rx):
            assert all(abs(v / dl - round(v / dl)) < 1e-7 for v in p[:2])
            assert min(p[0] - 3, domain[0] - 3 - p[0]) >= 10
            assert domain[1] - 3 - p[1] >= 2
            assert data[round(p[0] / dl), round(p[1] / dl), 0] == 0
        dt = 0.95 * dl / (C * math.sqrt(2))
        gx = np.linspace(station - 2, station + 2, 41)
        heights = np.interp(gx, poly[:, 0], poly[:, 1]) - bottom
        patch = [[round((station - 2 - xorigin) / dl) * dl,
                  math.floor((heights.min() - .15) / dl) * dl, 0.0],
                 [round((station + 2 - xorigin) / dl) * dl,
                  math.ceil((heights.max() + .15) / dl) * dl, dl]]
        obj = dict(id=cid, stage=stage, family=family, station_m=station,
                   geometry=str(geo.relative_to(root)).replace("\\", "/"), geometry_sha256=sha(geo),
                   averaging=averaging, grid_m=dl, domain_m=domain, tx_m=tx, rx_m=rx,
                   tx_global_x_elevation_m=[station - .65, tx[1] + bottom],
                   rx_global_x_elevation_m=[station + .65, rx[1] + bottom],
                   origin_global_x_elevation_m=[xorigin, bottom],
                   pml_cells=[round(3 / dl), round(3 / dl), 0, round(3 / dl), round(3 / dl), 0],
                   dt_s=dt, expected_samples=math.ceil(1100e-9 / dt) + 1,
                   geometry_only=geometry_only, surface_coefficient_export_bounds=patch,
                   baseline_trace_index=old["trace_index"] if "trace_index" in old else rec["trace_index"])
        obj["vertical_reference_events"] = base_info.get(station, {}).get("vertical_reference_events", [])
        cases.append(obj)
        return obj

    full_geos = {}
    for station in ANCHORS:
        t = trace_map[round(station, 4)]
        src = bscan / "gpu_candidate_geometries" / t["geometry_file"]
        assert sha(src) == t["geometry_sha256"]
        geo = geomdir / f"full_x{station:.1f}.h5"
        shutil.copy2(src, geo)
        with h5py.File(geo) as h:
            data = h["data"][:]
            attrs = dict(h.attrs)
        r = receipt_map[round(station, 4)]
        base_info[station] = {"vertical_reference_events": vertical_events(data, r["tx_local_m"], r["rx_local_m"], .025)}
        full_geos[station] = (geo, data, attrs)
        case(f"full_x{station:.1f}_y", station, geo, data, attrs, .025, "y", "A", "full_smoothing")
    geo, data, attrs = full_geos[174.0]
    case("full_x174.0_n_geometry_audit", 174.0, geo, data, attrs, .025, "n", "A", "geometry_audit", True)

    slopes = {s: abs(float((np.interp(s + 2, poly[:, 0], poly[:, 1]) -
                             np.interp(s - 2, poly[:, 0], poly[:, 1])) / 4)) for s in ANCHORS}
    surface_station = max(slopes, key=slopes.get)
    srcgeo, source_data, attrs = full_geos[surface_station]
    for dl in [.025, .0125]:
        nx = ny = round(50 / dl)
        origin = surface_station - 25
        x = origin + (np.arange(nx) + .5) * dl
        y = 410 + (np.arange(ny) + .5) * dl
        heights = np.interp(x, poly[:, 0], poly[:, 1])
        data = (y[None, :, None] < heights[:, None, None]).astype(np.int16)
        if dl == .025:
            assert np.array_equal(data != 0, source_data != 0), "surface reconstruction changed the coarse geometry"
        surfg = geomdir / f"surface_only_d{dl:g}.h5"
        write_geometry(surfg, data, attrs, dl, origin)
        for avg in ["n", "y"]:
            case(f"surface_d{dl:g}_{avg}", surface_station, surfg, data, attrs, dl, avg, "A", "surface_only")
        air = np.zeros_like(data)
        airg = geomdir / f"air_d{dl:g}.h5"
        write_geometry(airg, air, attrs, dl, origin)
        case(f"air_d{dl:g}", surface_station, airg, air, attrs, dl, "y", "A", "air_reference")

    # Larger-domain controls use the same parent geometry, source coordinates,
    # lower boundary and grid. Reduced domains are exact integer-cell slices.
    with h5py.File(a.mother) as parent:
        parent_data = parent["data"]
        for station in [150.0, 174.0, 198.0]:
            origin = station - 35
            data = parent_data[round(origin / .025):round((origin + 70) / .025), :, :]
            attrs = dict(parent.attrs)
            attrs["global_profile_x_origin_m"] = origin
            geo = geomdir / f"wide70_x{station:.1f}.h5"
            write_geometry(geo, data, attrs, .025, origin)
            case(f"wide70_x{station:.1f}", station, geo, data, attrs, .025, "y", "B", "expanded_domain")
    for station, (fullg, full, attrs0) in full_geos.items():
        for width, height, stage in [(50, 45, "B"), (40, 50, "B"), (40, 45, "B"), (30, 45, "C")]:
            ix0 = round((50 - width) / 2 / .025)
            ix1 = ix0 + round(width / .025)
            iy1 = round(height / .025)
            assert not np.any(full[:, iy1:, :] != 0), "vertical crop would remove geological material"
            data = full[ix0:ix1, :iy1, :].copy()
            attrs = dict(attrs0)
            origin = station - width / 2
            attrs["global_profile_x_origin_m"] = origin
            geo = geomdir / f"crop{width:g}x{height:g}_x{station:.1f}.h5"
            write_geometry(geo, data, attrs, .025, origin)
            case(f"crop{width:g}x{height:g}_x{station:.1f}", station, geo, data, attrs, .025, "y", stage, "reduced_domain")

    contract = dict(schema="line9-v31-smoothing-domain/1", id="line9_v31_smoothing_domain_20261013_r1",
                    authorization="User requested native dielectric smoothing followed by physics-preserving domain optimization",
                    gprMax="4.0.1", gpu=[0], gpu_precision="single", native_dtype="float32",
                    model_family="V3.1 candidate, 2D TMz, constant original materials, 5m nominal AGL",
                    record_ns=1100.0, stability_factor=.95, source="built-in unit impulse; ideal HertzianDipole z",
                    official_sfcw=dict(method="direct", f_start_hz=20_000_000, f_stop_hz=170_000_000,
                                       steps=501, window="gaussian", zero_pad=4, tail_taper=0, time_shift_s=0),
                    private_input_sha256={"mother_geometry": sha(a.mother), "source_zip": orig["archive_sha256"],
                                           "selection": sha(package / "config/selection_pdf_reconstructed.json"),
                                           "vectors": sha(package / "reference/line9_vectors_extracted.json")},
                    material_database_sha256=sha(geomdir / f"{DB}.json"), anchors_m=ANCHORS,
                    surface_diagnostic_station_m=surface_station, surface_4m_average_slopes=slopes,
                    surface_diagnostic="same continuous V3.1 surface, all subsurface material replaced by the original cover; fine-grid proxy is not analytic truth",
                    cross_grid_diagnostic="R=(H_surface-H_air)/H_air, computed from official complex responses; cancels grid-dependent 2D source scaling; not a replacement SFCW product",
                    max_new_solves=sum(not c["geometry_only"] for c in cases), max_wall_seconds=7200,
                    no_retry_after_case_failure=True,
                    domain_screen=dict(complex_band_relative_l2_max=.005,
                                       strong_bin_phase_rms_degrees_max=.5,
                                       time_window_relative_l2_max=.01, event_peak_relative_error_max=.01,
                                       event_peak_shift_ns_max=2.0, common_windows_ns=[[100, 400], [400, 1000]],
                                       event_halfwidth_ns=15,
                                       note="project-specific screening tolerances, not universal standards; requires all anchors and the 70m control; original 50m domain is not an absolute truth"),
                    stage_C_condition="Run 30x45 only if all 40x45 cases and 70m reference checks pass the frozen screen",
                    cases=cases)
    save(root / "contract.json", contract)
    shutil.copy2(Path(__file__), root / "study.py")
    print(json.dumps({"prepared_root": str(root), "case_count": len(cases), "max_solves": contract["max_new_solves"],
                      "surface_station_m": surface_station, "stage_counts": {s: sum(c["stage"] == s for c in cases) for s in ["A", "B", "C"]}}, ensure_ascii=False))


def scene_for(root, c, contract):
    import gprMax
    if gprMax.__version__ != "4.0.1":
        raise RuntimeError("gprMax 4.0.1 is required")
    geom = root / c["geometry"]
    assert sha(geom) == c["geometry_sha256"]
    assert sha(geom.parent / f"{DB}.json") == contract["material_database_sha256"]
    scene = gprMax.Scene()
    scene.add(gprMax.Title(name=c["id"]))
    scene.add(gprMax.Discretisation(p1=(c["grid_m"],) * 3))
    scene.add(gprMax.DomainMode(mode="TM"))
    scene.add(gprMax.Domain(p1=(*c["domain_m"], float("inf"))))
    scene.add(gprMax.TimeStepStabilityFactor(f=contract["stability_factor"]))
    scene.add(gprMax.TimeWindow(time=contract["record_ns"] * 1e-9))
    scene.add(gprMax.PMLFormulation(formulation="HORIPML"))
    scene.add(gprMax.PMLThickness(thickness=tuple(c["pml_cells"])))
    scene.add(gprMax.GeometryObjectsRead(p1=(0, 0, 0), geofile=str(geom), material_database=DB, averaging=c["averaging"]))
    scene.add(gprMax.Waveform(wave_type="impulse", amp=1.0, freq=1.0, id="line9_impulse"))
    scene.add(gprMax.HertzianDipole(p1=(*c["tx_m"][:2], float("inf")), polarisation="z", waveform_id="line9_impulse"))
    scene.add(gprMax.Rx(p1=(*c["rx_m"][:2], float("inf")), id="line9_rx", outputs=["Ez"]))
    p1, p2 = c["surface_coefficient_export_bounds"]
    scene.add(gprMax.GeometryObjectsWrite(p1=tuple(p1), p2=tuple(p2), filename=str(root / "cases" / c["id"] / "surface_yee")))
    return scene


def prepare_mid(a):
    """Bounded intermediate-width follow-up; preserves the original contract."""
    import copy
    parent = a.parent.resolve()
    root = a.root.resolve()
    if root.exists():
        raise FileExistsError(root)
    prior = load(parent / "contract.json")
    observed = load(a.screen_results)
    top = next(s for s in observed["summary"] if s["domain_m"] == [50.0, 45.0])
    narrow = next(s for s in observed["summary"] if s["domain_m"] == [40.0, 45.0])
    assert top["all_station_screens_pass"] and not narrow["all_station_screens_pass"]
    root.mkdir(parents=True)
    geodir = root / "geometries"
    geodir.mkdir()
    shutil.copy2(parent / "geometries" / f"{DB}.json", geodir / f"{DB}.json")
    cases=[]
    for station in prior["anchors_m"]:
        full=next(c for c in prior["cases"] if c["id"]==f"full_x{station:.1f}_y")
        with h5py.File(parent/full["geometry"]) as h:
            data=h["data"][100:1900,:1800,:]
            attrs=dict(h.attrs)
        origin=station-22.5
        geo=geodir/f"mid45x45_x{station:.1f}.h5"
        write_geometry(geo,data,attrs,.025,origin)
        case=copy.deepcopy(full)
        case.update(id=f"mid45x45_x{station:.1f}",stage="B",family="reduced_domain",domain_m=[45.0,45.0],
                    geometry=str(geo.relative_to(root)).replace("\\","/"),geometry_sha256=sha(geo),
                    origin_global_x_elevation_m=[origin,410.0])
        for key in ["tx_m","rx_m"]:
            case[key][0]-=2.5
            assert min(case[key][0]-3,45-3-case[key][0])>=10
            assert 45-3-case[key][1]>=2
        for p in case["surface_coefficient_export_bounds"]:p[0]-=2.5
        assert not np.any(data[:,round(max(case['tx_m'][1],case['rx_m'][1])/.025):,:])
        cases.append(case)
    contract=copy.deepcopy(prior)
    contract.update(id="line9_v31_domain_mid_20261013_r1",parent_contract_sha256=sha(parent/"contract.json"),
                    cases=cases,max_new_solves=7,max_wall_seconds=600,
                    stage_C_condition="not applicable; original 30m stage was stopped",
                    adaptive_extension=dict(reason="40m width changed target-window responses; test the 45m intermediate width",
                                             original_study_solves=37,total_study_solve_limit=44,
                                             frozen_tolerances_changed=False,
                                             assessment="Compare to full50 at all seven anchors and directly to wide70 at three anchors. Preserve the failed original guard result."))
    save(root/"contract.json",contract)
    shutil.copy2(Path(__file__),root/"study.py")
    print(json.dumps(dict(root=str(root),new_solves=7,total_study_max=44,parent_contract_sha256=contract['parent_contract_sha256'])))


def solve_case(a):
    import gprMax
    root = a.root.resolve()
    contract = load(root / "contract.json")
    case = next(c for c in contract["cases"] if c["id"] == a.case)
    dest = root / "cases" / case["id"]
    dest.mkdir(parents=True, exist_ok=True)
    output = dest / "impulse.h5"
    if output.exists():
        raise FileExistsError(output)
    scene = scene_for(root, case, contract)
    gprMax.run(scenes=[scene], n=1, outputfile=str(dest / "impulse"), gpu=[0], gpu_precision="single",
               geometry_only=case["geometry_only"], hide_progress_bars=True)


def audit_output(root, case):
    path = root / "cases" / case["id"] / "impulse.h5"
    with h5py.File(path) as h:
        ds = h["rxs/rx1/Ez"]
        data = ds[:]
        dt = float(ds.attrs["SampleInterval"])
        assert data.dtype == np.dtype("float32") and data.shape == (case["expected_samples"],)
        assert np.isfinite(data).all() and np.max(np.abs(data)) > 0
        assert math.isclose(dt, case["dt_s"], rel_tol=1e-12)
        assert np.allclose(h["rxs/rx1"].attrs["Position"], case["rx_m"], rtol=0, atol=1e-8)
        from gprMax.toolboxes.SFCW.processing import load_source
        source = load_source(path, "/srcs/src1").samples
        assert len(source) == len(data) and np.count_nonzero(source) == 1
        tail_n = math.ceil(.05 * len(data))
        tail_db = float(20 * np.log10(np.max(np.abs(data[-tail_n:])) / np.max(np.abs(data))))
    return dict(raw_sha256=sha(path), dtype=str(data.dtype), samples=len(data), dt_s=dt,
                raw_peak=float(np.max(np.abs(data))), tail_db=tail_db,
                source_sha256=hashlib.sha256(source.tobytes()).hexdigest())


def process_sfcw(raw, output, log, spec):
    """Run the official toolbox with the frozen study frequency definition."""
    if output.exists():
        raise FileExistsError(output)
    assert spec['tail_taper'] == 0 and spec['time_shift_s'] == 0
    cmd = [sys.executable, '-m', 'gprMax.toolboxes.SFCW', 'process', str(raw),
           '--source', '/srcs/src1', '--receiver', '/rxs/rx1', '--component', 'Ez',
           '--f-start', str(spec['f_start_hz']), '--f-stop', str(spec['f_stop_hz']),
           '--steps', str(spec['steps']), '--method', spec['method'],
           '--window', spec['window'], '--zero-pad', str(spec['zero_pad']), '--output', str(output)]
    with log.open('w', encoding='utf-8') as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=120)
    if proc.returncode:
        raise RuntimeError(f'Official SFCW failed: {raw}')
    with h5py.File(output) as h:
        assert h['response'].shape == (spec['steps'],) and np.isfinite(h['response'][:]).all()
        assert np.array_equal(h['frequency'][:], np.linspace(spec['f_start_hz'], spec['f_stop_hz'], spec['steps']))


def run_stage(a):
    import gprMax
    import pycuda.driver as cuda
    root = a.root.resolve()
    contract = load(root / "contract.json")
    assert gprMax.__version__ == "4.0.1"
    cuda.init()
    assert "4090" in cuda.Device(0).name()
    start_file = root / "execution_start.json"
    if not start_file.exists():
        save(start_file, dict(start_unix=time.time(), gpu=cuda.Device(0).name(), gprMax=gprMax.__version__,
                             script_sha256=sha(Path(__file__)), contract_sha256=sha(root / "contract.json")))
    started = load(start_file)["start_unix"]
    selected = [c for c in contract["cases"] if c["stage"] == a.stage]
    if a.stage == "C":
        authorization = load(root / "stage_C_screen.json")
        assert authorization["proceed"] is True
    receipt = []
    recovery_manifest = load(root / a.recovery_manifest) if a.recovery_manifest else None
    if recovery_manifest:
        assert recovery_manifest["corrected_worker_sha256"] == sha(Path(__file__))
    for case in selected:
        if time.time() - started > contract["max_wall_seconds"]:
            raise RuntimeError("Study wall limit reached")
        dest = root / "cases" / case["id"]
        dest.mkdir(parents=True, exist_ok=True)
        if recovery_manifest and (dest / "receipt.json").exists():
            previous = load(dest / "receipt.json")
            assert previous["status"] == "PASS"
            for file in previous["files"]:
                assert sha(root / file["path"]) == file["sha256"]
            receipt.append(previous)
            continue
        recover = (dest / "attempt.json").exists()
        if recover:
            if recovery_manifest:
                repair = next(r for r in recovery_manifest["entries"] if r["case"] == case["id"])
            else:
                if a.recover_case != case["id"] or (dest / "receipt.json").exists():
                    raise FileExistsError(f"Already attempted {case['id']}; no retry")
                repair = load(root / "audit_recovery.json")
            assert repair["case"] == case["id"]
            assert sha(dest / "impulse.h5") == repair["raw_sha256"]
            assert sha(dest / "solver.log") == repair["log_sha256"]
            assert sha(Path(__file__)) == (recovery_manifest or repair)["corrected_worker_sha256"]
        else:
            save(dest / "attempt.json", dict(time_unix=time.time(), case=case["id"], contract_sha256=sha(root / "contract.json")))
        save(root / "progress.json", dict(stage=a.stage, case=case["id"], completed=len(receipt), total=len(selected)))
        log = dest / "solver.log"
        begun = time.monotonic()
        if not recover:
            with log.open("w", encoding="utf-8") as f:
                cmd = [sys.executable, str(Path(__file__).resolve()), "solve-case", "--root", str(root), "--case", case["id"]]
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, timeout=1800)
            save(dest / "solver_exit.json", dict(exit_code=result.returncode, wall_seconds=time.monotonic()-begun,
                                                 raw_sha256=sha(dest/"impulse.h5") if (dest/"impulse.h5").exists() else None))
            if result.returncode:
                save(root / "failure.json", dict(stage=a.stage, case=case["id"], exit_code=result.returncode, solver_log_sha256=sha(log)))
                raise RuntimeError(f"Solver failed: {case['id']}")
        log_text = log.read_text(encoding="utf-8", errors="replace")
        completion_times = re.findall(r"Simulation completed in ([^\n=]+)", log_text)
        assert completion_times or case["geometry_only"]
        fdtd_times = re.findall(r"Time taken: ([^\n]+)", log_text)
        item = dict(case=case["id"], wall_seconds=(duration_seconds(completion_times[-1]) if recover else time.monotonic() - begun),
                    official_build_and_solver_seconds=duration_seconds(completion_times[-1]) if completion_times else time.monotonic() - begun,
                    official_fdtd_seconds=duration_seconds(fdtd_times[-1]) if fdtd_times else 0.0,
                    timing_source="official completion log; wrapper timing unavailable" if recover else "parent wall clock",
                    audit_recovered_without_resolve=recover, log_sha256=sha(log), status="PASS")
        if not case["geometry_only"]:
            item.update(audit_output(root, case))
            raw = dest / "impulse.h5"
            sf_log = dest / "sfcw.log"
            process_sfcw(raw, dest / 'sfcw.h5', sf_log, contract['official_sfcw'])
            item["sfcw_sha256"] = sha(dest / "sfcw.h5")
        item["files"] = [{"path": str(p.relative_to(root)).replace("\\", "/"), "sha256": sha(p), "bytes": p.stat().st_size}
                         for p in sorted(dest.glob("*")) if p.is_file() and p.name != "receipt.json"]
        save(dest / "receipt.json", item)
        receipt.append(item)
        save(root / f"stage_{a.stage}_receipts.json", receipt)
        print(f"{a.stage}: {case['id']} PASS {item['wall_seconds']:.1f}s", flush=True)
    reused = []
    for ref in contract.get('reuse_references', []):
        dest = root / 'references' / ref['id']
        raw = dest / 'impulse.h5'
        assert sha(raw) == ref['raw_sha256']
        assert sha(dest / 'sfcw_20_170.h5') == ref['wide_sfcw_sha256']
        process_sfcw(raw, dest / 'sfcw.h5', dest / 'sfcw.log', contract['official_sfcw'])
        item = dict(case=ref['id'], status='PASS', new_fdtd=False,
                    files=[dict(path=p.relative_to(root).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
                           for p in sorted(dest.glob('*')) if p.is_file() and p.name != 'receipt.json'])
        save(dest / 'receipt.json', item)
        reused.append(item)
    completion = dict(status='COMPLETE', cases=len(receipt), receipts=receipt,
                      elapsed_since_first_stage_s=time.time() - started)
    if reused:
        completion['reused_references'] = reused
    save(root / f"stage_{a.stage}_complete.json", completion)
    save(root / "progress.json", dict(stage=a.stage, status="COMPLETE", completed=len(receipt), total=len(selected)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--root", type=Path, required=True)
    prep.add_argument("--source-package", type=Path, required=True)
    prep.add_argument("--bscan-root", type=Path, required=True)
    prep.add_argument("--mother", type=Path, required=True)
    mid = sub.add_parser("prepare-mid")
    mid.add_argument("--root", type=Path, required=True)
    mid.add_argument("--parent", type=Path, required=True)
    mid.add_argument("--screen-results", type=Path, required=True)
    solve = sub.add_parser("solve-case")
    solve.add_argument("--root", type=Path, required=True)
    solve.add_argument("--case", required=True)
    run = sub.add_parser("run-stage")
    run.add_argument("--root", type=Path, required=True)
    run.add_argument("--stage", choices=["A", "B", "C"], required=True)
    run.add_argument("--recover-case", default=None, help="Audit an explicitly hash-bound completed output without re-solving")
    run.add_argument("--recovery-manifest", default=None, help="Hash-bound audit recovery manifest; never retries FDTD")
    a = p.parse_args()
    {"prepare": prepare, "prepare-mid": prepare_mid, "solve-case": solve_case, "run-stage": run_stage}[a.command](a)


if __name__ == "__main__":
    main()
