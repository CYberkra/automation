"""Prepare authorized practical coarse-grid comparisons; no solver here."""
from __future__ import annotations

import argparse
import copy
import math
from pathlib import Path
import shutil

import h5py
import numpy as np

from study_line9_v31_smoothing_domain import DB, load, save, sha, write_geometry


def raster(contacts, selection, station, spacing):
    n = math.ceil(50.0 / spacing - 1e-10)
    width = n * spacing
    origin = station - width / 2
    xs = origin + (np.arange(n) + .5) * spacing
    ys = 410.0 + (np.arange(n) + .5) * spacing
    lines = {k: np.interp(xs, np.asarray(v)[:, 0], np.asarray(v)[:, 1])
             for k, v in contacts.items()}
    surface, clay, upper, lower_top, lower_bottom = [lines[k] for k in
        ['surface', 'clay_base', 'upper_mudstone_base', 'lower_mudstone_top', 'lower_mudstone_base']]
    lo, hi = selection['upper_mudstone_extent_m']
    upper_on = (xs >= lo) & (xs < hi)
    y = ys[None, :]
    data = np.full((n, n), 3, dtype=np.int16)
    data[(y >= lower_bottom[:, None]) & (y < np.minimum(lower_top, clay)[:, None])] = 2
    data[(y >= upper[:, None]) & (y < clay[:, None]) & upper_on[:, None]] = 2
    data[(y >= clay[:, None]) & (y < surface[:, None])] = 1
    data[y >= surface[:, None]] = 0
    # Independent nested classifier prevents sequential-assignment mistakes.
    ref = np.where(y >= surface[:, None], 0,
          np.where(y >= clay[:, None], 1,
          np.where(upper_on[:, None] & (y >= upper[:, None]) & (y < clay[:, None]), 2,
          np.where((y >= lower_bottom[:, None]) & (y < np.minimum(lower_top, clay)[:, None]), 2, 3))))
    assert np.array_equal(data, ref)
    return data[:, :, None], origin, width


def prepare(args):
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    geodir = root / 'geometries'
    geodir.mkdir()
    package = args.source_package.resolve()
    parent = args.parent.resolve()
    inputs = args.parent_inputs.resolve()
    prior = load(parent / 'contract.json')
    cs = load(package / 'model/contacts_m.json')
    selpath = package / 'config/selection_pdf_reconstructed.json'
    sel = load(selpath)
    assert cs['selection_sha256'] == sha(selpath)
    assert cs['geometry_model_version'] == 'V3.1_VERTICAL_PDF_SCALE_FIXED'
    contacts = cs['contact_points_x_elevation_m']
    thirtieth = getattr(args, 'spacing_fraction', None) == '1/30'
    spacings = [1.0 / 30.0] if thirtieth else [.0325, .0375]
    db = inputs / 'geometries' / f'{DB}.json'
    assert sha(db) == prior['material_database_sha256']
    shutil.copy2(db, geodir / db.name)
    # Private source contacts are needed only for preparation, never sent to Git.
    shutil.copy2(package / 'model/contacts_m.json', root / 'contacts_source.json')
    cases, checks = [], []
    for station in [150., 162., 174.]:
        base = next(c for c in prior['cases'] if c['id'] == f'full_x{station:.1f}_y')
        geo = inputs / base['geometry']
        assert sha(geo) == base['geometry_sha256']
        with h5py.File(geo) as h:
            rebuilt, ox, width = raster(contacts, sel, station, .025)
            assert np.array_equal(rebuilt, h['data'][:])
            attrs = dict(h.attrs)
        rec = load(parent / 'cases' / base['id'] / 'receipt.json')
        assert sha(parent / 'cases' / base['id'] / 'impulse.h5') == rec['raw_sha256']
        assert sha(parent / 'cases' / base['id'] / 'sfcw.h5') == rec['sfcw_sha256']
        checks.append(dict(station_m=station,regenerated_2p5cm_voxels_identical=True,
                           reference_case=base['id'],reference_raw_sha256=rec['raw_sha256'],
                           reference_sfcw_sha256=rec['sfcw_sha256']))
        for spacing in spacings:
            # An odd 39-cell TX/RX separation needs a half-cell origin shift
            # to keep both absolute x positions on nodes in a 1500-cell domain.
            raster_station = station - spacing/2 if thirtieth else station
            data, origin, width = raster(contacts, sel, raster_station, spacing)
            case = copy.deepcopy(base)
            cid = f'grid1over30m_x{station:.1f}' if thirtieth else f'grid{spacing*100:g}cm_x{station:.1f}'
            out = geodir / f'{cid}.h5'
            write_geometry(out, data, attrs, spacing, origin)
            pml = math.ceil(3. / spacing - 1e-10)
            case.update(id=cid,stage='A',family='practical_coarse_grid',grid_m=spacing,
                        geometry=out.relative_to(root).as_posix(),geometry_sha256=sha(out),
                        domain_m=[width,width],origin_global_x_elevation_m=[origin,410.],
                        pml_cells=[pml,pml,0,pml,pml,0],pml_physical_m=pml*spacing,
                        averaging='y',baseline_case=base['id'])
            case['dt_s'] = .95 * spacing / (299792458. * math.sqrt(2))
            case['expected_samples'] = math.ceil(1100e-9 / case['dt_s']) + 1
            for kind in ['tx','rx']:
                intended = np.array(base[f'{kind}_global_x_elevation_m'])
                local = intended - [origin,410.]
                # Explicit positive tie breaking; never fit a translation to results.
                snapped = np.floor(local / spacing + .5 + 1e-9) * spacing
                global_xy = snapped + [origin,410.]
                case[f'{kind}_intended_global_x_elevation_m'] = intended.tolist()
                case[f'{kind}_global_x_elevation_m'] = global_xy.tolist()
                case[f'{kind}_coordinate_delta_m'] = (global_xy-intended).tolist()
                case[f'{kind}_m'] = [*snapped.tolist(),0.]
                assert np.max(np.abs(global_xy-intended)) <= spacing/2 + 1e-8
                assert min(snapped[0]-pml*spacing,width-pml*spacing-snapped[0]) >= 10
                assert width-pml*spacing-snapped[1] >= 2
            case['actual_horizontal_separation_m'] = case['rx_m'][0]-case['tx_m'][0]
            if thirtieth:
                case['grid_fraction_m'] = [1, 30]
                case['horizontal_domain_shift_m'] = origin-(station-25.)
                assert data.shape == (1500,1500,1) and width == 50. and pml == 90
                assert math.isclose(pml*spacing,3.,abs_tol=1e-14)
                assert all(abs(case[f'{kind}_coordinate_delta_m'][0]) < 1e-10 for kind in ['tx','rx'])
                assert math.isclose(case['actual_horizontal_separation_m'],1.3,abs_tol=1e-10)
            # Small coefficient export around the same physical surface neighborhood.
            bounds = base['surface_coefficient_export_bounds']
            oldorigin = np.asarray(base['origin_global_x_elevation_m'])
            lower = np.floor((np.array(bounds[0][:2])+oldorigin-[origin,410.])/spacing)*spacing
            upper = np.ceil((np.array(bounds[1][:2])+oldorigin-[origin,410.])/spacing)*spacing
            case['surface_coefficient_export_bounds'] = [[*lower.tolist(),0.],[*upper.tolist(),spacing]]
            centercol = data[data.shape[0]//2,:,0]
            e = np.r_[0,np.flatnonzero(np.diff(centercol))+1,len(centercol)]
            case['center_column_intervals'] = [dict(material_id=int(centercol[a]),
                thickness_m=float((b-a)*spacing),bottom_elevation_m=float(410+a*spacing),
                top_elevation_m=float(410+b*spacing)) for a,b in zip(e[:-1],e[1:])]
            cases.append(case)
    contract = copy.deepcopy(prior)
    contract.update(id='line9_v31_coarse_grid_20261013_r1',schema='line9-v31-coarse-grid/1',
        authorization='User approved 开始验证 after the 2.5 vs 3.25–3.75cm small-station proposal',
        anchors_m=[150.,162.,174.],cases=cases,max_new_solves=6,max_wall_seconds=1200,
        parent_contract_sha256=sha(parent/'contract.json'),prepare_script_sha256=sha(Path(__file__)),
        source_contacts_sha256=sha(package/'model/contacts_m.json'),source_selection_sha256=sha(selpath),
        reference_checks=checks,stage_C_condition='not applicable',
        comparison_scope='Practical replacement including declared coordinate/domain/PML/voxel quantization; not isolated bulk-grid dispersion',
        quantization_policy='ceil domain to >=50m; domain centered at nominal station; bottom410m; nearest source points with positive ties; ceil PML to >=3m; no response-fitted correction',
        acceptance='Use original frozen domain_screen tolerances as project replacement screen; reference itself is not absolutely converged',
        current_model_modified=False,source_or_sfcw_changed=False)
    if thirtieth:
        contract.update(id='line9_v31_grid_1over30_20261013_r1',
            authorization='User explicitly requested 网格尝试1/30 m，即 3.333… cm。',
            grid_fraction_m=[1,30],max_new_solves=3,max_wall_seconds=600,
            quantization_policy='spacing=1.0/30.0 without decimal truncation; exact50m/90-cell3mPML/39-cell1.3m separation; x origin shifts -1/60m to preserve TX/RX x; bottom410m fixed; nearest y positive ties; no fitted correction')
    save(root/'contract.json',contract)
    shutil.copy2(Path(__file__).with_name('study_line9_v31_smoothing_domain.py'),root/'study.py')
    save(root/'preparation_audit.json',dict(status='PASS',case_count=len(cases),reference_checks=checks,
         contract_sha256=sha(root/'contract.json'),worker_sha256=sha(root/'study.py')))
    print(f'Prepared {len(cases)} cases; regenerated 2.5cm geometry exactly matches all three references')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--source-package',type=Path,required=True)
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--parent-inputs',type=Path,required=True)
    p.add_argument('--spacing-fraction',choices=['1/30'],default=None)
    prepare(p.parse_args())
