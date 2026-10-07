"""CPU-only Fermat-path hypothesis on existing Line9 geometry; never calls FDTD."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import minimize

from review_line9_result_packages import FREQ, C0, indices, sha, save
from diagnose_line9_layer_kinematics import primary_paths
from audit_line9_postprocessing import inverse, weights
from trace_line9_time_origin import correlation


def extract_columns(data, dy):
    columns = []
    for col in data:
        cuts = np.flatnonzero(np.diff(col) != 0)+1
        columns.append([(int(col[y]), int(col[y-1]), float(y*dy)) for y in cuts[::-1]])
    return columns


def optimise_path(tx, rx, x, curves, materials, refractive, shifts=(-10., -5., 0., 5., 10.)):
    """Minimise phase optical length through an ordered interface branch.

    Interface heights are linearly interpolated. This real-index ray proxy
    excludes frequency-dependent saddle points, spreading and diffraction.
    """
    k = len(curves)-1
    contact_order = np.r_[np.arange(k+1), np.arange(k-1, -1, -1)]
    media = np.r_[materials, materials[::-1]]
    midpoint = float((tx[0]+rx[0])/2)
    def nodes(v):
        y = np.array([np.interp(xx, x, curves[j]) for xx, j in zip(v, contact_order)])
        return np.vstack([tx[:2], np.column_stack([v, y]), rx[:2]])
    def objective(v):
        lengths = np.linalg.norm(np.diff(nodes(v), axis=0), axis=1)
        return float(np.dot(lengths, refractive[media].real)/C0*1e9)
    trials = []
    for shift in shifts:
        start = np.full(len(contact_order), np.clip(midpoint+shift, x[0], x[-1]))
        opt = minimize(objective, start, method='L-BFGS-B', bounds=[(x[0], x[-1])]*len(start),
                       options={'ftol': 1e-12, 'gtol': 1e-6, 'maxiter': 600, 'maxls': 40})
        trials.append(opt)
    successful = [r for r in trials if r.success and np.isfinite(r.fun)]
    if not successful:
        raise ValueError('No converged candidate path')
    opt = min(successful, key=lambda r: r.fun)
    points = nodes(opt.x)
    per_material = np.bincount(media, weights=np.linalg.norm(np.diff(points, axis=0), axis=1), minlength=len(refractive))
    return dict(phase_time95_ns=float(opt.fun), nodes_xy_m=points.tolist(),
                reflection_offset_from_midpoint_m=float(opt.x[k]-midpoint),
                path_length_by_material_m=per_material.tolist(),
                segment_materials=media.tolist(), successful_starts=len(successful),
                best_minus_worst_converged_time_ns=float(min(r.fun for r in successful)-max(r.fun for r in successful)),
                crossing_x_interval_m=[float(x[0]), float(x[-1])])


def ray_on_geometry(tx, rx, columns, signature, dx, n95, stride, data):
    k = len(signature)-1
    valid = np.array([tuple((a,b) for a,b,y in c[:k+1]) == signature for c in columns])
    mid = round((tx[0]+rx[0])/2/dx)
    if not valid[mid]:
        raise ValueError('Midpoint topology differs from audited column')
    left = right = mid
    while left > 80 and valid[left-1]:
        left -= 1
    while right < len(valid)-81 and valid[right+1]:
        right += 1
    ix = np.unique(np.r_[left, np.arange(left, right+1, stride), right])
    xx = ix*dx
    curves = np.array([[columns[i][j][2] for i in ix] for j in range(k+1)])
    materials = np.array([a for a,b in signature])
    result = optimise_path(tx, rx, xx, curves, materials, n95)
    points = np.array(result['nodes_xy_m'])
    wrong = near = 0.
    # Do not certify rays that cross the wrong material away from a contact.
    # Two native cells is only a declared interpolation ambiguity band.
    for p, q, expected in zip(points[:-1], points[1:], result['segment_materials']):
        length = float(np.linalg.norm(q-p))
        steps = max(2, int(np.ceil(length/(dx/2))))
        u = (np.arange(steps)+.5)/steps
        sample = p[None,:]+u[:,None]*(q-p)[None,:]
        cell = np.floor(sample/dx).astype(int)
        cell = np.clip(cell, [0,0], np.array(data.shape)-1)
        observed = data[cell[:,0], cell[:,1]]
        distance_to_contact = np.min(abs(sample[:,1,None]-np.column_stack([
            np.interp(sample[:,0], xx, y) for y in curves])), axis=1)
        contact = distance_to_contact <= 2*dx
        near += length*np.mean(contact)
        wrong += length*np.mean((observed != expected)&~contact)
    result['wrong_material_length_outside_contact_band_m'] = wrong
    result['near_contact_unresolved_length_m'] = near
    result['geometry_check'] = 'NO_OUTSIDE_BAND_MISMATCH' if wrong == 0 else 'REJECT_WRONG_MATERIAL_SEGMENT'
    result['boundary_curve_sample_spacing_m'] = stride*dx
    return result


def checks():
    from review_line9_result_packages import ray_time
    x = np.linspace(-50., 50., 501)
    tx = np.array([0., 10.]);rx = np.array([1.3, 10.])
    surfaces = np.array([np.zeros(len(x)), np.full(len(x), -7.)])
    r = optimise_path(tx, rx, x, surfaces, np.array([0,1]), np.array([1.,3.]))
    expected = ray_time(np.array([10.,7.]), np.array([1.,3.]),1.3)
    flat_error = abs(r['phase_time95_ns']-expected)
    assert flat_error < 1e-4
    slope = .1
    r = optimise_path(tx, rx, x, np.array([slope*x]), np.array([0]), np.array([1.]))
    normal = np.array([-slope,1.]);normal /= np.linalg.norm(normal)
    reflected = tx-2*np.dot(tx,normal)*normal
    expected = np.linalg.norm(rx-reflected)/C0*1e9
    inclined_error = abs(r['phase_time95_ns']-expected)
    assert inclined_error < 1e-4
    return dict(flat_layer_phase_time_error_ns=flat_error, inclined_mirror_phase_time_error_ns=inclined_error)


def corrected_primary(g, materials, path):
    target = g['basal_sand']
    k = g['boundaries'].index(target)
    local = primary_paths(g, materials)[k]
    thickness = np.zeros(len(materials));thickness[0] = g['midpoint_agl_m']
    previous = 0.
    for b in g['boundaries'][:k+1]:
        thickness[b['above']] += b['depth_m']-previous
        previous = b['depth_m']
    n95 = indices(materials,95e6)
    offset = target['time95_ns']*1e-9-2*np.dot(thickness,n95.real)/C0
    n = np.array([indices(materials,f) for f in FREQ])
    delta = n@(np.array(path['path_length_by_material_m'])-2*thickness)/C0-offset
    return local, local*np.exp(-2j*np.pi*FREQ*delta)


def main(root, review, cache, out, stride):
    if out.exists():
        raise ValueError('Fresh ray diagnostic directory required')
    report = dict(calls_solver=False,calls_training=False,script_sha256=sha(__file__),
                  self_checks=checks(), boundary_curve_native_stride=stride, packages=[])
    out.mkdir(parents=True)
    for ap in sorted(review.glob('*_audit.json')):
        audit = json.loads(ap.read_text('utf-8'));p=root/audit['package'];rows=audit['records']
        mf=next((p/'geometries').glob('*.json'));gf=next((p/'geometries').glob('*.h5'))
        if sha(mf)!=audit['materials_sha256'] or sha(gf)!=audit['geometry_sha256']:
            raise ValueError('Model identity changed')
        materials=json.loads(mf.read_text('utf-8'))['materials'];n95=indices(materials,95e6)
        with h5py.File(gf) as h:
            data=h['data'][:,:,0];dx=float(h.attrs['dx_dy_dz'][0])
        columns=extract_columns(data,dx)
        # Every common station for the height comparison; fixed 9 stations in pkg1.
        selected=np.unique(np.linspace(0,len(rows)-1,9).round().astype(int)) if 'pkg1' in audit['package'] else np.arange(70)
        item=dict(package=audit['package'],audit_sha256=sha(ap),geometry_sha256=sha(gf),records=[])
        local_spectra=[];ray_spectra=[];accepted=[]
        for j in selected:
            row=rows[j];g=row['geometry'];k=g['boundaries'].index(g['basal_sand'])
            signature=tuple((b['above'],b['below']) for b in g['boundaries'][:k+1])
            raw=p/'cases'/row['id']/'profile.h5'
            if sha(raw)!=row['native_sha256']:
                raise ValueError('Native identity changed')
            with h5py.File(raw) as h:
                tx=np.array(h['srcs/src1'].attrs['Position']);rx=np.array(h['rxs/rx1'].attrs['Position'])
            path=ray_on_geometry(tx,rx,columns,signature,dx,n95,stride,data)
            path.update(station_index=int(j),id=row['id'],chainage_m=row['chainage_m'],
                        local_phase_time95_ns=g['basal_sand']['time95_ns'])
            item['records'].append(path)
            if path['geometry_check']=='NO_OUTSIDE_BAND_MISMATCH':
                old,new=corrected_primary(g,materials,path)
                local_spectra.append(old);ray_spectra.append(new);accepted.append(j)
            if len(item['records'])%10==0:
                print(f'{audit["package"]}: {len(item["records"])}/{len(selected)} rays inspected',flush=True)
        cp=cache/(audit['package']+'.npz')
        with np.load(cp) as z:
            response=z['response'][:,accepted].copy()
        item['accepted_station_indices']=[int(i) for i in accepted]
        item['private_cache_sha256']=sha(cp)
        item['comparison']={}
        if accepted:
            for window in ['hann','blackman']:
                actual,t=inverse(response,FREQ,weights(window,501))
                for name,spectra in [('local',local_spectra),('ray_proxy',ray_spectra)]:
                    model,_=inverse(np.column_stack(spectra),FREQ,weights(window,501))
                    centers=t[np.argmax(abs(model),axis=0)]*1e9
                    gate=abs(t[:,None]*1e9-centers[None,:])<=12
                    item['comparison'][window+'_'+name]=dict(
                        complex_correlation=correlation(model[gate],actual[gate]),
                        peak_ns=centers.tolist())
        report['packages'].append(item)
        print(json.dumps(dict(package=audit['package'],accepted=len(accepted),comparison=item['comparison']),ensure_ascii=False),flush=True)
    report['limits']=('Fermat phase-ray approximation on an interpolated, same-topology interface branch. '
        'Only the best of five converged starts, no proof of global minimum or unique physical path. '
        'No diffraction, geometric spreading, oblique Fresnel correction, full frequency-dependent ray bending, PML or numerical precision certificate. '
        'Wrong-material segments outside two native cells of contacts are rejected; contact-band geometry remains unresolved. '
        'No processing of actual data with the model; this is not an FDTD or migration result.')
    save(out/'refracted_paths.json',report)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['root','review','cache','out']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--stride',type=int,choices=[4,8],default=8)
    a=p.parse_args();main(a.root,a.review,a.cache,a.out,a.stride)
