"""Test, without deleting data, an alternative cover-reverberation path hypothesis."""
import argparse
import json
from pathlib import Path

import numpy as np

from diagnose_line9_layer_kinematics import primary_paths
from review_line9_result_packages import indices, C0, FREQ, sha, save, rel
from audit_line9_postprocessing import inverse, weights
from trace_line9_time_origin import correlation


def main(root, review, cache, out):
    if out.exists():
        raise ValueError('Fresh hypothesis report directory required')
    out.mkdir(parents=True)
    result = dict(calls_solver=False, calls_training=False, script_sha256=sha(__file__), packages=[])
    for ap in sorted(review.glob('*_audit.json')):
        audit = json.loads(ap.read_text('utf-8'))
        rows = audit['records']
        mf = next((root/audit['package']/'geometries').glob('*.json'))
        if sha(mf) != audit['materials_sha256']:
            raise ValueError('Materials changed')
        materials = json.loads(mf.read_text('utf-8'))['materials']
        n = np.array([indices(materials, f) for f in FREQ])
        spectra = []
        closures = []
        for r in rows:
            g = r['geometry']
            bottom = g['cover_base']
            k = g['boundaries'].index(bottom)
            primary = primary_paths(g, materials)[k]
            a, b = bottom['above'], bottom['below']
            if a != 1:
                raise ValueError('Expected one cover layer')
            r_bottom = (n[:, a]-n[:, b])/(n[:, a]+n[:, b])
            r_top = (n[:, a]-n[:, 0])/(n[:, a]+n[:, 0])
            q = r_bottom*r_top*np.exp(-4j*np.pi*FREQ*n[:, a]*bottom['depth_m']/C0)
            if np.max(abs(q)) >= 1:
                raise ValueError('Nonconvergent internal-bounce series')
            # Down/up through cover one extra time, without repeating the air path.
            spectra.append(primary*q)
            closures.append(rel(sum(primary*q**j for j in range(12)), primary/(1-q)))
        if max(closures) > 1e-10:
            raise ValueError('Geometric-series algebra check failed')
        cp = cache/(audit['package']+'.npz')
        with np.load(cp) as z:
            np.testing.assert_array_equal(z['ids'], [r['id'] for r in rows])
            response = z['response'].copy()
        multiple = np.column_stack(spectra)
        item = dict(package=audit['package'], audit_sha256=sha(ap), cache_sha256=sha(cp),
                    geometric_series_relative_L2_max=max(closures), windows={})
        for name in ['hann', 'blackman']:
            model, t = inverse(multiple, FREQ, weights(name, 501))
            actual, _ = inverse(response, FREQ, weights(name, 501))
            centers = t[np.argmax(abs(model), axis=0)]*1e9
            gate = abs(t[:, None]*1e9-centers[None, :])<=12
            actual_peaks = t[np.argmax(np.where(gate, abs(actual), 0), axis=0)]*1e9
            item['windows'][name] = dict(predicted_peak_ns_endpoints=centers[[0,-1]].tolist(),
                complex_correlation=correlation(model[gate], actual[gate]),
                actual_minus_predicted_peak_ns_percentile0_50_100=np.percentile(actual_peaks-centers, [0,50,100]).tolist())
        result['packages'].append(item)
        print(json.dumps(item, ensure_ascii=False), flush=True)
    result['conclusion'] = ('Cover second-bounce arrivals can approach apparent deep ridges, but complex correlations '
        'are only about0.41-0.49. Timing proximity alone does not certify this as the dominant cause. '
        'No envelope subtraction or demultiple operation applied. This is a local plane-wave hypothesis, '
        'not a matched physical ablation or PML attribution.')
    save(out/'cover_reverberation.json', result)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['root', 'review', 'cache', 'out']:
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.review, a.cache, a.out)
