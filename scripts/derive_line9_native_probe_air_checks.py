"""Exploratory native surface upward peak timing and flux direction; no fit."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    s=json.loads((a.results/'analysis.json').read_text('utf-8'))
    assert sha(a.raw)==s['native_sha256']
    rows=[]
    with h5py.File(a.raw) as h:
        dt=float(h.attrs['dt']);t=np.arange(20351)*dt*1e9;rx=h['rxs/rx1/Ez'][:-1]
        ix=np.flatnonzero((t>=300)&(t<=370));q=ix[np.argmax(abs(rx[ix]))];rx_peak=float(t[q])
        for p in m['probes']:
            if p['band']!='air_surface':continue
            r=[h[f"rxs/rx{v['receiver_index']}"] for v in p['anchors']]
            ez=r[0]['Ez'][:-1]
            hx=(r[0]['Hx'][:-1]+r[0]['Hx'][1:]+r[1]['Hx'][:-1]+r[1]['Hx'][1:])*.25
            hy=(r[0]['Hy'][:-1]+r[0]['Hy'][1:]+r[2]['Hy'][:-1]+r[2]['Hy'][1:])*.25
            sx=-ez*hy;sy=ez*hx
            ix=np.flatnonzero((t>=280)&(t<=350));j=ix[np.argmax(sy[ix])]
            gate=(t>=t[j]-4)&(t<=t[j]+4);angle=float(np.degrees(np.arctan2(np.sum(sy[gate]),np.sum(sx[gate]))))
            delay=float(np.hypot(170.65-p['x_m'],39.025-p['y_m'])/299792458*1e9)
            row=next(r for r in s['surface_delay_correlations'] if r['id']==p['id'])
            rows.append({'id':p['id'],'x_m':p['x_m'],'y_m':p['y_m'],'max_upward_flux_time_ns':float(t[j]),
                         'geometric_air_delay_ns':delay,'sum_arrival_ns':float(t[j]+delay),
                         'local8ns_flux_direction_deg':angle,
                         'geometric_direction_to_rx_deg':float(np.degrees(np.arctan2(39.025-p['y_m'],170.65-p['x_m']))),
                         'signed_waveform_cosine':row['signed_cosine']})
    out={'status':'EXPLORATORY_FULL41_SURFACE_TIMING_DIRECTION_CHECK','script_sha256':sha(__file__),
         'parent_analysis_sha256':sha(a.results/'analysis.json'),'native_sha256':s['native_sha256'],
         'rx_late_abs_Ez_peak_ns':rx_peak,'search_upward_flux_ns':[280,350],'direction_gate_relative_ns':[-4,4],
         'rows':rows,'limits':'Post-observation extrema and gates; flux and field peaks are different operators. Agreement of approximate air timing/direction does not certify one exit point or isolated ray.'}
    assert not a.out.exists();a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'rx_late_abs_Ez_peak_ns':rx_peak,'selected_examples':[r for r in rows if r['x_m'] in [160,162,164]]}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['package','results','raw','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
