"""Independent native reread: four-neighbor weights, gates and geometric delays."""
import argparse
import json
from pathlib import Path
import hashlib
import h5py
import numpy as np


def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()


def main(a):
    m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    s=json.loads((a.results/'analysis.json').read_text('utf-8'))
    assert digest(a.raw)==s['native_sha256'] and digest(a.arrays)==s['arrays_sha256']
    assert digest(Path(__file__).with_name('analyze_line9_native_probes.py'))==s['script_sha256']
    errors=[];directions=[];corr_errors=[]
    with h5py.File(a.raw) as h,h5py.File(a.package/'baseline/native_H0.h5') as old,np.load(a.arrays) as derived:
        stored={name:derived[name] for name in ['Ez','Hx','Hy','Sx','Sy']}
        for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
            assert h[key][:].tobytes()==old[key][:].tobytes()
        dt=float(h.attrs['dt']);t=np.arange(20351)*dt;ns=t*1e9;mask=(ns>=300)&(ns<=370)
        rx=h['rxs/rx1/Ez'][:-1];correlations={r['id']:r for r in s['surface_delay_correlations']}
        for k,p in enumerate(m['probes']):
            rr=[h[f"rxs/rx{q['receiver_index']}"] for q in p['anchors']]
            ez=rr[0]['Ez'][:-1]
            # Four values receive equal weight; independent of producer's staged averaging.
            hx=(rr[0]['Hx'][:-1]+rr[0]['Hx'][1:]+rr[1]['Hx'][:-1]+rr[1]['Hx'][1:])*.25
            hy=(rr[0]['Hy'][:-1]+rr[0]['Hy'][1:]+rr[2]['Hy'][:-1]+rr[2]['Hy'][1:])*.25
            for field,x in [('Ez',ez),('Hx',hx),('Hy',hy),('Sx',-ez*hy),('Sy',ez*hx)]:
                y=stored[field][k];err=float(np.linalg.norm(x-y)/max(np.linalg.norm(x),np.finfo(float).tiny))
                assert err<1e-14,(p['id'],field,err);errors.append(err)
            sx=-ez*hy;sy=ez*hx
            row=s['gates'][k];assert row['id']==p['id']
            for gate,bounds in m['gates_native_ns'].items():
                ids=np.flatnonzero((ns>=bounds[0])&(ns<=bounds[1]));z=row['gates'][gate]
                ij=[np.sum(sx[ids],dtype=np.longdouble)*dt,np.sum(sy[ids],dtype=np.longdouble)*dt]
                np.testing.assert_allclose(ij,[z['integrated_Sx_J_m2'],z['integrated_Sy_J_m2']],rtol=1e-10,atol=1e-24)
                angle=float(np.degrees(np.arctan2(*ij[::-1])));er=abs(angle-z['direction_deg'])
                assert er<1e-8;directions.append(er)
                for sign,fun in [('max',max),('min',min)]:
                    j=fun(ids,key=lambda q:sy[q]);assert z[sign+'_Sy_time_ns']==ns[j]
                    np.testing.assert_allclose(z[sign+'_Sy_W_m2'],sy[j],rtol=1e-14,atol=1e-25)
            if p['id'] in correlations:
                r=correlations[p['id']]
                delay=np.sqrt((170.65-p['x_m'])**2+(39.025-p['y_m'])**2)/299792458.
                u=(t[mask]-delay)/dt;j=np.floor(u).astype(int);alpha=u-j
                assert min(j)>=0 and max(j)+1<len(ez)
                shifted=ez[j]*(1-alpha)+ez[j+1]*alpha
                v=rx[mask];corr=float(np.sum(shifted*v)/np.sqrt(np.sum(shifted**2)*np.sum(v**2)))
                er=abs(corr-r['signed_cosine']);assert er<1e-10;corr_errors.append(er)
    out={'status':'PASS_INDEPENDENT287_NATIVE_PROBES_AND41_GEOMETRIC_CORRELATIONS',
         'script_sha256':digest(__file__),'analysis_sha256':digest(a.results/'analysis.json'),
         'field_comparisons':len(errors),'max_field_relative_L2':max(errors),'gate_directions':len(directions),
         'max_direction_error_deg':max(directions),'correlations':len(corr_errors),'max_correlation_error':max(corr_errors),
         'main_receiver_source_bitwise_equal':True,
         'limits':'Independent arithmetic/raw provenance audit, not a separated-path or site certificate.'}
    assert not a.out.exists();a.out.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8');print(json.dumps(out))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['raw','package','results','arrays','out']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
