"""Raw native DFT plus independent padded-IFFT audit of common-window metrics."""
import argparse,json
from pathlib import Path
import numpy as np
from audit_line9_v401_delivery import native_response,sha,read

def main(a):
    assert not a.out.exists();m=read(a.package/'manifest.json');p=read(a.public/'analysis.json');r=read(a.comparison/'common_gate_comparison.json')
    assert r['analysis_sha256']==sha(a.public/'analysis.json') and r['independent_audit_sha256']==sha(a.public/'independent_audit.json')
    identities={row['id']:row['native_sha256'] for row in p['native']};t=np.arange(4008)/(4008*300000.);f0=20e6;checks={};spectra={}
    for s in m['stations']:
        x=s['chainage_m'];wanted=[f'{v}_x{round(x*100):05d}_{role}' for v in ['high','low'] for role in ['H0','H1']]
        if all(sid in identities for sid in wanted):
            for sid in wanted:spectra[sid]=native_response(a.source/(sid+'.h5'),.025,identities[sid])
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w=w/w.mean();rows={q['chainage_m']:q for q in r['metrics'][window]};worst=0.;count=0;expected=[]
        for s in m['stations']:
            x=s['chainage_m'];prefixes={v:f'{v}_x{round(x*100):05d}_' for v in ['high','low']}
            if not all(prefixes[v]+role in spectra for v in prefixes for role in ['H0','H1']):continue
            expected.append(x);lo=min(s['templates'][v]['basal_gate_ns'][0] for v in prefixes);hi=max(s['templates'][v]['basal_gate_ns'][1] for v in prefixes);assert rows[x]['common_gate_ns']==[lo,hi];k=(t*1e9>=lo)&(t*1e9<=hi);values={}
            for v,prefix in prefixes.items():
                waves={}
                for role in ['H0','H1']:
                    padded=np.zeros(4008,complex);padded[:501]=spectra[prefix+role]*w
                    waves[role]=(np.fft.ifft(padded)*8*np.exp(2j*np.pi*f0*t))[k]
                b=waves['H0'];q=waves['H1'];d=q-b;norm=np.linalg.norm
                values[v]=dict(H0_norm=float(norm(b)),delta_norm=float(norm(d)),H1_norm=float(norm(q)),H0_over_delta=float(norm(b)/norm(d)),H1_vs_delta_correlation=float(abs(np.vdot(q,d))/(norm(q)*norm(d))))
                for key,val in values[v].items():
                    recorded=rows[x]['variants'][v][key];err=abs(val-recorded)/max(1,abs(recorded));assert err<1e-8;worst=max(worst,err);count+=1
            for key,val in [('low_over_high_delta_norm',values['low']['delta_norm']/values['high']['delta_norm']),('low_over_high_H0_norm',values['low']['H0_norm']/values['high']['H0_norm'])]:
                err=abs(val-rows[x][key])/max(1,abs(rows[x][key]));assert err<1e-8;worst=max(worst,err);count+=1
        assert set(rows)==set(expected);checks[window]=dict(scalar_checks=count,max_scaled_error=worst)
    result=dict(status='PASS_RAW_NATIVE_DFT_PADDED_IFFT_COMMON_GATES',auditor_sha256=sha(__file__),comparison_sha256=sha(a.comparison/'common_gate_comparison.json'),analysis_sha256=sha(a.public/'analysis.json'),checks=checks,limits='Independent arithmetic/source audit,not physical error budget or validation of the local template.')
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','public','comparison','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
