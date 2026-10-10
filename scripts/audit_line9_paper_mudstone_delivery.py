"""Recheck native bytes, retrieved files, official CLI and derived summaries."""
import argparse
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from line9_v401_version_controls import save


def main(a):
    manifest=json.loads((a.raw/'retrieval_manifest.json').read_text('utf-8-sig'))
    for r in manifest:
        p=a.raw/r['file']; assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes']
    c=json.loads((a.raw/'execution_contract.json').read_text('utf-8'))
    v=json.loads((a.raw/'completed_verification.json').read_text('utf-8'))
    analysis=json.loads((a.public/'analysis.json').read_text('utf-8'))
    prepared={g['id']:g for g in c['study_manifest']['groups']}
    assert analysis['contract_sha256']==v['contract_sha256']==sha(a.raw/'execution_contract.json')
    assert len(c['groups'])==len(v['groups'])==2
    with h5py.File(a.cli) as cli,h5py.File(a.numerical) as h:
        np.testing.assert_array_equal(cli['frequency'][:],h['frequency_Hz'][:])
        np.testing.assert_allclose(cli['response'][:]/.025,h['config2/response'][:,1],rtol=1e-12,atol=0)
        x=cli['time_response/complex_bandpass'][:]/.025; y=h['config2/hann/complex_bandpass'][:,1]
        cli_error=float(np.linalg.norm(x-y)/np.linalg.norm(y)); assert cli_error<1e-12
        np.testing.assert_allclose(cli['time_response/real_bandpass'][:]/.025,
                                   h['config2/hann/real_bandpass'][:,1],rtol=1e-12,atol=1e-12)
        t=h['time_s'][:]*1e9; mask=(t>=250)&(t<=550)
        summary={}
        for w in ['hann','blackman']:
            old=h[f'config0/{w}/complex_bandpass'][:][mask]
            new=h[f'config2/{w}/complex_bandpass'][:][mask]
            op=np.max(abs(old),axis=0); npk=np.max(abs(new),axis=0)
            summary[w]=dict(old_peak_abs=op.tolist(),new_peak_abs=npk.tolist(),
                            delta_peak_gain=float(npk[2]/op[2]),
                            delta_peak_gain_dB=float(20*np.log10(npk[2]/op[2])),
                            H0_peak_gain=float(npk[0]/op[0]),
                            old_delta_over_H0_peak=float(op[2]/op[0]),
                            new_delta_over_H0_peak=float(npk[2]/npk[0]),
                            delta_relative_to_H0_improvement=float((npk[2]/npk[0])/(op[2]/op[0])),
                            delta_relative_to_H0_improvement_dB=float(20*np.log10((npk[2]/npk[0])/(op[2]/op[0]))))
    for g,r in zip(c['groups'],v['groups']):
        p=a.raw/(g['id']+'.h5'); assert sha(p)==r['native_sha256']
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1' and h.attrs['Iterations']==40703
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1])
            assert h['rxs/rx1/Ez'].dtype==np.float64 and np.isfinite(h['rxs/rx1/Ez'][:]).all()
        for key in ['input','geometry','material']:
            assert sha(a.package/prepared[g['id']][key])==g[key+'_sha256']
    # Publish certificates and text logs only; private native/model H5 stay local.
    for r in manifest:
        src=a.raw/r['file']
        if src.suffix=='.h5':continue
        dst=a.public/r['file']; dst.parent.mkdir(parents=True,exist_ok=True)
        assert not dst.exists(); shutil.copyfile(src,dst)
    for name in ['manifest.json']:
        shutil.copyfile(a.package/name,a.public/'preparation_manifest.json')
    result=dict(status='PASS_TWO_NATIVE_AND_OFFICIAL_CLI_DELIVERY_NOT_FIELD_VALIDATION',
        audit_script_sha256=sha(__file__),retrieved_file_count=len(manifest),raw_precision='float64',
        official_CLI_complex_profile_relative_L2=cli_error,summary=summary,
        archive_sha256=sha(Path(str(a.raw)+'.zip')),numerical_sha256=sha(a.numerical),
        completed_FDTD_runs=2,failed_before_time_stepping=1,
        moving_Bscan=False,no_new_solves_in_audit=True,
        limitations='Delta and H0 compare known model contrasts, not SNR or field calibration; no moving-line or spatial convergence certification.')
    save(a.public/'delivery_verification.json',result)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['raw','public','package','cli','numerical']:p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
