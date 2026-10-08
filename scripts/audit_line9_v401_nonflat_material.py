"""Independent nonflat single-factor audit: native, voxels, materials, DFT, inverse."""
import argparse
import copy
import json
from pathlib import Path

import h5py
import numpy as np
from audit_line9_v401_delivery import native_response, read, sha


def main(a):
    assert not a.out.exists()
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    m=read(a.package/'manifest.json');pub=read(a.public/'analysis.json')
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')==pub['contract_sha256']
    assert sha(a.source/'sfcw.h5')==pub['numerical_sha256']
    group={g['id']:g for g in m['groups']}
    records={g['id']:g for g in v['groups']}
    arrays=[];errors=[]
    for name in ['baseline_H0','baseline_H1','low_H0','low_H1']:
        if name=='baseline_H0':
            g=m['reused_H0']['group'];digest=m['reused_H0']['native_sha256'];p=a.package/'reused_nonflat_H0.h5'
        else:
            g=group[name];digest=records[name]['native_sha256'];p=a.source/(name+'.h5')
            for key in ['input','material','geometry']:assert sha(a.package/g[key])==g[key+'_sha256']
        q=native_response(p,.025,digest);arrays.append(q)
        with h5py.File(p) as h:
            assert h.attrs['gprMax']=='4.0.1'
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[8400,1700,1])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            assert h.attrs['Iterations']==20352
            np.testing.assert_allclose(h.attrs['dt'],g['dt_s'],rtol=1e-14,atol=0)
            for key,pos in [('srcs/src1',g['tx_m']),('rxs/rx1',g['rx_m'])]:
                np.testing.assert_array_equal(h[key].attrs['GridPosition'],np.rint(np.array(pos)/.025))
                np.testing.assert_allclose(h[key].attrs['Position'],pos,rtol=0,atol=1e-12)
            source=h['srcs/src1/excitation/samples'][:]
            if name=='baseline_H0':reference_source=source
            else:np.testing.assert_array_equal(source,reference_source)
    z=np.column_stack(arrays);z=np.column_stack([z,z[:,1]-z[:,0],z[:,3]-z[:,2]])
    f=20e6+np.arange(501)*300000.
    with h5py.File(a.source/'sfcw.h5') as h:
        np.testing.assert_array_equal(h['frequency_Hz'][:],f)
        stored=h['response'][:]
        for j in range(6):
            error=float(np.linalg.norm(z[:,j]-stored[:,j])/np.linalg.norm(z[:,j]));assert error<1e-9;errors.append(error)
    assert group['baseline_H1']['geometry_sha256']==group['low_H1']['geometry_sha256']
    assert group['low_H0']['geometry_sha256']==m['reused_H0']['group']['geometry_sha256']
    baseline=read(a.package/group['baseline_H1']['material'])
    for name in ['low_H0','low_H1']:
        changed=read(a.package/group[name]['material']);restored=copy.deepcopy(changed)
        original=baseline['materials']['material_002_mudstone']
        assert original['base']['electric_conductivity_s_per_m']==.003
        assert changed['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m']==.0003
        restored['database']=copy.deepcopy(baseline['database'])
        restored['materials']['material_002_mudstone']['metadata']=copy.deepcopy(original['metadata'])
        restored['materials']['material_002_mudstone']['base']['electric_conductivity_s_per_m']=.003
        assert restored==baseline
    physical=[]
    for g in m['groups']:
        physical.append([line for line in (a.package/g['input']).read_text('utf-8').splitlines() if not line.startswith('#title:')])
    assert physical[0]==physical[1]==physical[2]
    # Independently verify H0 differs only in the bottom-connected sandstone.
    changed_count=0
    with h5py.File(a.package/group['low_H0']['geometry']) as h0,h5py.File(a.package/group['baseline_H1']['geometry']) as h1:
        before_all=h1['data'][:,:,0];after_all=h0['data'][:,:,0]
        for j in range(8400):
            before=before_all[j];after=after_all[j]
            depth=0
            while depth<len(before) and before[depth]==3:depth+=1
            expected=before.copy();expected[:depth]=2
            np.testing.assert_array_equal(expected,after);changed_count+=depth
    assert changed_count>0
    t=np.arange(4008)/(4008*300000.);metrics={}
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        metrics[window]={};w=w/w.mean()
        for gate in ['basal','deep']:
            lo,hi=m[gate+'_gate_ns'];mask=(t*1e9>=lo)&(t*1e9<=hi)
            x=np.exp(2j*np.pi*t[mask,None]*f)@(z*w[:,None])/501
            norm=lambda j:np.linalg.norm(x[:,j])
            corr=lambda j,k:float(abs(np.vdot(x[:,j],x[:,k]))/(norm(j)*norm(k)))
            values=dict(baseline_H0_over_delta=float(norm(0)/norm(4)),changed_H0_over_delta=float(norm(2)/norm(5)),
                        changed_delta_over_baseline_delta=float(norm(5)/norm(4)),changed_H0_over_baseline_H0=float(norm(2)/norm(0)),
                        changed_H1_vs_delta_correlation=corr(3,5),baseline_H1_vs_delta_correlation=corr(1,4))
            for key,value in values.items():assert abs(value-pub['metrics'][window][gate][key])<1e-9
            metrics[window][gate]=values
    result=dict(status='PASS_INDEPENDENT_NONFLAT_NATIVE_VOXELS_MATERIALS_DFT_INVERSE',
                script_sha256=sha(__file__),analysis_sha256=sha(a.public/'analysis.json'),
                contract_sha256=sha(a.source/'execution_contract.json'),
                independent_DFT_relative_L2=errors,bottom_connected_voxels_changed=changed_count,
                direct_inverse_metrics=metrics,physical_change_only='mudstone sigmaDC .003 -> .0003 S/m',
                limits='One nonflat station, no field calibration, entire-line or 3D certification.')
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','package','public','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
