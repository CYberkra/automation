"""Array/HDF5 contract checks only; never run FDTD or claim a field result."""
import argparse
import json
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf

from build_pdf_profile_geometry import digest, save_json
import run_line9_2d_first as run
import analyze_line9_2d_sfcw as analysis


def check(out):
    if out.exists(): raise ValueError('Fresh test evidence path required')
    out.mkdir(parents=True)
    task=json.loads(run.TASK.read_text('utf-8'))
    ids=task['acquisition']['preview_ids'];positions=np.array([int(n.rsplit('s',1)[1])*.5 for n in ids])
    np.testing.assert_array_equal(positions[:-1],np.arange(0,195,2))
    assert positions[-1]==195 and len(ids)==99 and len(set(ids))==99
    assert len(task['acquisition']['preview_pending_ids'])==97
    assert set(task['acquisition']['reuse_pilot_output'])=={'full2d_s0000','full2d_s0390'}
    np.testing.assert_array_equal(analysis.FREQUENCY,20e6+.3e6*np.arange(501))
    dt=.025/(299792458*np.sqrt(2));iterations=round(1.2e-6/dt)+1
    source_times=(np.arange(iterations)+.5)*dt
    a=(np.pi*100e6*(source_times-14.14e-9))**2
    source=40*(1-2*a)*np.exp(-a)
    shift=1697;receiver=np.zeros(iterations);receiver[shift:]=.4*source[:-shift]
    raw=out/'fixture.h5'
    with h5py.File(raw,'x') as h:
        h.attrs.update(gprMax='4.0.0',dt=dt,Iterations=iterations,dx_dy_dz=[.025]*3,nx_ny_nz=[16000,3000,1])
        sg=h.create_group('srcs/src1');sg.attrs['GridPosition']=[100,100,0]
        e=sg.create_group('excitation');e.create_dataset('samples',data=source)
        e.attrs.update(SampleInterval=dt,TimeSampleOffset=.5*dt,SourceType='HertzianDipole',
            SpatialScale=.025,Polarisation='z',WaveformAmplitude=40.,WaveformFrequency=100e6)
        rx=h.create_group('rxs/rx1');rx.attrs['GridPosition']=[152,100,0];rx.create_dataset('Ez',data=receiver)
    g=dict(id='fixture',input=str(raw.with_suffix('.in').resolve()),input_sha256='',profile_x_m=220.,
           acquisition_s_m=0.,tx_m=[2.5,2.5,.0125],rx_m=[3.8,2.5,.0125])
    raw.with_suffix('.in').write_text('Array fixture; not a gprMax input',encoding='utf-8')
    g['input_sha256']=digest(raw.with_suffix('.in'))
    contract=out/'execution_contract.json';save_json(contract,dict(task_sha256='fixture',stage='pilots',groups=[g]))
    with patch.object(run,'check_files'):
        assert run.audit(contract,True)['status']=='PASS'
        for kind in ['float32','bad_time','bad_polarization']:
            with h5py.File(raw,'r+') as h:
                if kind=='float32':
                    del h['rxs/rx1/Ez'];h['rxs/rx1'].create_dataset('Ez',data=receiver.astype('float32'))
                elif kind=='bad_time': h.attrs['dt']=dt*2
                else: h['srcs/src1/excitation'].attrs['Polarisation']='y'
            try: run.audit(contract,True)
            except ValueError: pass
            else: raise AssertionError(kind+' should be rejected')
            with h5py.File(raw,'r+') as h:
                if kind=='float32':
                    del h['rxs/rx1/Ez'];h['rxs/rx1'].create_dataset('Ez',data=receiver)
                elif kind=='bad_time': h.attrs['dt']=dt
                else: h['srcs/src1/excitation'].attrs['Polarisation']='z'
    product,_,record=analysis.process_trace(raw)
    expected=.4/.025*np.exp(-2j*np.pi*analysis.FREQUENCY*((shift-.5)*dt))
    phase_error=analysis.relative(product.response,expected)
    assert phase_error<1e-9
    inverse={}
    for window in ['rectangular','hann']:
        result=sf.reconstruct_time_response(product,window=window,zero_pad_factor=8,time_shift=0)
        take=np.arange(0,len(result.time),53)
        manual=np.exp(2j*np.pi*result.time[take,None]*analysis.FREQUENCY[None,:])@(result.weights*expected)/501
        inverse[window]=analysis.relative(result.complex_bandpass[take],manual)
        assert inverse[window]<1e-9
        np.testing.assert_array_equal(result.real_bandpass,2*result.complex_bandpass.real)
    save_json(out/'checks.json',dict(status='PASS_ARRAY_HDF5_CONTRACT_NOT_FDTD',calls_solver=False,
        checks=['99coarse endpoints/order/nonuniform final interval','97remaining +2reused',
                '501 exact tones retained','float32 rejected before conversion','wrong dt rejected',
                'wrong polarization rejected','half-step origin and current-element normalization',
                'two complex inverse sums','signed bandpass convention'],
        spectral_phase_relative_L2=phase_error,inverse_relative_L2=inverse,
        independent_DFT_relative_L2=record['independent_DFT_relative_L2'],
        code_sha256={n:digest(ROOT/'scripts'/n) for n in ['run_line9_2d_first.py','analyze_line9_2d_sfcw.py','check_line9_2d_first_run.py']}))
    print(f'Array/HDF5 checks PASS; phase L2={phase_error:.3g}; no FDTD')


ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    check(p.parse_args().out.resolve())
