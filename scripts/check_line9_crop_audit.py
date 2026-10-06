"""Synthetic native crop audit rejection checks; no FDTD and no field data."""
import argparse
import tempfile
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
from build_pdf_profile_geometry import digest, save_json
import line9_crop_controls as crop


def check(out):
    if out.exists():
        raise ValueError('Fresh check evidence required')
    checks=[]
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp);inp=root/'profile.in';inp.write_text('#title: array fixture\n',encoding='utf-8')
        ref=root/'reference.h5';raw=root/'profile.h5'
        for file,nx in [(ref,16000),(raw,8000)]:
            with h5py.File(file,'w') as h:
                h.attrs['gprMax']='4.0.0';h.attrs['dx_dy_dz']=[.025]*3
                h.attrs['nx_ny_nz']=[nx,3000,1];h.attrs['dt']=1e-9;h.attrs['Iterations']=4
                h['rxs/rx1/Ez']=np.array([0.,1.,-.3,.1])
                h['srcs/src1/excitation/samples']=np.array([0.,1.,0.,0.])
                h['srcs/src1'].attrs['GridPosition']=[4000,2000,1]
                h['rxs/rx1'].attrs['GridPosition']=[4052,2000,1]
                for k,v in dict(SpatialScale=.025,Polarisation='z',WaveformAmplitude=40.,WaveformFrequency=100e6).items():
                    h['srcs/src1/excitation'].attrs[k]=v
        contract=root/'execution_contract.json'
        save_json(contract,dict(reference_raw=str(ref),reference_raw_sha256=digest(ref),
            groups=[dict(id='fixture',input=str(inp),input_sha256=digest(inp),width_m=200,tx_m=[100.,50.,.0125],rx_m=[101.3,50.,.0125])]))
        with patch.object(crop.base,'check_files'):
            crop.audit(contract,True);checks.append('Native source/grid/time invariants accepted')
            mutations=[('shape','nx_ny_nz',[6400,3000,1],[8000,3000,1]),
                       ('time','dt',2e-9,1e-9)]
            for label,key,value,restore in mutations:
                with h5py.File(raw,'r+') as h:h.attrs[key]=value
                try:crop.audit(contract,True)
                except (AssertionError,ValueError):checks.append(label+' change rejected')
                else:raise AssertionError('Invalid crop accepted')
                with h5py.File(raw,'r+') as h:h.attrs[key]=restore
            with h5py.File(raw,'r+') as h:h['srcs/src1/excitation/samples'][1]=-1.
            try:crop.audit(contract,True)
            except AssertionError:checks.append('Source polarity change rejected')
            else:raise AssertionError('Source change accepted')
            with h5py.File(raw,'r+') as h:
                h['srcs/src1/excitation/samples'][1]=1.
                del h['rxs/rx1/Ez'];h['rxs/rx1/Ez']=np.array([0.,1.,-.3,.1],dtype=np.float32)
            try:crop.audit(contract,True)
            except ValueError:checks.append('Original FP32 trace rejected')
            else:raise AssertionError('Precision cast accepted')
    out.mkdir(parents=True)
    save_json(out/'checks.json',dict(status='PASS',calls_solver=False,checks=checks,
        scope='Synthetic native HDF5 rejection only; not domain equivalence or physical validation'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    check(p.parse_args().out.resolve())
