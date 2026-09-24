"""Paired ideal-interface/air diagnostic via official V4 SFCW; no solver."""
import argparse,csv,hashlib,json
from pathlib import Path
from dataclasses import replace
import numpy as np
import h5py
from scipy.constants import c
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ground',type=Path,required=True);p.add_argument('--air',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 sources=[load_source(x) for x in (a.ground,a.air)];rx=[load_receiver(x,receiver_path='name:measurement',component='Ex') for x in (a.ground,a.air)]
 for x in (a.ground,a.air):
  with h5py.File(x,'r') as f:
   assert f.attrs['gprMax']=='4.0.0' and f['rxs/rx1/Ex'].dtype==np.float64
   assert np.array_equal(f.attrs['dx_dy_dz'],[.05]*3)
   assert np.allclose(f['srcs/src1'].attrs['Position'],[12,11.3,19]) and np.allclose(f['rxs/rx1'].attrs['Position'],[12,12.6,19])
 assert rx[0].dt==rx[1].dt and rx[0].time_offset==rx[1].time_offset and len(rx[0].samples)==len(rx[1].samples)
 assert np.array_equal(sources[0].samples,sources[1].samples) and sources[0].dt==sources[1].dt and sources[0].time_offset==sources[1].time_offset
 assert all(np.all(np.isfinite(x.samples)) for x in rx)
 dt=rx[0].dt;n=min(len(rx[0].samples),int(np.floor(800e-9/dt))+1);freq=np.linspace(20e6,170e6,501)
 delta=rx[0].samples-rx[1].samples;times=rx[0].times;rows=[];conditions=[];spectra={}
 for taper_ns in (0,200,400):
  count=round(taper_ns*1e-9/dt);fraction=0 if not count else (count-.25)/n
  hs=[direct_frequency_response(s,replace(y,samples=y.samples[:n]),freq,tail_taper_fraction=fraction) for s,y in zip(sources,rx)]
  hd=direct_frequency_response(sources[0],replace(rx[0],samples=delta[:n]),freq,tail_taper_fraction=fraction)
  assert all(np.all(z.source_valid) and np.all(np.isfinite(z.response)) for z in hs+[hd])
  diff=hs[0].response-hs[1].response;err=float(max(abs(diff-hd.response))/max(abs(diff)))
  assert err<1e-9,'Official transform linearity mismatch'
  spectra[taper_ns]=diff
  conditions.append(dict(taper_ns=taper_ns,taper_samples=count,fraction=fraction,linearity_relative_max=err,max_abs_response=float(max(abs(diff)))))
  rows.extend(zip([taper_ns]*501,freq,hs[0].response.real,hs[0].response.imag,hs[1].response.real,hs[1].response.imag,diff.real,diff.imag))
 diagnostics=[]
 for lo,hi in ((0,80),(80,140),(140,400),(400,800)):
  mask=(times>=lo*1e-9)&(times<hi*1e-9);y=delta[mask];ix=np.argmax(abs(y));diagnostics.append(dict(start_ns=lo,stop_ns=hi,peak_abs_V_m=float(abs(y[ix])),peak_time_ns=float(times[mask][ix]*1e9),rms_V_m=float(np.sqrt(np.mean(y*y)))))
 result=dict(inputs={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in [('ground',a.ground),('air',a.air)]},samples=n,dt=dt,geometric_surface_path_ns=float(np.hypot(30,1.3)/c*1e9),conditions=conditions,raw_difference_windows=diagnostics,taper200_vs400_max_difference_normalized_by_peak=float(max(abs(spectra[200]-spectra[400]))/max(abs(spectra[400]))),solver_invoked=False,official_sfcw=True,convergence_certified=False,site_materials_measured=False)
 (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
 with (a.output/'frequency_response.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['taper_ns','frequency_Hz','ground_real','ground_imag','air_real','air_imag','difference_real','difference_imag']);w.writerows(rows)
 np.savez_compressed(a.output/'paired_raw.npz',time_s=times,ground=rx[0].samples,air=rx[1].samples,difference=delta)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained');ax[0].plot(times*1e9,delta,lw=.7);ax[0].axvline(result['geometric_surface_path_ns'],color='black',ls='--',label='Geometric surface path');ax[0].set(xlim=(0,400),xlabel='Time [ns]',ylabel='Ground minus air Ex [V/m]');ax[0].legend(fontsize=8)
 for taper,h in spectra.items():ax[1].semilogy(freq/1e6,abs(h),label=f'{taper} ns taper')
 ax[1].set(xlabel='Frequency [MHz]',ylabel='Abs difference response [V/(m A)]');ax[1].legend();fig.savefig(a.output/'interface.png',dpi=170);plt.close(fig)
 print(json.dumps(result))
if __name__=='__main__':main()
