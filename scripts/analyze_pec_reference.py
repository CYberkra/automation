"""Unfitted infinite-PEC image reference for paired official SFCW outputs."""
import argparse,csv,json,hashlib
from pathlib import Path
from dataclasses import replace
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response,reconstruct_time_response
from audit_official_sfcw import transverse_dipole

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--pec',type=Path,required=True);p.add_argument('--air',type=Path,required=True);p.add_argument('--dielectric',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 files={'PEC':a.pec,'air':a.air,'dielectric':a.dielectric};rx={k:load_receiver(v,receiver_path='name:measurement',component='Ex') for k,v in files.items()};src={k:load_source(v) for k,v in files.items()};base=rx['air'];freq=np.linspace(20e6,170e6,501);n=min(len(base.samples),int(np.floor(800e-9/base.dt))+1)
 for k in files:
  assert rx[k].dt==base.dt and rx[k].time_offset==base.time_offset and len(rx[k].samples)==len(base.samples)
  assert np.array_equal(src[k].samples,src['air'].samples) and src[k].spatial_scale==.05 and src[k].time_offset==src['air'].time_offset
 theory=-transverse_dipole(freq,np.hypot(30,1.3),.05);stats=[];rows=[];hs={};template=None
 for taper in (0,200,400):
  count=round(taper*1e-9/base.dt);fraction=0 if not count else (count-.25)/n
  spectra={}
  for k in files:
   fr=direct_frequency_response(src[k],replace(rx[k],samples=rx[k].samples[:n]),freq,tail_taper_fraction=fraction)
   assert np.all(fr.source_valid) and np.all(np.isfinite(fr.response));spectra[k]=fr.response;template=fr
  h=spectra['PEC']-spectra['air'];hs[taper]=h;error=abs(h-theory)/abs(theory)
  stats.append(dict(taper_ns=taper,max_relative_complex_error=float(max(error)),median_relative_complex_error=float(np.median(error)),max_absolute_amplitude_error_dB=float(max(abs(20*np.log10(abs(h)/abs(theory))))),max_absolute_phase_error_deg=float(max(abs(np.angle(h/theory,deg=True))))))
  rows.extend(zip([taper]*501,freq,h.real,h.imag,theory.real,theory.imag,error))
  if taper==400:dielectric=spectra['dielectric']-spectra['air']
 recon={k:reconstruct_time_response(replace(template,response=v),window='hann',zero_pad_factor=8,time_shift=0) for k,v in {'PEC':hs[400],'PEC theory':theory,'dielectric':dielectric}.items()}
 peaks={k:float(v.time[np.argmax(abs(v.complex_envelope))]*1e9) for k,v in recon.items()}
 raw=[]
 for k in ('PEC','dielectric'):
  y=rx[k].samples-base.samples;t=base.times;ff=np.fft.rfftfreq(len(y),base.dt);power=abs(np.fft.rfft((y-y.mean())*np.hanning(len(y))))**2
  raw.append(dict(case=k,peak_before80ns_V_m=float(max(abs(y[t<80e-9]))),dominant_MHz=float(ff[np.argmax(power)]/1e6)))
 result=dict(inputs={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in files.items()},reference='Negative transverse current element at mirror position z=-11m; r=hypot(30,1.3)m; no fitted parameters',conditions=stats,bandlimited_hann_zero_pad8_envelope_peak_ns=peaks,raw_diagnostics=raw,physical_accuracy_threshold=None,convergence_certified=False,solver_invoked=False,notes=['PEC validates only this reference geometry, not a dielectric halfspace.','Bandlimited sidelobes are not causality violations; zero padding adds no physical resolution.'])
 (a.output/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
 with (a.output/'frequency_comparison.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['taper_ns','frequency_Hz','PEC_reflected_real','PEC_reflected_imag','theory_real','theory_imag','relative_complex_error']);w.writerows(rows)
 np.savez_compressed(a.output/'reconstruction.npz',time_s=recon['PEC'].time,**{k.replace(' ','_'):v.complex_envelope for k,v in recon.items()})
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained')
 for taper,h in hs.items():ax[0].semilogy(freq/1e6,abs(h-theory)/abs(theory),label=f'{taper} ns taper')
 ax[0].set(xlabel='Frequency [MHz]',ylabel='PEC reflected relative complex error');ax[0].legend(fontsize=8)
 for k,v in recon.items():ax[1].plot(v.time*1e9,abs(v.complex_envelope),label=k)
 ax[1].set(xlim=(60,150),xlabel='Time [ns]',ylabel='Hann band-limited envelope [V/(m A)]');ax[1].legend(fontsize=8)
 for k in ('PEC','dielectric'):ax[2].plot(base.times*1e9,rx[k].samples-base.samples,label=k,lw=.6)
 ax[2].set(xlim=(55,145),xlabel='Time [ns]',ylabel='Raw air-subtracted Ex [V/m]');ax[2].legend(fontsize=8)
 fig.savefig(a.output/'PEC_reference.png',dpi=170);plt.close(fig);print(json.dumps(result))
if __name__=='__main__':main()
