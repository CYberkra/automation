"""Paired fixed/automatic PML mechanism diagnostic using official SFCW."""
import argparse,csv,hashlib,json
from pathlib import Path
from dataclasses import replace
import numpy as np
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response,reconstruct_time_response

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ('fixed','auto','air','output'):p.add_argument('--'+k,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);paths={k:getattr(a,k) for k in ('fixed','auto','air')}
 src={k:load_source(v) for k,v in paths.items()};rx={k:load_receiver(v,receiver_path='name:measurement',component='Ex') for k,v in paths.items()};base=rx['air'];freq=np.linspace(20e6,170e6,501);n=min(len(base.samples),int(np.floor(800e-9/base.dt))+1)
 for k in paths:
  assert rx[k].dt==base.dt and rx[k].time_offset==base.time_offset and len(rx[k].samples)==len(base.samples)
  assert src[k].dt==src['air'].dt and src[k].time_offset==src['air'].time_offset and src[k].spatial_scale==src['air'].spatial_scale
  assert np.array_equal(src[k].samples,src['air'].samples) and np.all(np.isfinite(rx[k].samples))
 t=base.times;diff={k:rx[k].samples-base.samples for k in ('fixed','auto')};windows=[]
 for lo,hi in ((0,80),(80,100),(100,140),(140,400),(400,800)):
  mask=(t>=lo*1e-9)&(t<hi*1e-9);row=dict(start_ns=lo,stop_ns=hi)
  for k,y in diff.items():
   v=y[mask];row[k]=dict(peak_abs_V_m=float(max(abs(v))),rms_V_m=float(np.sqrt(np.mean(v*v))))
  old=row['auto']['peak_abs_V_m'];new=row['fixed']['peak_abs_V_m'];row['fixed_over_auto_peak']=new/old if old else None;row['reduction_dB']=float(20*np.log10(old/new)) if new and old else None;row['fixed_exact_zero']=bool(np.all(diff['fixed'][mask]==0));windows.append(row)
 rows=[];conditions=[];responses={};template=None
 for taper in (0,200,400):
  count=round(taper*1e-9/base.dt);fraction=0 if not count else (count-.25)/n;h={}
  for k in paths:
   fr=direct_frequency_response(src[k],replace(rx[k],samples=rx[k].samples[:n]),freq,tail_taper_fraction=fraction)
   assert np.all(fr.source_valid) and np.all(np.isfinite(fr.response));h[k]=fr.response;template=fr
  fixed=h['fixed']-h['air'];auto=h['auto']-h['air'];delta=fixed-auto;assert np.all(abs(auto)>0)
  conditions.append(dict(taper_ns=taper,max_difference_over_auto_peak=float(max(abs(delta))/max(abs(auto))),max_pointwise_relative_difference=float(max(abs(delta)/abs(auto))),median_pointwise_relative_difference=float(np.median(abs(delta)/abs(auto)))))
  rows.extend(zip([taper]*501,freq,fixed.real,fixed.imag,auto.real,auto.imag,abs(delta)/abs(auto)))
  responses[taper]={'fixed':fixed,'auto':auto}
 recon={k:reconstruct_time_response(replace(template,response=v),window='hann',zero_pad_factor=8,time_shift=0) for k,v in responses[400].items()}
 result=dict(inputs={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in paths.items()},windows=windows,spectral_comparisons=conditions,solver_invoked=False,convergence_certified=False,causal_intervention='Explicit sigma on all six faces, not side-only. Identical air reference and source histories; unchanged physical materials and grid.',bandlimited_peak_ns={k:float(v.time[np.argmax(abs(v.complex_envelope))]*1e9) for k,v in recon.items()})
 (a.output/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
 with (a.output/'frequency_comparison.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['taper_ns','frequency_Hz','fixed_real','fixed_imag','auto_real','auto_imag','pointwise_relative_difference']);w.writerows(rows)
 np.savez_compressed(a.output/'raw_difference.npz',time_s=t,**diff)
 np.savez_compressed(a.output/'reconstruction.npz',time_s=recon['fixed'].time,**{k:v.complex_envelope for k,v in recon.items()})
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained')
 for k,y in diff.items():ax[0].plot(t*1e9,y,label=k,lw=.7)
 ax[0].set(xlim=(55,100),xlabel='Time [ns]',ylabel='Air-subtracted raw Ex [V/m]');ax[0].legend()
 for taper,data in responses.items():ax[1].semilogy(freq/1e6,abs(data['fixed']-data['auto'])/abs(data['auto']),label=f'{taper} ns taper')
 ax[1].set(xlabel='Frequency [MHz]',ylabel='Relative change from automatic PML');ax[1].legend(fontsize=8)
 for k,v in recon.items():ax[2].plot(v.time*1e9,abs(v.complex_envelope),label=k)
 ax[2].set(xlim=(60,150),xlabel='Time [ns]',ylabel='Hann band-limited envelope [V/(m A)]');ax[2].legend();fig.savefig(a.output/'fixedPML.png',dpi=170);plt.close(fig);print(json.dumps(result))
if __name__=='__main__':main()
