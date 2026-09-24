"""Native-time official SFCW target/background differences; no gain or detectability claim."""
import argparse,json,hashlib
from dataclasses import replace
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c
from gprMax.toolboxes.SFCW.processing import load_source,load_receiver,direct_frequency_response

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
names=['DEP_BG','DEP_05','DEP_10','DEP_20'];paths={k:Path(f'artifacts/research_checks/2026-09-25_{k}/{k}.h5') for k in names}
src={k:load_source(v) for k,v in paths.items()};rx={k:load_receiver(v,receiver_path='name:measurement',component='Ex') for k,v in paths.items()}
dt=src['DEP_BG'].dt;f=np.linspace(20e6,170e6,501);n=min(len(rx['DEP_BG'].samples),int(np.floor(1200e-9/dt))+1)
arrays={'frequency_Hz':f};checks={};responses={}
for k in names:
 with h5py.File(paths[k]) as h:
  assert np.array_equal(h.attrs['nx_ny_nz'],[1,1280,2000]) and h[rx[k].path].dtype==np.float64
  assert np.allclose(h[rx[k].path].parent.attrs['Position'][1:],[16.65,45],rtol=0,atol=1e-12)
 assert rx[k].dt==src[k].dt==dt and rx[k].time_offset==0 and src[k].time_offset==dt/2
 assert len(rx[k].samples)==len(rx['DEP_BG'].samples) and np.array_equal(src[k].samples,src['DEP_BG'].samples)
 assert np.all(np.isfinite(rx[k].samples));checks[k+'_metadata_finite_source']=True
 arrays[k+'_raw']=rx[k].samples
 for taper in (200,400):
  result=direct_frequency_response(src[k],replace(rx[k],samples=rx[k].samples[:n]),f,tail_taper_fraction=(round(taper*1e-9/dt)-.25)/n)
  assert np.all(result.source_valid) and np.all(np.isfinite(result.response));responses[k,taper]=result.response
arrays['background']=responses['DEP_BG',200];rows={};t=np.arange(8192)/(8192*300e3);arrays['envelope_time_s']=t;window=np.hanning(501)
for depth in (5,10,20):
 k=f'DEP_{depth:02d}';delta=responses[k,200]-responses['DEP_BG',200];other=responses[k,400]-responses['DEP_BG',400]
 arrays[k+'_difference']=delta;arrays[k+'_difference_taper400']=other
 envelope=np.abs(np.fft.ifft(delta*window,n=8192))*8192/window.sum();arrays[k+'_envelope']=envelope
 center=2*(15+3*4+(depth-3)*3)/c;gate=abs(t-center)<=80e-9;index=np.flatnonzero(gate)[np.argmax(envelope[gate])]
 norm=np.linalg.norm(delta);rawdiff=rx[k].samples-rx['DEP_BG'].samples
 rows[k]={'center_depth_m':depth,'vertical_lossless_center_time_ns':center*1e9,'gated_envelope_peak_time_ns':float(t[index]*1e9),'gated_envelope_peak':float(envelope[index]),'spectral_RMS':float(norm/np.sqrt(501)),
 'difference_to_background_norm_dB':float(20*np.log10(norm/np.linalg.norm(arrays['background']))),'tail_window_difference_over_target_norm':float(np.linalg.norm(other-delta)/norm),
 'raw_difference_pre100ns_peak_over_full_peak':float(max(abs(rawdiff[:int(100e-9/dt)]))/max(abs(rawdiff)))}
for row in rows.values():row['spectral_RMS_relative_5m_dB']=float(20*np.log10(row['spectral_RMS']/rows['DEP_05']['spectral_RMS']))
result={'checks':checks,'rows':rows,'spectrum_units':'(V/m)/A','envelope':'Hann weighted positive-band complex baseband envelope; IFFT zero-padding only interpolates. No physical noise floor, gain, or per-depth normalization.','source_frequency_count':501,'inputs':{str(v):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths.values()},'solver_invoked':False,'detectability_claim':False,'grid_and_boundary_validation_pending':True}
(a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');np.savez_compressed(a.output/'arrays.npz',**arrays)
np.savetxt(a.output/'spectra.csv',np.column_stack([f]+[v for k in names[1:] for v in (arrays[k+'_difference'].real,arrays[k+'_difference'].imag)]),delimiter=',',header='frequency_Hz,'+','.join(k+s for k in names[1:] for s in ('_real','_imag')),comments='')
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for depth in (5,10,20):
 k=f'DEP_{depth:02d}';ax[0].semilogy(f/1e6,abs(arrays[k+'_difference']),label=f'{depth} m');ax[1].semilogy(t*1e9,arrays[k+'_envelope'],label=f'{depth} m')
ax[0].set(xlabel='Frequency [MHz]',ylabel='Target - background magnitude [(V/m)/A]')
ax[1].set(xlabel='Time [ns]',ylabel='Hann-band envelope [(V/m)/A]',xlim=(100,700))
for axis in ax:axis.legend();axis.grid(alpha=.2)
fig.savefig(a.output/'depth_comparison.png',dpi=160);plt.close(fig)
print(json.dumps(result))
