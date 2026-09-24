"""Compare completed CUDA time-window/PML controls using official SFCW transforms.

Periodograms are diagnostics, not hardware spectra. The static-tail prediction is
an unfitted analytic mechanism check, not a correction offered for geological data.
"""
import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import epsilon_0
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response
from audit_official_sfcw import transverse_dipole


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',action='append',required=True,help='label=raw.h5')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    frequency=np.linspace(20e6,170e6,501)
    records={};spectra={};rows=[];band_rows=[];tails=[]
    for entry in a.case:
        label,path=entry.split('=',1);path=Path(path)
        source=load_source(path);rx=load_receiver(path,receiver_path='name:measurement',component='Ex')
        with h5py.File(path,'r') as f:
            if f.attrs['gprMax']!='4.0.0' or f['rxs/rx1/Ex'].dtype!=np.float64:
                raise ValueError('Expected V4 double output')
            if not np.array_equal(f.attrs['dx_dy_dz'],[.05,.05,.05]):raise ValueError('Unexpected grid')
            if not np.allclose(f['srcs/src1'].attrs['Position'],[12,11.3,19]) or not np.allclose(f['rxs/rx1'].attrs['Position'],[12,12.6,19]):raise ValueError('Unexpected geometry')
        if not np.all(np.isfinite(rx.samples)) or np.count_nonzero(source.samples)!=1:raise ValueError('Invalid raw data')
        theory=transverse_dipole(frequency,1.3,source.spatial_scale)
        static=-source.spatial_scale*np.sum(source.samples)*source.dt/(4*np.pi*epsilon_0*1.3**3)
        records[label]=dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),samples=len(rx.samples),dt=rx.dt,static_prediction_V_m=float(static))
        spectra[label]={}
        for start,stop in ((200,400),(600,800)):
            if rx.times[-1]*1e9<stop-.2:continue
            y=rx.samples[(rx.times>=start*1e-9)&(rx.times<stop*1e-9)]
            f=np.fft.rfftfreq(len(y),rx.dt);power=abs(np.fft.rfft((y-y.mean())*np.hanning(len(y))))**2
            tails.append(dict(case=label,start_ns=start,stop_ns=stop,mean_V_m=float(y.mean()),std_V_m=float(y.std()),peak_to_peak_V_m=float(np.ptp(y)),dominant_MHz=float(f[np.argmax(power)]/1e6),periodogram_fraction_20_170MHz=float(power[(f>=20e6)&(f<=170e6)].sum()/power.sum())))
        for window_ns in (200,400,600,800):
            if rx.times[-1]*1e9<window_ns-.2:continue
            n=min(len(rx.samples),int(np.floor(window_ns*1e-9/rx.dt))+1)
            selected=replace(rx,samples=rx.samples[:n])
            for taper_ns in (0,100,200,400):
                if taper_ns>window_ns/2:continue
                # Match taper sample count, rather than holding fraction fixed as N changes.
                count=int(round(taper_ns*1e-9/rx.dt))
                fraction=0 if count==0 else (count-.25)/n
                response=direct_frequency_response(source,selected,frequency,tail_taper_fraction=fraction)
                h=response.response
                if not np.all(response.source_valid) or not np.all(np.isfinite(h)):raise ValueError('Invalid official response')
                key=f'{window_ns}ns_{taper_ns}ns'
                spectra[label][key]=h
                error=abs(h-theory)/abs(theory)
                row=dict(case=label,window_ns=window_ns,actual_last_sample_ns=float(selected.times[-1]*1e9),taper_ns=taper_ns,taper_samples=count,official_taper_fraction=fraction,max_relative_complex_error=float(max(error)),median_relative_complex_error=float(np.median(error)))
                if taper_ns==0:
                    w=2*np.pi*frequency
                    missing_tail=(static/source.samples[0])*np.exp(-1j*w*(n*rx.dt-source.time_offset))/(-np.expm1(-1j*w*rx.dt))
                    residual=abs(h-(theory-missing_tail))/abs(theory)
                    row['max_residual_with_unfitted_static_tail_prediction']=float(max(residual))
                    row['median_residual_with_unfitted_static_tail_prediction']=float(np.median(residual))
                rows.append(row)
                band_rows.extend(zip([label]*501,[key]*501,frequency,h.real,h.imag,error))
    comparisons=[]
    baseline=spectra.get('baseline',{})
    for label,data in spectra.items():
        if label=='baseline':continue
        for key,h in data.items():
            if key in baseline:
                difference=abs(h-baseline[key])/abs(theory)
                comparisons.append(dict(case=label,condition=key,max_complex_difference_relative_to_theory=float(max(difference)),median_complex_difference_relative_to_theory=float(np.median(difference))))
    result=dict(cases=records,late_diagnostics=tails,response_conditions=rows,matched_comparisons=comparisons,solver_invoked=False,official_sfcw_transform=True,convergence_certified=False,static_model_fitted_parameters=0,limitations=['PML40 changes both thickness and inner interface position; not a pure thickness convergence test.','Periodogram removes segment mean and uses Hann; spectral fraction is a diagnostic, not physical energy.','No grid/CFL refinement or geological late-arrival preservation test.'])
    (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    with (a.output/'frequency_results.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.writer(f);wr.writerow(['case','condition','frequency_Hz','response_real','response_imag','relative_complex_error']);wr.writerows(band_rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for label,data in spectra.items():
        for key in ('400ns_200ns','800ns_200ns','800ns_400ns'):
            if key in data:
                axes[0].semilogy(frequency/1e6,abs(data[key]-theory)/abs(theory),label=f'{label}: {key}')
    axes[0].set(xlabel='Frequency [MHz]',ylabel='Relative complex error vs harmonic theory');axes[0].legend(fontsize=7)
    for label,path in ((v.split('=',1)) for v in a.case):
        rx=load_receiver(path,receiver_path='name:measurement',component='Ex')
        block=max(1,round(10e-9/rx.dt));n=len(rx.samples)//block
        blocks=rx.samples[:n*block].reshape(n,block)
        axes[1].plot((np.arange(n)+.5)*block*rx.dt*1e9,blocks.mean(axis=1),label=label)
    axes[1].axhline(static,color='black',ls='--',label='Predicted static limit')
    axes[1].set(xlabel='Time [ns]; nonoverlapping 10 ns block means',ylabel='Ex [V/m]',xlim=(50,800),ylim=(-.04,0))
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.25)
    fig.savefig(a.output/'controls.png',dpi=160);plt.close(fig)
    print(json.dumps({'cases':list(records),'conditions':len(rows),'comparisons':len(comparisons)}))


if __name__=='__main__':main()
