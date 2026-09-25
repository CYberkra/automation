"""Official paired SFCW versus the archived, independently checked halfspace reference."""
import argparse, csv, hashlib, json
from dataclasses import replace
from pathlib import Path
import h5py
import numpy as np
from scipy.constants import c
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response


from research_spectral_metrics import metrics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('air','medium','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--grid-z',type=float,choices=(.04,.05),default=.04)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    reference_path=Path('artifacts/research_checks/2026-09-24_halfspace_comparison/comparison.csv')
    old=np.genfromtxt(reference_path,delimiter=',',names=True)
    freq=old['frequency_Hz'];ref=old['reference_real']+1j*old['reference_imag']
    baseline=old['FDTD_real']+1j*old['FDTD_imag']
    assert np.array_equal(freq,np.linspace(20e6,170e6,501)) and np.all(abs(ref)>0)
    paths={'air':a.air,'medium':a.medium}
    src={k:load_source(v) for k,v in paths.items()}
    rx={k:load_receiver(v,receiver_path='name:measurement',component='Ex') for k,v in paths.items()}
    dt=rx['air'].dt
    assert np.isclose(dt,1/(c*np.sqrt(2/.05**2+1/a.grid_z**2)),rtol=1e-12,atol=0)
    for k,path in paths.items():
        with h5py.File(path,'r') as h:
            assert np.allclose(h.attrs['dx_dy_dz'],[.05,.05,a.grid_z],rtol=0,atol=1e-12)
            assert np.array_equal(h.attrs['nx_ny_nz'],[480,480,round(24/a.grid_z)])
        assert rx[k].dt==dt and rx[k].time_offset==0 and len(rx[k].samples)==len(rx['air'].samples)
        assert src[k].dt==dt and src[k].time_offset==dt/2 and src[k].spatial_scale==.05
        assert np.array_equal(src[k].samples,src['air'].samples)
        assert np.all(np.isfinite(rx[k].samples))
    n=min(len(rx['air'].samples),int(np.floor(800e-9/dt))+1)
    spectra={};conditions=[];rows=[]
    for taper in (0,200,400):
        count=round(taper*1e-9/dt);fraction=0 if not count else (count-.25)/n
        h={}
        for k in paths:
            response=direct_frequency_response(src[k],replace(rx[k],samples=rx[k].samples[:n]),freq,tail_taper_fraction=fraction)
            assert np.all(response.source_valid) and np.all(np.isfinite(response.response))
            h[k]=response.response
        spectra[taper]=h['medium']-h['air']
        conditions.append(dict(taper_ns=taper,**metrics(spectra[taper],ref)))
        rows.extend(zip([taper]*len(freq),freq,spectra[taper].real,spectra[taper].imag))
    measured=spectra[400];phase=np.angle(measured/ref,deg=True);amp=20*np.log10(abs(measured)/abs(ref))
    # Same unfitted axial/interface approximations recorded before this pair.
    u=np.sin(np.pi*freq*dt)/(c*dt/a.grid_z);q1=2*np.arcsin(u);q2=2*np.arcsin(3*u)
    predicted_phase=-(q1/a.grid_z-2*np.pi*freq/c)*30*180/np.pi
    gamma=(np.sin(q1)-np.sin(q2))/(np.sin(q1)+np.sin(q2))
    predicted_amp=20*np.log10(abs(gamma)/.5)
    diff=rx['medium'].samples-rx['air'].samples;early=diff[rx['air'].times<80e-9]
    current=metrics(measured,ref);prior=metrics(baseline,ref)
    result=dict(grid_z_m=a.grid_z,dt_s=dt,full_samples=len(diff),transformed_samples=n,
                inputs={k:dict(path=str(v),sha256=hashlib.sha256(v.read_bytes()).hexdigest()) for k,v in {**paths,'reference':reference_path}.items()},
                baseline_5cm=prior,refined=current,conditions=conditions,
                reduction_percent={k:100*(1-current[k]/prior[k]) for k in prior},
                early_0_80ns_peak_V_m=float(max(abs(early))),early_exact_zero=bool(np.all(early==0)),
                taper_200_vs400_max_relative_change=float(max(abs(spectra[200]-measured)/abs(measured))),
                phase_170MHz_deg=float(phase[-1]),amplitude_170MHz_dB=float(amp[-1]),
                max_phase_prediction_residual_deg=float(max(abs(phase-predicted_phase))),
                max_amplitude_prediction_residual_dB=float(max(abs(amp-predicted_amp))),
                grid_convergence_certified=False,physical_accuracy_threshold=None,solver_invoked=False,
                limitation='One vertical refinement with dt and PML sampling changes; lateral mesh unchanged; no fitted correction.')
    (a.output/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    with (a.output/'frequency_comparison.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['taper_ns','frequency_Hz','reflected_real','reflected_imag']);w.writerows(rows)
    np.savez_compressed(a.output/'comparison.npz',frequency_Hz=freq,reference=ref,baseline=baseline,refined=measured,time_s=rx['air'].times,raw_difference=diff,phase_prediction_deg=predicted_phase,amplitude_prediction_dB=predicted_amp)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained')
    for label,h in [('dz=5 cm',baseline),(f'dz={a.grid_z*100:g} cm',measured)]:
        ax[0].plot(freq/1e6,20*np.log10(abs(h)/abs(ref)),label=label)
        ax[1].plot(freq/1e6,np.angle(h/ref,deg=True),label=label)
        ax[2].semilogy(freq/1e6,abs(h-ref)/abs(ref),label=label)
    for axis,ylabel in zip(ax,('Amplitude error [dB]','Phase error [degrees]','Relative complex error')):
        axis.set(xlabel='Frequency [MHz]',ylabel=ylabel);axis.legend()
    fig.savefig(a.output/'refinement.png',dpi=170);plt.close(fig)
    print(json.dumps(result))


if __name__=='__main__':main()
