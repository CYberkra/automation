"""Validate independent halfspace quadrature and compare saved official SFCW CSV."""
import argparse,csv,json,time,hashlib,sys,importlib.metadata
from pathlib import Path
import numpy as np
import empymod
from halfspace_sommerfeld_reference import reflected,transverse

def metric(h,ref):
 e=abs(h-ref)/abs(ref)
 return dict(max_relative_complex_error=float(max(e)),median_relative_complex_error=float(np.median(e)),max_abs_amplitude_error_dB=float(max(abs(20*np.log10(abs(h)/abs(ref))))),max_abs_phase_error_deg=float(max(abs(np.angle(h/ref,deg=True)))))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);start=time.perf_counter()
 assert empymod.__version__=='1.10.6'
 f=np.linspace(20e6,170e6,501);data=np.genfromtxt(a.input,delimiter=',',names=True);data=data[data['taper_ns']==400];assert np.array_equal(data['frequency_Hz'],f)
 measured=data['fixed_real']+1j*data['fixed_imag']
 responses={n:np.array([reflected(fi,order=n) for fi in f]) for n in (128,256,512)};reference=responses[512];long_tail=np.array([reflected(fi,order=512,umax=3) for fi in f])
 convergence={str(n):float(max(abs(v-reference)/abs(reference))) for n,v in responses.items() if n!=512};convergence['umax3_vs2']=float(max(abs(long_tail-reference)/abs(reference)))
 fs=empymod.dipole(src=[0,0,-15],rec=[0,1.3,-15],depth=[],res=[np.inf],epermH=[1],freqtime=f,ab=11,xdirect=True,verb=1)*.05
 fs_error=float(max(abs(fs-transverse(f,1.3))/abs(transverse(f,1.3))))
 sample=f[[0,100,250,400,500]];no_interface=np.array([reflected(fi,er=1,order=256) for fi in sample]);pec_theory=-transverse(sample,np.hypot(30,1.3));pec_checks={}
 for sigma in (1e6,1e10):pec_checks[str(sigma)]=float(max(abs(np.array([reflected(fi,er=1,sigma_lower=sigma,order=512) for fi in sample])-pec_theory)/abs(pec_theory)))
 checks={'finite_all_references':bool(all(np.all(np.isfinite(v)) for v in responses.values()) and np.all(np.isfinite(long_tail))), 'integration_changes_below_1e_minus5':all(v<1e-5 for v in convergence.values()),'fullspace_units_phase_below_1e_minus7':fs_error<1e-7,'matched_material_exact_zero':bool(np.all(no_interface==0)),'conductive_limit_approaches_PEC':pec_checks['10000000000.0']<pec_checks['1000000.0'] and pec_checks['10000000000.0']<1e-5}
 if not all(checks.values()):
  (a.output/'validation_failure.json').write_text(json.dumps(dict(checks=checks,convergence=convergence,fullspace_error=fs_error,pec=pec_checks),indent=2));raise ValueError('Reference numerical checks failed')
 # This comparison is unfitted; do not replace the delivered spectrum by a correction.
 old_budget=json.loads(Path('artifacts/research_checks/2026-09-24_PEC_sources/phase_budget_before_result.json').read_text());pred=np.array(old_budget['predicted_response_phase_error_deg']);phase=np.angle(measured/reference,deg=True)
 result=dict(reference_package=empymod.__version__,reference_type='Full-wave spectral kernel with project branch-aware Gauss integration; exact zero conductivity, infinite planar interface',inputs=dict(path=str(a.input),sha256=hashlib.sha256(a.input.read_bytes()).hexdigest()),checks=checks,integration_relative_changes=convergence,fullspace_max_relative_error=fs_error,PEC_conductivity_limit_errors=pec_checks,fdtd_vs_reference=metric(measured,reference),phase_at_170MHz_deg=float(phase[-1]),max_phase_difference_from_prior_axial_budget_deg=float(max(abs(phase-pred))),reference_numerical_check_thresholds_only=True,fdtd_physical_accuracy_threshold=None,grid_convergence_certified=False,FDTD_solver_invoked=False,wall_s=time.perf_counter()-start)
 (a.output/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
 with (a.output/'comparison.csv').open('w',newline='') as out:
  w=csv.writer(out);w.writerow(['frequency_Hz','FDTD_real','FDTD_imag','reference_real','reference_imag','relative_complex_error','amplitude_error_dB','phase_error_deg','axial_phase_prediction_deg']);w.writerows(zip(f,measured.real,measured.imag,reference.real,reference.imag,abs(measured-reference)/abs(reference),20*np.log10(abs(measured)/abs(reference)),phase,pred))
 np.savez_compressed(a.output/'reference_arrays.npz',frequency_Hz=f,order128=responses[128],order256=responses[256],order512=reference,order512_umax3=long_tail)
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,ax=plt.subplots(1,3,figsize=(14,4),layout='constrained');ax[0].plot(f/1e6,20*np.log10(abs(measured)/abs(reference)));ax[0].set(xlabel='Frequency [MHz]',ylabel='FDTD/reference amplitude error [dB]');ax[1].plot(f/1e6,phase,label='Observed dielectric error');ax[1].plot(f/1e6,pred,'--',label='Prior axial Yee budget');ax[1].set(xlabel='Frequency [MHz]',ylabel='Phase error [degrees]');ax[1].legend(fontsize=8);ax[2].semilogy(f/1e6,abs(measured-reference)/abs(reference));ax[2].set(xlabel='Frequency [MHz]',ylabel='Relative complex error');fig.savefig(a.output/'halfspace_comparison.png',dpi=170);plt.close(fig)
 versions={n:importlib.metadata.version(n) for n in ('empymod','numpy','scipy','matplotlib')};identity={'executable':sys.executable,'python':sys.version,'versions':versions,'kernel_sha256':hashlib.sha256(Path(empymod.kernel.__file__).read_bytes()).hexdigest(),'integration_script_sha256':hashlib.sha256(Path('scripts/halfspace_sommerfeld_reference.py').read_bytes()).hexdigest()};(a.output/'runtime_identity.json').write_bytes((json.dumps(identity,indent=2)+'\n').encode());print(json.dumps(result))
if __name__=='__main__':main()
