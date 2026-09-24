"""Bind unresolved FDTD evidence to the existing evaluation contract, not labels."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from research_evaluation_contract import waveform_metrics

p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
root=Path('artifacts/research_checks/2026-09-25_deep_controls_analysis');evidence=json.loads((root/'results.json').read_text(encoding='utf-8'))
with np.load(root/'arrays.npz') as data:
 f=data['frequency_Hz'];t=data['time_s'];window=np.hanning(len(f));signals={}
 for name in ('baseline','FINE_200'):
  analytic=np.fft.ifft(data[name]*window,n=len(t))*len(t)/window.sum()
  signals[name]=(2*np.real(analytic*np.exp(2j*np.pi*f[0]*t)))[:,None]
mask=((t>=440e-9)&(t<=600e-9))[:,None]
# Fixed before any candidate execution; same real passband synthesis for both grids.
gates={name:waveform_metrics(s,s,mask,reference_kind='paired_contrast',state='numerically_unresolved',scope='contrast') for name,s in signals.items()}
assert all(not g['available'] and not g['absolute_preservation_eligible'] and g['metrics'] is None for g in gates.values())
observed=float(np.linalg.norm((signals['FINE_200']-signals['baseline'])[mask])/np.linalg.norm(signals['baseline'][mask]))
record={'reference_id':'depth20_paired_contrast_v1','reference_kind':'paired_contrast','reference_state':'numerically_unresolved','scope':'contrast',
 'physical_or_constructed_definition':'2D TMx hypothetical wet block minus identical background, 20m center depth. Common Hann-band real synthesis, no alignment or gain.',
 'mask_id':'fixed_440_600ns','mask_sha256':hashlib.sha256(mask.tobytes()).hexdigest(),'reference_sha256':hashlib.sha256((root/'arrays.npz').read_bytes()).hexdigest(),
 'axis_units':['s','single modeled receiver'],'amplitude_units':'(V/m)/A; Hann band synthesis','frozen_before_candidates':True,
 'common_protected_objects':['cover/rock interfaces remain protected, not pure clutter'],'numerical_credibility':{'fullband_grid_relative_L2':evidence['comparisons_to_original_2p5cm']['FINE']['all_relative_L2_difference'],'fixed_window_real_waveform_relative_L2':observed},
 'uncertainty_kind':'Observed grid sensitivity, NOT a rigorous error bound or a measured noise distribution',
 'contract_probe':gates,'probe_note':'Identity-input eligibility probes only; no processing candidates have been evaluated. Even exact self-match must not produce physical preservation labels.',
 'physical_label_eligible':False,'training_eligible':False,'unique_label':None,'gain_quality_label':None,'gain_note':'No independent unattenuated target response; depth ratios are not a compensation ground truth.',
 'next_use':'Explicit mechanism diagnostics may compare results on both grids, but must not promote either to physical truth.',
 'source_evidence_sha256':hashlib.sha256((root/'results.json').read_bytes()).hexdigest(),'solver_invoked':False}
(a.output/'reference_eligibility.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
np.savez_compressed(a.output/'reference_arrays.npz',time_s=t,mask=mask,**signals)
print(json.dumps({'eligibility':False,'identity_probe_rejected':True,'window_grid_difference':observed}))
