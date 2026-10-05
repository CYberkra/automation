"""Energy, passivity and known-delay checks for the independent plane reference."""
import argparse
import json
from pathlib import Path
import numpy as np
from analyze_uav_local_plane import C,F,OMEGA,reference
from hs_capsule_identity import sha256


def main(out):
    if out.exists():
        raise ValueError('fresh output required')
    checks={}
    k0=OMEGA/C
    for case,n in [('halfspace',3),('slab',2)]:
        r,t=reference(case)
        checks[case+'_energy_balance_max_error']=float(abs(abs(r)**2+n*abs(t)**2-1).max())
        if checks[case+'_energy_balance_max_error']>1e-12:
            raise ValueError('lossless flux balance failed')
    r,t=reference('halfspace')
    checks['halfspace_surface_polarity_error']=float(abs(r*np.exp(2j*k0*4)+.5).max())
    checks['halfspace_known_transmission_delay_error']=float(abs(t-.5*np.exp(-1j*OMEGA*8/C)).max())
    if max(checks['halfspace_surface_polarity_error'],checks['halfspace_known_transmission_delay_error'])>1e-12:
        raise ValueError('polarity/delay identity failed')
    for case in ('conductive','debye'):
        r,t=reference(case)
        flux=abs(r)**2+2*abs(t)**2
        checks[case+'_lossy_flux_minmax']=[float(flux.min()),float(flux.max())]
        if not ((flux>=0)&(flux<1)).all():
            raise ValueError('passivity failed')
    prediction={}
    for dl in (.05,.025,.0125):
        dt=dl/(C*np.sqrt(2))
        numerical_k=2*np.arcsin(3*np.sqrt(2)*np.sin(OMEGA*dt/2))/dl
        numerical_air_k=2*np.arcsin(np.sqrt(2)*np.sin(OMEGA*dt/2))/dl
        phase=-(numerical_k-numerical_air_k-2*k0)*4*180/np.pi
        prediction[str(dl)]={'four_metre_epsilon9_minus_air_phase_error_at170MHz_deg':float(phase[-1]),
                             'half_cell_air_reflection_location_phase_at170MHz_deg':float(k0[-1]*dl*180/np.pi)}
    out.mkdir(parents=True)
    result=dict(status='PASS',checks=checks,lossless_Yee_grid_predictions=prediction,
                scope='Analytic reference identities and precomputed grid-phase prediction, not native convergence proof',
                code_sha256=sha256(__file__),reference_code_sha256=sha256(Path(__file__).parent/'analyze_uav_local_plane.py'))
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    main(ap.parse_args().out)
