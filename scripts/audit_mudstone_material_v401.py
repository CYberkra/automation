"""Verify actual frozen material JSON against gprMax4.0.1 constitutive API."""
import argparse
import hashlib
import json
from pathlib import Path

import gprMax
from gprMax import config
from gprMax.materials import DispersiveMaterial
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(a):
    assert not a.out.exists() and gprMax.__version__=='4.0.1'
    f=20e6+np.arange(501)*300000.;rows=[]
    for path in a.material:
        db=json.loads(path.read_text('utf-8'))
        p=db['materials']['material_002_mudstone'];base=p['base']
        native=DispersiveMaterial(2,'mudstone');native.type='debye'
        native.er=base['relative_permittivity'];native.se=base['electric_conductivity_s_per_m']
        native.poles=len(p['poles'])
        native.deltaer=[x['relative_permittivity_difference'] for x in p['poles']]
        native.tau=[x['relaxation_time_s'] for x in p['poles']]
        pol=sum(x['relative_permittivity_difference']/(1+2j*np.pi*f*x['relaxation_time_s']) for x in p['poles'])
        dc=base['electric_conductivity_s_per_m']/(2*np.pi*f*config.e0)
        direct=base['relative_permittivity']+pol-1j*dc
        official=np.array([native.calculate_er(x) for x in f])
        error=float(np.linalg.norm(direct-official)/np.linalg.norm(direct));assert error<1e-14
        np.testing.assert_allclose(-official.imag,-pol.imag+dc,rtol=1e-14,atol=0)
        i=250;assert f[i]==95e6
        rows.append(dict(material_sha256=sha(path),sigma_DC_S_m=base['electric_conductivity_s_per_m'],
                         official_relative_L2=error,epsilon_real95=float(official[i].real),
                         polarization_loss95=float(-pol[i].imag),DC_loss95=float(dc[i]),
                         total_epsilon_loss95=float(-official[i].imag),
                         effective_sigma95_S_m=float(-official[i].imag*2*np.pi*f[i]*config.e0)))
    result=dict(status='PASS_FROZEN_JSON_OFFICIAL_V401_SPECTRUM_COMPONENTS',
                script_sha256=sha(__file__),materials_py_sha256=sha(Path(gprMax.__file__).parent/'materials.py'),
                version=gprMax.__version__,tones=501,rows=rows,
                scope='Checks actual formula includes exactly the specified DC plus Debye terms. Does not certify independent physical origin or site representativeness of their chosen values.')
    a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--material',type=Path,nargs='+',required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
