"""Read-only mechanism checks, distinct from unique-ray or convergence proof."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import h5py
import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from hs_capsule_identity import sha256 as sha

C=299792458.;EPS0=8.8541878128e-12


def read(p):
    return json.loads(p.read_text('utf-8'))


def optical(material,f):
    b=material['base'];sigma=b['electric_conductivity_s_per_m']
    eps=b['relative_permittivity']-1j*sigma/(2*np.pi*f*EPS0)
    derivative=1j*sigma/(2*np.pi*EPS0*f*f)
    for p in material.get('poles',[]):
        delta=p['relative_permittivity_difference'];tau=p['relaxation_time_s']
        eps+=delta/(1+2j*np.pi*f*tau)
        derivative+=-2j*np.pi*tau*delta/(1+2j*np.pi*f*tau)**2
    n=np.sqrt(eps)
    return n,float((n+f*derivative/(2*n)).real),float(-2*np.pi*f*n.imag/C)


def study(root,package,arrays,public):
    m=read(package/'manifest.json');materials=read(package/'geometries/line9_materials_v0_2.json')['materials']
    nc,ng,alpha=optical(materials['material_001_cover'],95e6)
    nm,_,alpha_mud=optical(materials['material_002_mudstone'],95e6)
    ns,_,_=optical(materials['material_003_sandstone'],95e6)
    rt=abs((nc-1)/(nc+1));rb=abs((nc-nm)/(nc+nm));rows=[]
    for x in [190.,170.]:
        g=next(g for g in m['groups'] if g['profile_x_m']==x and g['kind']=='local' and g['role']=='H0')
        p=next(p for p in g['probes'] if p['band']=='cover_middle' and p['x_m']==x-20)
        d=p['surface_y_m']-p['cover_y_m']
        with h5py.File(arrays/(g['id']+'_arrays.h5')) as h:
            t=h['sfcw_time_s'][:]*1e9;z=h['hann/main_complex_bandpass'][:];peaks=[]
            for lo,hi in [(100,220),(220,270) if x==170 else (300,350)]:
                k=np.flatnonzero((t>=lo)&(t<=hi));j=k[np.argmax(abs(z[k]))];peaks.append(dict(time_ns=float(t[j]),amplitude=float(abs(z[j]))))
        ratio=peaks[1]['amplitude']/peaks[0]['amplitude']
        prediction=float(rt*rb*np.exp(-2*alpha*d))
        mud_d=p['cover_y_m']-p['basal_y_m']
        transmission=abs(4*nc*nm/(nc+nm)**2)
        basal_prediction=float(transmission*abs((nm-ns)/(nm+ns))/rb*np.exp(-2*alpha_mud*mud_d))
        h1=next(g for g in m['groups'] if g['profile_x_m']==x and g['kind']=='local' and g['role']=='H1')
        with h5py.File(arrays/(h1['id']+'_arrays.h5')) as h:
            delta=h['hann/main_complex_bandpass'][:]-z;k=np.flatnonzero((t>=330)&(t<=430))
            basal_ratio=float(abs(delta[k]).max()/peaks[0]['amplitude'])
        rows.append(dict(profile_x_m=x,local_vertical_cover_thickness_m=d,H0_first_cover_peak=peaks[0],H0_upper_peak=peaks[1],
            observed_amplitude_ratio=ratio,observed_amplitude_ratio_db=float(20*np.log10(ratio)),
            normal_incidence_extra_roundtrip_ratio=prediction,normal_incidence_ratio_db=float(20*np.log10(prediction)),
            local_vertical_mud_thickness_m=mud_d,mud_double_pass_propagation_loss_db=float(-20/np.log(10)*2*alpha_mud*mud_d),
            observed_basal_delta_to_first_cover_ratio=basal_ratio,normal_incidence_basal_to_first_cover_ratio=basal_prediction))
    thickness_change=rows[0]['local_vertical_cover_thickness_m']-rows[1]['local_vertical_cover_thickness_m']
    trend=dict(thickness_change_m=thickness_change,
        predicted_first_peak_advance_ns=2*thickness_change*ng/C*1e9,
        predicted_extra_roundtrip_peak_advance_ns=4*thickness_change*ng/C*1e9,
        observed_first_peak_advance_ns=rows[0]['H0_first_cover_peak']['time_ns']-rows[1]['H0_first_cover_peak']['time_ns'],
        observed_upper_peak_advance_ns=rows[0]['H0_upper_peak']['time_ns']-rows[1]['H0_upper_peak']['time_ns'])
    g=next(g for g in m['groups'] if g['id']=='fine_x170_H0_local')
    raw=root/'native'/g['id']/'profile.h5';source=sf.load_source(raw);r=sf.load_receiver(raw,'/rxs/rx1','Ez')
    freq=20e6+np.arange(501)*300000.;fq=sf.direct_frequency_response(source,r,freq,tail_taper_fraction=0)
    full=sf.reconstruct_time_response(fq,window='hann',zero_pad_factor=8,time_shift=0);t=full.time*1e9;gate=(t>=220)&(t<=245);crop=[]
    for stop in [260,300,340,360,380,600,1800,2280]:
        z=r.samples.copy();z[r.times>=stop*1e-9]=0
        cut=sf.reconstruct_time_response(sf.direct_frequency_response(source,replace(r,samples=z),freq,tail_taper_fraction=0),window='hann',zero_pad_factor=8,time_shift=0)
        j=np.flatnonzero(gate)[np.argmax(abs(cut.complex_bandpass[gate]))]
        crop.append(dict(diagnostic_keep_before_ns=stop,upper_window_ns=[220,245],
            complex_relative_L2=float(np.linalg.norm((cut.complex_bandpass-full.complex_bandpass)[gate])/np.linalg.norm(full.complex_bandpass[gate])),
            peak_time_ns=float(t[j]),peak_abs=float(abs(cut.complex_bandpass[j])/.025)))
    delta=2*208-g['tx_m'][0]-g['rx_m'][0]
    right=dict(profile_x_m=170,right_PML_inner_x_m=208,
        vacuum_geometric_roundtrip_lower_bound_ns=float(np.hypot(delta,g['tx_m'][1]-g['rx_m'][1])/C*1e9),
        limits='Geometric bound is physical, not a numerical precursor certificate. Finite-band inverse has noncausal sidelobes. Truncation is diagnostic only, no formal processing change.',crop_contribution_checks=crop)
    result=dict(status='READ_ONLY_LAYER_REPEAT_AND_LATE_CONTRIBUTION_CHECKS',materials_sha256=sha(package/'geometries/line9_materials_v0_2.json'),
        source_script_sha256=sha(__file__),reference_frequency_Hz=95e6,cover_group_index=ng,cover_attenuation_Np_m=alpha,
        mud_attenuation_Np_m=alpha_mud,
        cover_internal_surface_reflection_abs=float(rt),cover_mud_reflection_abs=float(rb),paired_station_checks=rows,thickness_trend=trend,right_boundary=right,
        limits='Local vertical, normal-incidence planar approximation: neglects lateral geometry, angle, spreading and frequency weighting. Agreement in trend/order is supporting evidence, not fitted ray closure. Formal source, records, SFCW windows and models unchanged.')
    (public/'mechanism_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['root','package','arrays','public']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();print(json.dumps(study(a.root,a.package,a.arrays,a.public)))
