"""Independent native event reread and adaptive planar vector-field integration."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import jv
from hs_capsule_identity import sha256 as sha

C=299792458.;EPS=8.8541878128e-12;MU=1/(EPS*C*C)


def epsilon(m,f):
    b=m['base'];e=b['relative_permittivity']+b['electric_conductivity_s_per_m']/(2j*np.pi*f*EPS)
    for p in m.get('poles',[]):e+=p['relative_permittivity_difference']/(1+2j*np.pi*f*p['relaxation_time_s'])
    return e


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'));s=read(a.surface/'analysis.json');p=read(a.planar/'analysis.json');c=read(a.planar/'contract.json');m=read(a.package/'manifest.json')
    assert sha(a.raw)==s['native_sha256'] and sha(a.arrays)==p['arrays_sha256']
    assert sha(Path(__file__).with_name('analyze_line9_surface_rereflection.py'))==s['script_sha256']
    assert sha(Path(__file__).with_name('analyze_line9_planar_dimensional_returns.py'))==p['script_sha256']
    for name,d in c['helper_sha256'].items():assert sha(Path(__file__).with_name(name))==d
    db=read(a.package/'baseline/materials.json')['materials'];cover=db['material_001_cover'];mud=db['material_002_mudstone']
    f0=95e6;step=f0*1e-5;n0=np.sqrt(epsilon(cover,f0));ng=float((n0+f0*(np.sqrt(epsilon(cover,f0+step))-np.sqrt(epsilon(cover,f0-step)))/(2*step)).real)
    assert abs(ng-s['cover_group_index95'])<1e-8
    event_error=[];signs=0;lookup={(z['band'],z['x_m']):z for z in m['probes']}
    with h5py.File(a.raw) as h:
        dt=float(h.attrs['dt']);t=np.arange(20351)*dt*1e9
        for row in s['rows']:
            values={}
            for band in ['cover_upper','cover_middle','air_surface']:
                z=lookup[band,row['x_m']];rr=[h[f"rxs/rx{q['receiver_index']}"] for q in z['anchors']]
                e=rr[0]['Ez'][:-1];hx=(rr[0]['Hx'][:-1]+rr[0]['Hx'][1:]+rr[1]['Hx'][:-1]+rr[1]['Hx'][1:])/4
                values[band]=e*hx
            for key,band,bounds,sgn in [('upper_up','cover_upper',[160,210],1),
                ('upper_later_down','cover_upper',[row['upper_up']['time_ns']+2,row['upper_up']['time_ns']+15],-1),
                ('middle_first_up','cover_middle',[120,175],1),
                ('middle_later_down','cover_middle',[row['upper_up']['time_ns']+15,row['upper_up']['time_ns']+70],-1),
                ('air_first_up','air_surface',[row['upper_up']['time_ns'],row['upper_up']['time_ns']+15],1)]:
                v=values[band];ids=np.flatnonzero((t>=bounds[0])&(t<=bounds[1]));j=max(ids,key=lambda j:sgn*v[j]);ref=row[key]
                assert j==ref['index'] and t[j]==ref['time_ns']
                win=(t>=t[j]-2)&(t<=t[j]+2);area=float(np.sum(v[win],dtype=np.longdouble)*dt)
                np.testing.assert_allclose(area,ref['local4ns_signed_J_m2'],rtol=1e-12,atol=1e-24)
                event_error.append(float(abs(v[j]-ref['Sy_W_m2'])/abs(v[j])))
                if key!='air_first_up':assert sgn*area>0;signs+=1
    with np.load(a.arrays) as h:f=h['frequency_Hz'];response=h['response'];time=h['time_s']
    adaptive=[];admittance=[]
    for k,thickness in enumerate(c['cover_thicknesses_m']):
        for fi in [0,250,500]:
            freq=f[fi];k0=2*np.pi*freq/C;height=sum(c['air_source_rx_heights_m']);dx=c['offset_m'];e1=epsilon(cover,freq);e2=epsilon(mud,freq)
            def root(z):
                z=np.sqrt(z+0j);return -z if z.imag>0 else z
            def integrand(theta,evan=False):
                q=k0*np.cosh(theta) if evan else k0*np.sin(theta)
                ky0=-1j*k0*np.sinh(theta) if evan else k0*np.cos(theta)
                ky1=root(k0*k0*e1-q*q);ky2=root(k0*k0*e2-q*q)
                phase=np.exp(-2j*ky1*thickness);components=[]
                for y1,y2 in [(ky1,ky2),(ky1/e1,ky2/e2)]:
                    r=(ky0-y1)/(ky0+y1);b=(y1-y2)/(y1+y2);once=(1-r*r)*b*phase;loop=-r*b*phase
                    # Stable geometric remainder, independently from full-r-once-twice.
                    z=np.array([r,once,once*loop,once*loop**2/(1-loop)])
                    tangent=np.tan(ky1*thickness);Y=y1*(y2+1j*y1*tangent)/(y1+1j*y2*tangent)
                    refl=(ky0-Y)/(ky0+Y);admittance.append(float(abs(refl-z.sum())))
                    components.append(z)
                wave=np.exp(-1j*ky0*height)
                line=wave*(1j if evan else 1)*np.cos(q*dx)*components[0]*(-freq*MU)
                j0=jv(0,q*dx);j2=jv(2,q*dx)
                angular=(j0-j2)*components[0]-(ky0/k0)**2*(j0+j2)*components[1]
                measure=1j*k0*np.cosh(theta) if evan else k0*np.sin(theta)
                point=wave*measure*angular*(-2*np.pi*freq*MU/(8*np.pi))
                return np.stack([line,point])
            prop,_=quad_vec(integrand,0,np.pi/2,epsabs=1e-16,epsrel=1e-12)
            evan,_=quad_vec(lambda u:integrand(u,True),0,np.arcsinh(40/(k0*height)),epsabs=1e-16,epsrel=1e-12)
            exact=prop+evan;reference=response[fi,k];err=abs(exact-reference)/abs(exact)
            assert np.max(err)<1e-5,err;adaptive.append({'thickness_m':thickness,'frequency_Hz':freq,'relative_errors_dimension_component':err.tolist()})
    assert max(admittance)<1e-10
    # Direct complex summation on both named event windows (no FFT call).
    first=(time*1e9>=150)&(time*1e9<=250);late=(time*1e9>=300)&(time*1e9<=450);ids=np.flatnonzero(first|late)
    K=np.exp(2j*np.pi*time[ids,None]*f)/501;metric_errors=[];fullband_changes=[]
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();z=(K@(response.reshape(501,-1)*w[:,None])).reshape(len(ids),3,2,4)
        for k in range(3):
            fullratios=[]
            for dim,label in enumerate(['2D_TE_line','3D_transverse_point']):
                expected=np.linalg.norm(z[late[ids],k,dim,2])/np.linalg.norm(z[first[ids],k,dim,1])
                actual=p['metrics'][window][k]['dimensions'][label]['twice_late_over_once_first']
                error=float(abs(expected-actual)/abs(expected));assert error<1e-10;metric_errors.append(error)
                # Parseval removes any boundary clipping by named time windows.
                full=np.linalg.norm(response[:,k,dim,2]*w)/np.linalg.norm(response[:,k,dim,1]*w)
                actual=p['metrics'][window][k]['dimensions'][label]['twice_full_over_once_full']
                error=float(abs(full-actual)/abs(full));assert error<1e-10;metric_errors.append(error);fullratios.append(float(full))
            fullband_changes.append({'window':window,'thickness_m':c['cover_thicknesses_m'][k],
                'weighted_fullband_twice_over_once_by_dimension':fullratios,
                'point_over_line_relative_return_change_dB':float(20*np.log10(fullratios[1]/fullratios[0]))})
    result={'status':'PASS_NATIVE205_EVENTS_AND_ADAPTIVE_PLANAR9_CASES','script_sha256':sha(__file__),
            'surface_analysis_sha256':sha(a.surface/'analysis.json'),'planar_analysis_sha256':sha(a.planar/'analysis.json'),
            'native_events':len(event_error),'maximum_native_peak_relative_error':max(event_error),'signed_windows_checked':signs,
            'group_index_finite_difference_error':abs(ng-s['cover_group_index95']),
            'adaptive_case_count':len(adaptive),'adaptive_cases':adaptive,'maximum_adaptive_component_error':max(np.max(r['relative_errors_dimension_component']) for r in adaptive),
            'maximum_reflection_admittance_error':max(admittance),'relative_return_metric_checks':len(metric_errors),'maximum_metric_relative_error':max(metric_errors),
            'fullband_parseval_return_changes':fullband_changes,
            'limits':'Arithmetic/native provenance/planar integrals only, no nonflat3D/site/unique path certification.'}
    assert not a.out.exists();a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k not in ['adaptive_cases','fullband_parseval_return_changes']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['surface','planar','package','raw','arrays','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
