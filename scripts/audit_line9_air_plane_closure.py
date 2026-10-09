"""Raw reread and explicit spatial/time matrix sums; no FFT transport reuse."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha
from gprMax.toolboxes.SFCW import processing as sf

C=299792458.;EPS=8.8541878128e-12;MU=1/(EPS*C*C)


def taper(n,alpha):
    if alpha==0:return np.ones(n)
    x=np.arange(n)/(n-1);v=np.ones(n)
    lo=x<alpha/2;hi=x>1-alpha/2
    v[lo]=(1+np.cos(np.pi*(2*x[lo]/alpha-1)))/2
    v[hi]=(1+np.cos(np.pi*(2*x[hi]/alpha-2/alpha+1)))/2
    return v


def angular(spec,f,x,sector=.9,roll=.25):
    n=512;loc=(n-41)//2+np.arange(41);m=np.arange(n);signed=np.where(m<n//2,m,m-n)
    # Direct forward spatial sums from 41 samples; no FFT or zero-pad array.
    forward=np.exp(-2j*np.pi*loc[:,None]*m[None,:]/n)
    backward=np.conj(forward).T/n
    q=2*np.pi*signed/(n*.5);k=2*np.pi*f[:,None]/C;valid=abs(q)[None,:]<=sector*k
    ky=np.sqrt(np.maximum(k*k-q[None,:]**2,0));imp=np.zeros_like(ky)
    imp[valid]=(2*np.pi*f[:,None]*MU*np.ones((1,n)))[valid]/ky[valid]
    ratio=abs(q)[None,:]/k
    if roll==0:weight=valid.astype(float)
    else:weight=(1+np.cos(np.pi*np.clip((ratio-sector+roll)/roll,0,1)))/2
    sp=np.einsum('fci,ij->fcj',spec*taper(41,.25),forward)
    up1=(sp[:,0]+imp*sp[:,1])/2*weight;down1=(sp[:,0]-imp*sp[:,1])/2*weight
    up2=(sp[:,2]+imp*sp[:,3])/2*weight;down2=(sp[:,2]-imp*sp[:,3])/2*weight
    predicted_up=up1*np.exp(-4j*ky);predicted_down=down1*np.exp(4j*ky)
    vectors=[sp[:,2]*weight,predicted_up+predicted_down,predicted_up,predicted_down,up2,down2]
    return np.stack([v@backward for v in vectors]+[spec[:,2]*taper(41,.25)],axis=-1)


def metric(z):
    norm=np.linalg.norm;target,pred,pu,pd,au,ad,raw=[z[:,:,i] for i in range(7)]
    cosine=lambda a,b:float(abs(np.vdot(a,b))/(norm(a)*norm(b)))
    return {'total_transport_relative_error':float(norm(pred-target)/norm(target)),
            'up_only_vs_total_relative_error':float(norm(pu-target)/norm(target)),
            'up_to_up_transport_relative_error':float(norm(pu-au)/norm(au)),
            'up_to_up_complex_cosine':cosine(pu,au),'up_to_total_complex_cosine':cosine(pu,target),
            'upper_down_over_up_norm':float(norm(ad)/norm(au)),'predicted_down_over_up_norm':float(norm(pd)/norm(pu)),
            'sector_over_raw_relative_change':float(norm(target-raw)/norm(raw))}


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    c=read(a.smooth/'contract.json');r=read(a.smooth/'analysis.json');hard=read(a.hard/'analysis.json');parent=read(a.parent/'analysis.json');m=read(a.package/'manifest.json')
    assert sha(a.smooth/'contract.json')==r['contract_sha256'] and sha(a.profiles)==r['arrays_sha256']
    assert sha(a.raw)==c['raw_sha256']==parent['native_sha256'] and sha(a.package/'manifest.json')==parent['manifest_sha256']
    assert sha(a.hard/'analysis.json')==c['hard_analysis_sha256']
    for n,d in c['helper_sha256'].items():assert sha(Path(__file__).with_name(n))==d
    assert sha(Path(__file__).with_name('analyze_line9_air_smooth_closure.py'))==r['script_sha256']
    source=sf.load_source(a.raw);f=np.array(c['frequency_Hz']);dt=parent['dt_s'];t=np.arange(20351)*dt
    with np.load(a.profiles) as h:
        x=h['x_m'];times=h['time_s'];spectra={name:h['spectra_'+name] for name,_ in c['variant_native_gates']}
        profiles={name+'_'+window:h[name+'_'+window] for name in spectra for window in ['hann','blackman']}
    gates={name:np.ones(len(t)) for name in spectra}
    for name,bounds in c['variant_native_gates'][1:]:
        mask=np.flatnonzero((t*1e9>=bounds[0])&(t*1e9<=bounds[1]));gate=np.zeros(len(t));gate[mask]=taper(len(mask),.25);gates[name]=gate
    rawerrors=[];event_errors=[];event_directions={b:[] for b in ['cover_lower','cover_middle','cover_upper','air_surface']}
    lookup={(p['band'],p['x_m']):p for p in m['probes']};tone=[0,250,500]
    K=np.exp(-2j*np.pi*f[tone,None]*t)*dt
    src=source.dt*(np.exp(-2j*np.pi*f[tone,None]*source.times)@source.samples)*.025
    with h5py.File(a.raw) as h:
        assert h['srcs/src1'].attrs['Position'][1]>c['plane_y_m'][1]
        for j,b in enumerate(['air_34','air_38']):
            for k,xx in enumerate(x):
                p=lookup[b,xx];rr=[h[f"rxs/rx{v['receiver_index']}"] for v in p['anchors']]
                e=rr[0]['Ez'][:-1];hx=(rr[0]['Hx'][:-1]+rr[0]['Hx'][1:]+rr[1]['Hx'][:-1]+rr[1]['Hx'][1:])/4
                for name,gate in gates.items():
                    for ci,v in [(j*2,e),(j*2+1,hx)]:
                        value=(K@(v*gate))/src;reference=spectra[name][tone,ci,k]
                        error=float(np.linalg.norm(value-reference)/np.linalg.norm(reference));assert error<1e-10;rawerrors.append(error)
        for row in hard['native_late_events']:
            for band,event in row['bands'].items():
                p=lookup[band,row['x_m']];rr=[h[f"rxs/rx{v['receiver_index']}"] for v in p['anchors']]
                e=rr[0]['Ez'][:-1];hx=(rr[0]['Hx'][:-1]+rr[0]['Hx'][1:]+rr[1]['Hx'][:-1]+rr[1]['Hx'][1:])/4
                hy=(rr[0]['Hy'][:-1]+rr[0]['Hy'][1:]+rr[2]['Hy'][:-1]+rr[2]['Hy'][1:])/4
                sy=e*hx;sx=-e*hy;lo,hi=read(a.hard/'contract.json')['native_exploratory_up_windows_ns'][band]
                ix=np.flatnonzero((t*1e9>=lo)&(t*1e9<=hi));j=ix[np.argmax(sy[ix])];assert j==event['index']
                mask=(t*1e9>=t[j]*1e9-4)&(t*1e9<=t[j]*1e9+4)
                sums=np.array([np.sum(sx[mask],dtype=np.longdouble),np.sum(sy[mask],dtype=np.longdouble)])*dt
                ref=np.array([event['signed8ns_Sx_J_m2'],event['signed8ns_Sy_J_m2']]);err=float(np.linalg.norm(sums-ref)/np.linalg.norm(ref));assert err<1e-12;event_errors.append(err)
                event_directions[band].append({'x_m':row['x_m'],'Sy_integral_positive':bool(sums[1]>0),'Sx_integral_positive':bool(sums[0]>0)})
    # Verify both planes and all points between them are native homogeneous air.
    geometry=a.package/'baseline/geometry.h5';assert sha(geometry)==m['baseline_sha256']['geometry.h5']
    with h5py.File(geometry) as h:
        data=h['data'];j0=round(c['plane_y_m'][0]/.025);j1=round(c['plane_y_m'][1]/.025)
        for xx in x:assert (data[round(xx/.025),j0:j1+1,0]==0).all()
    metric_errors=[];profile_errors=[];hard_gate_checks=[];core=(x>=160)&(x<=168)
    choose=(times*1e9>=250)&(times*1e9<=400);tt=times[choose]
    for name,sp in spectra.items():
        values=angular(sp,f,x)
        for window in ['hann','blackman']:
            # Explicit cosine window and direct inverse sum, no FFT/IFFT.
            u=np.arange(501)/500;w=.5-.5*np.cos(2*np.pi*u) if window=='hann' else .42-.5*np.cos(2*np.pi*u)+.08*np.cos(4*np.pi*u);w/=w.mean()
            zz=(np.exp(2j*np.pi*tt[:,None]*f)@(values.reshape(501,-1)*w[:,None])/501).reshape(len(tt),41,7)
            actual=profiles[name+'_'+window][choose];error=float(np.linalg.norm(zz-actual)/np.linalg.norm(actual));assert error<1e-10;profile_errors.append(error)
            row=next(z for z in r['sensitivity_metrics'] if z['variant']==name and z['window']==window and all(z[k]==v for k,v in c['primary'].items()))
            for gate,bounds in c['gates_sfcw_ns'].items():
                select=(tt*1e9>=bounds[0])&(tt*1e9<=bounds[1]);recomputed=metric(zz[select][:,core]);ref=row['gates'][gate]
                error=max(abs(recomputed[k]-ref[k])/max(1,abs(ref[k])) for k in recomputed);assert error<1e-10;metric_errors.append(error)
            hv=angular(sp,f,x,roll=0)
            hz=(np.exp(2j*np.pi*tt[:,None]*f)@(hv.reshape(501,-1)*w[:,None])/501).reshape(len(tt),41,7)
            mm=metric(hz[:,core]);hard_gate_checks.append({'variant':name,'window':window,**mm})
    sign_counts={b:{'positive_vertical_windows':sum(z['Sy_integral_positive'] for z in vv),'positive_horizontal_windows':sum(z['Sx_integral_positive'] for z in vv),'total':len(vv)} for b,vv in event_directions.items()}
    out={'status':'PASS_RAW_PLANES_AND_EXPLICIT_SPATIAL_TEMPORAL_TRANSPORT_SUMS','script_sha256':sha(__file__),
         'smooth_analysis_sha256':sha(a.smooth/'analysis.json'),'hard_analysis_sha256':sha(a.hard/'analysis.json'),
         'native_three_tone_field_groups':len(rawerrors),'native_three_tone_max_relative_L2':max(rawerrors),
         'native_late_events':len(event_errors),'native_event_max_relative_error':max(event_errors),
         'native_direction_sign_counts':sign_counts,'explicit_primary_profiles':len(profile_errors),
         'explicit_profile_max_relative_L2':max(profile_errors),'explicit_primary_metric_rows':len(metric_errors),
         'explicit_metric_max_scaled_error':max(metric_errors),'hard_cut_timegate_checks':hard_gate_checks,
         'homogeneous_air_between_planes_verified':True,
         'limits':'Arithmetic/native aperture/air transport reference only; no exact full-space/ray/bounce/energy/site certification'}
    assert not a.out.exists();save=lambda p,v:p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');save(a.out,out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['hard_cut_timegate_checks','native_direction_sign_counts']}))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['smooth','profiles','raw','package','parent','hard','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
