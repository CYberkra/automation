"""Independent raw E/H reread and explicit spatial/temporal sums."""
import argparse,json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def taper(n,alpha):
    x=np.arange(n)/(n-1);v=np.ones(n)
    if alpha:
        low=x<alpha/2;high=x>1-alpha/2
        v[low]=(1+np.cos(np.pi*(2*x[low]/alpha-1)))/2
        v[high]=(1+np.cos(np.pi*(2*x[high]/alpha-2/alpha+1)))/2
    return v


def main(a):
    read=lambda p:json.loads(Path(p).read_text('utf-8'))
    result=read(a.analysis/'analysis.json');c=read(a.analysis/'contract.json');m=read(a.package/'manifest.json')
    assert result['arrays_sha256']==sha(a.arrays) and result['contract_sha256']==sha(a.analysis/'contract.json')
    assert c['raw_sha256']==sha(a.raw) and c['manifest_sha256']==sha(a.package/'manifest.json')
    with np.load(a.arrays) as h:data={k:h[k] for k in h.files}
    dt=m['dt_s'];t=data['time_s'];f=20e6+np.arange(501)*300000.;np.testing.assert_array_equal(f,data['frequency_Hz'])
    checked=0;previous=0
    with h5py.File(a.raw) as raw,h5py.File(a.reference) as old:
        src=raw['srcs/src1/excitation/samples'][:];assert src.tobytes()==old['srcs/src1/excitation/samples'][:].tobytes()
        source_offset=float(raw['srcs/src1/excitation'].attrs['TimeSampleOffset'])
        assert source_offset==.5*dt
        for name,r in old['rxs'].items():
            for field in r:assert raw['rxs/'+name+'/'+field][:].tobytes()==r[field][:].tobytes();previous+=1
        for i,p in enumerate(m['probes'][287:]):
            r0,r1,r2=[raw[f"rxs/rx{q['receiver_index']}"] for q in p['anchors']]
            e=r0['Ez'][:];h0=r0['Hx'][:];hb=r1['Hx'][:];v0=r0['Hy'][:];vl=r2['Hy'][:]
            expected=[e[:-1],(h0[:-1]+hb[:-1]+h0[1:]+hb[1:])/4,(v0[:-1]+vl[:-1]+v0[1:]+vl[1:])/4]
            for key,v in zip(['Ez','Hx','Hy'],expected):np.testing.assert_allclose(data[key][i],v,rtol=1e-14,atol=1e-16);checked+=1
    assert previous==1436 and checked==369
    native=np.stack([data[key][i*41:(i+1)*41] for i in range(3) for key in ['Ez','Hx']])
    pick=np.array([0,83,250,417,500]);freq=f[pick];ss=dt*(np.exp(-2j*np.pi*freq[:,None]*(source_offset+np.arange(len(src))*dt))@src)
    for name,bounds in c['native_variants']:
        win=np.ones(len(t))
        if bounds:
            win[:]=0;ix=np.flatnonzero((t*1e9>=bounds[0])&(t*1e9<=bounds[1]));win[ix]=taper(len(ix),.25)
        z=np.einsum('ft,cxt->fcx',np.exp(-2j*np.pi*freq[:,None]*t),native*win)*dt/ss[:,None,None]/.025
        np.testing.assert_allclose(z,data['spectra_'+name][pick],rtol=2e-10,atol=1e-12)
    # Independent complex material expression and explicit DFT matrices.
    omega=2*np.pi*f;eps0=8.8541878128e-12;light=299792458.;mu=1/(eps0*light*light)
    wt=omega*6.4567e-9;eps=11+.5*(1-1j*wt)/(1+wt*wt)-1j*.001/(omega*eps0)
    n=512;pos=(n-41)//2+np.arange(41);modes=np.arange(n);signed=np.where(modes<n//2,modes,modes-n)
    forward=np.exp(-2j*np.pi*pos[:,None]*modes[None,:]/n);backward=np.conj(forward).T/n
    q=2*np.pi*signed/(n*.5);k0=omega[:,None]/light;ky=np.sqrt(k0*k0*eps[:,None]-q[None,:]**2)
    weight=(1+np.cos(np.pi*np.clip((abs(q)[None,:]/k0-.9+.25)/.25,0,1)))/2
    valid=abs(q)[None,:]<=.9*k0;imp=np.zeros_like(ky);imp[valid]=(omega[:,None]*mu*np.ones((1,n)))[valid]/ky[valid]
    samples=np.unique(np.r_[np.arange(0,4008,127),np.flatnonzero((data['sfcw_time_s']*1e9>=220)&(data['sfcw_time_s']*1e9<=330))[::9]])
    temporal=np.exp(2j*np.pi*data['sfcw_time_s'][samples,None]*f)
    profile_checks=[]
    for variant,pair in [('full',[0,1]),('full',[1,2]),('full',[0,2]),('native_220_380',[0,2])]:
        spec=data['spectra_'+variant];sp=(spec*taper(41,.25))@forward;i,j=pair;dy=(j-i)*1.5
        up=(sp[:,2*i]+imp*sp[:,2*i+1])*weight/2;down=(sp[:,2*i]-imp*sp[:,2*i+1])*weight/2
        au=(sp[:,2*j]+imp*sp[:,2*j+1])*weight/2;ad=(sp[:,2*j]-imp*sp[:,2*j+1])*weight/2
        phase=np.where(valid,ky,0)*dy;pu=up*np.exp(-1j*phase);pd=down*np.exp(1j*phase)
        vals=[sp[:,2*j]*weight,pu+pd,pu,pd,au,ad]
        explicit=np.stack([v@backward for v in vals]+[spec[:,2*j]*taper(41,.25)],axis=-1)
        for window in ['hann','blackman']:
            ind=np.arange(501);w=.5-.5*np.cos(2*np.pi*ind/500) if window=='hann' else .42-.5*np.cos(2*np.pi*ind/500)+.08*np.cos(4*np.pi*ind/500);w/=w.mean()
            v=(temporal@(explicit.reshape(501,-1)*w[:,None])/501).reshape(len(samples),41,7)
            actual=data[variant+'_'+str(i)+str(j)+'_'+window][samples]
            error=float(np.linalg.norm(v-actual)/np.linalg.norm(v));assert error<2e-10
            profile_checks.append({'variant':variant,'pair':pair,'window':window,'relative_error':error,'time_samples':len(samples)})
    core=(data['x_m']>=160)&(data['x_m']<=168);metric_count=0
    for row in result['metrics']:
        if any(row[key]!=value for key,value in c['primary'].items()):continue
        z=data[row['variant']+'_'+str(row['pair'][0])+str(row['pair'][1])+'_'+row['window']]
        for gate,bounds in c['gates_sfcw_ns'].items():
            mask=(data['sfcw_time_s']*1e9>=bounds[0])&(data['sfcw_time_s']*1e9<=bounds[1]);v=z[mask][:,core]
            target,pred,pu,pd,au,ad,unfiltered=[v[:,:,k] for k in range(7)]
            norm=lambda a:np.sqrt(np.sum(abs(a)**2))
            own={'total_relative_error':norm(pred-target)/norm(target),'up_relative_error':norm(pu-au)/norm(au),
                 'down_relative_error':norm(pd-ad)/norm(ad),'up_complex_cosine':abs(np.sum(pu.conj()*au))/(norm(pu)*norm(au)),
                 'down_complex_cosine':abs(np.sum(pd.conj()*ad))/(norm(pd)*norm(ad)),
                 'destination_down_up_norm_ratio':norm(ad)/norm(au),'sector_raw_relative_change':norm(target-unfiltered)/norm(unfiltered)}
            for key,value in own.items():np.testing.assert_allclose(value,row['gates'][gate][key],rtol=1e-12,atol=1e-13)
            metric_count+=1
    v={'status':'PASS_RAW369_FIELDS_PREVIOUS1436_CHANNELS_EXPLICIT_COVER_SUMS','analysis_sha256':sha(a.analysis/'analysis.json'),
       'auditor_sha256':sha(__file__),'native_collocated_fields':checked,'previous_bitwise_channels':previous,
       'five_tone_field_groups':4*6*41,'explicit_spatial_temporal_profiles':profile_checks,'primary_gate_metric_rows':metric_count,
       'source_time_offset_s':source_offset,'source_bitwise_equal':True,'new_solver_runs':0,'limits':'Independent formula and finite-array computations; not all-q sampling/physical/site/bounce certification.'}
    assert not a.out.exists();a.out.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8');print(json.dumps(v))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['analysis','arrays','package','raw','reference','out']:p.add_argument('--'+n,type=Path,required=True)
    main(p.parse_args())
