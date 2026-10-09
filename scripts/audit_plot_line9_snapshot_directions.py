"""Independently reconstruct all Rx-aperture direction estimates from raw H5; plot diagnostics."""
import argparse
import json
from pathlib import Path, PureWindowsPath
import hashlib
import h5py
import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main(a):
    read = lambda p: json.loads(Path(p).read_text('utf-8'))
    analysis = read(a.out/'analysis.json')
    assert analysis['script_sha256'] == sha(Path(__file__).with_name('analyze_line9_snapshot_directions.py'))
    c = read(a.source/'execution_contract.json');v = read(a.source/'completed_verification.json')
    assert sha(a.source/'execution_contract.json') == analysis['contract_sha256'] == v['contract_sha256']
    assert sha(a.source/'completed_verification.json') == analysis['certificate_sha256']
    assert analysis['solver_runs'] == 0
    for row in analysis['runtime_snapshot_source_identities']:
        assert sha(a.runtime_sources/row['name']) == row['sha256']
    cache = a.cache/'aligned_probes_and_maps.npz'
    assert sha(cache) == analysis['cache_sha256']
    with np.load(cache, allow_pickle=False) as h:
        z = {k:h[k] for k in h.files}
    xx, yy, times = z['x_m'],z['y_m'],z['time_ns']
    rx = np.array(analysis['receiver_m'])
    xi = np.flatnonzero(abs(xx-rx[0]) <= 1+1e-12)
    yi = np.flatnonzero(abs(yy-rx[1]) <= 1+1e-12)
    records = sorted([r for r in v['groups'][0]['snapshots'] if r['roi']=='local_view'], key=lambda r:r['iteration'])
    fields, te, th = [],[],[]
    for row in records:
        p = a.source/'profile_snaps'/PureWindowsPath(row['file']).name
        assert sha(p) == row['sha256']
        with h5py.File(p,'r') as h:
            fields.append([h[k][xi[0]:xi[-1]+1,yi[0]:yi[-1]+1,0] for k in ('Ez','Hx','Hy')])
            te.append(h.attrs['time']);th.append(h.attrs['magnetic_time'])
    fields = np.asarray(fields)
    np.testing.assert_allclose(np.array(te)[1:-1]*1e9, times, rtol=0, atol=1e-12)
    hxhy = fields[:,1:]
    # Independent interpolation: numpy interp per sample, or barycentric polynomial.
    target = np.array(te)[1:-1]
    flat = hxhy.reshape(250,-1)
    interpolated = np.stack([np.interp(target,th,flat[:,j]) for j in range(flat.shape[1])],axis=1).reshape((248,*hxhy.shape[1:]))
    alpha = (np.array(te)[1:-1]-np.array(th)[1:-1])/(np.array(th)[2:]-np.array(th)[1:-1])
    aa = alpha[:,None,None,None]
    polynomial = aa*(aa-1)/2*hxhy[:-2] + (1-aa**2)*hxhy[1:-1] + aa*(aa+1)/2*hxhy[2:]
    backward = (1+aa)*hxhy[1:-1]-aa*hxhy[:-2]
    by_method = {'forward_linear':interpolated,'quadratic':polynomial,'backward_linear':backward}
    ez = fields[1:-1,0]
    directions = []
    for row in analysis['aperture_gate_directions']:
        rad = row['radius_m'];lo,hi = row['gate_ns']
        mask = (abs(xx[xi,None]-rx[0]) <= rad+1e-12) & (abs(yy[None,yi]-rx[1]) <= rad+1e-12)
        selected = (times>=lo)&(times<=hi)
        h = by_method[row['alignment']]
        sx,sy = (-ez*h[:,1]),(ez*h[:,0])
        actual = np.array([sx[selected][:,mask].mean(),sy[selected][:,mask].mean()])
        expected = [row['mean_Sx_W_m2'],row['mean_Sy_W_m2']]
        np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=1e-17)
        angle = float(np.degrees(np.arctan2(actual[1],actual[0])))
        assert abs(angle-row['angle_from_positive_x_deg']) < 1e-9
        assert mask.sum() == row['spatial_samples'] and selected.sum()==row['time_samples']
        # A separate segment-line solve for the extrapolated contact.
        contacts=[]
        if np.all(actual>0):
            slope=actual[1]/actual[0]
            for j in range(len(xx)-1):
                ground_slope=(z['surface_m'][j+1]-z['surface_m'][j])/(xx[j+1]-xx[j])
                if ground_slope==slope:continue
                contact=(z['surface_m'][j]-ground_slope*xx[j]-rx[1]+slope*rx[0])/(slope-ground_slope)
                if xx[j] <= contact <= xx[j+1] and contact < rx[0]:contacts.append(contact)
        np.testing.assert_allclose(contacts,row['straight_air_backprojection_surface_x_m'],rtol=0,atol=1e-9)
        directions.append(row)
    assert len(directions)==27
    # All saved strip samples checked directly against the native arrays at sparse times;
    # scalar-loop indices avoid sharing the producer's vectorized indexing.
    audited_strip_steps=[]
    for k in [49,99,149,169,179,199,229]:
        row=records[k+1];name=PureWindowsPath(row['file']).name
        following=PureWindowsPath(records[k+2]['file']).name
        with h5py.File(a.source/'profile_snaps'/name) as h,h5py.File(a.source/'profile_snaps'/following) as n:
            for strip in range(8):
                for i in range(0,350,17):
                    y=int(z['iy'][strip,i]);e=h['Ez'][i,y,0]
                    weight=(h.attrs['time']-h.attrs['magnetic_time'])/(n.attrs['magnetic_time']-h.attrs['magnetic_time'])
                    hx=(1-weight)*h['Hx'][i,y,0]+weight*n['Hx'][i,y,0]
                    hy=(1-weight)*h['Hy'][i,y,0]+weight*n['Hy'][i,y,0]
                    np.testing.assert_allclose(z['fields'][k,:,strip,i],[e,hx,hy,-e*hy,e*hx],rtol=2e-11,atol=1e-17)
        audited_strip_steps.append(row['iteration'])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm,LogNorm
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    names=['空气中y=32.05m','空气中y=34.05m','空气中y=36.05m','空气中y=38.05m',
           '地表空气侧约+0.35m','覆盖层内距底约1/4厚度','覆盖层内距底约1/2厚度','覆盖层内距底约3/4厚度']
    fig,axes=plt.subplots(4,2,figsize=(15,14),sharex=True,sharey=True,layout='constrained')
    for k,ax in enumerate(axes.ravel()):
        im=ax.pcolormesh(xx,times,z['fields'][:,4,k],cmap='gray',norm=SymLogNorm(1e-6,vmin=-.01,vmax=.01),shading='nearest')
        ax.axvline(rx[0],color='red',ls='--',lw=.8)
        ax.set(title=names[k],ylim=[400,100],xlabel='模型x / m',ylabel='原生时间 / ns')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='总场上向分量 Sy估计 / W/m²；正=向上，负=向下；超限截色',shrink=.65)
    fig.suptitle('190m高损耗、无底砂对照：原生波场时空观察\n共同物理尺度；磁场前向线性对齐；不是移动天线B-scan、不是各阶波能量分解',fontsize=14)
    fig.savefig(a.out/'native_signed_vertical_transport.png',dpi=140);plt.close(fig)
    # Air only: avoid material-interface interpolation in direction arrows.
    air=(yy[None,:]>z['surface_m'][:,None]+.3)&(yy[None,:]<40.3)
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained',sharex=True,sharey=True)
    late=[r for r in directions if r['gate_ns']==[330,350] and r['alignment']=='forward_linear' and r['radius_m']==.5][0]
    contact=late['straight_air_backprojection_surface_x_m'][0]
    for j,ax in enumerate(axes.ravel()):
        sx,sy=z['maps'][j];magnitude=np.hypot(sx,sy)
        im=ax.pcolormesh(xx,yy,np.where(air,magnitude,np.nan).T,cmap='gray_r',norm=LogNorm(vmin=1e-10,vmax=1e-4),shading='nearest')
        ax.plot(xx,z['surface_m'],color='green',label='真实地表');ax.plot(xx,z['bottom_m'],color='cyan',label='覆盖层底')
        qi=np.arange(5,350,10);qj=np.arange(5,315,10)
        mag=magnitude[np.ix_(qi,qj)]
        usable=air[np.ix_(qi,qj)]&(mag>1e-8)
        u=np.where(usable,sx[np.ix_(qi,qj)]/np.maximum(mag,1e-100),np.nan)
        w=np.where(usable,sy[np.ix_(qi,qj)]/np.maximum(mag,1e-100),np.nan)
        ax.quiver(xx[qi],yy[qj],u.T,w.T,color='tab:orange',angles='xy',scale_units='xy',scale=1.3,width=.003)
        ax.scatter([c['groups'][0]['tx_m'][0],rx[0]],[c['groups'][0]['tx_m'][1],rx[1]],c=['orange','red'],s=25)
        ax.axhspan(40.5,41.5,color='red',alpha=.15)
        if j==5:
            ax.plot([contact,rx[0]],[np.interp(contact,xx,z['surface_m']),rx[1]],'b--',lw=1.5,label='方向外推假设')
        ax.set(xlim=[150,176],ylim=[24,41.5],xlabel='模型x / m',ylabel='模型y / m',title=f'原生t={int(z["map_centres_ns"][j])}±4ns平均；{int(z["map_counts"][j])}帧')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='空气中总场 |平均E×H|估计 / W/m²；共同尺度',shrink=.65)
    fig.suptitle('原生晚波传播方向：磁场对齐后，橙色箭头仅表示总场方向\n蓝虚线仅为空气直线外推假设；没有分离反射次数，灰度饱和不等于能量占比',fontsize=14)
    fig.savefig(a.out/'native_air_direction_sequence.png',dpi=140);plt.close(fig)
    with h5py.File(a.source/'profile.h5') as h:
        native=h['rxs/rx1/Ez'][z['iterations']]
    fig,axes=plt.subplots(3,1,figsize=(13,10),layout='constrained')
    axes[0].plot(times,native,label='原生主Rx Ez');axes[0].plot(times,z['rx_traces'][:,0],'--',label='最近快照中心Ez（y差0.025m）')
    axes[0].set(xlim=[280,380],ylabel='Ez / V/m',title='独立原生Rx与最近快照点，未经源归一化或时间拟合');axes[0].legend()
    axes[1].plot(times,z['rx_traces'][:,3],label='Sx：正=向右');axes[1].plot(times,z['rx_traces'][:,4],label='Sy：正=向上')
    axes[1].axhline(0,color='black',lw=.5);axes[1].set(xlim=[280,380],ylabel='总场E×H估计 / W/m²',xlabel='原生时间 / ns');axes[1].legend()
    for method,marker in [('forward_linear','o'),('backward_linear','x'),('quadratic','s')]:
        for gate,color in [([150,230],'gray'),([300,370],'tab:green'),([330,350],'tab:blue')]:
            rows=[r for r in directions if r['alignment']==method and r['gate_ns']==gate]
            axes[2].plot([r['radius_m'] for r in rows],[r['angle_from_positive_x_deg'] for r in rows],marker=marker,color=color,label=f'{method} / {gate[0]}–{gate[1]}ns',lw=.7)
    axes[2].set(xlabel='Rx周围方形半宽 / m',ylabel='方向相对+x / °',title='三种时间对齐 × 三种孔径 × 三个事后时间窗；属于敏感性诊断')
    axes[2].legend(ncol=3,fontsize=8)
    fig.suptitle('190m高损耗H0：原生晚波的局部方向核查，无新求解',fontsize=14)
    fig.savefig(a.out/'native_rx_direction_sensitivity.png',dpi=140);plt.close(fig)
    proof={'status':'PASS_INDEPENDENT250_APERTURE_READS_27_DIRECTION_ESTIMATES_AND_STRIP_SAMPLES',
           'script_sha256':sha(__file__),'analysis_sha256':sha(a.out/'analysis.json'),
           'source_snapshots_independently_read':250,'direction_estimates':27,'strip_checked_iterations':audited_strip_steps,
           'limits':'Arithmetic and source identity audit only; no energy closure, unique path, physical threshold or site certification.'}
    (a.out/'independent_audit.json').write_text(json.dumps(proof,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(proof))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('source','runtime-sources','cache','out'):p.add_argument('--'+name,type=Path,required=True)
    main(p.parse_args())
