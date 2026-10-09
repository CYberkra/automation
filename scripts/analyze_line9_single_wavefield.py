"""Independent native observer check and exact-band single-station SFCW report."""
import argparse
import json
from pathlib import Path
import h5py
import numpy as np
from audit_line9_v401_delivery import native_response, sha


def main(a):
    read=lambda p:json.loads(p.read_text('utf-8'))
    assert not a.out.exists()
    c=read(a.source/'execution_contract.json');v=read(a.source/'completed_verification.json')
    assert v['completed'] and v['contract_sha256']==sha(a.source/'execution_contract.json')
    row=v['groups'][0];assert row['observer_receiver_and_source_bitwise_equal'] and row['snapshot_count']==500
    raw=a.source/'observer_H0.h5';assert sha(raw)==row['native_sha256']
    old=a.baseline/'high_x19000_H0.h5';target=a.baseline/'high_x19000_H1.h5'
    with h5py.File(raw) as h,h5py.File(old) as b:
        for key in ['rxs/rx1/Ez','srcs/src1/excitation/samples']:
            assert h[key].dtype==b[key].dtype==np.float64 and h[key][:].tobytes()==b[key][:].tobytes()
        rx=h['rxs/rx1/Ez'][:];dt=h.attrs['dt'];assert len(rx)==20352
        times=np.arange(len(rx))*dt+h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
    p=read(a.previous/'analysis.json');pa=read(a.previous/'independent_audit.json')
    assert pa['analysis_sha256']==sha(a.previous/'analysis.json') and pa['status'].startswith('PASS')
    receipts={r['id']:r['native_sha256'] for r in p['native']}
    z0=native_response(raw,.025,row['native_sha256']);zold=native_response(old,.025,receipts['high_x19000_H0'])
    z1=native_response(target,.025,receipts['high_x19000_H1'])
    assert z0.tobytes()==zold.tobytes()
    spectra=np.column_stack([z1,z0,z1-z0]);f=20e6+np.arange(501)*300000.
    t=np.arange(4008)/(4008*300000.);k=(t*1e9>=300)&(t*1e9<=450);metrics={}
    a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    for window,w in [('hann',np.hanning(501)),('blackman',np.blackman(501))]:
        w/=w.mean();profiles=np.exp(2j*np.pi*t[k,None]*f)@(spectra*w[:,None])/501
        padded=np.zeros((4008,3),complex);padded[:501]=spectra*w[:,None]
        check=np.fft.ifft(padded,axis=0)*8*np.exp(2j*np.pi*20e6*t[:,None])
        err=float(np.linalg.norm(check[k]-profiles)/np.linalg.norm(profiles));assert err<1e-9
        peaks=[float(t[k][abs(profiles[:,j]).argmax()]*1e9) for j in range(3)]
        metrics[window]=dict(H1_H0_delta_peaks_ns=peaks,independent_inverse_relative_L2=err,observer_spectrum_bitwise_equal=True)
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        kk=(times*1e9>=250)&(times*1e9<=450)
        axes[0,0].plot(times[kk]*1e9,rx[kk],lw=.8,color='#333333')
        axes[0,0].set(title='190m H0原生A-scan：与旧结果逐位一致',xlabel='原生时间 / ns（含源时间参考）',ylabel='Ez / V/m')
        limit=float(abs(profiles.real).max())
        im=axes[0,1].imshow(profiles.real,cmap='gray',vmin=-limit,vmax=limit,aspect='auto',interpolation='nearest',extent=[-.5,2.5,t[k][-1]*1e9,t[k][0]*1e9])
        axes[0,1].set_xticks([0,1,2],['既有H1总场','新H0＋快照','底砂差场H1−H0'])
        axes[0,1].set(title='同一站位的场量列；共同物理灰度，非空间测线',ylabel='SFCW时间 / ns')
        fig.colorbar(im,ax=axes[0,1],label='(V/m)/(A·m)')
        for j,label in enumerate(['既有H1总场','新H0（底砂已移除）','底砂对比差场']):
            axes[1,0].plot(t[k]*1e9,abs(profiles[:,j]),label=label)
            axes[1,1].plot(t[k]*1e9,profiles[:,j].real,label=label)
        for ax,title,ylabel in [(axes[1,0],'SFCW复包络（先复数差分再取模）','复包络 / (V/m)/(A·m)'),(axes[1,1],'SFCW带符号实部','实部 / (V/m)/(A·m)')]:
            ax.set(title=title,xlabel='SFCW时间 / ns',ylabel=ylabel);ax.legend(fontsize=8);ax.grid(alpha=.15)
            ax.axvspan(345.6180971390552,374.6081170991351,color='#d7e5dd',alpha=.25)
        fig.suptitle(f'单站快照批接收验收 / 190m / AGL约8m / 4.0.1原生FP64 / 精确501点 / {window}\n无AGC、插值、幅相拟合；H0不是噪声，差场不是纯一次波。只有一站，不代表完整B-scan。')
        fig.savefig(a.out/f'single_station_{window}.png',dpi=130);plt.close(fig)
    result=dict(status='PASS_LOCAL_NATIVE_BITWISE_DFT_AND_DIRECT_INVERSE',script_sha256=sha(__file__),
        contract_sha256=sha(a.source/'execution_contract.json'),verification_sha256=sha(a.source/'completed_verification.json'),
        old_native_sha256=sha(old),new_native_sha256=sha(raw),observer_native_and_spectrum_bitwise_equal=True,
        source_previous_analysis_sha256=sha(a.previous/'analysis.json'),metrics=metrics,
        limits='One station native-observer passivity and CPU arithmetic only; no unique path attribution or full spatial scan.')
    (a.out/'analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['source','baseline','previous','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
