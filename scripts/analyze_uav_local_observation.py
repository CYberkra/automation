"""Reuse verified air fields for declared source-weight sensitivity; no hardware fit."""
import argparse
import json
from pathlib import Path
import numpy as np
from analyze_uav_local_free_space import FREQ, transfer
from hs_capsule_identity import sha256


def main(study,out):
    if out.exists():
        raise ValueError('fresh output required')
    c=json.loads((study/'execution_contract.json').read_text('utf-8'))
    audit=json.loads((study/'completed_verification.json').read_text('utf-8'))
    if audit['status']!='PASS' or audit['contract_sha256']!=sha256(study/'execution_contract.json'):
        raise ValueError('verified native capsule required')
    g=next(g for g in c['groups'] if g['polarisation']=='y' and g['direction']=='forward')
    raw=Path(g['input']).with_suffix('.h5')
    v=next(v for v in audit['groups'] if v['id']==g['id'])
    if sha256(raw)!=v['raw_sha256']:
        raise ValueError('raw identity changed')
    # Frozen mechanism choices: equal current moment or equal free-space dipole
    # radiated power, identical 95 MHz anchor. Neither is actual accepted port power.
    weights={'current_moment':np.ones(501),'radiated_power':95e6/FREQ}
    time=np.linspace(0,120e-9,1201)
    window=np.hanning(501)
    kernel=np.exp(2j*np.pi*time[:,None]*FREQ[None,:])
    spectra={}; profiles={}
    for i,name in ((1,'broadside'),(2,'axial')):
        response,_,_=transfer(raw,i,'Ey',g['spacing_m'])
        spectra[name]=response
        for branch,w in weights.items():
            profiles[branch,name]=2*.3e6*np.real(kernel@(response*window*w))
    # One shared reference for every panel; no per-trace peak normalization.
    scale=max(np.max(abs(profiles['current_moment',name])) for name in spectra)
    metrics={}
    for branch in weights:
        metrics[branch]={'axial_to_broadside_band_energy_amplitude_dB':float(20*np.log10(
            np.linalg.norm(spectra['axial']*window*weights[branch])/
            np.linalg.norm(spectra['broadside']*window*weights[branch]))),
            'relative_peaks_shared_current_reference':{name:float(np.max(abs(profiles[branch,name]))/scale) for name in spectra}}
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei']
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
    axes[0].plot(FREQ/1e6,20*np.log10(abs(spectra['axial']/spectra['broadside'])))
    axes[0].axhline(0,color='gray',lw=.8)
    axes[0].set(title='同频轴向/侧向耦合；近场含全部项',xlabel='频率(MHz)',ylabel='幅比(dB)')
    for ax,branch,title in zip(axes[1:],weights,['等电流元、Hann频带','等自由空间辐射功率、Hann频带']):
        for name,label in [('broadside','侧向基线'),('axial','轴向基线')]:
            ax.plot(time*1e9,profiles[branch,name]/scale,label=label)
        ax.set(title=title,xlabel='重建时间(ns)',ylabel='同一固定参考下的相对幅值',xlim=(0,50),ylim=(-1.8,1.8))
        ax.legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('接收/频谱口径辅助：仅复用空气偶极；不是设备功率、端口S21或地下回波')
    fig.savefig(out/'observation_sensitivity.png',dpi=140);plt.close(fig)
    result=dict(status='COMPLETED_DECLARED_MECHANISM_SENSITIVITY',metrics=metrics,
                native_raw_sha256=sha256(raw),analysis_code_sha256=sha256(__file__),
                transfer_code_sha256=sha256(Path(__file__).parent/'analyze_uav_local_free_space.py'),
                source_weights={'current_moment':'1','radiated_power':'95MHz/f'},
                spectral_window='Hann501;20-170MHz;df0.3MHz',same_shared_display_scale=True,
                calls_solver=False,field_fitted=False,
                scope='Air-only current-moment/radiated-power hypothesis sensitivity; not port model, relief, deep-layer or field validation')
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();main(a.study,a.out)
