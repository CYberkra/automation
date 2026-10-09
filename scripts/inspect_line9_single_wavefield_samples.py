"""Read selected native snapshots locally; show late weak fields on one fixed scale."""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

    read = lambda p: json.loads(p.read_text('utf-8'))
    cert = read(a.public/'completed_verification.json')
    contract = read(a.public/'execution_contract.json')
    assert cert['completed'] and cert['contract_sha256'] == sha(a.public/'execution_contract.json')
    records = {Path(r['file']).name:r for r in cert['groups'][0]['snapshots']}
    samples = sorted(a.samples.glob('*.h5'))
    assert len(samples) == 20
    dt = contract['study_manifest']['dt_s']
    arrays, checked = {}, []
    for path in samples:
        r = records[path.name]
        assert sha(path) == r['sha256']
        roi = contract['study_manifest']['ROIs'][r['roi']]
        with h5py.File(path) as h:
            assert h.attrs['gprMax'] == '4.0.1' and h.attrs['iteration'] == r['iteration']
            assert abs(h.attrs['time']-r['iteration']*dt) < 1e-20
            assert abs(h.attrs['magnetic_time']-(r['iteration']-.5)*dt) < 1e-20
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],roi['shape'])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],roi['spacing'])
            np.testing.assert_allclose(h.attrs['origin'],roi['bounds'][:3],atol=1e-12,rtol=0)
            assert set(h) == {'Ex','Ey','Ez','Hx','Hy','Hz'}
            for field in h:
                x=h[field][:]
                assert x.dtype==np.float64 and list(x.shape)==roi['shape'] and np.isfinite(x).all()
                if field=='Ez':arrays[(r['roi'],r['iteration'])]=x[:,:,0].T
        checked.append(dict(file=path.name,sha256=r['sha256'],roi=r['roi'],iteration=r['iteration']))
    idle=read(a.samples/'remote_idle.json')
    assert idle['status']=='PASS_OWNED_PYTHON_EXITED' and not idle['owned_pids']
    with h5py.File(a.geometry) as h: data=h['data'][:,:,0]
    expected=contract['study_manifest']['groups'][0]['geometry_sha256']
    assert sha(a.geometry)==expected
    xx=(np.arange(len(data))+.5)*.025
    surface=np.max(np.where(data!=0,np.arange(data.shape[1])[None,:]+1,0),axis=1)*.025
    mudtop=np.max(np.where(data==2,np.arange(data.shape[1])[None,:]+1,0),axis=1)*.025
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    steps=[4216,4760,5100,5372,5508,5644,5780,5916,6052]
    # Fixed post hoc display scale, chosen to show the native receiver packet (~0.05 V/m).
    # Stronger fields clip; no frame normalization and no energy attribution.
    scale=.1;threshold=.0001
    fig,axes=plt.subplots(3,3,figsize=(15,10),layout='constrained',sharex=True,sharey=True)
    g=contract['groups'][0];keep=(xx>=145)&(xx<=180)
    for ax,j in zip(axes.ravel(),steps):
        im=ax.imshow(arrays[('local_view',j)],origin='lower',extent=[145,180,10,41.5],
            cmap='gray',norm=SymLogNorm(threshold,vmin=-scale,vmax=scale),aspect='auto')
        ax.plot(xx[keep],surface[keep],color='lime',lw=.8)
        ax.plot(xx[keep],mudtop[keep],color='cyan',lw=.8)
        ax.scatter([g['tx_m'][0],g['rx_m'][0]],[g['tx_m'][1],g['rx_m'][1]],c=['orange','red'],s=15)
        ax.axhspan(40.5,41.5,color='red',alpha=.16)
        ax.set(xlim=[145,180],ylim=[23.5,41.5],title=f'原生t={j*dt*1e9:.3f}ns')
    fig.supxlabel('模型local x / m（剖面里程=x+20）');fig.supylabel('模型y / m')
    fig.suptitle('190m高损耗H0：空气中弱场追踪 / 同一固定灰度 / 原生Ricker时标\n'
        '绿色实际地表，青色覆盖层底，Tx橙／Rx红，红色阴影为顶部PML。±0.1 V/m之外饱和；仅显示诊断，不是能量占比。')
    fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.65,label='Ez / V/m（固定对称对数；线性阈值0.0001 V/m）')
    image=a.public/'weak_air_return_sequence_r2.png';assert not image.exists()
    fig.savefig(image,dpi=125);plt.close(fig)
    result=dict(status='PASS_LOCAL20_NATIVE_SNAPSHOTS_SIX_FIELDS_HASH_AND_TIMING',
        script_sha256=sha(__file__),contract_sha256=sha(a.public/'execution_contract.json'),
        selected_transport=read(a.transport),native_files_locally_read=20,all500_verified_on_remote=True,
        fields=['Ex','Ey','Ez','Hx','Hy','Hz'],checked=checked,remote_idle=idle,
        weak_display=dict(file=image.name,sha256=sha(image),fixed_limit_V_m=scale,linear_threshold_V_m=threshold,
            displayed_iterations=steps,post_hoc_display_only=True),
        limits='Selected20 files independently read locally; all500 audit ran on target. No Poynting flux, unique path or physical closure.')
    (a.public/'local_snapshot_audit_r2.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (a.public/'remote_idle.json').write_text(json.dumps(idle,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],native_files_locally_read=20,owned_remote_pids=idle['owned_pids'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['samples','transport','geometry','public']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
