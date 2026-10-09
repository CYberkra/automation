"""Render verified dual-ROI Ricker snapshots and export small native receipts."""
import argparse
import json
from pathlib import Path
import shutil
import zipfile

import h5py
import numpy as np
from hs_capsule_identity import sha256 as sha


def main(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from PIL import Image
    read = lambda p: json.loads(p.read_text('utf-8'))
    assert not a.out.exists()
    c = read(a.execution/'execution_contract.json'); v = read(a.execution/'completed_verification.json')
    events = [json.loads(x) for x in (a.execution/'execution.jsonl').read_text('utf-8').splitlines()]
    assert events[-1]['status'] == 'COMPLETED' and events[-1]['traces'] == 1
    assert events[-1]['verification_sha256'] == sha(a.execution/'completed_verification.json')
    assert v['completed'] and v['contract_sha256'] == sha(a.execution/'execution_contract.json')
    row = v['groups'][0]; assert row['observer_receiver_and_source_bitwise_equal'] and row['snapshot_count'] == 500
    m = c['study_manifest']; g = c['groups'][0]; package = Path(c['package'])
    dt = m['dt_s']; steps = m['snapshot_iterations'][::2]
    records = {(r['roi'], r['iteration']):r for r in row['snapshots']}
    raw = Path(row['native_path']); assert sha(raw) == row['native_sha256']
    with h5py.File(raw) as h:
        rx = h['rxs/rx1/Ez'][:]; source = h['srcs/src1/excitation/samples'][:]
        native_times = np.arange(len(rx))*dt+h['rxs/rx1/Ez'].attrs['TimeSampleOffset']
        source_peak = (source.argmax()*dt+h['srcs/src1/excitation'].attrs['TimeSampleOffset'])*1e9
    with h5py.File(package/m['groups'][0]['geometry']) as h: data = h['data'][:, :, 0]
    xx = (np.arange(data.shape[0])+.5)*.025
    surface = np.max(np.where(data!=0, np.arange(data.shape[1])[None, :]+1, 0), axis=1)*.025
    mudtop = np.max(np.where(data==2, np.arange(data.shape[1])[None, :]+1, 0), axis=1)*.025
    values = {}; limits = {}; hashes = []
    for name, roi in m['ROIs'].items():
        arrays = []
        for j in steps:
            r = records[(name,j)]; path = Path(r['file']); assert sha(path) == r['sha256']
            with h5py.File(path) as h:
                assert h.attrs['iteration'] == j and abs(h.attrs['time']-j*dt)<1e-20
                x = h['Ez'][:, :, 0]; assert x.dtype == np.float64 and list(x.shape) == roi['shape'][:2]
                assert np.isfinite(x).all()
            arrays.append(x.T)
            hashes.append(dict(roi=name, iteration=j, sha256=r['sha256']))
        values[name] = arrays
        limits[name] = max(float(abs(x).max()) for j,x in zip(steps,arrays) if j*dt>=100e-9)
    a.out.mkdir(parents=True)
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'DejaVu Sans']; plt.rcParams['axes.unicode_minus'] = False
    fig = plt.figure(figsize=(14,8), dpi=85, layout='constrained')
    grid = fig.add_gridspec(2,2,height_ratios=[3,1]); images = []
    names = list(m['ROIs'])
    for col,name in enumerate(names):
        roi = m['ROIs'][name]; bounds = roi['bounds']; lim = limits[name]
        ax = fig.add_subplot(grid[0,col])
        im = ax.imshow(values[name][0], origin='lower', extent=[bounds[0],bounds[3],bounds[1],bounds[4]],
            cmap='gray', aspect='auto', norm=SymLogNorm(lim*1e-4,vmin=-lim,vmax=lim,base=10))
        images.append(im)
        keep = (xx>=bounds[0]) & (xx<=bounds[3])
        ax.plot(xx[keep],surface[keep],color='lime',lw=.9,label='实际地表')
        ax.plot(xx[keep],mudtop[keep],color='cyan',lw=.9,label='实际覆盖层底／泥岩顶')
        ax.scatter([g['tx_m'][0],g['rx_m'][0]],[g['tx_m'][1],g['rx_m'][1]],c=['orange','red'],s=18,label='Tx橙／Rx红')
        for lo,hi in [(0,2),(208,210)]:ax.axvspan(lo,hi,color='red',alpha=.18)
        for lo,hi in [(0,2),(40.5,42.5)]:ax.axhspan(lo,hi,color='red',alpha=.18)
        if name=='global_view':ax.add_patch(plt.Rectangle((145,10),35,31.5,fill=False,ec='#c040dc',lw=1,label='局部观察区域'))
        ax.set(xlabel='模型local x / m（剖面里程=x+20）',ylabel='模型y / m',xlim=[bounds[0],bounds[3]],ylim=[bounds[1],bounds[4]],
            title=('全域粗看0.25m（含边界）' if col==0 else '局部细看0.1m')+' / H0：已移除底砂')
        ax.legend(fontsize=7,loc='lower left')
        fig.colorbar(im,ax=ax,shrink=.7,label='原生Ez / V/m；固定对称对数灰度')
    ax = fig.add_subplot(grid[1,:]); tn = native_times*1e9
    ax.plot(tn,rx,color='#222222',lw=.8,label='原生A-scan：与旧H0逐位一致')
    keep=(tn>=250)&(tn<=450); lim=max(abs(rx[keep]))
    ax.set(xlim=[250,450],ylim=[-lim*1.05,lim*1.05],xlabel='原生时间 / ns（含源时间参考）',ylabel='接收Ez / V/m')
    ax.legend(fontsize=8); cursor=ax.axvline(0,color='red',lw=1); title=fig.suptitle('')
    selected={min(steps,key=lambda j:abs(j*dt*1e9-t)) for t in [100,150,200,250,280,300,320,340,360,400,450]}
    frames=[]; static=[]
    for i,j in enumerate(steps):
        for im,name in zip(images,names):im.set_data(values[name][i])
        cursor.set_xdata([j*dt*1e9]*2)
        title.set_text(f'190m高损耗H0原生Ricker波场：t={j*dt*1e9:.2f}ns，源峰{source_peak:.3f}ns\n'
            '非SFCW时标；红色阴影为PML，绿色地表，青色覆盖层底。各视野固定灰度由100–500ns设定，早期可能饱和。')
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()).convert('P',palette=Image.Palette.ADAPTIVE,colors=128))
        if j in selected:
            p=a.out/f'wavefield_{j:05d}.png';fig.savefig(p,dpi=110)
            static.append(dict(file=p.name,iteration=j,time_ns=j*dt*1e9,sha256=sha(p)))
    gif=a.out/'single_wavefield.gif';frames[0].save(gif,save_all=True,append_images=frames[1:],duration=90,loop=0,optimize=True)
    plt.close(fig)
    for name in ['execution_contract.json','execution.jsonl','completed_verification.json','preflight_verification.json','input_audit.json']:
        shutil.copyfile(a.execution/name,a.out/name)
    shutil.copyfile(raw,a.out/'observer_H0.h5')
    receipt=dict(status='PASS_VERIFIED_NATIVE_DUAL_ROI_RENDER_NO_UNIQUE_PATH_CLAIM',renderer_sha256=sha(__file__),
        contract_sha256=sha(a.execution/'execution_contract.json'),verification_sha256=sha(a.execution/'completed_verification.json'),
        frames=len(frames),snapshot_count=500,source_peak_ns=source_peak,scales=limits,symlog_threshold_fraction=1e-4,
        selected_snapshot_hashes=hashes,static_frames=static,gif_sha256=sha(gif),gif_bytes=gif.stat().st_size,
        limits='NativeRicker visualization,not SFCW. H0 has no basal sandstone. Fixed scales can clip early frames; no physical energy fractions or unique path certified.')
    (a.out/'render_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    with zipfile.ZipFile(a.out.with_suffix('.zip'),'x',zipfile.ZIP_DEFLATED) as z:
        for p in a.out.iterdir():z.write(p,p.name)
    print(json.dumps(dict(zip=str(a.out.with_suffix('.zip')),sha256=sha(a.out.with_suffix('.zip')),gif_bytes=gif.stat().st_size)))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--execution',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
