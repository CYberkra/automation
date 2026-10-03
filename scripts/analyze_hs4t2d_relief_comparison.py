"""Compare locally rerun baseline, 0.8 m relief and matched planar reference."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import h5py
import numpy as np

from hs_capsule_identity import sha256
from plot_hs4_bscan_v0_2 import diagnostic_arrays
from sfcw_official_loader_v0_2 import FREQ, verify_official_runtime, OFFICIAL_PROCESSING_SHA256
from gprMax.toolboxes.SFCW.processing import load_source, load_receiver, direct_frequency_response, reconstruct_time_response


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--run-directory',type=Path,help='relocated raw archive; original contract and output hashes remain mandatory')
    args=ap.parse_args()
    if args.out.exists():
        raise ValueError('new output directory required')
    verify_official_runtime()
    contract=json.loads(args.contract.read_text('utf-8'))
    run=args.run_directory or Path(contract['run_directory'])
    events=[json.loads(s) for s in (run/'series_execution.jsonl').read_text('utf-8').splitlines()]
    if events[-1].get('status')!='COMPLETED' or events[-1].get('traces')!=363:
        raise ValueError('complete 363-trace run required')
    if events[0]['contract_sha256']!=sha256(args.contract):
        raise ValueError('contract identity changed')
    responses, products, midpoints, raw_summary = {}, {}, None, {}
    reference_source, reference_rx = None, None
    for name in ['baseline','relief2p0','flat']:
        folder=run/name
        audit=json.loads((folder/'audit.json').read_text('utf-8'))
        traces=[]
        for k,row in enumerate(audit,1):
            path=folder/row['file']
            if row['trace']!=k or sha256(path)!=row['sha256']:
                raise ValueError('raw output identity differs')
            src=load_source(path)
            rx=load_receiver(path,receiver_path='/rxs/rx1',component='Ey')
            if reference_source is None:
                reference_source,reference_rx=src,rx
            if not (src.dt==reference_source.dt and src.time_offset==reference_source.time_offset
                    and np.array_equal(src.samples,reference_source.samples)
                    and rx.dt==reference_rx.dt and rx.time_offset==reference_rx.time_offset
                    and rx.samples.shape==reference_rx.samples.shape):
                raise ValueError('source or receiver time axis differs')
            with h5py.File(path,'r') as h:
                if h['rxs/rx1/Ey'].dtype!=np.float64 or str(h.attrs['gprMax'])!='4.0.0':
                    raise ValueError('version/dtype audit failed')
            traces.append(rx.samples)
        samples=np.stack(traces,axis=1)
        taper=(round(200e-9/reference_rx.dt)-.25)/len(samples)
        r=direct_frequency_response(reference_source,replace(reference_rx,samples=samples),FREQ,tail_taper_fraction=taper)
        responses[name]=r
        products[name]=reconstruct_time_response(r,zero_pad_factor=8,window='hann')
        current=np.asarray([(row['source_m'][0]+row['receiver_m'][0])/2 for row in audit])
        if midpoints is not None and not np.array_equal(midpoints,current):
            raise ValueError('stations differ across comparison')
        midpoints=current
        raw_summary[name]={'traces':len(traces),'dtype':'float64','solver':'4.0.0',
                           'audit_sha256':sha256(folder/'audit.json'),
                           'raw_peak':float(np.max(np.abs(samples)))}
    time_ns=products['baseline'].time*1e9
    for product in products.values():
        if not np.array_equal(time_ns,product.time*1e9):
            raise ValueError('reconstructed axes differ')
    diffs={}
    for name in ['baseline','relief2p0']:
        difference=replace(responses[name],response=responses[name].response-responses['flat'].response,
                           receiver_spectrum=responses[name].receiver_spectrum-responses['flat'].receiver_spectrum)
        diffs[name]=reconstruct_time_response(difference,zero_pad_factor=8,window='hann')
    args.out.mkdir(parents=True)
    tm=time_ns[time_ns<=250]
    residual={name:diagnostic_arrays(time_ns,products[name].real_bandpass)[3] for name in products}
    saved_window=time_ns<=250
    payload={'time_ns':tm,'midpoint_m':midpoints,'residual_time_ns':tm,'frequency_Hz':FREQ,
             'full_reconstruction_time_ns':time_ns}
    for name in products:
        payload[name+'_raw_signed']=products[name].real_bandpass[saved_window]
        payload[name+'_complex_envelope']=products[name].complex_envelope[saved_window]
        payload[name+'_frequency_complex']=responses[name].response
        payload[name+'_rank1_residual']=residual[name]
    for name in diffs:
        payload[name+'_flat_difference_signed']=diffs[name].real_bandpass[saved_window]
        payload[name+'_flat_difference_complex_envelope']=diffs[name].complex_envelope[saved_window]
    np.savez_compressed(args.out/'comparison_arrays.npz',**payload)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family']='Microsoft YaHei'
    modes=[('原始带符号响应', {name:products[name].real_bandpass[time_ns<=250] for name in ['baseline','relief2p0']}),
           ('去道内均值 + 去1阶残差', {name:residual[name] for name in ['baseline','relief2p0']}),
           ('与同条件平界面作复频谱差分', {name:diffs[name].real_bandpass[time_ns<=250] for name in ['baseline','relief2p0']})]
    limits={title:float(np.quantile(np.abs(np.stack(list(values.values()))),.999)) or 1 for title,values in modes}
    for suffix,window in [('full',(0,250)),('underground',(160,220))]:
        fig,axes=plt.subplots(3,2,figsize=(13,11),constrained_layout=True)
        mask=(tm>=window[0])&(tm<=window[1])
        for i,(title,values) in enumerate(modes):
            for j,name in enumerate(['baseline','relief2p0']):
                ax=axes[i,j]
                im=ax.pcolormesh(midpoints,tm[mask],values[name][mask],shading='nearest',cmap='gray',
                                 vmin=-limits[title],vmax=limits[title],rasterized=True)
                ax.set_ylim(window[1],window[0])
                ax.set_xlabel('实际收发中点 X（m）')
                ax.set_ylabel('重建时间（ns）')
                ax.set_title(('原版：0.4 m 起伏' if name=='baseline' else '新版：0.8 m 起伏')+'\n'+title)
                fig.colorbar(im,ax=ax,label='场/源响应（同一行共用色标）')
        fig.suptitle('gprMax V4.0.0 / CUDA FP64 / 121道；官方20–170 MHz / 501频点 / Hann\n'
                     '无AGC；差分表示偏离平界面的响应，不是地下干净真值；SVD可能删除共同地层')
        fig.savefig(args.out/f'bscan_comparison_{suffix}.png',dpi=145)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,5),constrained_layout=True)
    env_limit=float(np.quantile(np.abs(np.stack([2*diffs[n].complex_envelope for n in diffs])),.999)) or 1
    mask=(time_ns>=160)&(time_ns<=220)
    for ax,name in zip(axes,['baseline','relief2p0']):
        im=ax.pcolormesh(midpoints,time_ns[mask],2*np.abs(diffs[name].complex_envelope[mask]),
                         shading='nearest',cmap='viridis',vmin=0,vmax=env_limit,rasterized=True)
        ax.set_ylim(220,160)
        ax.set_xlabel('实际收发中点 X（m）')
        ax.set_ylabel('重建时间（ns）')
        ax.set_title(('原版0.4 m' if name=='baseline' else '新版0.8 m')+'：平界面差分包络')
        fig.colorbar(im,ax=ax,label='2|差分复包络|（共享色标）')
    fig.suptitle('先作复频谱差分，再重建取模；不将亮度/纹理当作形态保真验收')
    fig.savefig(args.out/'flat_difference_envelope.png',dpi=145)
    plt.close(fig)
    metrics={}
    ground=(time_ns>=85)&(time_ns<=120)
    underground=(time_ns>=160)&(time_ns<=220)
    for name in diffs:
        matrix=diffs[name].real_bandpass
        metrics[name]={'flat_difference_ground_energy':float(np.sum(matrix[ground]**2)),
                       'flat_difference_underground_energy':float(np.sum(matrix[underground]**2)),
                       'flat_difference_underground_peak_abs':float(np.max(np.abs(matrix[underground]))),
                       'underground_envelope_max_time_ns_by_trace':time_ns[underground][np.argmax(np.abs(diffs[name].complex_envelope[underground]),axis=0)].tolist()}
    report={'status':'COMPLETED','raw_audit':raw_summary,'contract_sha256':sha256(args.contract),
            'solver_version':'4.0.0','backend':'CUDA','precision':'float64','traces':363,
            'official_processing_sha256':OFFICIAL_PROCESSING_SHA256,'plot_color_limits':limits,
            'comparison_metrics':metrics,'processing':'20-170 MHz/501 tones/0.3MHz/Hann/8x pad; same 200ns tail taper; no AGC',
            'interpretation':'Development geometry contrast, not physical acceptance. Numerical/PML attribution unresolved.',
            'analysis_code_sha256':sha256(__file__)}
    report['saved_time_domain_ns']=[0,250]
    report['complete_frequency_complex_products_saved']=True
    (args.out/'results.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    for name in products:
        for row in json.loads((run/name/'audit.json').read_text('utf-8')):
            if sha256(run/name/row['file'])!=row['sha256']:
                raise ValueError('raw output changed during analysis')
    (args.out/'manifest.json').write_text(json.dumps([{'file':p.name,'bytes':p.stat().st_size,'sha256':sha256(p)}
        for p in sorted(args.out.iterdir())],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(metrics,ensure_ascii=False)[:1500])


if __name__=='__main__':
    main()
