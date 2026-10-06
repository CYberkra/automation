"""Independent raw-HDF5/DFT/inverse/window checks and readable weak-return plots."""
import argparse
import json
import math
from pathlib import Path
import h5py
import numpy as np
from build_pdf_profile_geometry import digest, save_json


def main(raws, arrays, report, out):
    if out.exists():
        raise ValueError('Fresh independent evidence required')
    expected=json.loads(report.read_text('utf-8'))
    if digest(arrays)!=expected['arrays_sha256']:
        raise ValueError('Transferred arrays differ from recorded report')
    with np.load(arrays) as z:
        data={k:z[k].copy() for k in z.files}
    frequencies=20e6+np.arange(501)*.3e6
    np.testing.assert_array_equal(data['frequency_Hz'],frequencies)
    source_ref=None;errors=[]
    for i,(raw,nx) in enumerate(zip(raws,[16000,8000,6400])):
        with h5py.File(raw) as h:
            y=h['rxs/rx1/Ez'];src=h['srcs/src1/excitation/samples']
            attrs=h['srcs/src1/excitation'].attrs
            if str(h.attrs['gprMax'])!='4.0.0' or y.dtype!=np.float64 or src.dtype!=np.float64:
                raise ValueError('Native version/precision differs')
            np.testing.assert_array_equal(h.attrs['nx_ny_nz'],[nx,3000,1])
            np.testing.assert_array_equal(h.attrs['dx_dy_dz'],[.025]*3)
            np.testing.assert_array_equal(y[:],data['native_Ez_V_m'][:,i])
            x=src[:];yv=y[:]
            if source_ref is None:source_ref=x.copy()
            np.testing.assert_array_equal(x,source_ref)
            tx=np.arange(len(x))*float(attrs['SampleInterval'])+float(attrs['TimeSampleOffset'])
            ty=np.arange(len(yv))*float(y.attrs['SampleInterval'])+float(y.attrs['TimeSampleOffset'])
            np.testing.assert_array_equal(ty,data['native_time_s'])
            pieces=[]
            for first in range(0,501,61):
                f=frequencies[first:first+61,None]
                xf=float(attrs['SampleInterval'])*(np.exp(-2j*np.pi*f*tx)@x)
                yf=float(y.attrs['SampleInterval'])*(np.exp(-2j*np.pi*f*ty)@yv)
                pieces.append(yf/xf/float(attrs['SpatialScale']))
            independent=np.concatenate(pieces)
        official=data['response_complex'][:,i]
        error=float(np.linalg.norm(independent-official)/np.linalg.norm(official))
        if error>1e-10:
            raise ValueError('Independent501tone complex phase/amplitude mismatch')
        errors.append(error)
    weights=np.hanning(501);weights/=np.mean(weights)
    take=np.arange(0,len(data['sfcw_time_s']),53)
    inverse=2*np.real(np.exp(2j*np.pi*data['sfcw_time_s'][take,None]*frequencies[None,:])@(weights[:,None]*data['response_complex'])/501)
    np.testing.assert_allclose(inverse,data['Hann_signed'][take],rtol=1e-10,atol=1e-8)
    for row in expected['windows']:
        lo,hi=np.asarray(row['window_ns'])*1e-9
        for item,i in zip(row['comparisons'],[1,2]):
            for domain,time,key in [('native','native_time_s','native_Ez_V_m'),('Hann_SFCW','sfcw_time_s','Hann_signed')]:
                mask=(data[time]>=lo)&(data[time]<hi)
                ref=data[key][mask,0];delta=data[key][mask,i]-ref
                independent=math.sqrt(math.fsum(float(v)*float(v) for v in delta)/math.fsum(float(v)*float(v) for v in ref))
                if not math.isclose(independent,item[domain]['relative_L2'],rel_tol=1e-11,abs_tol=1e-14):
                    raise ValueError('Independent fsum window metric mismatch')
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(4,2,figsize=(14,13),layout='constrained')
    names=['完整400m','切片200m','切片160m']
    for offset,(time,key,label,unit) in enumerate([('native_time_s','native_Ez_V_m','原始Ricker总场','Ez / (V/m)'),('sfcw_time_s','Hann_signed','SFCW Hann总场','源归一化场响应（非S21）')]):
        t=data[time]*1e9;values=data[key]
        for col,(lo,hi) in enumerate([(300,600),(600,1200)]):
            mask=(t>=lo)&(t<hi)
            ax=axes[2*offset,col];residual=axes[2*offset+1,col]
            for i in range(3):ax.plot(t[mask],values[mask,i],label=names[i],lw=.8)
            for i in [1,2]:residual.plot(t[mask],values[mask,i]-values[mask,0],label=names[i]+'减完整域',lw=.8)
            for a,title in [(ax,label),(residual,label+'差值')]:
                a.set(title=f'{title}；{lo}–{hi}ns',xlabel='时间 / ns',ylabel=unit)
                a.legend(fontsize=8);a.grid(alpha=.15)
    fig.suptitle('X220同站位：专看弱回波与有符号差值；每窗独立物理纵轴，无逐道归一化/增益\n2.5cm/FP64/材料/15m/1200ns保留；20–170MHz/501复数频点；一站结果不代表整线等价')
    fig.savefig(out/'weak_windows_and_residuals.png',dpi=140);plt.close(fig)
    save_json(out/'independent_checks.json',dict(status='PASS_INDEPENDENT_RAW_PHASE_INVERSE_AND_METRICS_NOT_PHYSICAL_ACCEPTANCE',calls_solver=False,
        raw_sha256=[digest(p) for p in raws],arrays_sha256=digest(arrays),report_sha256=digest(report),
        independent501_DFT_relative_L2=errors,independent_inverse_checked_times=len(take),
        all_window_L2_fsum_recomputed=True,script_sha256=digest(Path(__file__))))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raws',type=Path,nargs=3,required=True)
    p.add_argument('--arrays',type=Path,required=True);p.add_argument('--report',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();main(a.raws,a.arrays,a.report,a.out)
