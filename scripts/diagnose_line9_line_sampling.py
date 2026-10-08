"""Local-template sampling warning; no solver or actual-wavefield bandlimit claim."""
import argparse,json
from pathlib import Path
import numpy as np
from hs_capsule_identity import sha256 as sha

def main(a):
    assert not a.out.exists();m=json.loads((a.package/'manifest.json').read_text('utf-8'))
    s=[r for r in m['stations'] if r['regular_line_station']];rows=[]
    for left,right in zip(s,s[1:]):
        spacing=left['chainage_m']-right['chainage_m'];delay=right['geometry']['basal_sand']['time95_ns']-left['geometry']['basal_sand']['time95_ns']
        phase=2*np.pi*95e6*delay*1e-9
        rows.append(dict(x0_m=left['chainage_m'],x1_m=right['chainage_m'],spacing_m=spacing,local_phase95_delay_change_ns=delay,
            local_phase95_change_rad=phase,exceeds_pi_at95=bool(abs(phase)>np.pi),local_constant_slope_max_step_at95_m=None if delay==0 else float(spacing/(2*95e6*abs(delay)*1e-9))))
    worst=max(rows,key=lambda r:abs(r['local_phase95_change_rad']));a.out.mkdir(parents=True)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(2,1,figsize=(12,8),layout='constrained')
    axes[0].plot([r['chainage_m'] for r in s],[r['geometry']['basal_sand']['time95_ns'] for r in s],'.-',label='局部一次底砂路径95MHz相位传播时间')
    axes[0].set(ylabel='相位传播时间 / ns',title='几何辅助局部模板，不是实际FDTD回波');axes[0].legend()
    mid=[(r['x0_m']+r['x1_m'])/2 for r in rows];phase=[abs(r['local_phase95_change_rad']) for r in rows]
    axes[1].plot(mid,phase,'.-',label='相邻2m站位的95MHz相位变化绝对值');axes[1].axhline(np.pi,c='red',ls='--',label='π：局部恒斜率Nyquist预警线')
    axes[1].set(ylabel='相位变化 / rad',xlabel='剖面里程 / m',title='超过π意味着该局部模板的相位条纹不能由2m采样可靠跟踪');axes[1].legend(fontsize=9)
    for ax in axes:ax.set_xlim(220,140);ax.grid(alpha=.2)
    fig.suptitle('批前空间采样诊断：2m用于趋势粗筛，不能认证细斜纹消失\n模型辅助预警不等于真实波场带限；未加入侧向/多次响应，不改变任何雷达数据')
    figure=a.out/'local_template_sampling_warning.png';fig.savefig(figure,dpi=140);plt.close(fig)
    result=dict(status='LOCAL_TEMPLATE_SAMPLING_WARNING_NOT_ACTUAL_WAVEFIELD_NYQUIST_CERTIFICATE',script_sha256=sha(__file__),manifest_sha256=sha(a.package/'manifest.json'),new_solves=0,
        frequency_Hz=95e6,rows=rows,worst=worst,exceeding_intervals=sum(r['exceeds_pi_at95'] for r in rows),figure_sha256=sha(figure),
        limits='Only local phase95 propagation estimates frozen before solves, not actual observed peak slopes or all501-tone/nonlocal/multiple spatial bandwidth. Nyquist line uses local constant-slope reasoning. Dense190->180m proposed because it overlaps original high-loss strong branch, not because it resolves all80m geometry.')
    (a.out/'sampling_warning.json').write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(dict(status=result['status'],worst=worst,exceeding_intervals=result['exceeding_intervals'])))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['package','out']:p.add_argument('--'+key,type=Path,required=True)
    main(p.parse_args())
