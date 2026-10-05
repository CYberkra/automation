"""Analytical review of the user's low-loss Debye candidate; no FDTD calls."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import gprMax
from gprMax import config
from gprMax.materials import DispersiveMaterial

C = 299792458.
EPS0 = config.e0
TAU = 6.4567e-9
DEPTH = 18.5
NEW = dict(epsilon_infinity=11., sigma_S_m=.001, delta_epsilon=.5, tau_s=TAU)
OLD = dict(epsilon_infinity=18.017, sigma_S_m=.003, delta_epsilon=7.878, tau_s=TAU)


def dielectric(f, p):
    w = 2*np.pi*f
    return p['epsilon_infinity'] + p['delta_epsilon']/(1+1j*w*p['tau_s']) - 1j*p['sigma_S_m']/(w*EPS0)


def group_delay(f, p):
    w = 2*np.pi*f
    e = dielectric(f, p)
    derivative = (-1j*p['tau_s']*p['delta_epsilon']/(1+1j*w*p['tau_s'])**2
                  + 1j*p['sigma_S_m']/(EPS0*w*w))
    return 2*DEPTH/C*np.real(np.sqrt(e)+w*derivative/(2*np.sqrt(e)))


def metrics(f):
    w = 2*np.pi*f
    e = dielectric(f, NEW)
    old_e = dielectric(f, OLD)
    n, old_n = np.sqrt(e), np.sqrt(old_e)
    k, old_k = w/C*n, w/C*old_n
    alpha, old_alpha = -k.imag, -old_k.imag
    rock_n = np.sqrt(9-1j*.001/(w*EPS0))
    reflection = (n-rock_n)/(n+rock_n)
    old_reflection = (old_n-rock_n)/(old_n+rock_n)
    surface, old_surface = (1-n)/(1+n), (1-old_n)/(1+old_n)
    gain_absorption = 20/np.log(10)*2*DEPTH*(old_alpha-alpha)
    gain_reflection = 20*np.log10(abs(reflection/old_reflection))
    return dict(frequency_Hz=f, epsilon_real=e.real, epsilon_loss=-e.imag,
                phase_refractive_index=n.real, wavelength_m=2*np.pi/k.real,
                amplitude_1_e_length_m=1/alpha,
                two_way_absorption_dB=-20/np.log(10)*2*DEPTH*alpha,
                old_two_way_absorption_dB=-20/np.log(10)*2*DEPTH*old_alpha,
                normal_interface_abs_reflection=abs(reflection),
                normal_interface_reflection_dB=20*np.log10(abs(reflection)),
                old_normal_interface_reflection_dB=20*np.log10(abs(old_reflection)),
                normal_surface_reflection_dB=20*np.log10(abs(surface)),
                old_normal_surface_reflection_dB=20*np.log10(abs(old_surface)),
                absorption_improvement_dB=gain_absorption,
                interface_change_dB=gain_reflection,
                absorption_plus_interface_gain_dB=gain_absorption+gain_reflection,
                bulk_two_way_group_delay_ns=group_delay(f, NEW)*1e9)


def main(out):
    if out.exists():
        raise ValueError('refuse to overwrite review')
    if gprMax.__version__ != '4.0.0':
        raise ValueError('reviewed local V4.0.0 implementation required')
    f = np.linspace(20e6, 170e6, 501)
    result = metrics(f)
    native = DispersiveMaterial(2, 'cover_review')
    native.type='debye'; native.er=11.; native.se=.001
    native.poles=1; native.deltaer=[.5]; native.tau=[TAU]
    expected = dielectric(f, NEW)
    official = np.array([native.calculate_er(v) for v in f])
    native_error = np.linalg.norm(expected-official)/np.linalg.norm(official)
    assert native_error < 1e-14
    # Separate real/imaginary Debye + conductivity formula, independent of sqrt.
    x = 2*np.pi*f*TAU
    assert np.allclose(expected.real, 11+.5/(1+x*x), rtol=1e-14)
    assert np.allclose(-expected.imag, .5*x/(1+x*x)+.001/(2*np.pi*f*EPS0), rtol=1e-14)
    # Check analytical group delay against phase-wavenumber finite differences.
    df = 100.
    derivative_fd = (np.real(2*np.pi*(f+df)/C*np.sqrt(dielectric(f+df, NEW)))
                     -np.real(2*np.pi*(f-df)/C*np.sqrt(dielectric(f-df, NEW))))/(2*df)
    error_group = np.max(abs(2*DEPTH*derivative_fd/(2*np.pi)/group_delay(f, NEW)-1))
    assert error_group < 1e-8
    assert np.all(-expected.imag > 0) and np.all(result['amplitude_1_e_length_m'] > 0)
    # Initial Yee estimate only: lossless continuum n(real), axial 2D wave,
    # maximum air CFL time step. Does NOT model gprMax's discrete Debye update.
    mesh=[]
    for dl in (.025, .05):
        dt=dl/(C*np.sqrt(2))
        n=np.sqrt(dielectric(f, NEW)).real
        k=2*np.pi*f/C*n
        k_grid=2/dl*np.arcsin(n*np.sqrt(2)*np.sin(2*np.pi*f*dt/2))
        mesh.append(dict(spacing_m=dl, air_CFL_2D_dt_s=dt,
                         cells_per_wavelength_at170=float(result['wavelength_m'][-1]/dl),
                         estimated_unwrapped_37m_phase_error_at170_deg=float((k_grid[-1]-k[-1])*2*DEPTH*180/np.pi),
                         scope='Lossless axial Yee illustration; ignores discrete Debye/loss, interfaces, angles; not measured FDTD error'))
    exact = metrics(np.array([20,40,95,170])*1e6)
    rows=[{key:float(value[i]) for key,value in exact.items()} for i in range(4)]
    out.mkdir(parents=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei']
    plt.rcParams['axes.unicode_minus']=False
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for key,label in [('old_two_way_absorption_dB','归档覆盖层'),('two_way_absorption_dB','建议覆盖层')]:
        axes[0].plot(f/1e6,result[key],label=label)
    axes[0].set(title='18.5m地下双程：仅介质吸收',xlabel='频率 (MHz)',ylabel='场幅损耗 (dB)')
    axes[0].legend()
    axes[1].plot(f/1e6,result['bulk_two_way_group_delay_ns'])
    axes[1].set(title='建议覆盖层：双程体传播群时延',xlabel='频率 (MHz)',ylabel='群时延 (ns)，未含空气路径')
    for ax in axes:ax.grid(alpha=.25)
    fig.suptitle('低损耗配方解析复核：非完整接收回波、非实测标定、未求解')
    fig.savefig(out/'material_propagation_review.png',dpi=150)
    plt.close(fig)
    summary=dict(status='ANALYTICALLY_REVIEWED_CANDIDATE_NOT_EXECUTION_CONTRACT',
                 calls_solver=False,calls_training=False,field_fitted=False,
                 recommended_role='Low-loss weak-dispersion mechanism control; not measured field material',
                 candidate=NEW, historical=OLD, bedrock=dict(epsilon_r=9,sigma_S_m=.001),
                 epsilon_static_dielectric_part=11.5, debye_only_loss_peak_MHz=1/(2*np.pi*TAU)/1e6,
                 table=rows, full_band_group_delay_range_ns=[float(np.min(result['bulk_two_way_group_delay_ns'])),float(np.max(result['bulk_two_way_group_delay_ns']))],
                 two_way_absorption_spectral_spread_dB=float(np.ptp(result['two_way_absorption_dB'])),
                 mesh_initial_diagnostics=mesh,
                 checks=dict(native_V4_calculate_er_relative_L2=float(native_error),group_delay_finite_difference_max_relative_error=float(error_group)),
                 gain_scope='Frequency-specific bulk absorption plus single normal-incidence interface; no transmission/spreading/antenna/direct/background/processing gains',
                 sign_convention='exp(+j omega t), epsilon=real-j*loss, field propagation=exp(-j k L)',
                 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args().out)
