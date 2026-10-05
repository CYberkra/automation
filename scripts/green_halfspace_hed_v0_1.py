"""Independent half-space reference: x-directed Hertzian dipole over a lossy dispersive half-space.

Spectral (Sommerfeld) formulation, e^{+i omega t} convention consistent with
the validated free-space dipole_field (spatial factor exp(-i k r), lossy
epsilon = eps' - i eps''). Reflected co-polar Ex at receiver azimuth phi=0:

  Ex^R = K ∫ (k_rho/k_z0) exp(-i k_z0 (z+z'))
         [R_TE (J0+J2) - R_TM (k_z0^2/k0^2) (J0-J2)] dk_rho

with beta = k_rho*rho inside J0/J2. The constant K = -omega*mu0/(8*pi) is the
analytic Weyl normalization (validated against the FDTD-verified dipole_field);
branch choices and R signs are pinned by the PEC image check.

Self-checks (action `selftest`):
1. eps_r=1, sigma=0, no Debye  -> reflected field ~ 0.
2. PEC limit (sigma huge)      -> reflected = -image via dipole_field.
3. Direct spectral integral    -> matches dipole_field across band and heights.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import jv


def j0(x):
    return jv(0, x)


def j2(x):
    return jv(2, x)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from analyze_uav_local_free_space import dipole_field, FREQ

C_SI = 299792458.
MU0 = 4*np.pi*1e-7
EPS0 = 1/(MU0*C_SI*C_SI)

N_PROP = 16384  # t-grid for k_rho in [0, k0], k_rho = k0 sin t
N_EVAN = 8192   # s-grid for k_rho = sqrt(k0^2+s^2), decaying tail
S_MAX_NP = 32.0  # integrate evanescent decay to exp(-S_MAX_NP)


def epsr_cover(freqs, eps_inf=18.017, sigma=0.003, d_eps=7.878, tau=6.4567e-9):
    """Archived cover material complex relative permittivity (eps' - i eps'')."""
    w = 2*np.pi*freqs
    wt = w*tau
    er1 = eps_inf + d_eps/(1+wt**2)
    er2 = d_eps*wt/(1+wt**2) + sigma/(w*EPS0)
    return er1 - 1j*er2


def _spectral_integrals(k0, k1sq, rho, zsum, use_r=True):
    """Return the two k_rho integrals I_TE, I_TM such that
    Ex^R = K [I_TE + I_TM] with R_TE/R_TM inside; plus direct dyadic integrals
    when use_r=False (then R=1 forced and dyadic = free-space xx)."""
    k1sq = np.complex128(k1sq)
    # propagating region: k_rho = k0 sin t, measure (k_rho/k_z0) dk_rho = k0 sin t dt
    t = (np.arange(N_PROP)+0.5)/N_PROP*(np.pi/2)
    kz0 = k0*np.cos(t)
    kr = k0*np.sin(t)
    # evanescent region: k_rho = sqrt(k0^2 + s^2), k_z0 = -i s, measure = i ds
    smax = S_MAX_NP/zsum
    s = (np.arange(N_EVAN)+0.5)/N_EVAN*smax
    kz0 = np.concatenate([kz0, -1j*s])
    kr = np.concatenate([kr, np.sqrt(k0**2+s**2)])
    w_prop = k0*np.sin(t)*(np.pi/2)/N_PROP
    w_evan = np.full(N_EVAN, smax/N_EVAN, dtype=complex)  # multiplied by i later
    beta = np.outer(kr, [rho])[:, 0] if np.isscalar(rho) else None
    b = kr*rho
    J0m, J2m = j0(b), j2(b)
    W = np.exp(-1j*kz0*zsum)
    kz1 = np.sqrt(k1sq - kr**2)
    kz1 = np.where(kz1.imag > 0, -kz1, kz1)  # decay into ground: Im(k_z1) <= 0
    k1 = np.sqrt(k1sq)
    eps1r = k1sq/k0**2
    R_TE = (kz0-kz1)/(kz0+kz1)
    R_TM = (eps1r*kz0-kz1)/(eps1r*kz0+kz1)
    if not use_r:
        R_TE = np.ones_like(kz0); R_TM = np.ones_like(kz0)
    integ_te = R_TE*(J0m+J2m)*W
    integ_tm = -R_TM*(kz0**2/k0**2)*(J0m-J2m)*W
    val_te = np.sum(integ_te[:N_PROP]*w_prop) + 1j*np.sum(integ_te[N_PROP:])*w_evan[0]
    val_tm = np.sum(integ_tm[:N_PROP]*w_prop) + 1j*np.sum(integ_tm[N_PROP:])*w_evan[0]
    # direct free-space xx dyadic: 2J0 - (k_rho^2/k0^2)(J0-J2)
    integ_dir = (2*J0m - (kr**2/k0**2)*(J0m-J2m))*W
    val_dir = np.sum(integ_dir[:N_PROP]*w_prop) + 1j*np.sum(integ_dir[N_PROP:])*w_evan[0]
    return val_te + val_tm, val_dir


def k_norm_of(freq):
    """Per-frequency Weyl constant: Ex = -omega*mu0/(8*pi) * (spectral integral)."""
    return -MU0*2*np.pi*freq/(8*np.pi)


def reflected_ex(freqs, rho, zsum, eps1r):
    """Reflected co-pol Ex of x-HED, receiver at (rho, phi=0), both at height zsum/2."""
    freqs = np.asarray(freqs, dtype=float)
    out = np.empty_like(freqs, dtype=complex)
    for i, f in enumerate(freqs):
        k0 = 2*np.pi*f/C_SI
        te_tm, _ = _spectral_integrals(k0, k0**2*eps1r[i], rho, zsum)
        out[i] = k_norm_of(f)*te_tm
    return out


def direct_ex_spectral(freqs, rho, zdiff):
    """Direct (free-space) co-pol Ex via the same spectral machinery."""
    freqs = np.asarray(freqs, dtype=float)
    out = np.empty_like(freqs, dtype=complex)
    for i, f in enumerate(freqs):
        k0 = 2*np.pi*f/C_SI
        _, vd = _spectral_integrals(k0, k0**2*1.0, rho, zdiff)
        out[i] = k_norm_of(f)*vd
    return out


def image_ex(freqs, rho, zsum):
    """PEC image: x-dipole of moment -I*dl at mirror position, x-component."""
    exact = dipole_field(np.array([rho, 0., zsum]), 'x')
    return -exact[:, 0]


def selftest(out):
    rho, zsum = 1.3, 16.0
    # check 3 (direct machinery): spectral direct matches dipole_field over band and heights
    direct_err = {}
    for f in (20e6, 95e6, 170e6):
        for z in (16.0, 8.0, 3.0):
            vd = direct_ex_spectral(np.array([f]), rho, z)[0]
            idx = int(np.argmin(np.abs(FREQ-f)))
            ez = dipole_field(np.array([rho, 0., z]), 'x')[idx, 0]
            direct_err[f'f{f/1e6:g}MHz_z{z:g}m'] = float(abs(vd-ez)/abs(ez))
    # check 1: eps_r=1 -> reflected ~ 0 over band
    ones = np.ones(len(FREQ), dtype=complex)
    r1 = reflected_ex(FREQ, rho, zsum, ones)
    exact_dir = dipole_field(np.array([rho, 0., 0.]), 'x')[:, 0]
    rel1 = float(np.max(np.abs(r1))/np.max(np.abs(exact_dir)))
    # check 2: PEC (sigma=1e12) -> -image
    er_pec = 1.0 - 1j*1e12/(2*np.pi*FREQ*EPS0)
    r2 = reflected_ex(FREQ, rho, zsum, er_pec)
    img = image_ex(FREQ, rho, zsum)
    rel2 = float(np.linalg.norm(r2-img)/np.linalg.norm(img))
    result = {'status': 'PASS' if (max(direct_err.values()) < 1e-6 and rel1 < 1e-6 and rel2 < 1e-3) else 'FAIL',
              'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'k_norm_form': '-omega*mu0/(8*pi), per frequency (analytic Weyl normalization)',
              'direct_spectral_vs_dipole_field_rel': direct_err,
              'eps1_reflected_over_direct_max': rel1,
              'pec_vs_negative_image_relL2': rel2,
              'tolerances': {'direct': 1e-6, 'eps1': 1e-6, 'pec': 1e-3},
              'convention': 'e^{+i wt}; spatial exp(-i k r); lossy eps = eps1 - i*eps2',
              'scope': 'Independent quasi-analytic reference for x-HED over half-space, co-pol Ex at phi=0.'}
    out.mkdir(parents=True)
    (out/'selftest.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action', choices=('selftest',))
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    selftest(a.out)


if __name__ == '__main__':
    main()
