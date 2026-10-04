"""Input-derived materials/timing; no event identification from correlation."""
from dataclasses import dataclass
import json
from pathlib import Path
import h5py
import numpy as np
from hs_capsule_identity import sha256

C0 = 299792458.0
EPS0 = 8.8541878128e-12
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CAPSULE = ROOT / 'artifacts/research_checks/2026-10-02_hs4t2d_transect'
FREQUENCIES = np.linspace(20e6, 170e6, 501)


@dataclass(frozen=True)
class Material:
    epsilon_inf: float
    conductivity: float
    mu_r: float
    poles: tuple = ()

    @property
    def static_dielectric_permittivity(self):
        return self.epsilon_inf + sum(delta for delta, _ in self.poles)

    def epsilon(self, frequency):
        f = np.asarray(frequency, dtype=float)
        if not np.isfinite(f).all() or np.any(f <= 0):
            raise ValueError('positive finite frequencies required; DC conduction is separate')
        omega = 2 * np.pi * f
        value = self.epsilon_inf + self.conductivity / (1j * omega * EPS0)
        for delta, tau in self.poles:
            value = value + delta / (1 + 1j * omega * tau)
        return value

    def index(self, frequency):
        return np.sqrt(self.mu_r * self.epsilon(frequency))

    def group_index(self, frequency):
        f = np.asarray(frequency, dtype=float)
        step = np.minimum(1e3, f / 1000)
        return self.index(f).real + f * (self.index(f + step).real - self.index(f - step).real) / (2 * step)


def read_model(path):
    lines = path.read_text('utf-8').splitlines()
    def values(command):
        return next(s.split()[1:] for s in lines if s.startswith(command))
    definitions = {}
    for s in lines:
        if s.startswith('#material:'):
            tokens = s.split()
            definitions[tokens[-1]] = Material(*map(float, tokens[1:4]))
    for s in lines:
        if s.startswith('#add_dispersion_debye:'):
            tokens = s.split()
            count = int(tokens[1])
            material = definitions[tokens[-1]]
            poles = tuple((float(tokens[2 + 2*k]), float(tokens[3 + 2*k])) for k in range(count))
            definitions[tokens[-1]] = Material(material.epsilon_inf, material.conductivity, material.mu_r, poles)
    boxes = [(tuple(map(float, s.split()[1:7])), s.split()[-1]) for s in lines if s.startswith('#box:')]
    cover = sorted([box for box, name in boxes if name == 'cover'], key=lambda b: b[0])
    if not cover or any(b[5] != cover[0][5] for b in cover):
        raise ValueError('flat ground and cover interface required')
    return {'domain_m': list(map(float, values('#domain:'))), 'spacing_m': list(map(float, values('#dx_dy_dz:'))),
            'pml_cells': list(map(int, values('#pml_cells:'))),
            'source_m': list(map(float, values('#hertzian_dipole:')[1:4])),
            'receiver_m': list(map(float, values('#rx:')[:3])),
            'materials': definitions, 'cover_boxes': cover, 'ground_z_m': cover[0][5]}


def profile_at_x(model, x):
    x = np.asarray(x, dtype=float)
    boxes = model['cover_boxes']
    edges = np.array([b[0] for b in boxes] + [boxes[-1][3]])
    if np.any(x < edges[0]) or np.any(x >= edges[-1]):
        raise ValueError('station outside declared profile')
    return np.array([b[2] for b in boxes])[np.searchsorted(edges, x, side='right') - 1]


def vertical_roundtrip_ns(height, cover_depth, rock_depth, cover_index, rock_index):
    return 2 * (height + np.asarray(cover_depth) * cover_index + np.asarray(rock_depth) * rock_index) / C0 * 1e9


def analyse_capsule(capsule):
    path = capsule / 'hs4t2d_t01.in'
    model = read_model(path)
    identity = {r['file']: r['sha256'] for r in json.loads((capsule / 'manifest.json').read_text('utf-8'))}
    npz = capsule / 'hs4t2d_bscan_official.npz'
    if any(sha256(p) != identity[p.name] for p in (path, npz)):
        raise ValueError('capsule identity differs')
    with np.load(npz) as data:
        time, signed = data['t'].copy(), data['S'].copy()
    if signed.ndim != 2 or not np.isfinite(signed).all() or len(time) != len(signed):
        raise ValueError('invalid saved signed B-scan')
    midpoints, raw_identities = [], []
    dt = None
    for k in range(1, signed.shape[1] + 1):
        raw = capsule / f'hs4t2d_t{k:02d}.h5'
        digest = sha256(raw)
        if digest != identity[raw.name]:
            raise ValueError('raw identity differs')
        with h5py.File(raw) as h:
            values = h['rxs/rx1/Ey'][:]
            if str(h.attrs['gprMax']) != '4.0.0' or values.dtype != np.float64 or not np.isfinite(values).all():
                raise ValueError('raw validity differs')
            tx, rx = (np.asarray(h[name].attrs['Position']) for name in ('srcs/src1', 'rxs/rx1'))
            if not np.isclose(rx[0] - tx[0], 1.3) or tx[2] != 27 or rx[2] != 27:
                raise ValueError('unexpected acquisition geometry')
            midpoints.append((tx[0] + rx[0]) / 2)
            current_dt = float(h.attrs['dt'])
            if dt is not None and current_dt != dt:
                raise ValueError('raw timestep changes')
            dt = current_dt
        raw_identities.append({'file': raw.name, 'sha256': digest})
    interface_z = profile_at_x(model, midpoints)
    cover, rock = (model['materials'][name] for name in ('cover', 'rock'))
    height = model['source_m'][2] - model['ground_z_m']
    depth = model['ground_z_m'] - interface_z
    phase = vertical_roundtrip_ns(height, depth[:, None], 0, cover.index(FREQUENCIES).real, 0)
    group = vertical_roundtrip_ns(height, depth[:, None], 0, cover.group_index(FREQUENCIES), 0)
    bottom = vertical_roundtrip_ns(height, depth[:, None], interface_z[:, None], cover.group_index(FREQUENCIES), rock.group_index(FREQUENCIES))
    spacing = np.asarray(model['spacing_m'])
    active = np.rint(np.asarray(model['domain_m']) / spacing) > 1
    cfl = 1 / (C0 * np.sqrt(np.sum(1 / spacing[active]**2)))
    window = (time >= 165) & (time <= 200)
    peaks = time[window][np.argmax(np.abs(signed[window]), axis=0)]
    vertical95 = vertical_roundtrip_ns(height, depth, 0, cover.group_index(95e6), 0)
    a, b = peaks - peaks.mean(), vertical95 - vertical95.mean()
    denominator = np.linalg.norm(a) * np.linalg.norm(b)
    correlation = float(np.dot(a, b) / denominator) if denominator else None
    cppw = C0 / FREQUENCIES / cover.index(FREQUENCIES).real / max(spacing[active])
    return {'status': 'DIAGNOSTIC_NOT_EVENT_IDENTIFICATION',
            'model': {k: v for k, v in model.items() if k not in ('materials', 'cover_boxes')},
            'identities': {'input_sha256': sha256(path), 'saved_bscan_sha256': sha256(npz), 'raw_h5': raw_identities},
            'material': {'epsilon_inf': cover.epsilon_inf, 'delta_epsilon_tau': cover.poles,
                         'static_dielectric_permittivity': cover.static_dielectric_permittivity,
                         'conductivity_S_per_m': cover.conductivity,
                         'epsilon_real_20_170MHz': cover.epsilon(np.array([20e6, 170e6])).real.tolist()},
            'numerics': {'active_axes': np.array(list('xyz'))[active].tolist(), 'raw_dt_s': dt,
                         'vacuum_CFL_dt_s': float(cfl), 'raw_to_CFL_ratio': float(dt / cfl),
                         'minimum_cover_cells_per_phase_wavelength_in_band': float(cppw.min()),
                         'scope': 'dimensional vacuum CFL/wavelength diagnostic, not full dispersive stability or convergence'},
            'vertical_ray_diagnostic': {'phase_roundtrip_interface_ns_range': [float(phase.min()), float(phase.max())],
                'group_roundtrip_interface_ns_range': [float(group.min()), float(group.max())],
                'group_roundtrip_bottom_Z0_ns_range': [float(bottom.min()), float(bottom.max())],
                'scope': 'vertical midpoint proxy, both legs; excludes bistatic obliquity/reflection phase/peak interference; not event identity'},
            'window_peak_diagnostic': {'window_ns': [165, 200], 'peak_time_ns': peaks.tolist(),
                'vertical95_group_proxy_correlation': correlation, 'scope': 'signed-window maxima can switch lobes; low correlation cannot identify bottom/PML'},
            'event_183_8ns_verdict': 'UNRESOLVED: cannot identify from local-depth correlation alone',
            'deprecated_inference': 'Historical182.57ns bottom estimate was one-way; old Debye formula misread epsilon_inf and pole increment.',
            'geometry_3d_comparison': 'NOT_CHECKED by this 2D diagnostic'}


def write_analysis(capsule, output):
    if output.exists():
        raise ValueError('new output required; never overwrite historical evidence')
    report = analyse_capsule(capsule)
    report['code_identities'] = {'hs4t2d_physics_diagnostic.py': sha256(__file__)}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('material', 'numerics', 'vertical_ray_diagnostic', 'event_183_8ns_verdict')}, indent=2))
