"""GSSI-400-like antenna scaled x4 to the ~100 MHz class (declared assumed design).

Wraps gprMax.toolboxes.GPRAntennaModels.GSSI.antenna_like_GSSI_400 (2 mm,
300x300x178 mm, Stadler et al. 2022 optimisation) and applies exact EM
similitude scaling s=4: all geometry coordinates x4 (8 mm class), absorber
conductivity /4, lumped resistances unchanged, excitation replaced by an
impulse voltage source so the official SFCW actual-source chain applies.
The scaled antenna is an ASSUMED design, not the validated 400 MHz model.
"""
import gprMax
from gprMax.toolboxes.GPRAntennaModels.GSSI import antenna_like_GSSI_400
from gprMax.user_objects.user_objects import GeometryUserObject
from gprMax.user_objects.cmds_multiuse import (Material as MaterialCmd,
                                               VoltageSource as VoltageSourceCmd,
                                               Rx as RxCmd,
                                               Waveform as WaveformCmd,
                                               ExcitationFile as ExcitationFileCmd)

SCALE = 4.0
WAVEFORM_ID = 'hs4_impulse'
# Original (unscaled) preset values, for audit records
ORIGINAL = {'resolution_m': 0.002, 'casesize_m': (0.3, 0.3, 0.178),
            'sourceresistance_ohm': 257.97407389585214,
            'receiverresistance_ohm': 288.92728542970417,
            'absorber_er': 2.42966922703319, 'absorber_sigma_Sm': 0.03839822151712033}


def scaled_antenna(x, y, z):
    """Return scene objects for the x4-scaled antenna.

    (x, y, z) is the antenna geometric centre in x/y and skid bottom in z,
    matching the convention of antenna_like_GSSI_400.
    """
    objs = antenna_like_GSSI_400(0.0, 0.0, 0.0, resolution=0.002)
    out = []
    tx_scaled = None
    rx_id = None
    for o in objs:
        if isinstance(o, (ExcitationFileCmd, WaveformCmd)):
            continue  # excitation replaced below
        if isinstance(o, MaterialCmd):
            if o.kwargs['id'] == 'absorber':
                o.kwargs['se'] = o.kwargs['se'] / SCALE
            out.append(o)
            continue
        if isinstance(o, VoltageSourceCmd):
            tx_scaled = tuple(SCALE * v for v in o.kwargs['p1'])
            continue  # replaced by impulse-fed source below
        if isinstance(o, RxCmd):
            p1_scaled = tuple(SCALE * v + t for v, t in zip(o.kwargs['p1'], (x, y, z)))
            o.kwargs['p1'] = p1_scaled
            o.point = p1_scaled  # Rx.build uses self.point, not kwargs
            rx_id = o.kwargs['id']
            out.append(o)
            continue
        if isinstance(o, GeometryUserObject):
            for attr in ('p1', 'p2', 'p3'):
                if attr in o.kwargs:
                    o.kwargs[attr] = tuple(SCALE * v + t for v, t in zip(o.kwargs[attr], (x, y, z)))
            out.append(o)
            continue
        raise ValueError(f'unhandled antenna object type: {type(o).__name__}')
    if tx_scaled is None or rx_id is None:
        raise ValueError('scaled antenna lacks source/receiver anchors')
    wave = gprMax.Waveform(wave_type='impulse', amp=1.0, freq=100e6, id=WAVEFORM_ID)
    vs = gprMax.VoltageSource(polarisation='y',
                              p1=tuple(v + t for v, t in zip(tx_scaled, (x, y, z))),
                              resistance=ORIGINAL['sourceresistance_ohm'],
                              waveform_id=WAVEFORM_ID)
    out.extend((wave, vs))
    return out


CASE_SIZE = tuple(SCALE * v for v in ORIGINAL['casesize_m'])  # 1.2 x 1.2 x 0.712 m
GRID_HINT_M = SCALE * ORIGINAL['resolution_m']  # 8 mm


def coarse_adapt(objs, dl):
    """Adapt the scaled antenna to a coarse main grid (V4 subgrids are CPU-only).

    Declared modifications, each recorded in the returned log:
    - boxes with a sub-cell axis: pcb sheets (8 mm) dropped; hdpe skid and the
      pec shield front wall thickened to exactly one cell;
    - feed-pin edges (8 mm long) dropped; feed is the voltage source itself.
    Plates/triangles are zero-thickness surface objects and survive unchanged.
    """
    out, log = [], []
    for o in objs:
        t = type(o).__name__
        if t == 'Edge':
            log.append({'action': 'drop', 'type': t, 'material': o.kwargs.get('material_id'),
                        'p1': o.kwargs['p1'], 'p2': o.kwargs['p2'],
                        'reason': 'sub-cell feed pin; feed is the voltage source point'})
            continue
        if t == 'Box':
            p1, p2 = list(o.kwargs['p1']), list(o.kwargs['p2'])
            cells = [(b - a) / dl for a, b in zip(p1, p2)]
            thin = [i for i, c in enumerate(cells) if c < 0.5]
            if thin:
                mat = o.kwargs.get('material_id')
                if mat == 'pcb':
                    log.append({'action': 'drop', 'type': t, 'material': mat,
                                'p1': o.kwargs['p1'], 'p2': o.kwargs['p2'],
                                'reason': 'sub-cell pcb sheet (8 mm) not representable'})
                    continue
                for i in thin:
                    mid = 0.5 * (p1[i] + p2[i])
                    p1[i], p2[i] = mid - dl / 2, mid + dl / 2
                log.append({'action': 'thicken', 'type': t, 'material': mat,
                            'p1': o.kwargs['p1'], 'p2': o.kwargs['p2'],
                            'axes': thin, 'new_p1': p1, 'new_p2': p2,
                            'reason': 'sub-cell axis thickened to one cell'})
                o.kwargs['p1'], o.kwargs['p2'] = tuple(p1), tuple(p2)
        out.append(o)
    return out, log
