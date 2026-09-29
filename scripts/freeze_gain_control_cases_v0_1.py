"""Freeze the gain-class control-case contract v0.1 (design user-approved).

User approved the gain-class controlled-case design ("批准", 2026-09-29;
docs/research/2026-09-29_gain_class_control_case_design.md). This contract
pre-registers everything the gain ladder needs so the study cannot drift:
attenuation levels, window structures, NC dual state, the original-event
reference declaration, and the clipping model boundary. The ONLY quantity
left for later calibration is the clip tolerance itself (design §4 step 2-4:
derived from effects-table evidence, user-confirmed, never reused from other
classes).

Freeze-script assertion discipline: design-doc presence and key pre-registered
values asserted before writing; build runs twice in-process, byte-identical.
"""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'configs/research/gain_control_cases_v0.1.json'
DESIGN = ROOT / 'docs/research/2026-09-29_gain_class_control_case_design.md'

LEVELS_DB = [-3.0, -6.0, -10.0, -20.0, -40.0]
STRUCTURES = ['event_only', 'event_band', 'time_global']
NC_STATES = ['nc_untouched', 'nc_attenuated']


def sha256(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build():
    text = DESIGN.read_text(encoding='utf-8')
    for token in ('−3 / −6 / −10 / −20 / −40', 'event_band', 'nc_attenuated',
                  '原始未衰减事件', 'simulated_full_scale'):
        assert token in text, f'design doc missing pre-registered token: {token}'

    doc = {
        'schema': 'gain-control-cases/0.1',
        'date': '2026-09-29',
        'status': 'frozen',
        'frozen_by': 'user "批准" (2026-09-29) of '
                     'docs/research/2026-09-29_gain_class_control_case_design.md',
        'carrier': 't3 development families BG CO33 signed radargrams (S2X/S2TZX/'
                   'C3mX; official SFCW chain, Hann window, 95 MHz; identical '
                   'loading chain to study_t3_damage_ladder)',
        'attenuation_levels_db': LEVELS_DB,
        'linear_factors': {str(db): round(10.0 ** (db / 20.0), 9) for db in LEVELS_DB},
        'structures': {
            'event_only': 'multiplicative factor k inside the Fermat event window '
                          '(er=18, +/-10 ns) only',
            'event_band': 'factor k inside the event window; raised-cosine taper '
                          'to 1 over +/-20 ns around the window (local lossy '
                          'body edge model)',
            'time_global': 'factor k on all samples from the earliest event-window '
                           'edge onward (extra path-loss model)',
        },
        'nc_states': {
            'nc_untouched': 'NC window (250-400 ns) not attenuated; gain must not '
                            'inflate NC energy',
            'nc_attenuated': 'NC window carries the same factor k; genuine '
                             'compensation restores it equally (distinguishes true '
                             'compensation from event-directed amplification)',
        },
        'reference': {
            'type': 'original_unattenuated_event',
            'declaration': 'constructed_reference: the attenuation is imposed by '
                           'us, so the original event is known; quality labels '
                           'measure recovery OF the known original, the opposite '
                           'direction of the background-class ladder (whose '
                           'reference is the damaged input itself)',
        },
        'clipping_model': {
            'normalization': 'radargram scaled so global max|x| = 1 before gain',
            'full_scale': 1.0,
            'status': 'simulated_full_scale_array_construction_not_device_adc',
            'clip_tolerance': None,
            'clip_tolerance_rule': 'to be derived from effects-table evidence '
                                   '(design step 3-4) and user-confirmed; no '
                                   'existing number may be reused',
        },
        'identity_anchor_expectations': {
            'event_only': 'measured a on the event window equals the linear '
                          'factor exactly (float tolerance 1e-12 relative)',
            'nc_untouched': 'NC energy ratio vs original equals 1.0 exactly',
            'nc_attenuated': 'NC energy ratio vs original equals k^2 exactly',
            'control_unattenuated': 'identity/no-attenuation control: a=1, D_e=0',
        },
        'metric_battery': ['a (event amplitude factor vs original)',
                           'D_e/A/H/rho (event window vs original)',
                           'nc_ratio_vs_original (mean-square ratio)',
                           'peak_abs diagnostic'],
        'scope_limits': [
            'values are properties of the t3 development families plus the '
            'dispersion material assumption; not field attenuation models',
            'background-class frozen contracts untouched; this contract adds the '
            'gain class only',
            'test families {C5,C8} untouched; no solver runs; no training',
            'clip conclusions are internal-scale comparisons; not device claims',
        ],
        'evidence': {
            'design_doc': {'path': 'docs/research/2026-09-29_gain_class_control_case_design.md',
                           'sha256': sha256(DESIGN)},
        },
    }
    return doc


def main():
    text = json.dumps(build(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    again = json.dumps(build(), ensure_ascii=False, indent=1, sort_keys=True) + '\n'
    assert text == again, 'two build runs disagree'
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('anchors verified | wrote', OUT.relative_to(ROOT))
    print('sha256', hashlib.sha256(text.encode('utf-8')).hexdigest())


if __name__ == '__main__':
    main()
