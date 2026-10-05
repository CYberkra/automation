"""Small counterexamples for the read-only observation audit; no FDTD run."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from gprMax.toolboxes.SFCW import processing as sf
from audit_sim2real_observation_v0_1 import F, main, raw_tail, relative, response
from hs_capsule_identity import sha256


def check(out):
    if out.exists():
        raise ValueError('refuse to overwrite checks')
    results = []

    def require(name, condition):
        if not condition:
            raise AssertionError(name)
        results.append(dict(name=name, status='PASS'))

    # Native Yee clock origins and polarized current-element length both matter.
    samples = np.zeros(1000)
    samples[20] = 1
    source = sf.SampledSignal('toy-source', samples, 1e-9, .5e-9,
                              spatial_scale=.05)
    receiver = sf.SampledSignal('toy-receiver', 2*samples, 1e-9, 1e-9)
    product = response(source, receiver, 0)
    expected = 40*np.exp(-2j*np.pi*F*.5e-9)
    require('actual_clock_half_step_phase', relative(product.response, expected)<1e-12)
    require('native_y_length_not_x_spacing', np.allclose(abs(product.response), 40))

    # A strong early background cannot certify completeness of a weak late event.
    background = np.zeros(651)
    background[20] = 1
    late = np.zeros(651)
    late[590:] = 1e-9
    t = np.arange(651, dtype=float)*1e-9
    b = dict(receiver=SimpleNamespace(samples=background, times=t))
    r = dict(receiver=SimpleNamespace(samples=background+late, times=t))
    tail = raw_tail(r, b, 18.5)
    require('target_tail_uses_own_peak', tail['last_sample_over_isolated_peak']==1)
    require('tiny_total_tail_does_not_hide_open_target_tail',
            np.sqrt(np.mean((background+late)[585:]**2))<1e-8
            and tail['last_tenth_rms_over_isolated_peak']>.9)

    # Omitting a pure propagation delay preserves amplitude but corrupts complex
    # comparison. This is an exact toy, not validation of the spreading model.
    reference = np.exp(-F/1e8).astype(complex)
    delta = 7
    k = 2*np.pi*F/299792458 - .1j
    measured = reference*np.exp(-2j*k*delta)
    attenuation_only = reference*np.exp(2*k.imag*delta)
    require('attenuation_only_can_match_magnitude_but_fail_phase',
            np.allclose(abs(measured), abs(attenuation_only))
            and relative(attenuation_only, measured)>.5)
    require('complex_propagation_includes_delay',
            relative(reference*np.exp(-2j*k*delta), measured)==0)

    receiver_late = sf.SampledSignal('late-toy', late, 1e-9, 0)
    impulse = np.zeros(651)
    impulse[20] = 1
    source_late = sf.SampledSignal('source-toy', impulse, 1e-9, 0, spatial_scale=.05)
    require('receiver_taper_changes_real_late_response',
            relative(response(source_late, receiver_late, 20).response,
                     response(source_late, receiver_late, 0).response)>.01)

    # Existing output must be rejected before loading or changing any evidence.
    out.mkdir(parents=True)
    try:
        main(out)
    except ValueError as error:
        require('existing_audit_output_refused', 'new output' in str(error))
    else:
        raise AssertionError('existing output was not refused')
    result = dict(status='PASS_ARRAY_COUNTEREXAMPLES_NOT_FDTD_VALIDATION',
                  checks=results, calls_solver=False, calls_training=False,
                  code_sha256=sha256(__file__))
    (out/'checks.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    check(parser.parse_args().out)
