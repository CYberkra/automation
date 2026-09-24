"""Prepare inputs/argument records for the OFFICIAL V4 SFCW tool; no solver import.

The archived v1 geometry is reused; no signal processing is implemented here.
"""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, required=True)
    args = parser.parse_args()
    old = ROOT/'configs/research/v4_calibration_packet_v1'
    manifest = json.loads((old/'manifest.json').read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    for run in manifest['runs']:
        original = old/run['input']
        assert hashlib.sha256(original.read_bytes()).hexdigest() == run['input_sha256']
        text = original.read_text().replace('V4 calibration', 'V4 official SFCW calibration')
        text = text.replace('#waveform: ricker 1 80e6 pulse80', '#waveform: impulse 1 1 impulse')
        text = text.replace(' inf pulse80', ' inf impulse')
        assert 'ricker' not in text and 'pulse80' not in text
        target = args.output/run['input']
        target.write_text(text, encoding='utf-8', newline='\n')
        run['input_sha256'] = hashlib.sha256(target.read_bytes()).hexdigest()
    source_files = ['gprMax/_version.py', 'gprMax/toolboxes/SFCW/README.rst',
                    'gprMax/toolboxes/SFCW/cli.py', 'gprMax/toolboxes/SFCW/processing.py',
                    'gprMax/toolboxes/SFCW/examples/cylinder_sfcw_2D.in']
    manifest.update(packet_id='V4-CAL-02-SFCW', status='review_only_official_sfcw_not_runtime_validated',
                    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    parent_packet='V4-CAL-01 geometry only; Ricker excitation superseded',
                    source_timing='built-in impulse, amplitude 1, start zero; actual source Yee-time metadata required',
                    official_source_sha256={p:hashlib.sha256((args.source_root/p).read_bytes()).hexdigest() for p in source_files},
                    sfcw=dict(implementation='gprMax.toolboxes.SFCW', method='direct',
                              frequency_start_Hz=20000000, frequency_stop_Hz=170000000,
                              frequency_step_Hz=300000, frequency_count=501,
                              receiver='name:measurement', component='Ez',
                              source_selection='resolve actual unique source from official inspect output',
                              window='rectangular', zero_pad=1, time_shift_s=0,
                              tail_taper=0, source_floor_db=-100,
                              source_floor_note='official default; numerical guard, not physical accuracy certificate',
                              independent_check='official homodyne on M00_base only; no additional FDTD',
                              primary_output='complex frequency response with source/receiver spectra and validity mask',
                              hardware_export_equivalence=False),
                    proposed_run_order=['M00_base','M01_base','M00_fine','M01_fine',
                                        'M00_boundary','M01_boundary','M00_long','M01_long'])
    manifest['proposed_limits']['total_official_postprocess_wall_minutes'] = 30
    (args.output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'packet':manifest['packet_id'],'inputs':len(manifest['runs']),
                      'official_source_files_hashed':len(source_files),'solver_executed':False}))


if __name__ == '__main__':
    main()
