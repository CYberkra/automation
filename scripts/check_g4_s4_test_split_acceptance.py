"""Read-only acceptance for S4 test-family ladder and capability exports."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
R1 = ROOT / 'artifacts/research_checks/2026-09-27_g4_mission_capability_test_r1'
R2 = ROOT / 'artifacts/research_checks/2026-09-27_g4_mission_capability_test_r2'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def main():
    results = [(d / 'results.json').read_bytes() for d in (R1, R2)]
    assert results[0] == results[1], 'test capability exports are not byte-identical'
    print('test capability r1/r2 byte-identical:', sha256(R1 / 'results.json'))

    data = json.loads(results[0].decode('utf-8'))
    assert data['split'] == 'test' and data['families'] == ['c5', 'c8']
    assert data['ladder_records_merged'] == 12768
    assert data['hard_limits']['G4'].startswith('NOT relieved')
    assert data['hard_limits']['training_eligible'] is False
    assert data['hard_limits']['physical_acceptance_threshold'] is None
    assert {r['geometry'] for r in data['capability']} == {'mt'}
    assert {r['family'] for r in data['capability']} == {'c5', 'c8'}
    candidates = {r['candidate_id'] for r in data['capability']} | {
        r['candidate_id'] for r in data['negative_control']}
    assert len(candidates) == 57, f'expected 57 candidates, got {len(candidates)}'
    assert all(r['geometry'] == 'mt' and r['family'] in ('c5', 'c8')
               for r in data['negative_control'])
    print('capability:', data['n_capability_rows'], 'rows;',
          data['n_negative_control_rows'], 'negative-control rows; 57 candidates; MT only')

    for run, directory in (('r1', R1), ('r2', R2)):
        manifest = read_json(directory / 'run_manifest.json')
        assert manifest['split'] == 'test' and manifest['ladder_run'] == run
        assert manifest['exporter_sha256'] == sha256(ROOT / 'scripts/run_g4_mission_capability.py')
        assert manifest['results_sha256'] == sha256(directory / 'results.json')
        assert manifest['mission_tolerance_sha256'] == sha256(
            ROOT / 'configs/research/g4_mission_tolerance_v0.2.json')
        assert manifest['event_table_sha256'] == sha256(
            ROOT / 'configs/research/batch2d_v1_event_table_v0.1.json')
        assert len(manifest['ladder_inputs']) == 8
        for k, item in enumerate(manifest['ladder_inputs']):
            chunk = ROOT / f'artifacts/research_checks/2026-09-27_damage_ladder_test_{run}_c{k:02d}'
            assert item['chunk'] == k
            assert item['records_path'].replace('\\', '/') == (
                f'artifacts/research_checks/2026-09-27_damage_ladder_test_{run}_c{k:02d}/records.json')
            assert item['manifest_path'].replace('\\', '/') == (
                f'artifacts/research_checks/2026-09-27_damage_ladder_test_{run}_c{k:02d}/run_manifest.json')
            assert item['records_sha256'] == sha256(chunk / 'records.json')
            assert item['manifest_sha256'] == sha256(chunk / 'run_manifest.json')
            ladder_manifest = read_json(chunk / 'run_manifest.json')
            assert ladder_manifest['split'] == 'test'
            assert ladder_manifest['families_used'] == ['c5', 'c8']
            assert ladder_manifest['gates']['test_groups_only'] is True
            assert ladder_manifest['gates']['dev_groups_only'] is False
            assert ladder_manifest['hard_limits']['fdtd_solver_invoked'] is False
        print(run, 'source manifest: 8 actual chunks and hashes verified')

    print('ALL G4 S4 TEST-SPLIT ACCEPTANCE CHECKS PASSED')


if __name__ == '__main__':
    main()
