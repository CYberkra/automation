"""Bounded17->paired audit->four wavefield controls->audited figures, no retries."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import zipfile


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def archive(folder, target):
    assert folder.is_dir() and not target.exists()
    with zipfile.ZipFile(target, 'x', zipfile.ZIP_DEFLATED) as z:
        for p in sorted(folder.rglob('*')):
            assert not p.is_symlink()
            if p.is_file():
                z.write(p, p.relative_to(folder).as_posix())
    print(json.dumps(dict(path=str(target), sha256=sha(target))))


def main(config_path):
    c = json.loads(config_path.read_text('utf-8'))
    work = Path(c['work']); work.mkdir(exist_ok=False)
    connection = json.loads(Path(c['connection']).read_text('utf-8'))
    ssh = ['ssh', *connection['options'], connection['host']]
    deadline = time.monotonic() + 10800
    journal = work / 'pipeline.jsonl'
    cancel = Path(c['cancel'])

    def event(status, **details):
        with journal.open('a', encoding='utf-8') as f:
            f.write(json.dumps(dict(status=status, unix_s=time.time(), **details)) + '\n')

    def guard():
        assert not cancel.exists(), 'USER_STOP'
        assert time.monotonic() < deadline, 'Pipeline10800s deadline'
        for p, digest in c['source_identities'].items():
            assert sha(p) == digest, 'Pipeline source changed; preserve attempts'

    def command(args, tag, timeout=900):
        guard()
        with (work / (tag + '.stdout.log')).open('xb') as out, (work / (tag + '.stderr.log')).open('xb') as err:
            r = subprocess.run(args, stdout=out, stderr=err, timeout=timeout)
        assert r.returncode == 0, (tag, r.returncode)
        event('STEP_PASSED', step=tag)

    def ps_literal(value):
        return "'" + str(value).replace("'", "''") + "'"

    def remote_ps(code, tag, timeout=900):
        encoded = base64.b64encode(("$ErrorActionPreference='Stop';$ProgressPreference='SilentlyContinue';" + code).encode('utf-16-le')).decode()
        command([*ssh, 'powershell', '-NoProfile', '-EncodedCommand', encoded], tag, timeout)
        return (work / (tag + '.stdout.log')).read_text('utf-8')

    def remote_python(args, tag, timeout=900):
        return remote_ps('& ' + ps_literal(c['remote_python']) + ' ' + ' '.join(ps_literal(x) for x in args) + ";if($LASTEXITCODE -ne 0){throw 'Python step failed'}", tag, timeout)

    def last_json(text):
        for line in reversed(text.splitlines()):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                pass
        raise ValueError('No structured output')

    def retrieve(record, destination, tag):
        target = Path(destination); assert not target.exists()
        zip_path = work / (tag + '.zip')
        command(['scp', *connection['options'], connection['host'] + ':' + record['path'].replace('\\', '/'), str(zip_path)], tag + '_scp')
        assert sha(zip_path) == record['sha256']
        target.mkdir(parents=True)
        with zipfile.ZipFile(zip_path) as z:
            for item in z.infolist():
                path = Path(item.filename)
                assert not item.is_dir() and not path.is_absolute() and '..' not in path.parts
                dst = target / path; assert not dst.exists()
                assert dst.resolve().is_relative_to(target.resolve())
                dst.parent.mkdir(parents=True, exist_ok=True); dst.write_bytes(z.read(item))
        event('RETRIEVED_HASH_VERIFIED', tag=tag, source_zip_sha256=record['sha256'], directory=str(target))

    def local_script(name, args, tag, timeout=900):
        command([sys.executable, str(Path(c['repo']) / 'scripts' / name), *map(str, args)], tag, timeout)

    event('PIPELINE_STARTED', config_sha256=sha(config_path),
        scope='Wait existing17,join45+2 with one missing,independent exact-tone analysis,then4new wavefield controls. No old attempt retry,training or material replacement.')
    try:
        count = 0
        while True:
            guard()
            root = ps_literal(c['continuation_root'] + '/execution')
            output = remote_ps('$r=' + root + ";$e=@(Get-Content -LiteralPath ($r+'/execution.jsonl')|ForEach-Object{ConvertFrom-Json $_});[pscustomobject]@{last=$e[-1].status;terminal_traces=$e[-1].traces;completed=@($e|Where-Object{$_.status -eq 'COMPLETED' -and $_.group}).Count}|ConvertTo-Json -Compress", f'poll_{count:03d}', 45)
            state = last_json(output)
            event('CONTINUATION_OBSERVED', **state)
            assert state['last'] != 'FAILED', 'Continuation failed;do not retry'
            if state['last'] == 'COMPLETED' and state['completed'] == 17 and state.get('terminal_traces') == 17:
                # Terminal event precedes runner exit by a few instructions.
                time.sleep(2)
                break
            count += 1
            time.sleep(45)
        script = c['continuation_root'] + '/repo/scripts/line9_v401_dense_continuation.py'
        record = last_json(remote_python([script, 'combine', '--execution', c['continuation_root'] + '/execution', '--out', c['combined_remote']], 'combine'))
        retrieve(record, c['combined_local'], 'combined')
        source, public, numerical = c['combined_local'], c['dense_public'], c['dense_numerical']
        local_script('analyze_line9_dense_loss_line.py', ['--source', source, '--package', c['dense_package'], '--out', public, '--numerical', numerical], 'dense_analysis')
        local_script('audit_line9_dense_loss_results.py', ['--source', source, '--package', c['dense_package'], '--public', public, '--numerical', numerical, '--out', public + '/independent_audit.json'], 'dense_independent_audit')
        local_script('compare_line9_dense_common_gates.py', ['--package', c['dense_package'], '--public', public, '--numerical', numerical, '--out', public + '/basal_zoom'], 'dense_common_gates')
        local_script('audit_line9_dense_common_gates.py', ['--source', source, '--package', c['dense_package'], '--public', public, '--comparison', public + '/basal_zoom', '--out', public + '/basal_zoom/independent_audit.json'], 'dense_common_audit')
        for name in ['execution.jsonl', 'continuation_execution.jsonl', 'snapshot_manifest.json', 'combined_certificate.json']:
            shutil.copyfile(Path(source) / name, Path(public) / name)
        for local, remote_name in [(Path(public) / 'analysis.json', 'analysis.json'), (Path(public) / 'independent_audit.json', 'sfcw_independent_audit.json')]:
            command(['scp', *connection['options'], str(local), connection['host'] + ':' + c['combined_remote'] + '/' + remote_name], 'upload_' + remote_name)
        event('COMBINED45_PLUS2_INDEPENDENTLY_AUDITED_ONE_MISSING')
        freeze_output = remote_ps('& ' + ps_literal(c['wave_root'] + '/freeze.cmd') + ";if($LASTEXITCODE -ne 0){throw 'Wave freeze failed'}", 'wave_freeze')
        frozen = last_json(freeze_output)
        assert frozen['status'] == 'FROZEN_APPROVED' and frozen['groups'] == 4
        event('FOUR_WAVEFIELD_CONTRACT_FROZEN', frozen=frozen)
        local_script('run_line9_ssh_lease_controller.py', ['--connection', c['connection'], '--log', str(work / 'wave_controller.jsonl'), '--cancel', c['cancel'], '--remote-execution', c['wave_root'] + '/execution', '--contract-sha256', frozen['contract_sha256'], '--remote-run-cmd', c['wave_root'].replace('/', '\\') + '\\run.cmd', '--max-seconds', str(min(7200, int(deadline - time.monotonic()) - 900))], 'wave_run', 7900)
        script = c['wave_root'] + '/repo/scripts/line9_v401_cover_wavefield_run.py'
        remote_python([script, 'collect', '--out', c['wave_root'] + '/execution', '--export-out', c['wave_root'] + '/native_export_r1'], 'wave_collect')
        archive_script = c['wave_root'] + '/repo/scripts/run_line9_bounded_continuation_pipeline.py'
        record = last_json(remote_python([archive_script, 'archive', '--folder', c['wave_root'] + '/native_export_r1', '--zip', c['wave_root'] + '/native_export_r1.zip'], 'wave_archive'))
        retrieve(record, c['wave_local'], 'wave_native')
        local_script('analyze_line9_v401_cover_wavefield.py', ['--source', c['wave_local'], '--out', c['wave_public'], '--numerical', c['wave_numerical']], 'wave_analysis')
        local_script('audit_line9_v401_cover_analysis.py', ['--source', c['wave_local'], '--public', c['wave_public'], '--numerical', c['wave_numerical'], '--out', c['wave_public'] + '/independent_audit.json'], 'wave_independent_audit')
        remote_python([c['wave_root'] + '/repo/scripts/render_line9_v401_cover_wavefield.py', '--execution', c['wave_root'] + '/execution', '--out', c['wave_root'] + '/render_r1'], 'wave_render', 1200)
        record = last_json(remote_python([archive_script, 'archive', '--folder', c['wave_root'] + '/render_r1', '--zip', c['wave_root'] + '/render_r1.zip'], 'wave_movie_archive'))
        retrieve(record, c['wave_public'] + '/wavefield', 'wave_movie')
        for name in ['execution_contract.json', 'preflight_verification.json', 'completed_verification.json', 'execution.jsonl', 'independent_input_audit.json', 'export_receipt.json']:
            shutil.copyfile(Path(c['wave_local']) / name, Path(c['wave_public']) / name)
        for directory in [Path(c['dense_public']), Path(c['wave_public'])]:
            (directory / 'delivery_manifest.json').write_text(json.dumps(dict(status='AUDITED_NUMERICAL_DELIVERY_PHYSICAL_INTERPRETATION_PENDING',
                files={p.relative_to(directory).as_posix(): sha(p) for p in sorted(directory.rglob('*')) if p.is_file()}), indent=2) + '\n', encoding='utf-8')
        event('PIPELINE_COMPLETED_AUDITED_FIGURES_AND_MOVIE', dense_public=c['dense_public'], wave_public=c['wave_public'],
            limits='No unique path/site/real-antenna conclusion from automatic arithmetic audit. Human scientific review and Git delivery remain separate.')
        if c.get('git_branch'):
            repo = Path(c['repo'])
            branch = subprocess.check_output(['git', '-C', str(repo), 'branch', '--show-current'], text=True).strip()
            assert branch == c['git_branch'], 'Branch changed;do not stage'
            assert subprocess.run(['git', '-C', str(repo), 'diff', '--cached', '--quiet']).returncode == 0, 'Other staged work;do not stage'
            report = repo / c['completion_report']; assert not report.exists()
            report.write_text('# 有界接续与波场批次机器验收完成\n\n'
                '17项未运行输入已一次完成；合并45新增+2复用，low_x18675_H1为既有中断缺道。原46项不记为完整完成。\n\n'
                '合并结果已在本机从原生H5独立重算精确501点DFT、逆变换和共同窗。随后四项波场控制一次完成，原几何观察器接收/source逐位一致，500帧六分量/时间/网格在ROG审计；完整场史留ROG，本机有动画、收据和哈希，不冒称本机重读全部场史。\n\n'
                f"- 密采样灰度/审核：`{Path(c['dense_public']).relative_to(repo).as_posix()}`\n"
                f"- 四组灰度/审核/波场GIF：`{Path(c['wave_public']).relative_to(repo).as_posix()}`\n\n"
                '材料仍为研究假设。本文件仅记录自动验收通过，强波路径及三维/真实天线/实测差异仍待科学审查；未训练。\n', encoding='utf-8')
            owned = [str(Path(c['dense_public']).relative_to(repo)), str(Path(c['wave_public']).relative_to(repo)), c['completion_report']]
            command(['git', '-C', str(repo), 'add', '--', *owned], 'git_stage_owned_deliveries')
            staged = subprocess.check_output(['git', '-C', str(repo), 'diff', '--cached', '--name-only', '-z']).decode().split('\0')
            allowed = [p.replace('\\', '/') for p in owned]
            assert all(not p or any(p == a or p.startswith(a + '/') for a in allowed) for p in staged)
            command(['git', '-C', str(repo), 'diff', '--cached', '--check'], 'git_staged_check')
            command(['git', '-C', str(repo), 'commit', '-m', 'research: audit bounded dense continuation and four wavefield controls'], 'git_commit_results')
            command(['git', '-C', str(repo), '-c', 'credential.helper=', '-c', 'credential.helper=!gh auth git-credential', '-c', 'http.version=HTTP/1.1', '-c', 'http.lowSpeedLimit=1000', '-c', 'http.lowSpeedTime=30', 'push', 'origin', c['git_branch']], 'git_push_results', 180)
            head = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()
            remote = subprocess.check_output(['git', '-C', str(repo), 'ls-remote', 'origin', 'refs/heads/' + c['git_branch']], text=True).split()[0]
            assert head == remote
            event('GIT_RESULTS_PUSH_VERIFIED', commit=head, branch=c['git_branch'])
    except BaseException as exc:
        if not cancel.exists():
            cancel.write_text('Owned pipeline failed;stop its active lease controller\n', encoding='utf-8')
        event('PIPELINE_FAILED_PRESERVE_ATTEMPTS_NO_RETRY', error=str(exc))
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='action', required=True)
    run = sub.add_parser('run'); run.add_argument('--config', type=Path, required=True)
    arc = sub.add_parser('archive'); arc.add_argument('--folder', type=Path, required=True); arc.add_argument('--zip', type=Path, required=True)
    a = p.parse_args()
    main(a.config.resolve()) if a.action == 'run' else archive(a.folder, a.zip)
