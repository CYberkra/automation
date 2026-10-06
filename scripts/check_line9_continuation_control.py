"""Check cancellation before/during owned execution without starting a solver."""
import json
import os
from pathlib import Path
import sys
import tempfile
from unittest.mock import Mock, patch

import hs4_station_grid_controls as control


def check():
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory)
        for mode in ['cancel_before','lease_missing','lease_expired','cancel_during']:
            out=root/mode;out.mkdir()
            cancel=out/'USER_STOP';lease=out/'session_heartbeat'
            if mode=='cancel_before': cancel.touch()
            if mode!='lease_missing': lease.touch()
            if mode=='lease_expired': os.utime(lease,(1,1))
            c=dict(status='FROZEN_APPROVED',python=sys.executable,groups=[dict(id='fixture',input=str(out/'profile.in'))],
                max_runs=1,code_identities={},source_identities={},cancel_file=str(cancel),lease_file=str(lease),
                max_lease_age_s=120,gpu_lock=str(out/'lock'))
            path=out/'execution_contract.json';path.write_text(json.dumps(c),encoding='utf-8')
            process=Mock();process.poll.return_value=None
            def launch(*args,**kwargs):
                cancel.touch();return process
            with patch.object(control,'audit'),patch.object(control,'live_resources',return_value={}),\
                 patch.object(control.subprocess,'Popen',side_effect=launch) as popen,\
                 patch.object(control,'terminate_owned_tree') as terminate:
                try: control.run(path)
                except RuntimeError as error:
                    expected='USER_CANCEL_REQUESTED' if mode.startswith('cancel') else 'SESSION_LEASE_EXPIRED'
                    assert str(error)==expected,(mode,str(error))
                else: raise AssertionError('Cancellation ignored')
                if mode=='cancel_during':
                    popen.assert_called_once();terminate.assert_called_once_with(process)
                    events=[json.loads(line) for line in (out/'execution.jsonl').read_text('utf-8').splitlines()]
                    assert events[-1]['status']=='FAILED' and events[-1]['error']=='USER_CANCEL_REQUESTED'
                else:
                    popen.assert_not_called();terminate.assert_not_called()
    print('PASS:4 cancellation/lease checks; no solver started')


if __name__=='__main__': check()
