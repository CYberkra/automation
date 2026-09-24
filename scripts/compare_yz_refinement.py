"""Compare the two archived YZ cover grid levels; no solver invocation."""
import json
from pathlib import Path
import numpy as np
from analyze_vertical_refinement import metrics

root=Path('artifacts/research_checks/2026-09-25_yz_fine_analysis')
old=Path('artifacts/research_checks/2026-09-25_yz_layers_analysis')
coarse=json.loads((old/'results.json').read_text(encoding='utf-8'))
fine=json.loads((root/'fine/results.json').read_text(encoding='utf-8'))
with np.load(old/'arrays.npz') as a,np.load(root/'fine/arrays.npz') as b:
 assert np.array_equal(a['frequency_Hz'],b['frequency_Hz']) and np.array_equal(a['covered'],b['covered'])
 difference=metrics(b['COVER_200'],a['COVER_200'])
first=coarse['comparisons']['COVER_200'];second=fine['comparisons']['COVER_200']
result={'coarse_2p5cm':first,'fine_1p25cm':second,'coarse_error_over_fine_error':{k:first[k]/second[k] for k in first},'fine_vs_coarse':difference,
 'claim_limit':'Two levels support or contradict refinement improvement; no fitted order, no deep-target accuracy claim, no correction applied.'}
(root/'comparison.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result))
