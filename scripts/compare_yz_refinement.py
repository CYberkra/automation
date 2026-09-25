"""Compare archived YZ grids into a new directory; no solver invocation."""
import argparse
import json
from pathlib import Path
import numpy as np
from research_spectral_metrics import metrics

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New directory; must not exist')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = ROOT/'artifacts/research_checks/2026-09-25_yz_fine_analysis'
    old = ROOT/'artifacts/research_checks/2026-09-25_yz_layers_analysis'
    coarse = json.loads((old/'results.json').read_text(encoding='utf-8'))
    fine = json.loads((root/'fine/results.json').read_text(encoding='utf-8'))
    with np.load(old/'arrays.npz') as a, np.load(root/'fine/arrays.npz') as b:
        assert np.array_equal(a['frequency_Hz'], b['frequency_Hz']) and np.array_equal(a['covered'], b['covered'])
        difference = metrics(b['COVER_200'], a['COVER_200'])
    first, second = coarse['comparisons']['COVER_200'], fine['comparisons']['COVER_200']
    result = {'coarse_2p5cm': first, 'fine_1p25cm': second,
              'coarse_error_over_fine_error': {k: first[k]/second[k] for k in first},
              'fine_vs_coarse': difference,
              'claim_limit': 'Two levels support or contradict refinement improvement; no fitted order, no deep-target accuracy claim, no correction applied.'}
    with (args.output/'comparison.json').open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
