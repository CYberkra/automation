"""Timing diagnostics; correlation cannot identify an interface or PML echo."""
import argparse
from pathlib import Path
from hs4t2d_physics_diagnostic import DEFAULT_CAPSULE, write_analysis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capsule', type=Path, default=DEFAULT_CAPSULE)
    parser.add_argument('--out', type=Path, required=True, help='new JSON; never overwrite historical evidence')
    args = parser.parse_args()
    write_analysis(args.capsule, args.out)


if __name__ == '__main__':
    main()
