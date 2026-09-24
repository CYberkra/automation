"""Replay the archived air-only recipe through the official V4 SFCW API; no solve."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np
from gprMax.toolboxes.SFCW.processing import (load_source,load_receiver,
    direct_frequency_response,reconstruct_time_response,write_sfcw_output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('recipe',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=json.loads(a.recipe.read_text(encoding='utf-8'))
    if a.output.exists():raise SystemExit('Refusing overwrite')
    root=Path(__file__).resolve().parents[1];raw=root/r['raw_path']
    if hashlib.sha256(raw.read_bytes()).hexdigest()!=r['raw_sha256']:raise ValueError('Raw hash mismatch')
    source=load_source(raw);receiver=load_receiver(raw,receiver_path='name:measurement',component='Ex')
    receiver=replace(receiver,samples=receiver.samples[:r['receiver_samples_used']])
    response=direct_frequency_response(source,receiver,np.linspace(r['frequency_start_Hz'],r['frequency_stop_Hz'],r['tones']),tail_taper_fraction=r['tail_taper_fraction'])
    if not np.all(response.source_valid):raise ValueError('Invalid source bins')
    time=reconstruct_time_response(response,window=r['window'],zero_pad_factor=r['zero_pad'],time_shift=r['time_shift'])
    write_sfcw_output(a.output,response,time)


if __name__=='__main__':main()
