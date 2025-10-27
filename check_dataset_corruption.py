#!/usr/bin/env python3
"""Check if dataset files are corrupted."""

import os

files = ['Suzuki.jsonl', 'C_N_Coupling.jsonl', 'C_O_Coupling.jsonl', 'Amide_formation.jsonl']

for fname in files:
    fpath = os.path.join('data/reaction_dataset', fname)
    if os.path.exists(fpath):
        with open(fpath, 'rb') as f:
            first_bytes = f.read(100)
            print(f"{fname}:")
            print(f"  First 100 bytes: {first_bytes}")
            print(f"  Starts with '{{'? {first_bytes.startswith(b'{')}")
            print()
