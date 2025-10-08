from __future__ import annotations

import argparse
import os
from typing import List
import sys

# Ensure project root on path when running from source
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


def _gather_reaction_smiles(source: str) -> List[str]:
    # Source options: 'auto' (default) uses chemtools.precedent._load
    if source == "auto":
        from chemtools import precedent as _prec
        rows = _prec._load()
        rs = []
        for r in rows:
            s = r.get("reaction_smiles")
            if s:
                rs.append(str(s))
        # Deduplicate, keep order
        seen = set()
        out: List[str] = []
        for s in rs:
            if s not in seen:
                seen.add(s)
                out.append(s)
        return out
    # If a path is provided, attempt to read JSONL with 'reaction_smiles'
    path = source
    rsmi: List[str] = []
    if os.path.isdir(path):
        for name in os.listdir(path):
            if not name.lower().endswith('.jsonl'):
                continue
            rsmi.extend(_read_jsonl(os.path.join(path, name)))
    elif os.path.isfile(path):
        rsmi.extend(_read_jsonl(path))
    else:
        raise SystemExit(f"Source not found: {source}")
    # Deduplicate
    seen = set()
    out: List[str] = []
    for s in rsmi:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _read_jsonl(path: str) -> List[str]:
    rs: List[str] = []
    import json
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except Exception:
                continue
            # Try common fields
            s = rec.get('reaction_smiles')
            if not s:
                smiles = rec.get('smiles') or {}
                r = smiles.get('reactants') or ''
                p = smiles.get('products') or ''
                s = f"{r}>>{p}" if (r or p) else ''
            if s:
                rs.append(str(s))
    return rs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Precompute DRFP fingerprints and save as NPZ")
    ap.add_argument('--source', default='auto', help="'auto' to use chemtools.precedent dataset, or path to JSONL file/dir")
    ap.add_argument('--out', default=os.path.join('artifacts', 'drfp_index.npz'), help='Output NPZ path')
    ap.add_argument('--n-bits', type=int, default=4096)
    ap.add_argument('--radius', type=int, default=3)
    args = ap.parse_args(argv)

    from chemtools.reaction_similarity import drfp_available, encode_drfp
    import numpy as np

    if not drfp_available():
        raise SystemExit("DRFP is not installed. Please pip install drfp")

    rs = _gather_reaction_smiles(args.source)
    if not rs:
        raise SystemExit("No reaction SMILES found from source")

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)

    fps = []
    keys = []
    failed = 0
    for i, s in enumerate(rs, 1):
        fp = encode_drfp(s, n_bits=args.n_bits, radius=args.radius)
        if fp is None:
            failed += 1
            continue
        fps.append(fp)
        keys.append(s)
        if i % 1000 == 0:
            print(f"encoded {i}/{len(rs)}...", flush=True)

    if not fps:
        raise SystemExit("All DRFP encodings failed; aborting")

    fps_arr = np.vstack([np.array(x, dtype='uint8') for x in fps])
    keys_arr = np.array(keys, dtype=object)

    np.savez_compressed(args.out, fps=fps_arr, keys=keys_arr, n_bits=np.array(args.n_bits), radius=np.array(args.radius))
    print(f"Saved {len(keys)} fingerprints to {args.out} (failed: {failed})")
    return 0


if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
