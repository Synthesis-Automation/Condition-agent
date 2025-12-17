#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POC classifier for Suzuki + Buchwald using:
- atomic SMARTS features
- derived boolean features
- simple reaction-type rules

Requires RDKit.
Usage:
  python classify_suzuki_buchwald_poc.py --reactant "Brc1ccccc1" --reactant "OB(O)c1ccccc1"
"""
import argparse
import json
from pathlib import Path

def load_features(path_atomic: str, path_derived: str):
    atomic = json.loads(Path(path_atomic).read_text(encoding="utf-8"))
    derived = json.loads(Path(path_derived).read_text(encoding="utf-8"))
    return atomic, derived

def compile_smarts(atomic):
    from rdkit import Chem
    compiled = []
    for f in atomic:
        pats = []
        for s in f["detect"]["smarts_any"]:
            m = Chem.MolFromSmarts(s)
            if m is None:
                raise ValueError(f"Bad SMARTS for {f['token']}: {s}")
            pats.append(m)
        compiled.append((f["token"], pats))
    return compiled

def eval_atomic(mol, compiled):
    out = {}
    for token, pats in compiled:
        out[token] = any(mol.HasSubstructMatch(p) for p in pats)
    return out

def eval_derived(values, derived_list):
    # simple topo via repeated passes (POC)
    pending = {d["token"]: d for d in derived_list}
    changed = True
    while changed and pending:
        changed = False
        for token in list(pending.keys()):
            d = pending[token]
            expr = d["derive"]
            def get(t): return bool(values.get(t, False))
            ok = None
            if "all_of" in expr:
                ok = all(get(t) for t in expr["all_of"])
            if "any_of" in expr:
                ok = (ok if ok is not None else True) and any(get(t) for t in expr["any_of"])
            if "none_of" in expr:
                ok = (ok if ok is not None else True) and (not any(get(t) for t in expr["none_of"]))
            if ok is None:
                raise ValueError(f"Unsupported derive expr in {token}: {expr}")
            # If all dependencies are already known, we can finalize (POC heuristic)
            deps = set(expr.get("all_of", []) + expr.get("any_of", []) + expr.get("none_of", []))
            if deps.issubset(values.keys()):
                values[token] = ok
                pending.pop(token)
                changed = True
    # Any remaining pending tokens: best-effort evaluate with missing treated as False
    for token, d in pending.items():
        expr = d["derive"]
        def get(t): return bool(values.get(t, False))
        ok = True
        if "all_of" in expr: ok = ok and all(get(t) for t in expr["all_of"])
        if "any_of" in expr: ok = ok and any(get(t) for t in expr["any_of"])
        if "none_of" in expr: ok = ok and (not any(get(t) for t in expr["none_of"]))
        values[token] = ok
    return values

def mixture_or(feature_dicts):
    out = {}
    for fd in feature_dicts:
        for k, v in fd.items():
            out[k] = out.get(k, False) or bool(v)
    return out

def match_reaction_types(mix_features, reaction_types):
    hits = []
    for rt in reaction_types:
        when = rt["when"]
        ok = True
        for t in when.get("all_of", []):
            ok = ok and bool(mix_features.get(t, False))
        for t in when.get("any_of", []):
            ok = ok and bool(mix_features.get(t, False))
        for t in when.get("none_of", []):
            ok = ok and (not bool(mix_features.get(t, False)))
        if ok:
            hits.append(rt["id"])
    return hits

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reactant", action="append", required=True, help="Reactant SMILES (repeatable)")
    ap.add_argument("--atomic", default="calculable_features.atomic.suzuki_buchwald.poc.json")
    ap.add_argument("--derived", default="calculable_features.derived.suzuki_buchwald.poc.json")
    ap.add_argument("--reaction-types", default="reaction_types.suzuki_buchwald.poc.json")
    args = ap.parse_args()

    from rdkit import Chem

    atomic, derived = load_features(args.atomic, args.derived)
    compiled = compile_smarts(atomic)

    per = []
    for smi in args.reactant:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Bad SMILES: {smi}")
        feats = eval_atomic(mol, compiled)
        feats = eval_derived(feats, derived)
        per.append(feats)

    mix = mixture_or(per)

    rtypes = json.loads(Path(args.reaction_types).read_text(encoding="utf-8"))
    hits = match_reaction_types(mix, rtypes)

    print("Matched reaction types:", hits)
    # Print a few key tokens
    keys = ["sp2_electrophile_present","suzuki_boron_partner_present","buchwald_amine_partner_present"]
    print("Key mixture tokens:", {k: mix.get(k, False) for k in keys})

if __name__ == "__main__":
    main()
