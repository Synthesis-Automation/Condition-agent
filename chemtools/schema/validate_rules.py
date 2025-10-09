#!/usr/bin/env python3
"""
validate_rules.py  —  tiny stdlib-only validator for "reagents + conditions" files.

Usage:
  python validate_rules.py FILE1.json [FILE2.json ...]

Checks:
  • roles ∈ {metal_source, ligand, base, solvent, additive}
  • every reagent has "amount" (string; "" allowed)
  • amount format:
      - solvent: "xM" or "x–yM" (en-dash or hyphen), or ""
      - others : "x%" or "x–y%" (en-dash or hyphen), or ""
  • conditions.temperature_C/time_h: [number, number]
  • Supports both layouts:
      - entries[] (CN/Suzuki): each entry has reagents/conditions
      - selectors[].then (+ paths) (Amide): each then-block has reagents/conditions
"""
import json, sys, re, os, math

ALLOWED = {"metal_source","ligand","base","solvent","additive"}
DASH = r"[–-]"  # en-dash or hyphen

RE_PCT_SINGLE = re.compile(rf"^\s*\d+(?:\.\d+)?\s*%\s*$", re.I)
RE_PCT_RANGE  = re.compile(rf"^\s*\d+(?:\.\d+)?\s*{DASH}\s*\d+(?:\.\d+)?\s*%\s*$", re.I)
RE_M_SINGLE   = re.compile(rf"^\s*\d+(?:\.\d+)?\s*M\s*$", re.I)
RE_M_RANGE    = re.compile(rf"^\s*\d+(?:\.\d+)?\s*{DASH}\s*\d+(?:\.\d+)?\s*M\s*$", re.I)

def ok_amount(role, s):
    if s is None:
        return False, "amount missing (must be string; '' allowed)"
    if not isinstance(s, str):
        return False, "amount must be a string"
    t = s.strip()
    if role == "solvent":
        if t == "": return True, ""
        if RE_M_SINGLE.match(t) or RE_M_RANGE.match(t): return True, ""
        return False, f"solvent amount must be 'xM' or 'x{DASH}yM' (got: {s!r})"
    # non-solvent
    if t == "": return True, ""
    if RE_PCT_SINGLE.match(t) or RE_PCT_RANGE.match(t): return True, ""
    return False, f"amount must be 'x%' or 'x{DASH}y%' (got: {s!r})"

def ok_range2(v):
    if v is None: return True, ""  # optional
    if not isinstance(v, list) or len(v)!=2: return False, "must be [min,max]"
    if not all(isinstance(x,(int,float)) for x in v): return False, "both values must be numbers"
    return True, ""

def validate_reagents(reagents, ctx, errs):
    if not isinstance(reagents, list):
        errs.append(f"{ctx}: reagents must be an array")
        return
    for i, r in enumerate(reagents):
        path = f"{ctx}.reagents[{i}]"
        if not isinstance(r, dict):
            errs.append(f"{path}: each reagent must be an object")
            continue
        name = r.get("name","")
        role = r.get("role", None)
        amt  = r.get("amount", "")
        if not isinstance(name, str):
            errs.append(f"{path}.name: must be string")
        if role not in ALLOWED:
            errs.append(f"{path}.role: {role!r} not in {sorted(ALLOWED)}")
        ok, msg = ok_amount(role, amt)
        if not ok:
            errs.append(f"{path}.amount: {msg}")

def validate_conditions(conditions, ctx, errs):
    if not isinstance(conditions, dict):
        errs.append(f"{ctx}: conditions must be an object")
        return
    for key in ("temperature_C", "time_h"):
        if key in conditions:
            ok, msg = ok_range2(conditions[key])
            if not ok:
                errs.append(f"{ctx}.{key}: {msg}")

def walk_then_block(tb, ctx, errs):
    if "reagents" in tb:
        validate_reagents(tb["reagents"], f"{ctx}", errs)
    else:
        errs.append(f"{ctx}: missing 'reagents'")
    if "conditions" in tb:
        validate_conditions(tb["conditions"], f"{ctx}.conditions", errs)
    else:
        errs.append(f"{ctx}: missing 'conditions'")
    # nested paths
    if "paths" in tb and isinstance(tb["paths"], list):
        for j, p in enumerate(tb["paths"]):
            if isinstance(p, dict):
                walk_then_block(p, f"{ctx}.paths[{j}]", errs)

def validate_file(path):
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    errs = []

    if "entries" in data and isinstance(data["entries"], list):
        for idx, e in enumerate(data["entries"]):
            if not isinstance(e, dict):
                errs.append(f"entries[{idx}]: must be object")
                continue
            if "reagents" not in e: errs.append(f"entries[{idx}]: missing 'reagents'")
            else: validate_reagents(e["reagents"], f"entries[{idx}]", errs)
            if "conditions" not in e: errs.append(f"entries[{idx}]: missing 'conditions'")
            else: validate_conditions(e["conditions"], f"entries[{idx}].conditions", errs)

    if "selectors" in data and isinstance(data["selectors"], list):
        for si, s in enumerate(data["selectors"]):
            tb = s.get("then")
            if isinstance(tb, dict):
                walk_then_block(tb, f"selectors[{si}].then", errs)
            else:
                errs.append(f"selectors[{si}]: missing 'then' object")

    return errs

def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    exit_code = 0
    for path in argv[1:]:
        try:
            errs = validate_file(path)
            if errs:
                print(f"[FAIL] {path}")
                for e in errs:
                    print("  -", e)
                exit_code = 1
            else:
                print(f"[OK]   {path}")
        except Exception as ex:
            print(f"[ERROR] {path}: {ex}")
            exit_code = 1
    return exit_code

if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
