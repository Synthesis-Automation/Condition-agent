#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
POC builder: SMARTS templating with attachment points for Suzuki + Buchwald feature families.

This script regenerates the "templated" atomic SMARTS features:
- aromatic/vinyl halides
- aromatic/vinyl sulfonate pseudohalides (OTf/OTs/OMs)
- aromatic/vinyl organoboron partners (boronic acid, Bpin)

Outputs a JSON list of atomic feature objects (token/type/description/detect).
"""
import json
from pathlib import Path

FRAG = {
  "AromaticCore": "[a:1]",
  "VinylCore": "[C;X3:1]=[C;X3]",
  "Cl": "[Cl]",
  "Br": "[Br]",
  "I": "[I]",
  "OTfTail": "S(=O)(=O)C(F)(F)F",
  "OTsTail": "S(=O)(=O)c1ccc(C)cc1",
  "OMsTail": "S(=O)(=O)C",
  "BoronicAcid": "[B](O)O",
  "Bpin": "[B]1OC(C)(C)C(C)(C)O1"
}

TPL = {
  "core_halide": "{core}{lg}",
  "core_sulfonate": "{core}O{tail}",
  "core_boron": "{core}{boron}"
}

def feat(token, desc, smarts, template=None, frags=None):
  out = {"token": token, "type": "bool", "description": desc, "detect": {"smarts_any": [smarts]}}
  if template:
    out["generated_from"] = {"template": template, "fragments": frags or []}
  return out

def main(out_path: str):
  atomic = []
  for core_name, prefix in [("AromaticCore","aromatic"), ("VinylCore","vinyl")]:
    core = FRAG[core_name]
    for lg_name, lg_label in [("Cl","chloride"), ("Br","bromide"), ("I","iodide")]:
      smarts = TPL["core_halide"].format(core=core, lg=FRAG[lg_name])
      atomic.append(feat(f"{prefix}_{lg_label}_present", f"{prefix} sp2 {lg_label}", smarts, "core_halide", [core_name, lg_name]))
    for tail_name, tail_label in [("OTfTail","triflate"), ("OTsTail","tosylate"), ("OMsTail","mesylate")]:
      smarts = TPL["core_sulfonate"].format(core=core, tail=FRAG[tail_name])
      atomic.append(feat(f"{prefix}_{tail_label}_present", f"{prefix} sp2 {tail_label} pseudohalide", smarts, "core_sulfonate", [core_name, tail_name]))
    for bor_name, bor_label in [("BoronicAcid","boronic_acid"), ("Bpin","bpin")]:
      smarts = TPL["core_boron"].format(core=core, boron=FRAG[bor_name])
      atomic.append(feat(f"{prefix}_{bor_label}_present", f"{prefix} organoboron partner: {bor_label}", smarts, "core_boron", [core_name, bor_name]))
  Path(out_path).write_text(json.dumps(atomic, ensure_ascii=False, indent=2), encoding="utf-8")
  print(f"Wrote {len(atomic)} atomic features to {out_path}")

if __name__ == "__main__":
  main("atomic_features.templated.poc.json")
