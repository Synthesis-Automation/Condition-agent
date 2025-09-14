# Role-Aware Featurization – Implementation Plan

## Scope & Goals
- Build a **role-aware** featurizer for amines, alcohols, and aryl halides that outputs a **single consistent vector** per reactant with presence masks.
- Stack three layers: **global descriptors**, **role-specific features**, **reactive-site centered ECFP**.

## Inputs / Outputs
- **Input:** RDKit `Mol` (optionally reaction SMILES with atom-maps).
- **Output:** `{ "vector": np.ndarray, "fields": List[str], "masks": Dict[str, int], "meta": Dict }`.

## Directory Structure
```
chem_feats/
  __init__.py
  registry.yaml          # feature registry & defaults
  smarts.py              # role SMARTS & Finders
  global_feats.py        # role-agnostic RDKit descriptors
  role_feats/
    amine.py
    alcohol.py
    aryl_halide.py
  fingerprints.py        # centered ECFP helpers
  assemble.py            # concat, order, masks, sentinel fill
  tests/
    test_find_roles.py
    test_vectors_shape.py
    test_reference_cases.py
```

## Feature Registry (minimal example)
```yaml
global:
  - {name: global.MW, type: float, default: 0.0}
  - {name: global.logP, type: float, default: 0.0}
amine:
  - {name: amine.present, type: int, default: 0}
  - {name: amine.class_ps3, type: int, default: -1}
alcohol:
  - {name: alcohol.present, type: int, default: 0}
aryl_halide:
  - {name: aryl_halide.present, type: int, default: 0}
fingerprints:
  amine: {bits: 512, radius: 2}
  alcohol: {bits: 512, radius: 2}
  aryl_halide: {bits: 512, radius: 2}
```

## Core Steps
1. **Role detection**
   - Prefer atom-maps; else match SMARTS per role (amine N, alcohol O, aryl ipso C + X).
2. **Global features**
   - MW, logP, TPSA, RB, HBD/HBA, aromatic rings, heteroatom count.
3. **Role-specific features**
   - **Amines:** degree (1°/2°/3°), aniline flag, α-branching, formal charge.
   - **Alcohols:** 1°/2°/3°, benzylic/allylic flags.
   - **Aryl halides:** halide one-hot (F/Cl/Br/I), ortho-block (0–2), ipso degree.
4. **Reactive-site fingerprints**
   - Centered Morgan (radius 2–3) using `fromAtoms=[center_idxs]`.
5. **Assemble vector**
   - Fixed field order from `registry.yaml`; fill missing with sentinels; emit `has_amine/alcohol/aryl_halide` masks.
6. **Validation**
   - Shape consistency across roles; deterministic outputs; speed baseline.

## Public API (Python)
```python
from chem_feats import featurize_mol, featurize_reaction
vec = featurize_mol(mol, roles=["amine","aryl_halide"])
rxn_vecs = featurize_reaction(rxn_smiles)  # returns dict per role + masks
```

## Reference Test Cases
- Aryl halides: p-chloroanisole (ortho_block=0), 2,6-di-tBu-bromobenzene (ortho_block=2).
- Amines: aniline vs tert-butylamine (class + aniline flag).
- Alcohols: benzyl alcohol vs t-BuOH (benzylic/tertiary flags).
- Mixed: `Ar–Cl + aniline` and `Ar–Br + t-BuNH2`.

## Milestones
- **M1 (Day 1–2):** SMARTS finders + global features + registry scaffold.
- **M2 (Day 3–4):** Role-specific features + centered ECFP + assembler.
- **M3 (Day 5):** Tests, SHAP sanity on small model, integration into ConditionCore pipeline.

## Acceptance Criteria
- Fixed-length vectors for any input; masks correctly reflect roles.
- All tests pass; <5 ms per molecule on average; fields documented from registry.
