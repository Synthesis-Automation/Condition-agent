# Attachment-Point SMARTS Templating — Codex Evaluation Instructions

This document describes how to **evaluate** the attachment-point templating approach used in the Suzuki/Buchwald POC.
Goal: confirm that template expansion is **correct**, **stable**, and **useful** for maintaining cross-coupling SMARTS families.

## Inputs (from the ZIP)

You should have these files available:

- `smarts_templates.suzuki_buchwald.poc.json`  
  Contains **fragments** (with attachment labels) + **templates** + an expansion plan.
- `calculable_features.atomic.suzuki_buchwald.poc.json`  
  Contains the **generated atomic SMARTS** used at runtime.

Optional:

- `build_suzuki_buchwald_templates_poc.py`  
  Demonstration script of how expansion was performed (may be adapted).

---

## 1) What “attachment-point templating” means (expected behavior)

### Fragments

A fragment is a SMARTS snippet with one or more labeled attachment atoms, e.g.:

- `Ar = "[c:1]"`
- `Vinyl = "[C;X3:1]=[C;X3]"`

The `:1` label indicates the intended attachment site.

### Templates

A template is a formatting rule that combines fragments into a full SMARTS:

- `core_halide: "{core}{lg}"`
- `core_sulfonate: "{core}O{tail}"`
- `core_boron: "{core}{boron}"`

### Expansion

Expansion enumerates combinations:

- `(core ∈ {Ar, Vinyl}) × (lg ∈ {Cl, Br, I})` for electrophiles
- `(core ∈ {Ar, Vinyl}) × (boron ∈ {B(OH)2, Bpin})` for Suzuki partners
- plus pseudohalides `(tail ∈ {OTf, OTs, OMs})` if present

**Key property:** templating runs at **build-time** only.
At runtime, the system uses the compiled atomic SMARTS in `calculable_features.atomic...json`.

---

## 2) Implement an evaluator (what to build)

Create a Python module, e.g. `template_eval.py`, that implements:

### A) Parse templates file

Load `smarts_templates...json` and extract:

- `fragments`: dict[name -> smarts]
- `templates`: dict[name -> format_str]
- `expansions`: list of expansion specifications (or infer from structure)

### B) Expand templates deterministically

Implement:

- `expand_templates(fragments, templates, expansions) -> list[GeneratedPattern]`

Where `GeneratedPattern` includes:

- `feature_token`
- `smarts`
- metadata: which template + which fragments

**Determinism requirement:** repeated runs produce identical outputs (order and content).

### C) Compile SMARTS with RDKit

For each generated SMARTS:

- `Chem.MolFromSmarts(smarts)` must not return `None`
- report compile failures with token + smarts + reason

### D) Cross-check against atomic features

Load `calculable_features.atomic...json` and verify:

1) Every generated feature token exists in atomic features
2) Every generated SMARTS string matches the atomic feature’s `detect.smarts_any` (exact string match preferred)
3) No generated token appears twice
4) No generated SMARTS string appears under two different tokens (unless explicitly allowed)

Output a summary report.

---

## 3) Functional validation with small molecule tests

Create a table of **test molecules** and expected matches.

### Electrophile tests

- Bromobenzene: `Brc1ccccc1` should match `aryl_bromide_present`
- Chlorobenzene: `Clc1ccccc1` should match `aryl_chloride_present`
- Iodobenzene: `Ic1ccccc1` should match `aryl_iodide_present`

### Organoboron tests (Suzuki)

- Phenylboronic acid: `OB(O)c1ccccc1` should match `aryl_boronic_acid_present`
- Phenyl Bpin: `B1OC(C)(C)C(C)(C)O1c2ccccc2` (or an equivalent Bpin SMILES) should match `aryl_bpin_present`

### Vinyl tests (optional; can be tricky)

- Vinyl bromide: `C=CBr` should match `vinyl_bromide_present` if vinyl patterns are included

Implement:

- `test_smiles_matches(smiles, expected_tokens)`
- evaluate using RDKit `HasSubstructMatch` on the generated atomic SMARTS

**Tip:** allow a small set of tests to be marked `xfail` if representation variability is expected (e.g., pseudohalides or boronate SMILES variants).

---

## 4) Metrics to report (Codex should produce a short evaluation report)

Produce a Markdown report, e.g. `TEMPLATE_EVAL_REPORT.md`, including:

### Build integrity

- # fragments, # templates, # expansions

- # generated patterns

- # SMARTS compile successes / failures

### Consistency with atomic features

- % generated tokens found in atomic features
- % generated SMARTS strings found in atomic features
- duplicates detected (tokens and SMARTS)

### Functional match tests

- pass/fail summary per test molecule
- list of unexpected matches (false positives) if any

---

## 5) Common pitfalls and what to check

1) **Aromatic vs aliphatic representation**

- `[c]` matches aromatic carbon; `C` matches aliphatic carbon
- ensure your “Ar–X” uses aromatic carbon patterns (e.g. `c[Br]`-like)

2) **Overly generic patterns**

- avoid patterns that match unintended motifs (e.g. `[#6][Br]` matches alkyl bromides too)

3) **SMARTS string normalization**

- keep the generated SMARTS strings stable (no unnecessary whitespace)
- if you must canonicalize, do it consistently (but exact string matching is simplest)

4) **Attachment label collisions**

- if templates combine fragments with multiple `:1` labels, ensure labels are unique or intentionally shared
- for simple two-piece concatenation (`{core}{lg}`) this is usually safe

5) **RDKit SMARTS compilation**

- compilation failures typically indicate malformed brackets or invalid valence expressions

---

## 6) Deliverables (Codex should output)

1) `template_eval.py` (or a small package) implementing expansion + validation  
2) `pytest` tests:
   - compilation tests
   - atomic-feature cross-check tests
   - functional molecule match tests
3) `TEMPLATE_EVAL_REPORT.md` summarizing results

Everything should be UTF-8.
