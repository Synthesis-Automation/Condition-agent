# Codex-ready checklist + pseudocode: reaction typing via Motifs + ΔHandles (no RDKit reaction SMARTS)

This design classifies reactions (e.g., “Suzuki”, “amide formation”) using:

1) **motif typing** for each molecule (from `organic_groups.v1.3.json` + `organic_compounds.v1.3.json`), and  
2) **simple handle deltas** between reactants and products (e.g., B decreases, halide decreases, carbonyl decreases).

It does **not** execute reactions, does **not** require atom mapping, and does **not** use RDKit reaction SMARTS.

---

## 0) Inputs and outputs

### Input

- A list of reaction strings, typically like:
  - `"<reactants_smiles> >> <products_smiles> (optional label)"`
  - where reactants and products are “dot-separated” molecules.
- Example (from your sample list):
  - `Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Simple Ph-Ph)`

### Output per reaction

```json
{
  "rxn": "Brc1...>>c1...",
  "ok": true,
  "predicted": "Suzuki",
  "confidence": 0.86,
  "evidence": {
    "reactant_hits": ["Ar-Br", "Ar-B(OH)2"],
    "product_hits": ["biaryl_like"],
    "delta": { "B": -1, "sp2_leaving_group": -1, "C(=O)O": 0, "C(=O)N": 0 }
  }
}
```

---

## 1) Checklist (implementation tasks)

### 1.1 Molecule motif typing (already needed for vendor catalog typing)

- [ ] Load `organic_groups.v1.3.json` (UTF-8)
- [ ] Load `organic_compounds.v1.3.json` (UTF-8)
- [ ] Compile **detect SMARTS** for each compound motif:
  - template-based: `single_bond = A_smarts + B_smarts`, `via_oxygen = A_smarts + "O" + B_smarts`
  - direct SMARTS: use `compound.smarts`
- [ ] Strip atom-maps from SMARTS for detection (`:1`, `:2`, …)
- [ ] Precompile to `RDKit Mol` queries once; cache in memory

### 1.2 Reaction parsing

- [ ] Parse strings of the form: `left>>right (optional label)`
- [ ] Ignore any trailing `(...)` label if present
- [ ] Reactant molecules: `left.split('.')`
- [ ] Product molecules: `right.split('.')`
- [ ] Convert each SMILES to RDKit Mol; keep invalid SMILES as errors

### 1.3 Feature extraction per reaction

- [ ] For each reactant/product molecule:
  - motif hits: list of compound IDs that match (`Ar-Br`, `Ar-B(OH)2`, etc.)
- [ ] Aggregate motif counts on each side:
  - `reactant_tag_counts`, `product_tag_counts`
- [ ] Compute **handle counters** on each side:
  - element counts: B, Sn, Zn, Mg, Li, halogens (F/Cl/Br/I), etc.
  - functional-group counts (optional but very useful):
    - amide: `C(=O)N`
    - acid: `C(=O)O` (acid/acid salt patterns as you prefer)
    - aldehyde: `[#6X3H1](=O)` etc.
- [ ] Compute deltas: `Δ = product - reactant` for each handle

### 1.4 Rule engine (reaction-type classifier)

- [ ] Implement a prioritized rule list (first-match wins, or score-based)
- [ ] Each rule uses only:
  - required reactant motifs/tags
  - optional required product motifs/tags
  - optional delta constraints (ΔB < 0, ΔHalide < 0, etc.)
- [ ] Return predicted class + confidence + evidence

### 1.5 Evaluation (optional)

- [ ] If labels exist in the sample data, compute accuracy by prefix (e.g., `Suzuki`, `C-N`, `Reductive Amination`)

---

## 2) Core idea: “Motif presence + ΔHandles” beats reactants-only

Reactants-only can be ambiguous (e.g., `RNH2 + RCOOH` could be salt, amidation, etc.).  
Adding **product confirmation** (e.g., product contains amide `C(=O)N`) and **handle depletion** (acid decreases) makes classification robust.

This is the practical interpretation of shorthand like:

- `RNH2 + RCOOH >> RNHCOR`  
as constraints:
- product must contain `C(=O)N` (amide)
- acid-like motif should decrease (or at least not persist as acid)

---

## 3) Pseudocode (Python, RDKit)

### 3.1 Utilities

```python
import json, re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
from rdkit import Chem

RE_MAP = re.compile(r":\d+(?=\])")

def strip_atom_maps(smarts: str) -> str:
    return RE_MAP.sub("", smarts)

def split_label(rxn_line: str) -> str:
    # Remove trailing " (....)" label if present
    # Keep simple: remove last parenthesized suffix only
    return re.sub(r"\s*\([^)]*\)\s*$", "", rxn_line).strip()

def parse_rxn(rxn_line: str) -> Tuple[str, str]:
    s = split_label(rxn_line)
    left, right = s.split(">>", 1)
    return left.strip(), right.strip()

def split_mols(side_smiles: str) -> List[str]:
    side_smiles = side_smiles.strip()
    if not side_smiles:
        return []
    return [x for x in side_smiles.split(".") if x]
```

### 3.2 Compile motif queries from groups + compounds

```python
def load_groups(groups_path: str) -> Dict[str, dict]:
    data = json.load(open(groups_path, "r", encoding="utf-8"))
    return {g["id"]: g for g in data["groups"]}

def load_compounds(compounds_path: str) -> List[dict]:
    data = json.load(open(compounds_path, "r", encoding="utf-8"))
    return data["compounds"]

def compile_compound_smarts(comp: dict, groups: Dict[str, dict]) -> List[str]:
    if "smarts" in comp:
        s = comp["smarts"]
        return s if isinstance(s, list) else [s]

    A = groups[comp["A"]]["smarts"]
    B = groups[comp["B"]]["smarts"]
    tmpl = comp["template"]

    if tmpl == "single_bond":
        return [f"{A}{B}"]
    if tmpl == "via_oxygen":
        return [f"{A}O{B}"]

    raise ValueError(f"Unknown template: {tmpl} for {comp.get('id')}")

def build_query_index(groups_path: str, compounds_path: str) -> Dict[str, List[Chem.Mol]]:
    groups = load_groups(groups_path)
    compounds = load_compounds(compounds_path)

    index: Dict[str, List[Chem.Mol]] = {}
    for comp in compounds:
        smarts_list = compile_compound_smarts(comp, groups)
        qlist = []
        for s in smarts_list:
            detect = strip_atom_maps(s)
            q = Chem.MolFromSmarts(detect)
            if q is None:
                raise ValueError(f"Bad SMARTS for {comp['id']}: {detect}")
            qlist.append(q)
        index[comp["id"]] = qlist
    return index
```

### 3.3 Motif typing for one molecule

```python
def type_molecule(mol: Chem.Mol, query_index: Dict[str, List[Chem.Mol]]) -> List[str]:
    hits = []
    for cid, qlist in query_index.items():
        if any(mol.HasSubstructMatch(q) for q in qlist):
            hits.append(cid)
    return hits
```

---

## 4) Handle counters (tiny, stable list)

Keep this list **small** and stable. Add only when needed.

### 4.1 Element counts

```python
def count_elements(mol: Chem.Mol) -> Dict[str, int]:
    out = {}
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        out[sym] = out.get(sym, 0) + 1
    return out
```

### 4.2 Functional-group counts (optional but high leverage)

Precompile a small set of SMARTS:

- amide: `C(=O)N`
- acid: `C(=O)O` (you may refine to acid vs ester if needed)
- aldehyde: `[#6X3H1](=O)`
- ketone: `[#6X3](=O)[#6X4,#6X3]` (example; tune if needed)

```python
FG_SMARTS = {
    "amide": "C(=O)N",
    "acid_or_ester": "C(=O)O",  # refine later if needed
    "aldehyde": "[#6X3H1](=O)",
}

FG_QUERIES = {k: Chem.MolFromSmarts(v) for k, v in FG_SMARTS.items()}

def count_functional_groups(mol: Chem.Mol) -> Dict[str, int]:
    out = {}
    for k, q in FG_QUERIES.items():
        if q is None:
            continue
        # count matches; set uniquify=True to avoid duplicates
        out[k] = len(mol.GetSubstructMatches(q, uniquify=True))
    return out
```

### 4.3 Aggregate counts per side

```python
def aggregate_side(mols: List[Chem.Mol], query_index) -> dict:
    motif_hits = []
    elem = {}
    fg = {}

    for mol in mols:
        motif_hits.extend(type_molecule(mol, query_index))

        ce = count_elements(mol)
        for k, v in ce.items():
            elem[k] = elem.get(k, 0) + v

        cf = count_functional_groups(mol)
        for k, v in cf.items():
            fg[k] = fg.get(k, 0) + v

    # motif counts
    mcnt = {}
    for h in motif_hits:
        mcnt[h] = mcnt.get(h, 0) + 1

    return {"motifs": mcnt, "elements": elem, "fg": fg}
```

### 4.4 Deltas

```python
def delta_counts(prod: dict, reac: dict) -> dict:
    out = {}
    for domain in ("elements", "fg"):
        d = {}
        keys = set(reac[domain]) | set(prod[domain])
        for k in keys:
            d[k] = prod[domain].get(k, 0) - reac[domain].get(k, 0)
        out[domain] = d
    return out
```

---

## 5) Rule engine (simple, explainable)

### 5.1 Define “motif buckets” used by rules

These are **lists of compound IDs** (from your compound typing) or **group tags** if you add them later.

Start with motif-ID buckets (less code complexity):

```python
BUCKETS = {
  "sp2_electrophile": {"Ar-Br","Ar-Cl","Ar-I","Ar-OTf","Ar-OTs","Ar-OMs",
                      "Vinyl-Br","Vinyl-Cl","Vinyl-I"},
  "organoboron": {"Ar-B(OH)2","Ar-Bpin","Vinyl-B(OH)2","Vinyl-Bpin"},
  "amine": {"R-NH2","R-NHR","Ar-NH2","Ar-NHR"},  # adapt to your real motif IDs
  "acid_or_acyl": {"R-CO2H","Ar-CO2H","R-COCl","Ar-COCl"},  # adapt
}
```

> If you don’t have `R-NH2` etc. as compounds yet, use functional-group SMARTS counters for “amine present” until you add them.

### 5.2 Helper predicates

```python
def has_any(motif_counts: dict, bucket: set) -> bool:
    return any(motif_counts.get(k, 0) > 0 for k in bucket)

def count_any(motif_counts: dict, bucket: set) -> int:
    return sum(motif_counts.get(k, 0) for k in bucket)
```

### 5.3 Rules (prioritized)

Order matters. Put very distinctive ones first.

```python
@dataclass
class Rule:
    name: str
    priority: int
    def match(self, reac: dict, prod: dict, delta: dict) -> Optional[float]:
        ...

class SuzukiRule(Rule):
    def match(self, reac, prod, delta):
        if not has_any(reac["motifs"], BUCKETS["sp2_electrophile"]):
            return None
        if not has_any(reac["motifs"], BUCKETS["organoboron"]):
            return None

        # optional confirmations
        dB = delta["elements"].get("B", 0)
        dHal = sum(delta["elements"].get(x, 0) for x in ["Cl","Br","I","F"])
        score = 0.75
        if dB < 0: score += 0.10
        if dHal < 0: score += 0.05
        return min(score, 0.95)

class AmideFormationRule(Rule):
    def match(self, reac, prod, delta):
        # reactant signature
        if not has_any(reac["motifs"], BUCKETS["amine"]):
            return None
        if not has_any(reac["motifs"], BUCKETS["acid_or_acyl"]):
            return None

        # product confirmation: amide increased/present
        amide_prod = prod["fg"].get("amide", 0)
        if amide_prod <= 0:
            return None

        # optional: acid/ester decreases
        score = 0.75
        if delta["fg"].get("acid_or_ester", 0) < 0:
            score += 0.10
        return min(score, 0.95)
```

Create a list:

```python
RULES = sorted([
  SuzukiRule(name="Suzuki", priority=10),
  AmideFormationRule(name="Amide formation", priority=20),
  # add Negishi/Kumada/Sonogashira/C-N/Reductive amination/SN2 similarly
], key=lambda r: r.priority)
```

Classification:

```python
def classify_reaction(reactant_mols, product_mols, query_index):
    reac = aggregate_side(reactant_mols, query_index)
    prod = aggregate_side(product_mols, query_index)
    delt = delta_counts(prod, reac)

    best = None
    for rule in RULES:
        conf = rule.match(reac, prod, delt)
        if conf is not None:
            best = (rule.name, conf, reac, prod, delt)
            break

    if best is None:
        return {"predicted": "Unknown", "confidence": 0.0, "evidence": {"reactant":reac,"product":prod,"delta":delt}}
    name, conf, reac, prod, delt = best
    return {"predicted": name, "confidence": conf, "evidence": {"reactant":reac,"product":prod,"delta":delt}}
```

---

## 6) Add more reaction families (still “simple”)

Use the same pattern: **reactant signature + optional product confirmation + Δhandles**.

### 6.1 Negishi / Kumada / Stille (very easy)

- Negishi: reactants contain `Zn` (element count) + `sp2_electrophile`, confirm ΔZn < 0 and ΔHalide < 0 (optional)
- Kumada: reactants contain `Mg`, same
- Stille: reactants contain `Sn`, same

### 6.2 Sonogashira

- reactants: `sp2_electrophile` + terminal alkyne (SMARTS `C#C[H]` or an `alkyne_terminal` motif)
- confirm: product still contains alkyne (`C#C`) but terminal H gone (optional), halide decreases

### 6.3 C–N coupling (Buchwald/Ullmann bucket)

- reactants: `sp2_electrophile` + `amine`
- confirm: product contains aryl–N motif (if you have `Ar-N` compound types) OR at least ΔHalide < 0 and amide not formed

### 6.4 Reductive amination

- reactants: carbonyl present (`aldehyde` or `ketone`) + amine present
- confirm: carbonyl decreases (Δaldehyde < 0) AND product contains amine (or N count unchanged but carbonyl decreased)

### 6.5 SN2 alkylation

- reactants: alkyl halide (`R-Br` etc.) + nucleophile (amine/alcohol/thiol)
- confirm: halide decreases and product contains corresponding substituted motif (optional)

---

## 7) Notes specific to your JSON v1.3

- Group SMARTS include atom maps for attachpoints (e.g., `[c:1]`, `[Br:2]`), but:
  - **for motif detection you should strip atom maps** before compiling SMARTS queries.
- `organic_compounds.v1.3.json` omits `anchors` when it is the default:
  - If `anchors` is absent, treat it as `core=A`, `fg=B`.

---

## 8) Minimal “smarts_templates” decision

For this classifier, you only need:

- `single_bond`
- `via_oxygen`

So you can:

- hardcode these 2 templates in code, and
- ignore/remove `smarts_templates.v1.json` entirely for reaction typing.

---

## 9) Deliverable expectation (what Codex should implement)

Create a module like:

- `motif_index.py`
  - `build_query_index(groups_path, compounds_path) -> query_index`
  - `type_molecule(mol, query_index) -> list[motif_id]`

- `reaction_parser.py`
  - `parse_rxn(line) -> (left, right)`
  - `split_mols(side_smiles) -> list[str]`
  - `mols_from_smiles_list(list[str]) -> list[Chem.Mol]`

- `features.py`
  - `aggregate_side(mols, query_index) -> dict`
  - `delta_counts(prod, reac) -> dict`

- `rules.py`
  - `RULES` list and rule classes
  - `classify_reaction(reactant_mols, product_mols, query_index)`

- `cli.py` (optional)
  - reads vendor rxns / dataset, writes JSONL results

---

If you want, I can also generate a small `rules.yaml` format (declarative) so you can update/add rules without touching Python code.
