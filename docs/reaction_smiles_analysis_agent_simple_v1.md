# Reaction-SMILES Analysis Agent — Simple Implementation (v1)

Goal: Given a **reaction SMILES** (`reactants>>products`), produce a **reliable** analysis (roles, bond changes, coarse reaction class, short mechanistic sketch) using:

- **Deterministic cheminformatics** for “ground truth” (RDKit + atom-mapping + bond diff)
- **GPT-5.2** for interpretation, constrained by **Structured Outputs**

Design rule: **Tools compute facts; GPT explains facts.**  
If tools can’t compute reliable facts → **return warnings / low confidence**, don’t invent.

---

## 1) Minimal Architecture

### Dataflow

1. **Parse & pre-check** reaction SMILES (RDKit)
2. **Conservative cleanup** (canonicalize, light standardize; keep raw)
3. **Detect spectators** (e.g., `Cl`, `[Na+]`) and optionally drop from mapping
4. **Atom mapping** (`rxnmapper`)
5. **Bond-change extraction** (formed/broken/order-change) from mapped reaction
6. **GPT call** with **structured JSON output** referencing computed bond changes
7. **Validation + QC gating** (if mapping/bond diff weak → degrade gracefully)

### Core outputs

- `rxn_smiles_raw`
- `rxn_smiles_clean`
- `mapped_rxn_smiles` (optional)
- `bond_changes[]` (tool-derived)
- `interpretation{...}` (GPT-derived, schema constrained)
- `warnings[]` + `confidence`

---

## 2) Dependencies (simple set)

- Python 3.10+
- `rdkit`
- `rxnmapper` (for atom mapping)
- `openai` (API client)
- Optional: `pydantic` (nice validation + schema generation)

> Keep everything UTF-8 (default in Python 3; explicitly set when writing files).

---

## 3) Reliability Policies (important)

### 3.1 Conservative cleanup only

Do **not** aggressively:

- tautomer-canonicalize
- strip stereochem
- neutralize everything

Instead:

- keep `rxn_smiles_raw`
- produce `rxn_smiles_clean` with **canonical SMILES** and minimal normalization
- store `standardization_actions[]` and `parse_warnings[]`

### 3.2 Spectator handling

Many datasets include salts/counterions (e.g., `.Cl` for amine·HCl).  
For mapping and bond diff:

- Detect simple ions as **spectators**
- Prefer mapping on “organic fragments”
- Keep spectators listed for transparency

Example input:
`CNC.Cl.Clc1cc(Cl)...>>CN(C)c1cc(N(C)C)...`

Likely spectators:

- `Cl` (counterion)
Likely nucleophile:
- amine fragment (`CNC` canonicalizes to `CN(C)C`)

### 3.3 QC gating

Never “force” a confident mechanism if the tool facts are weak.

Recommended QC checks:

- RDKit parse success for all major components
- Mapping yields high heavy-atom coverage (product atoms mapped back)
- Bond changes are not absurd (e.g., dozens of random changes)

If QC fails:

- return `interpretation.confidence <= 0.3`
- add warnings: `"mapping_failed"`, `"bond_diff_unreliable"`
- optionally skip mechanism steps

---

## 4) Output Schema (keep it thin)

### 4.1 Why thin schema

A huge schema doesn’t reduce model intelligence, but it increases:

- token overhead
- forced filling of uncertain fields
- truncation risk

So v1 uses:

- stable top-level structure
- small enums
- an `events[]` list for single-step or tandem reactions

### 4.2 v1 JSON structure (conceptual)

```json
{
  "schema_version": "reaction_analysis.v1",
  "input": {
    "rxn_smiles_raw": "...",
    "rxn_smiles_clean": "...",
    "parse_warnings": [],
    "standardization_actions": [],
    "spectators": []
  },
  "tool_facts": {
    "mapped_rxn_smiles": null,
    "mapping_qc": {
      "ok": true,
      "notes": []
    },
    "bond_changes": [
      { "id": "BC1", "change": "broken", "a1": 12, "a2": 57, "bond": "single" },
      { "id": "BC2", "change": "formed", "a1": 12, "a2": 89, "bond": "single" }
    ],
    "reaction_center_atoms": [12, 57, 89]
  },
  "interpretation": {
    "overall_class": "nucleophilic_substitution",
    "tags": ["SNAr"],
    "roles": {
      "electrophile": "…",
      "nucleophile": "…",
      "leaving_group": "…"
    },
    "events": [
      {
        "event_id": "E1",
        "event_type": "substitution",
        "bond_change_refs": ["BC1", "BC2"],
        "short_rationale": "C–Cl broken, C–N formed at aromatic carbon",
        "confidence": 0.78
      }
    ],
    "mechanism_summary": [
      "…",
      "…"
    ],
    "warnings": [],
    "confidence": 0.72
  }
}
```

### 4.3 Minimal enums

`overall_class` enum (coarse):

- `cross_coupling`
- `nucleophilic_substitution`
- `electrophilic_addition`
- `nucleophilic_addition`
- `elimination`
- `condensation`
- `cycloaddition`
- `rearrangement`
- `oxidation`
- `reduction`
- `protection_deprotection`
- `other`

`event_type` enum (small set):

- `bond_formation`
- `bond_cleavage`
- `bond_order_change`
- `substitution`
- `addition`
- `elimination`
- `rearrangement`
- `redox`
- `other`

---

## 5) Implementation Plan (functions)

### 5.1 `clean_reaction_smiles(rxn_smiles) -> CleanReport`

Responsibilities:

- Validate exactly one `>>`
- Split components on `.`
- RDKit parse each component
- Canonicalize SMILES
- Detect spectators (simple ions)
- Return:
  - `rxn_smiles_raw`
  - `rxn_smiles_clean`
  - `spectators[]`
  - `parse_warnings[]`
  - `standardization_actions[]`

### 5.2 `map_reaction(rxn_smiles_clean) -> MappingReport`

Responsibilities:

- Run rxnmapper
- Return `mapped_rxn_smiles`
- Compute simple QC: coverage, uniqueness of map numbers, etc.

### 5.3 `extract_bond_changes(mapped_rxn_smiles) -> BondChangeReport`

Responsibilities:

- Compare mapped reactants vs mapped products
- Identify:
  - formed bonds
  - broken bonds
  - bond order changes
- Derive `reaction_center_atoms`

### 5.4 `call_gpt_interpretation(payload) -> Interpretation`

Responsibilities:

- Call GPT-5.2 with `reasoning.effort` set based on complexity
- Use **Structured Outputs** (strict JSON schema) to enforce shape
- Constrain GPT:
  - “Use only tool_facts; do not invent reagents/conditions”
  - “If insufficient evidence, emit warnings + low confidence”

### 5.5 `finalize() -> ReactionAnalysis`

Responsibilities:

- Merge inputs + tool_facts + interpretation
- Apply gating:
  - if mapping_qc not ok → confidence cap + warnings

---

## 6) Suggested “Thinking Level” Strategy

Keep it simple:

- Default: `reasoning.effort = "low"` (fast + cheap)
- Escalate to `"high"` when:
  - >1 bond change event (likely tandem)
  - ambiguous class (many changes, rearrangements/redox)
  - mapping QC ok but interpretation uncertain

Do **not** use `"xhigh"` by default; reserve for rare hard cases.

---

## 7) Code Skeleton (pseudo-real)

### 7.1 Cleaning + spectator detection (RDKit)

```python
from dataclasses import dataclass
from typing import List, Optional, Dict
from rdkit import Chem

SPECTATORS = {"Cl", "Br", "I", "[Cl-]", "[Br-]", "[I-]", "[Na+]", "[K+]", "[Li+]"}

@dataclass
class CleanReport:
    rxn_smiles_raw: str
    rxn_smiles_clean: str
    spectators: List[str]
    parse_warnings: List[str]
    standardization_actions: List[str]

def _canon(smiles: str) -> Optional[str]:
    m = Chem.MolFromSmiles(smiles, sanitize=True)
    if m is None:
        return None
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)

def clean_reaction_smiles(rxn: str, drop_spectators: bool = True) -> CleanReport:
    warnings, actions, spectators = [], [], []

    if rxn.count(">>") != 1:
        raise ValueError("Invalid reaction SMILES: expected exactly one '>>'.")

    lhs, rhs = rxn.split(">>")
    lhs_parts = [p for p in lhs.split(".") if p]
    rhs_parts = [p for p in rhs.split(".") if p]

    def process(parts: List[str], side: str) -> List[str]:
        out = []
        for p in parts:
            c = _canon(p)
            if c is None:
                warnings.append(f"{side}: failed_to_parse: {p}")
                continue
            if drop_spectators and c in SPECTATORS:
                spectators.append(c)
                continue
            out.append(c)
        return out

    lhs_clean = process(lhs_parts, "LHS")
    rhs_clean = process(rhs_parts, "RHS")

    rxn_clean = ".".join(lhs_clean) + ">>" + ".".join(rhs_clean)
    return CleanReport(rxn, rxn_clean, sorted(set(spectators)), warnings, actions)
```

### 7.2 Mapping + bond diff (outline)

```python
def map_reaction(rxn_smiles_clean: str) -> dict:
    # Use rxnmapper here; return mapped_rxn_smiles + qc
    # qc fields: ok(bool), notes(list[str]), coverage(float)
    ...

def extract_bond_changes(mapped_rxn_smiles: str) -> dict:
    # Parse mapped reactants/products, compare bond sets by atom-map ids
    # Output: bond_changes[{id, change, a1, a2, bond}], reaction_center_atoms[]
    ...
```

### 7.3 GPT call (Structured Outputs)

```python
import json
from openai import OpenAI
client = OpenAI()

REACTION_SCHEMA = {
  "type": "object",
  "additionalProperties": False,
  "properties": {
    "overall_class": {"type": "string"},
    "tags": {"type": "array", "items": {"type": "string"}},
    "roles": {
      "type": "object",
      "additionalProperties": False,
      "properties": {
        "electrophile": {"type": "string"},
        "nucleophile": {"type": "string"},
        "leaving_group": {"type": "string"}
      },
      "required": ["electrophile", "nucleophile", "leaving_group"]
    },
    "events": {
      "type": "array",
      "items": {
        "type": "object",
        "additionalProperties": False,
        "properties": {
          "event_id": {"type": "string"},
          "event_type": {"type": "string"},
          "bond_change_refs": {"type": "array", "items": {"type": "string"}},
          "short_rationale": {"type": "string"},
          "confidence": {"type": "number"}
        },
        "required": ["event_id", "event_type", "bond_change_refs", "short_rationale", "confidence"]
      }
    },
    "mechanism_summary": {"type": "array", "items": {"type": "string"}},
    "warnings": {"type": "array", "items": {"type": "string"}},
    "confidence": {"type": "number"}
  },
  "required": ["overall_class", "tags", "roles", "events", "mechanism_summary", "warnings", "confidence"]
}

def call_gpt_interpretation(payload: dict, effort: str = "low") -> dict:
    resp = client.responses.create(
        model="gpt-5.2",
        reasoning={"effort": effort},
        input=[
            {"role": "system", "content": (
                "You are a reaction analysis engine. "
                "Use ONLY the provided tool_facts (bond changes, mapping_qc). "
                "If evidence is insufficient, add warnings and lower confidence. "
                "Do NOT invent reagents, catalysts, or conditions."
            )},
            {"role": "user", "content": json.dumps(payload, ensure_ascii=False)}
        ],
        text={
            "format": {
                "type": "json_schema",
                "name": "reaction_interpretation_v1",
                "strict": True,
                "schema": REACTION_SCHEMA
            }
        }
    )
    return json.loads(resp.output_text)
```

---

## 8) Simple Orchestrator

```python
def analyze_reaction(rxn_smiles: str) -> dict:
    clean = clean_reaction_smiles(rxn_smiles, drop_spectators=True)

    mapping = map_reaction(clean.rxn_smiles_clean)
    tool_facts = {"mapped_rxn_smiles": mapping.get("mapped_rxn_smiles"),
                 "mapping_qc": mapping.get("mapping_qc", {"ok": False, "notes": ["no_mapping"]})}

    if tool_facts["mapping_qc"]["ok"]:
        tool_facts.update(extract_bond_changes(tool_facts["mapped_rxn_smiles"]))
    else:
        tool_facts["bond_changes"] = []
        tool_facts["reaction_center_atoms"] = []

    payload = {
        "input": {
            "rxn_smiles_raw": clean.rxn_smiles_raw,
            "rxn_smiles_clean": clean.rxn_smiles_clean,
            "parse_warnings": clean.parse_warnings,
            "standardization_actions": clean.standardization_actions,
            "spectators": clean.spectators
        },
        "tool_facts": tool_facts
    }

    # Effort heuristic (simple):
    effort = "high" if len(tool_facts.get("bond_changes", [])) > 3 else "low"
    interp = call_gpt_interpretation(payload, effort=effort)

    # QC gating
    if not tool_facts["mapping_qc"]["ok"]:
        interp["warnings"].append("mapping_failed")
        interp["confidence"] = min(interp["confidence"], 0.3)

    return {
        "schema_version": "reaction_analysis.v1",
        "input": payload["input"],
        "tool_facts": tool_facts,
        "interpretation": interp
    }
```

---

## 9) Notes for Next Iteration (optional, still simple)

- Add a “fallback mode” that returns only `overall_class="other"` and warnings when:
  - parse fails
  - mapping fails
  - bond diff yields too many changes
- Add small spectator dictionary for common counterions (OTf, TsO, BF4, PF6)
- Add `stages` only if you regularly see tandem reactions (otherwise keep `events[]` flat)

---

## 10) Quick Checklist (operational)

- [ ] Always store raw + clean SMILES
- [ ] Log what you changed
- [ ] Drop spectators for mapping
- [ ] Anchor GPT to bond_changes IDs
- [ ] Enforce strict schema
- [ ] QC gate confidence (never overclaim)
