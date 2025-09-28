# Condition_Rule_Library (CRL) — Implementation Guide

A concise, **implementation‑ready** plan to build an evolving rule/recipe library your Manager/Dev/Critic loop can read from and write back to.

---

## 1) What CRL Is

A versioned, JSON‑backed library that stores, per **reaction family**:

- **Playbooks** – reusable condition recipes keyed by substrate context
- **Guards** – hard constraints / safety & compatibility checks
- **Priors** – soft boosts/penalties that bias ranking
- **Repairs** – structured fallback steps for retries/diagnostics
- **Metrics** – `n_runs`, `median_yield`, success rate, etc., to *promote* or *deprecate* entries

> CRL replaces “Template Library”. Keep the public API thin and stable so your agents remain decoupled.

---

## 2) Repository Layout

```
condition_agent/
  config.py
  condition_rule_library/
    schema.json           # CRL JSON Schema
    crl.json              # The library (by family)
    api.py                # Loader + query/eval APIs
    evolver.py            # Trace -> CRL updater (promotion/deprecation)
    bandit.py             # Thompson Sampling over playbooks (explore/exploit)
tests/
  test_crl_api.py
  test_evolver.py
```

Update `config.py`:

```python
class RunConfig:
    crl_path = "condition_agent/condition_rule_library/crl.json"
    crl_schema = "condition_agent/condition_rule_library/schema.json"
    runs_dir = "runs"
```

---

## 3) Minimal JSON Schema

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title": "Condition_Rule_Library",
  "type": "object",
  "required": ["families"],
  "properties": {
    "version": {"type": "string"},
    "schema_version": {"type": "string"},
    "defaults": {"type": "object"},
    "families": {
      "type": "object",
      "patternProperties": {
        "^[A-Za-z0-9._-]+$": {
          "type": "object",
          "properties": {
            "playbooks": {"type": "array"},
            "guards": {"type": "array"},
            "priors": {"type": "array"},
            "repairs": {"type": "array"}
          }
        }
      }
    }
  }
}
```

---

## 4) Example `crl.json` (Buchwald–C–N excerpt)

```json
{
  "version": "2025-09-28",
  "schema_version": "1.1",
  "defaults": {
    "concentration_M_range": [0.05, 0.25],
    "temperature_C_range": [80, 120]
  },
  "families": {
    "Buchwald-CN": {
      "playbooks": [
        {
          "id": "buchwald_cn_steric_aniline_arcl_v1",
          "name": "steric_aniline_ArCl",
          "status": "candidate",
          "when": {"substrate_class": "steric_aniline_ArCl"},
          "recipe": {
            "catalyst": "Pd/tBuBrettPhos",
            "ligand": "tBuBrettPhos",
            "base": "K3PO4",
            "solvent": "1,4-dioxane",
            "temperature_C": 100,
            "time_h": 12,
            "loadings": {"Pd_mol%": 1.5, "ligand_mol%": 3.0, "base_eq": 2.0}
          },
          "fallbacks": [{"ligand": "XPhos"}, {"base": "Cs2CO3"}],
          "metrics": {"n_runs": 0, "median_yield": null}
        }
      ],
      "guards": [
        {
          "id": "GUARD-ORTHO-STERICS-BULKY-LIGAND",
          "if": {"electrophile": {"class": "aryl chloride", "ortho_sub_count": {">=": 1}},
                 "nucleophile": {"class": "aniline"}},
          "action": "enforce",
          "params": {"require_ligand_in": ["tBuBrettPhos","XPhos","RuPhos"]}
        }
      ],
      "priors": [
        {
          "feature": "electrophile=aryl chloride & ortho_sub_count>=1",
          "effect": "boost",
          "targets": ["Pd/tBuBrettPhos + K3PO4 + 1,4-dioxane"]
        }
      ],
      "repairs": [
        {
          "id": "BH-REPAIR-LOW-CONVERSION",
          "if": {"outcome": "low_conversion"},
          "then": {"steps": ["+10C","base hop K2CO3→K3PO4→Cs2CO3","ligand hop SPhos↔XPhos↔tBuBrettPhos","Pd→2 mol%","solvent hop toluene↔dioxane↔DMA"]}
        }
      ]
    }
  }
}
```

---

## 5) Python API (thin & stable)

`condition_rule_library/api.py` (stubs):

```python
from typing import Dict, Any, List, Tuple
import json, pathlib
from jsonschema import validate

def load_crl(path: str, schema_path: str) -> Dict[str, Any]:
    data = json.loads(pathlib.Path(path).read_text(encoding="utf-8"))
    schema = json.loads(pathlib.Path(schema_path).read_text(encoding="utf-8"))
    validate(data, schema)
    return data

def _match_when(when: Dict[str, Any], features: Dict[str, Any]) -> bool:
    sc = when.get("substrate_class")
    if sc and sc == features.get("substrate_class"):
        return True
    pred = when.get("feature_predicate", {})
    for k, v in pred.items():
        if isinstance(v, dict) and ">=" in v:
            if not (features.get(k, 0) >= v[">="]): return False
        elif features.get(k) != v:
            return False
    return True

def select_playbooks(crl: Dict[str, Any], family: str, features: Dict[str, Any]) -> List[Dict[str, Any]]:
    fam = crl.get("families", {}).get(family, {})
    pbs = fam.get("playbooks", [])
    hits = [p for p in pbs if _match_when(p.get("when", {}), features)]
    def key(p):
        status_rank = {"stable": 0, "candidate": 1, "deprecated": 2}.get(p.get("status","candidate"), 1)
        n_runs = (p.get("metrics", {}) or {}).get("n_runs", 0)
        return (status_rank, -n_runs)
    return sorted(hits or pbs, key=key)

def apply_guards(crl: Dict[str, Any], family: str, candidate: Dict[str, Any], job_ctx: Dict[str, Any]) -> Dict[str, Any]:
    fam = crl.get("families", {}).get(family, {})
    msgs, patches, ok = [], {}, True
    for g in fam.get("guards", []):
        # TODO: implement real evaluator using g["if"]
        passed = True
        if not passed:
            act = g.get("action", "warn")
            if act == "enforce": ok = False
            msgs.append({"id": g.get("id"), "action": act})
    return {"ok": ok, "messages": msgs, "patch": patches}

def score_with_priors(crl: Dict[str, Any], family: str, candidate: Dict[str, Any], features: Dict[str, Any]) -> float:
    fam = crl.get("families", {}).get(family, {})
    bump = 0.0
    for pr in fam.get("priors", []):
        # TODO: map pr["feature"] to feature checks and add bumps
        pass
    return bump
```

---

## 6) Wire CRL into Manager / Dev / Critic

**Manager**
- `crl = load_crl(cfg.crl_path, cfg.crl_schema)`
- `playbooks = select_playbooks(crl, family, features)`
- Pass top‑N playbooks to Dev to seed candidates

**Dev**
- Expand playbook `recipe` → concrete candidates (+ cheap analogs)
- Fill missing params from `defaults`

**Critic**
- `apply_guards(...)` to veto or propose patches
- `score_with_priors(...)` → add to heuristic score before ranking
- Emit patches (base/ligand/solvent/temp) that Dev understands

---

## 7) Evolution Pipeline (Self‑Improving CRL)

1. **Trace** per run → `runs/<id>/trace.jsonl`
2. **Append** HTE result: `{"event":"lab_outcome","observed_yield":0.72}`
3. **Evolver** buckets by `(family, substrate_class)` → modal ligand/base/solvent; robust medians for temperature/time; common loadings
4. **Update metrics**: `n_runs`, `median_yield`, `p75_yield`, success rate
5. **Promote/Deprecate**:
   - Promote to `stable` if `n_runs ≥ 4` and `median_yield ≥ 0.55` (tune thresholds)
   - Deprecate supported underperformers
6. **Write back** to `crl.json`, append `CRL_CHANGELOG.md`

CLI:
```bash
python -m condition_agent.condition_rule_library.evolver --runs_dir runs \
  --crl condition_agent/condition_rule_library/crl.json --min_support 4 --promote_yield 0.55
```

---

## 8) Explore/Exploit with Bandit

- Arm key: `family::substrate_class::playbook_id`
- **Thompson Sampling** with success = `observed_yield ≥ target_yield`
- Update after each HTE result to balance exploration with exploitation

---

## 9) Validation, Tests, Governance

- JSON Schema validation in CI (pre‑commit)
- Unit tests for:
  - playbook selection & ordering
  - guard evaluation & patch shaping
  - evolver thresholds (promotion/deprecation)
- `CRL_CHANGELOG.md` updated on every evolution
- `status: candidate` tried only under explore mode

---

## 10) Migration from Template Library

- Move `templates/playbooks.json` → `condition_rule_library/crl.json` under `families.<family>.playbooks`
- Rename imports to CRL API
- Convert existing rules to **guards/priors/repairs**

---

## 11) Minimal “Done” Checklist

- [ ] `crl.json` + `schema.json` in repo  
- [ ] `api.py` wired to Manager/Dev/Critic  
- [ ] Evolver harvesting traces → CRL updates  
- [ ] Bandit for explore/exploit  
- [ ] Tests passing (selection, guards, evolver)  
- [ ] Changelog appended on updates
