# Lightweight Reaction Similarity & Condition Recommendation — Step‑by‑Step Build Plan

This plan delivers a **fully local, CPU‑friendly** pipeline for finding similar reactions from reaction SMILES and recommending transferable conditions using **DRFP** + **cheap chemistry filters** + **kNN voting**.

---

## 0) Goals & Success Criteria
**Goal:** Given a query reaction SMILES, retrieve the most similar reactions and recommend **ConditionCore** and **Tail** (base, solvent, T, time) suitable for transfer.  
**Win conditions:**
- High **precision@N** of retrieved neighbors within the same reaction family.
- Top‑1 **ConditionCore** match rate competitive with expert heuristics.
- Latency \< 300 ms per query on ≤100k reactions (CPU, no GPU).

---

## 1) Minimal Tech Stack
- **Python ≥3.9**
- Libraries: `drfp`, `rdkit-pypi`, `pandas`, `numpy`, `pyyaml`  
  *(Optional large-scale)*: `nmslib` for ANN (HNSW on Jaccard)  
  *(Optional extras)*: `faiss-cpu`, `rxnfp`, `torch`

```bash
pip install drfp rdkit-pypi pandas numpy pyyaml
# optional for >100k reactions
pip install nmslib
```

---

## 2) Repository Layout
```
react-sim-lite/
├─ data/
│  ├─ reactions.jsonl         # your dataset
│  ├─ lookup/                 # enumerations (LG types, roles…)
├─ features/
│  ├─ rules/                  # SMARTS and role rules
├─ src/
│  ├─ normalize.py            # SMILES canonicalization & ordering
│  ├─ features.py             # leaving group, nucleophile class, ortho-sub count
│  ├─ fps_drfp.py             # DRFP computation utilities
│  ├─ index_exact.py          # exact Tanimoto scan (≤100k)
│  ├─ index_hnsw.py           # NMSLIB Jaccard HNSW (optional)
│  ├─ retrieve.py             # end-to-end retrieval + filters
│  ├─ recommend.py            # kNN vote for ConditionCore/Tail
│  ├─ evaluate.py             # metrics (Recall@K, MAP@K, Core-Top1)
│  └─ io_schema.py            # read/write JSONL, schema checks
├─ configs/
│  └─ default.yaml            # knobs for radius, nBits, K, filters
├─ cli/
│  ├─ build_index.py
│  ├─ query.py
│  └─ eval.py
└─ README.md
```

---

## 3) Data Schema (JSONL)
Each line = one reaction record.
```json
{
  "rxn_id": "BHA-001234",
  "reaction_smiles": "reactants.reagents>>products",
  "class": "Buchwald",                     // normalized class label
  "condition_core": "Pd/XPhos",            // ConditionCore
  "tail": { "base": "K3PO4", "solvent": "toluene" },
  "temperature_c": 110,
  "time_h": 12.0,
  "yield_pct": 78.0,
  "metadata": { "lg": "ArCl", "nuc_class": "primary_aliphatic", "ortho_subs": 0 }
}
```
*If your dataset lacks `metadata`, you will compute it in Step 5.*

---

## 4) Normalization (src/normalize.py)
- Canonicalize molecule SMILES (RDKit `MolStandardize`).
- Deterministic ordering of sides: `reactants . reagents >> products`.
- Drop obvious spectators (salts, trivial solvents) unless they are mechanistically essential.
- Enforce atom counts reasonability checks.

**Deliverable:** a cleaned `reactions.jsonl` and a **hash** of normalized `reaction_smiles` for caching.

---

## 5) Fast Chemistry Features (src/features.py)
Implement lightweight, rule‑based extractors:
- **Leaving group (LG)**: Detect Ar–Cl/Br/I/OTf/OMs/… via SMARTS on reactant side.
- **Nucleophile class**: primary/secondary; aliphatic vs aniline; heteroaryl amines.
- **Ortho‑substitution count**: count ortho substituents vs the electrophilic carbon.
- **Heteroaryl flag**, **ring strain flag** (optional).

**Output:** augment each record with `metadata` fields used for **filters**.

---

## 6) Fingerprints (src/fps_drfp.py)
- Use **DRFP** with `n_bits=4096`, `r=3` (ngrams 0–3).
- Store both dense `np.uint8` and **sparse bit indices** (for NMSLIB).

**Artifacts:**  
- `artifacts/drfp.bits.npy` (packed bits for exact scoring)  
- `artifacts/drfp.sparse.pkl` (list of set bit indices for ANN)  
- `artifacts/rxn_ids.npy`

Config in `configs/default.yaml`:
```yaml
fingerprint:
  type: drfp
  n_bits: 4096
  radius: 3
retrieval:
  K: 200
  use_ann: false         # true for >100k
filters:
  by_class: true
  by_lg: true
  by_nuc_class: true
  max_ortho_subs: 2
recommendation:
  temp_agg: median
  time_agg: median
```

---

## 7) Indexing
### Option A — Exact (≤100k)
- Scan all vectors and compute **Tanimoto** on 0/1 arrays.
- Return top‑K.

### Option B — ANN (≥100k)
- Build **NMSLIB HNSW** with `space='jaccard_sparse'` from sparse bit indices.
- Query for top‑K, then **re-score** these with exact Tanimoto for final ranking.

**Deliverables:**  
- `artifacts/index_hnsw.bin` (optional)  
- `cli/build_index.py` to (re)build on demand.

---

## 8) Retrieval API (src/retrieve.py)
**Input:** query reaction SMILES.  
**Steps:**
1. Normalize and compute DRFP.
2. Hard **filters**: `class`, `lg`, `nuc_class`, `ortho_subs` (tunable).
3. Retrieve nearest neighbors (exact or ANN+re-score).
4. Return top‑K with scores + metadata required for recommendation.

**CLI:**
```bash
python -m cli.query \
  --config configs/default.yaml \
  --query "Clc1ccc(Cl)cc1.NCC>>Clc1ccc(NCC)cc1" \
  --topk 50 --filters class lg nuc_class
```

---

## 9) Recommendation Logic (src/recommend.py)
On the retrieved **top‑K** (after filters):
- **ConditionCore**: majority vote with **Laplace smoothing**.
- **Base**/**Solvent** (**Tail**): conditional frequency given Core and filters.
- **Temperature/Time**: median (or Tukey‑trimmed mean) of the subset that shares the chosen Core.
- **Guardrails**: cap temperature/time to empirical percentiles; preserve solvent/base compatibility (e.g., KOtBu not in protic solvent).

**Output:**
```json
{
  "recommended": {
    "condition_core": "Pd/XPhos",
    "base": "K3PO4",
    "solvent": "toluene",
    "temperature_c": 110,
    "time_h": 12
  },
  "neighbors": [ { "rxn_id": "...", "sim": 0.84 }, ... ]
}
```

---

## 10) Evaluation (src/evaluate.py)
- **Retrieval**: Recall@K, MAP@K within the same reaction class and matching LG/nucleophile class.
- **Recommendation**: Top‑1 **Core accuracy**, Tail **top‑1/top‑3** match rate; **MAE** for T/time vs held‑out ground truth.
- **Splits**: stratify by reaction class and leaving group; keep leakage out (don’t share identical product scaffolds across train/test).

**CLI:**
```bash
python -m cli.eval --config configs/default.yaml --k 50 200
```

---

## 11) Interfaces & I/O
- **JSONL in/out** (stable schema).
- Optional **REST** wrapper (FastAPI) or **CLI‑only** for desktop/offline.
- Provide a **deterministic seed** and exact version stamps for reproducibility.

---

## 12) Milestones (sequential)
1. **Data normalization & schema checks** → verified `reactions.jsonl`.
2. **Feature extraction** (LG, nuc class, ortho count) → metadata augmented.
3. **DRFP computation** and **exact retrieval** path working.
4. **Filters wired** → qualitative inspection of top‑K.
5. **Recommendation** (kNN voting) → first end‑to‑end result.
6. **Evaluation harness** → numbers for Retrieval/Recommendation.
7. **ANN option** (NMSLIB) for large datasets.
8. **Packaging**: configs, CLI, README with examples.
9. **(Optional)**: caching, simple GUI/Gradio, batch mode.

---

## 13) Minimal Code Sketches

**DRFP encode**
```python
from drfp import DrfpEncoder

def drfp_bits(rsmi: str, n_bits=4096, r=3):
    return DrfpEncoder.encode(rsmi, ngram_range=(0, r), n_bits=n_bits).astype('uint8')
```

**Exact Tanimoto**
```python
import numpy as np

def tanimoto(a: np.ndarray, b: np.ndarray):
    inter = np.sum((a & b) == 1)
    denom = np.sum(a==1) + np.sum(b==1) - inter
    return 0.0 if denom == 0 else inter / denom
```

**kNN voting (smoothed)**
```python
from collections import Counter, defaultdict

def vote_core(neighbors, alpha=1.0):
    c = Counter(n['condition_core'] for n in neighbors)
    labels = list(c.keys())
    tot = sum(c.values())
    scores = {L: (c[L] + alpha) / (tot + alpha*len(labels)) for L in labels}
    return max(scores, key=scores.get), scores
```

---

## 14) Extensions (when needed)
- **Reaction‑center FP** (RDKit difference FP) to re‑rank top‑K.
- **RXNFP** re‑ranker (Transformer) if you later want extra precision.
- **Learning‑to‑Rank** on top‑K with simple features (class, LG, nuc class, ortho count, DRFP sim).

---

## 15) Ready‑to‑Use Defaults
- DRFP: `n_bits=4096`, `radius=3`
- Retrieval: `K=200` (re-score top‑K if using ANN)
- Filters: `class`, `lg`, `nuc_class`, `max_ortho_subs=2`
- Aggregation: `median` for temperature & time
- ANN (if enabled): HNSW `M=32`, `efConstruction=200`, query `ef=200`

---

**That’s the full build plan.** Start with Steps 1–5 to get an end‑to‑end pipeline, then harden with Step 10 (evaluation) and add ANN if your dataset size demands it.
