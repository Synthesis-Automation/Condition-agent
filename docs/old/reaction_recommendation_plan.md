# Condition Recommendation — Solid Plan (Project Summary)

This plan distills everything we discussed into a concrete, buildable roadmap. It focuses on **Ullmann C–N** first (to prove reliability), then generalizes.

---

## 0) Goals & Success Criteria

**Primary goal:** Ship a **reliable, explainable, constraints-aware** condition recommender for Ullmann C–N and other reaction (later) that returns:  
`ConditionCore + Base + Solvent (+ Additive) + T/t`, with top-N alternates, precedents, confidence, and plate fallback.

**Success by V1 (Ullmann-only):**
- ≥80% **hit@3** (one of top-3 achieves ≥X% yield on holdout)  
- ≥0.7 **calibration** (80% confidence → ~80% empirical success)  
- ≤2s median latency for `recommend()`  
- 100% picks map to **registry UIDs**, pass inventory/green constraints

---

## 1) Architecture (who does what)

- **Local App = Orchestrator (recommended)**  
  Deterministic planning, tool execution, constraints, scoring, logging.
- **LLM = Interface + Explainer**  
  Parses user intent, calls coarse APIs (`recommend`, `design_plate`, `explain`), writes human-readable rationales.

**Why:** deterministic, auditable core + flexible UX; no hallucinated reagents; secure & reproducible.

---

## 2) Data Foundations

### 2.1 Compound Registry (single source of truth)
- **Fields:**  
  `uid`, `aliases{cas, token, name}`, `role ∈ {CATALYST, LIGAND, BASE, SOLVENT, ADDITIVE}`, optional `smiles`, `props`  
  (bases: pKa; solvents: KT/DN/AN/Hansen; ligands: %V_bur/TEP)
- **Targets:**  
  *≥95%* of referenced chemicals in reactions resolve to a `uid`.

### 2.2 Reaction Datasets (Ullmann first)
- **Fields:**  
  `reaction_id, rxn_type, condition_core, participants[{uid,role}], yield_value, yield_type, T_C, time_h, substrates{electrophile_smiles,nucleophile_smiles}, dataset_id, source(doi/lab), year`
- **ConditionCore decomposition:** `{metal_source, ligand, L:M, precatalyst?}`  
- **Quality:** normalize units, dedupe, keep negatives/low yields.

**Definition of Done (data):**  
- 100% rows with `rxn_type`; ≥90% with `yield_value`; ≥80% with `T_C` & `time_h`.  
- Zero free-text roles/units; **UID** resolution ≥95%.

---

## 3) Chemistry Tools (deterministic)

Expose these as local APIs the LLM can call:

1) `router.detect_family(reactants)` → `"Ullmann_CN" | "Buchwald_CN" | "SNAr" | ...", confidence`
2) `smiles.normalize(smiles)` → largest fragment, neutralize, aromatize  
3) `featurize.substrates(family, electrophile, nucleophile)` →  
   `{"LG":"I/Br/Cl/OTf","nuc_class":"aniline/indole/1°/2°","ortho_count":0–2,"para_EWG":true/false,"heteroaryl":bool,...}`
4) `precedent.knn(family, features, k=50)` → `{prototype_id, support, precedents[]}` (yield-weighted mean vector)
5) `cores.retrieve(prototype_id, constraints, top=20)` → list of cores (cosine, mean_yield, n, shrinkage)
6) `tails.retrieve(prototype_id, core, role, top=50)` → bases/solvents/additives shortlist
7) `ranker.score_packs(family, features, candidates, constraints)` → final scores & expected yields
8) `predict.tt(family, features, pack)` → `[T25,T50,T75]`, `[t25,t50,t75]`
9) `constraints.filter(candidates, rules)` → remove non-stock/hazard/green violations
10) `explain.precedents(pack, features)` → 2–3 closest with yields & citations
11) `log.write_event(query, picks, outcome?)` → telemetry for learning loop

> **SMARTS starter (Ullmann featurizer):**  
> Aryl–X: `c[Cl,Br,I]` ; Triflate: `OS(=O)(=O)C(F)(F)F` ;  
> Aniline N: `[$([NX3;H1][c]),$([NX3]([c])[H])]` ; Indole N–H: `[nH]1cccc2ccccc12` ; Phenol: `c[OH]`.

---

## 4) Modeling & Retrieval

### 4.1 Embeddings (offline, refreshed weekly)
- **Item vectors (per UID):** role-weighted co-usage + rxn_type + yield buckets → **PPMI → SVD(64)**  
  *(Optionally concat ECFP/physchem and props; z-score by block.)*
- **ConditionCore vectors:** core × (co-used tails + rxn_type + yield) → **PPMI → SVD(64)**  
- Build **per-role ANN** (base/solvent/additive) and **core ANN** (FAISS/HNSW).

### 4.2 Online recommendation (Ullmann v1)
1) **Normalize SMILES → Featurize** (LG, nuc_class, ortho/EWG, …)  
2) **Prototype** = yield-weighted mean of **kNN precedents** in same bin (or relaxed)  
3) **Core shortlist**: top-M by cosine(core_vec, prototype), with **shrinkage**  
   `μ_shrunk = (μ·n)/(n+5)`; `score_core = 0.6·cos + 0.3·μ_shrunk + 0.1·prior/constraints`
4) **Tail selection per core**: ANN retrieve per role → **MMR diversity** → **GBDT re-rank** using:  
   `cos_to_proto`, bin features, role rules (e.g., non-nucleophilic base), constraints (stock/green/Tmax)
5) **T/t prediction**: separate quantile regressors  
6) **Confidence & gating**: if support < threshold or confidence low → **24/96-well screen** (diverse cores/tails)

**Always return precedents + reasons.**

---

## 5) Evaluation (honest & actionable)

- **Splits:** by substrate scaffold **and** by publication (to prevent leakage).
- **Metrics:**  
  - *Topline*: **hit@k** (≥X% yield), NDCG@k, median top-1 yield  
  - *Calibration*: coverage vs stated confidence  
  - *Ops*: constraint satisfaction rate, latency
- **Slices:** LG (Cl/Br/I), nucleophile class (aniline/indole/1°/2°), data-rich vs data-poor bins.
- **Ablations:** no-precedent / no-structure / no-shrinkage to prove each component’s value.

---

## 6) Deliverables & APIs

- `POST /recommend_ullmann`  
  **In:** `{electrophile, nucleophile, constraints{inventory, green, Tmax, glovebox}}`  
  **Out:** `{bin, recommendations[{core, base_uid, solvent_uid, additive_uids?, T_C, time_h, confidence, precedents[]}], alternatives[], if_low_confidence{plate_csv}}`

- `POST /design_plate` → diversified 24/96-well screen CSV (UIDs, mol%, volumes)

- `GET /explain` → reasons + precedent citations for a chosen pack

- **Artifacts versioning:** `chemvec_v1`, `corevec_v1`, `ranker_v1`, `tt_v1`

---

## 7) Roadmap & Milestones

**Week 1–2 — Data & Tools**
- Finish **ETL cleanup** (roles, cores, units, SMILES)  
- Implement **router** & **Ullmann featurizer** (SMARTS)  
- Build **embeddings** (items + cores) & ANN indices

**Week 3 — Ullmann MVP**
- Implement prototype from **kNN precedents**  
- Implement **core shortlist** + **tail retrieve** + **GBDT re-rank**  
- Add **predict.tt**; **constraints.filter**; **precedents**  
- Ship `recommend_ullmann` API + MD/CSV output; basic UI

**Week 4 — Reliability & Eval**
- Build evaluation harness (scaffold + publication splits)  
- Calibrate confidence; add **gating → plate** fallback  
- Add **logging** & weekly refresh pipeline

**Week 5–6 — Hardening & Expansion Prep**
- Add **inventory/green** connectors; ops dashboards  
- Prepare **Buchwald** featurizer & embeddings (no release yet)  
- Optional: **core-graph** neighbors (cross-tech substitutes)

---

## 8) Risks & Mitigations

- **Messy identities/roles →** enforce UID resolution & role enums in ETL; coverage alerts.  
- **Sparse bins/overfitting →** shrinkage scoring; confidence gating; screen fallback.  
- **Publication/popularity bias →** inverse-propensity debias in re-ranker; keep negatives.  
- **Latency →** cache `prototype_id`; precompute ANN; lazy precedents fetch.  
- **Explainability →** require 2–3 precedents for every returned pack.

---

## 9) Extensions (post-V1)

- **Two-tower model** (substrate ↔ core) to bridge Pd↔Cu and OOD substrates.  
- **Core substitution graph** from substrate-binned profile correlations.  
- **More families**: Buchwald, SNAr, Suzuki, Amide coupling (add featurizer plugins).  
- **Ops-aware ranking**: cost, hazard, glovebox, greener-solvent scoring.  
- **Active learning**: uncertainty-driven plate designs; prioritize data collection.

---

## 10) Summary

- Start **Ullmann-only**, nail **determinism, constraints, and explanations**.  
- Use **registry-first data**, **SMARTS features**, **precedent-based prototypes**, **retrieve → re-rank**, and **confidence gating**.  
- Keep the **LLM as the face**, **local app as the brain**.  
- Log everything; iterate weekly.  
- Once reliable, expand to new families by adding featurizer plugins and reusing the same pipeline.




