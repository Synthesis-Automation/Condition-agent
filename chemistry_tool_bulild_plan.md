#

Here’s a crisp, **step-by-step plan** to build the chemistry tools your recommender needs. It’s pragmatic, RDKit-first, and ships Ullmann C–N quickly while leaving clean hooks for other families.

## 0) Ground rules (once)

- **Language/stack:** Python 3.11 + RDKit + pydantic (schemas) + FastAPI (HTTP tools).
- **IDs:** `uid` for registry items; `reaction_id` for rows; `core` as normalized string (e.g., `"Cu/Phenanthroline"`).
- **Error codes:** `INVALID_SMILES`, `NO_REACTIVE_FAMILY`, `NO_PRECEDENTS`, `CONSTRAINTS_EMPTY_SET`.
- **Testing:** pytest + golden cases per tool; ≥90% coverage on SMARTS functions.
- **Perf budgets:** each tool < 300 ms p95 (except kNN retrieval < 800 ms). Add simple LRU caching.

---

## 1) Define tool contracts (schemas & stubs)

Create pydantic models for I/O and implement no-op stubs so the API surface is stable.

**Tools to implement:**

- `smiles.normalize`
- `router.detect_family`
- `featurize.substrates`
- `condition_core.parse` (from reagent list → normalized core)
- `properties.lookup` (bases/solvents/ligands)
- `precedent.knn`
- `constraints.filter`
- `explain.precedents`

**Deliverables:**

- `chemtools/contracts.py` (pydantic models)
- `chemtools/api.py` (FastAPI routes, stubs return `NOT_IMPLEMENTED`)
- Unit test skeletons

---

## 2) SMILES normalizer (robust parsing)

**Function:** `smiles.normalize(smiles: str) -> {smiles_norm, fragments, largest_smiles}`

**Steps:**

1. RDKit parse; if fails → `INVALID_SMILES`.
2. Split salts/solvents on `.`; keep **largest organic fragment**.
3. Neutralize common forms (carboxylate ↔ acid, ammonium ↔ amine).
4. Kekulize/aromatize; canonical SMILES out.

**Tests:**

- Quats, HCl salts, hydrates, multiple fragments, tautomers.
- Benchmark p95 < 10 ms.

---

## 3) Registry resolver

**Function:** `registry.resolve(query:str|cas|uid) -> {uid, role, name, aliases, props}`

**Steps:**

1. Build alias map `{CAS, token, common names} → uid`.
2. Strict role enum (`CATALYST|LIGAND|BASE|SOLVENT|ADDITIVE`).
3. Return minimal props (pK_a for bases; KT/DN/AN/Hansen for solvents; %V_bur/TEP for ligands if available).

**Tests:**

- Ambiguous names (e.g., “soda ash” vs Na2CO3) → tie-break by priority rules.
- Missing entry → `NOT_FOUND`.

---

## 4) Reaction-family router (deterministic first)

**Function:** `router.detect_family(reactants:[smiles]) -> {family:str, confidence:float, hits:[rule_ids]}`

**Tier A – SMARTS rules (fast):**

- Ullmann/Buchwald C–N: aryl/vinyl **halide** or **OTf** + N-nucleophile.
- Suzuki: aryl halide + boron partner (`B(OH)2`, `Bpin`).
- SNAr: electron-poor halo-aryl (≥1 strong EWG ortho/para).
- Amide coupling: carboxylic acid (or activated) + amine.
- Sonogashira: aryl/vinyl halide + terminal alkyne.

**Tier B – tie-breaker classifier (optional, v2):**

- ECFP(reactants set, radius 2–3, 1024) → logistic regression over families.

**Tests:**

- 10–20 curated examples per family + 20 adversarial “near misses”.
- Output confidence: 0.9 when single rule match; 0.6 when multiple.

---

## 5) Featurizer plugins (start with Ullmann C–N)

**Function:** `featurize.substrates(family, electrophile, nucleophile) -> dict`

**Ullmann v1 fields (SMARTS + ring neighborhood):**

- `LG ∈ {I, Br, Cl, OTf}`
- `elec_class ∈ {aryl, vinyl, alkyl}`
- `ortho_count ∈ {0,1,2+}` (around LG ipso)
- `para_EWG: bool` (NO2, CF3, CN, C=O)
- `heteroaryl: bool` (electrophile ring)
- `nuc_class ∈ {aniline, indole, amine_primary, amine_secondary, phenol, amide_deactivated}`
- `n_basicity ∈ {aromatic_primary, aliphatic_primary, secondary, deactivated}`
- `steric_alpha ∈ {low, med, high}` near N/O

**SMARTS starters (illustrative):**

- Aryl-X: `[$(c-[Cl,Br,I]),$(c[Cl,Br,I])]`
- OTf: `OS(=O)(=O)C(F)(F)F` (or anionic variant)
- Aniline N: `[$([NX3;H1][c]),$([NX3]([c])[H])]`
- Indole N–H: `[nH]1cccc2ccccc12`
- Phenol: `c[OH]`
- Nitro: `[N+](=O)[O-]` ; CF3: `C(F)(F)F` ; CN: `C#N`

**Tests:**

- Deterministic outputs on golden SMILES (e.g., Brc1ccc(F)cc1 + aniline).
- Timing p95 < 25 ms.

---

## 6) ConditionCore normalizer

**Function:** `condition_core.parse(reagents:[{uid, role}], text?:str) -> {core:str, metal_source_uid, ligand_uid?, ratio?:float, precatalyst?:bool}`

**Steps:**

1. From participants, pick `metal_source` (CuX, Pd(OAc)2…), `ligand` (phen, XPhos, DMEDA).
2. Normalize to `"{MetalSymbol}/{LigandCanonical}"` (e.g., `"Cu/Phenanthroline"`).
3. If a named precatalyst (e.g., “XPhos Pd G3”) is used, set `precatalyst=true`, also resolve to same canonical core string.

**Tests:**

- Multiple ligands present; missing ligand (Cu only).
- Ambiguous text hints → prefer participants list.

---

## 7) Properties service

**Function:** `properties.lookup(uid|token) -> props`

**Scope:**

- **Bases**: pK_a (DMSO if possible), hardness flags (carbonate vs hydroxide), nucleophilicity risk.
- **Solvents**: KT (Kamlet–Taft α/β/π\*), DN/AN, dielectric, bp, safety flags (chlorinated, peroxide-former).
- **Ligands**: %V_bur, TEP (if library available), denticity family (bidentate diamine, phenanthroline).

**Use:**

- Featurizer augment (e.g., polar solvent flag), ranker features, constraints filter (green, hazard).

---

## 8) Precedent finder (kNN)

**Function:** `precedent.knn(family, features, k=50, relax:dict) -> {prototype_id, support, precedents[]}`

**Steps:**

1. Build coarse **bin key** (e.g., `LG:Br|NUC:aniline|ortho:0|para_EWG:1`) to retrieve candidates fast.
2. Compute similarity in feature space (Hamming + small numeric distances).
3. Weight neighbors by yield and recency; output:
   - `prototype_id` (cache key for averaged condition vector)
   - `support` (n)
   - top 10 precedents `{reaction_id, yield, core, base_uid, solvent_uid, T, t}`

**Tests:**

- Returns empty with `NO_PRECEDENTS` if none; otherwise consistent ranking.

---

## 9) Constraints filter

**Function:** `constraints.filter(candidates, rules) -> {allowed, blocked[{id, reason}]}`

**Rules:**

- Inventory whitelist/blacklist (UIDs)
- Environmental/green flags (no HMPA, no chlorinated, bp > X, etc.)
- Operational: `max_T`, glovebox/no-glovebox, aqueous-only.

Run this both **before** and **after** tail retrieval.

---

## 10) Explain precedents

**Function:** `explain.precedents(pack, features) -> {reasons[], precedents[]}`

**Compose:**

- 2–3 closest matches (substrates, conditions, yields, dataset_id/DOI).
- Short reasons: “Br–aniline; ortho=0; para-EWG present; core & KOH/H2O frequent high yield.”

---

## 11) Packaging as services

- Wire all tools behind FastAPI with clear routes.
- Add request/response logging (PII-free), latency histograms, error rate dashboards.
- Add simple **LRU cache** on `featurize` and `knn`.

---

## 12) QA & golden tests

- Build a **golden set** (≈100 reactions) across families with expected outputs:
  - Family classification
  - Featurizer dict
  - Core normalization
- Unit tests for each SMARTS; fuzz tests with random substitutions.
- Regression tests: lock JSON outputs with approvals (e.g., `pytest-approvaltests`).

---

## 13) Performance & robustness

- Benchmarks on 10k reactions: total ETL + featurize throughput ≥ 500 rps on a single node.
- Guardrails:
  - If **router** low-confidence → return top-2 families.
  - If **knn support** < threshold → tell caller to switch to **plate design**.
  - Always return **machine-parsable error codes**.

---

## 14) Documentation & examples

- Markdown docs with:
  - API examples (curl + Python snippets)
  - SMARTS reference
  - Error codes and recovery
  - “How to add a new family plugin” checklist

---

## 15) Roadmap to more families (after Ullmann)

**Order of addition (each is a ~1-page featurizer):**

1. **Buchwald C–N**: same LG features; ligand families; base (K3PO4, Cs2CO3); solvent (toluene, dioxane).
2. **SNAr**: EWG density score; phenoxide/amine class; polar aprotic solvents.
3. **Suzuki**: boron partner type; base family (carbonate/phosphate); water content; Pd/XPhos/Pd/SPhos cores.
4. **Amide coupling**: coupling agent vs in situ activation; base/solvent windows.

---

## Milestones (suggested)

- **Week 1:** `smiles.normalize`, `registry.resolve`, router v1, Ullmann featurizer v1, unit tests.
- **Week 2:** core normalizer, properties service, precedent.kNN (bin + features), constraints.filter.
- **Week 3:** explain.precedents, API wiring & caching, golden tests, perf pass.
- **Week 4:** expand tests, add SNAr featurizer, docs; handoff to recommender pipeline.
 
