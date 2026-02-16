# zScore-App–style Condition Recommender (v1, simple)

This note describes a **minimal, robust implementation** of a “best-seller” condition recommender using the **structure-free** Roche HTE dataset (e.g., `z-Score Peaks with FG.csv`) and the ranking idea from the associated paper/app. citeturn0search0turn0search1turn0search2

The recommender answers:

- **Given** a *reaction class* + *reacting functional-group pair* (e.g., `Suzuki-Miyaura` + `ArCl` × `ArB(OH)2`)
- **Return** a single “most reliable” **condition core** (Catalyst/Ligand/Base/Solvent) **with evidence**
- (Optional) build a **diversified HTE plate** (e.g., 24 wells) for hard/undocumented cases

---

## 0) Key constraints (dataset reality)

- No reactant/product structures; only coarse labels like `ArCl`, `RNH2`, etc.
- Many reaction types have **very few FG-pair buckets** (sometimes 1), so recommendations are **priors**, not structure-specific predictions.
- Use **support thresholds** + **back-off** to avoid overconfidence.

---

## 1) Inputs / Outputs

### Input A: dataset CSV

Expected columns (exact names from the repo’s CSV):

- `ELN_ID` (transformation ID)
- `Reaction Type`
- `FG A`, `FG B`, `FG_sorted`
- `z-Score` (string or number; may use comma decimal in some exports)
- condition columns: `Catalyst`, `Ligand`, `Base`, `Solvent`, `Additive`, `Coupling Reagent`, `Secondary Solvent`, `Tertiary Solvent`

### Input B: query

Minimum query key (order-invariant):

- `reaction_type`: str (e.g., `Suzuki-Miyaura`)
- `fg1`, `fg2`: str tokens from the dataset vocabulary (e.g., `ArCl`, `ArB(OH)2`)

### Output (Top-1 recommendation)

A JSON-like dict:

- `query`: {reaction_type, fg1, fg2, backoff_level}
- `condition_core`: {Catalyst, Ligand, Base, Solvent}
- `score`: float (best-seller score)
- `support`: {n_eln, n_rows}
- `evidence`: list of top precedent `ELN_ID`s with their best z-scores
- `notes`: warnings if slice is coarse / support low

---

## 2) Core definitions (keep simple)

### 2.1 Canonical FG pair

Always use **order-invariant** FG pair:

- If `FG_sorted` exists: parse it into (`FG1`, `FG2`)
- Else: sort (`FG A`, `FG B`) lexicographically to build `FG_sorted`

> Treat `FG A` vs `FG B` as dataset ordering, not electrophile/nucleophile roles.

### 2.2 Condition key (recommendation granularity)

Use a **condition core** to ensure repetition:

- `ConditionCore = (Catalyst, Ligand, Base, Solvent)`

Avoid full bundles (additive/coupling reagent/secondary solvent) at v1, because they fragment support.

**Normalization**:

- Convert NaN/empty strings to `None`
- Strip whitespace
- Optionally collapse “No Ligand” / blank → `None` (configurable)

---

## 3) Preprocessing (must-do)

### 3.1 Load and sanitize numeric z-scores

- Read with `pandas.read_csv(..., low_memory=False)`
- Convert `z-Score` to float:
  - replace comma decimal `,` with `.` if needed
  - coerce errors to NaN and drop

### 3.2 Deduplicate within a transformation (optimum-seeking)

Within each `(ELN_ID, ConditionCore)` keep the **maximum** `z-Score`.

This mimics the app/paper practice of focusing on **optima** rather than penalizing a condition for multiple mediocre wells.

---

## 4) “Best-seller” ranking (single recommendation)

### 4.1 Candidate slice (exact bucket)

Filter rows where:

- `Reaction Type == reaction_type`
- `FG_sorted == sort_pair(fg1, fg2)` (or `FG1==fg1 & FG2==fg2` orderless)

Compute:

- `n_eln_bucket = nunique(ELN_ID)` in this slice

### 4.2 Back-off rules (avoid brittle slices)

If `n_eln_bucket < MIN_ELN` (default 5), back off:

1. **Exact**: `(Reaction Type, FG_sorted)`
2. **Contains**: `(Reaction Type, FG contains fg1 OR fg2)`  
   (i.e., any partner as long as one FG matches)
3. **ReactionType only**: `(Reaction Type)`
4. **Global**: (all rows) — last resort

Return `backoff_level` in output.

---

## 5) Scoring function (simple, robust)

For each `ConditionCore` in the (possibly backed-off) slice:

1. Build `per_eln_best = max(z)` for each `ELN_ID` (after dedup)
2. Let `Z = sorted(per_eln_best, reverse=True)`
3. Support: `n_eln = len(Z)`
4. If `n_eln < MIN_ELN_PER_CONDITION` → skip
5. Score:
   - `top = Z[:TOP_N]` (default `TOP_N = 5`)
   - `score = median(top)`

**Tie-breakers** (in order):

1. higher `n_eln`
2. higher `score`
3. fewer missing fields in (Catalyst/Ligand/Base/Solvent)

This is the “Amazon best seller” analogue:

- reward high performance **and** repeated evidence
- median(top-N) resists outlier wells

---

## 6) Returning one recommendation (Top-1)

Select the condition with the maximum `(score, n_eln, completeness)`.

Include evidence:

- top 3–10 `ELN_ID`s with their best z-scores for that condition
- the distribution summary: min/median/max of `per_eln_best`

Also include warnings:

- `if backoff_level > 1`: “exact bucket sparse; backed off”
- `if reaction_type has <=2 FG pairs overall`: “coarse class; substrate effects not resolved”

---

These defaults match the spirit of the paper/app: prioritize **robust, repeated top-end performance** over one-off winners. citeturn0search0turn0search1turn0search2
