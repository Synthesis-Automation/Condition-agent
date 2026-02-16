# there will be two applications, one for a simple/straightforward and well documented reaction, find a proven conditions usually work well. The other is for difficult or even undocumented reactions, give a list of conditions (e.g. diversified HTE-24 plate) based on literature reports of close reactions and mechanistic understanding

Below is a clean way to treat your two applications as **two “modes” of the same recommender**, with different **precision/recall targets, evidence thresholds, and diversity constraints**.

---

## Shared backbone (both apps)

### 0) Normalize + featurize the input reaction

* **Standardize reactants/products** (RDKit cleanup, salts/solvents stripping, normalize charges/tautomers).
* **(Optional but high value)** atom-map (RxnMapper) → extract:

  * changed bonds / reacting centers
  * coupling partners / leaving groups
  * role assignment (electrophile, nucleophile, catalyst-like, base-like)
  * Note: the RxnMapper may often fails in complex reactions (e.g. tandem processes), a llm thinking model may help in these cases. llm may also be better get mechanistic insights.

* Compute **reaction fingerprint** (DRFP / rxnfp-style) + **substrate features**

  * electrophile class (Ar–X, alkyl–X, sulfonate, etc.)
  * nucleophile class (amine/boronate/thiol/enolate…)
  * sterics/electronics flags (heteroaryl, EWG-rich, ortho-substitution, etc.)

### 1) Retrieve literature precedents

Store each precedent as:

* `rxn_embedding`, `rxn_template` (or bond-change signature)
* `ConditionCore` (catalyst system) + `ConditionTail` (base/solvent/additives/T/time)
* outcome label (yield / success / notes), plus metadata (scale, substrate scope)

Retrieval strategy:

* **Exact-ish**: same bond-change signature / template bucket
* **Similarity**: nearest neighbors by reaction fingerprint
* **Feature-match**: filter by substrate class (e.g., heteroaryl chloride + secondary amine) (we used motifs labels such as Ar-Br, Alkyl-B(OH)2 in this matching process)

### 2) Rank candidate conditions

Score each candidate with a blend of:

* **evidence strength** (count, yield distribution, recency, lab/source reliability if you track it)
* **match quality** (reaction similarity + substrate feature overlap)
* **constraint fit** (availability, glovebox, max temp, solvent restrictions, etc.)

### 3) Re-rank under constraints

Hard constraints (must satisfy): inventory, temp limit, banned reagents, safety rules.
Soft constraints (optimize): cost, convenience, green solvents, robustness.

---

# App 1 — “Proven conditions for straightforward/documented reactions”

**Goal:** high precision, low exploration. Think “Amazon *best seller* with lots of reviews.”
          Ranking system in one sentence:
          Rank each condition bundle by (1) how similar the precedent reactions are to your query and (2) how often that condition worked, but don’t trust tiny sample sizes.

## Policy differences vs the shared backbone

### Retrieval (tight)

* Require **high similarity** (same transformation + close substrate class).
* Prefer **exact template/bond-change match** over pure embedding similarity.

### Ranking (evidence-heavy)

* Score dominated by **precedent count + success rate**:

  * median/trimmed-mean yield (or success probability)
  * number of independent precedents
  * robustness across substrate scope (variance penalty)

### Output (deterministic + boring, by design)

Return:

1. **Primary condition** (most supported, most robust)
2. **Two close variants** (swap base or solvent only)
3. “Why” in quantifiable terms:

   * *N* close precedents, yield range, most common catalyst/base/solvent

## Practical guardrails

* If evidence is thin (e.g., <3 good precedents), **auto-fallback to App 2 behavior** (diversified HTE) rather than pretending confidence.

---

# App 2 — “Difficult/undocumented → diversified HTE-24 plate”

**Goal:** maximize probability of *finding at least one hit* fast, while learning mechanism. Think “YouTube exploration” or “bandit sampling,” but constrained by chemistry.

## Policy differences vs the shared backbone

### Retrieval (broad)

* Allow broader similarity neighborhoods:

  * same bond-change but different substrate class
  * same substrate class but adjacent transformation
* Pull **condition families**, not just exact recipes:

  * catalyst families (Pd/BrettPhos vs Pd/XPhos; Ni/bpy vs Ni/phen; Cu/diamine vs Cu/proline, etc.)
  * base strength families
  * solvent polarity families

### Candidate generation (mechanism-aware “menus”)

Build a **shortlist**:

* 6–8 **ConditionCores** (orthogonal mechanistic bets)
* 4–6 **ConditionTails** (solvent/base/additive/temperature “envelopes”)

Mechanistic heuristics usually improve coverage:

* oxidative addition challenged? → stronger ligands / Ni vs Pd / additives
* nucleophile problematic? → base/solvent swap, additive to suppress deactivation
* solubility / phase issues? → polar aprotic vs alcohol vs mixed solvent
* catalyst poisoning suspected? → scavengers / alternate metals / ligand class switch

### Diversification objective (this is the key)

Select 24 wells to maximize:

* **Predicted success** (from your ranker/ML model)
* **Diversity** across *meaningful axes* (core and tail), e.g.
  `metal`, `ligand family`, `base pKa band`, `solvent class`, `additive type`, `temperature band`

A simple, effective selection rule:

* Start with top ~200 candidates (core×tail combos)
* Cluster them in “condition space”
* Pick 24 by **max-min** distance (diverse) with a **score floor** (not stupid)

---

## A robust HTE-24 construction pattern (6 cores × 4 tails)

This pattern is easy to implement and easy to interpret.

### Step A — choose 6 orthogonal **ConditionCores**

Pick cores to cover *different failure modes*. Example structure (generic placeholders):

* Core-1: “Most proven family for this class” (best prior)
* Core-2: Same metal, different ligand family (sterics/electronics shift)
* Core-3: Alternative metal family (e.g., Ni vs Pd; Cu vs Pd)
* Core-4: More “forcing” catalyst system (reactivity-up)
* Core-5: More “forgiving” system (functional-group tolerant)
* Core-6: Radical/photoredox or electro variant **if** plausibly relevant (optional)

### Step B — choose 4 **ConditionTails** (rows)

Each row is a coherent “environment”:

* Tail-A: mild base + medium-polar solvent + moderate temp
* Tail-B: stronger base + polar aprotic + higher temp
* Tail-C: inorganic base + mixed solvent / water-tolerant + additive
* Tail-D: solubility rescue (different solvent class) + additive or lower temp

### Step C — assemble a 4×6 plate

* Columns = cores
* Rows = tails
* Interpretation becomes: “which *core* works?” and “which *environment* works?”

**Output schema (what your agent should emit)**

* `plate_id`
* For each well `A1..D6`:

  * catalyst system (metal source, ligand)
  * base
  * solvent
  * additive(s)
  * temperature/time (or a temp band)
  * notes/rationale tag (e.g., “OA-boost”, “solubility”, “poisoning-mitigation”)

---

## Implementation note: when to switch apps automatically

Define a **confidence gate** after retrieval:

* If you have **enough close precedents** (e.g., ≥5 within tight similarity + decent success), run **App 1**.
* Otherwise run **App 2** (and explicitly label it as exploratory).

---

## If you want, I can turn this into a concrete spec

I can write a compact “v1 design” that includes:

* your condition schema (`ConditionCore`/`ConditionTail`)
* scoring function (evidence + similarity + constraints)
* diversification selector (cluster + max-min)
* the exact JSON output for **App 1** and **HTE-24** in **App 2**.

If you share one example reaction SMILES for each case (easy vs hard), I’ll show exactly how the two modes behave end-to-end.
