# Lightweight, Evolvable Upgrades for the Condition Recommendation System

This document refines your current pipeline:

- **Input:** reaction SMILES (e.g., `Brc1ccccc1.Nc1ccccc1>>`)
- **Normalize:** `chemtools.smiles.normalize_reaction`
- **Family detection:** `chemtools.router.detect_family`
- **Featurize:** `chemtools.featurizers.molecular.featurize` (+ coarse `bin`), drop `role_aware` for caching
- **Precedents:** `chemtools.precedent.knn` (+ optional DRFP re‑ranking: `relax.use_drfp`, `drfp_weight`, `drfp_n_bits`, `drfp_radius`)
- **Core selection:** Laplace‑smoothed vote over precedents + simple confidence
- **Base/Solvent:** most frequent for chosen core, filtered by `chemtools.constraints.apply_filter` (e.g., `no_chlorinated`, `no_HMPA`, `aqueous_only`, `min_bp_C`)
- **T/t:** medians from precedents, fallback to all
- **Explain:** `chemtools.explain.for_pack`
- **Output:** `{recommendation, alternatives, precedent_pack, reasons, filters, formatted}`

---

## P0 (1–2 days): High‑impact, low effort

### 1) Calibrated confidence + abstain
- Convert vote share → **calibrated probability** (isotonic/Platt) on a scaffold‑split holdout.
- **Abstain** when support is thin: `n_precedents < 20` or `p_top < 0.45`. Return 2 cores + a 4‑well HTE probe.

### 2) Family smoothing for Core votes (fix rare labels)
Shrink each core’s vote toward its **family** (e.g., biaryl‑P = XPhos/SPhos/RuPhos/BrettPhos):
```
shrink = n / (n + τ)        # τ≈20
score_core = shrink * vote_core + (1 - shrink) * mean(vote_family_excl_core)
```

### 3) Typed guardrails for tails
- **Suzuki:** require `H2O ≥ 5%`; allow `KF` add‑on for tough Bpin; avoid strong alkoxides.
- **C–O:** forbid `H2O > 1%`, prefer alkoxides.
- **Borylation:** prefer **KOAc**; ether solvents.

### 4) Robust temperature/time
- **Trimmed medians** (drop top/bottom 10%).  
- Clip `T` by solvent BP: `T = min(median_T, bp(solvent) − 10 °C, 120 °C)`.

### 5) Density‑aware *k* for kNN
- Start at `k=200`; if sim@rank200 < 0.2, reduce to 100; if > 0.4, increase to 300.

---

## P1 (1–2 weeks): Accuracy & scope, still lightweight

### 6) Dual‑index retrieval
- Keep **substrate ECFP kNN** (current) and add a **reaction fingerprint** index (DRFP or RXNFP‑lite).
- Retrieve from both; **blend** (e.g., `0.7*ECFP + 0.3*DRFP`) before voting.

### 7) Similar‑core expansion
Maintain a small **core similarity map** (JSON: families + descriptor neighbors such as `%V_bur`, Tolman angle). After scoring, append nearest neighbors of the top core(s) with a ×0.95 multiplier so they can compete in the tail step.

### 8) Tiny candidate ranker (ordering boost)
Train a **LightGBM ranker** to order (Core, Base, Solvent) tuples using existing features:
- vote share, support, DRFP score
- family‑smoothed vote
- practicality (air‑stable, solvent BP, green score)
- rule‑fit flags (e.g., `Suzuki_has_water`)

**Labels:** success (≥Y₀) or yield. Expect **+5–15 nDCG@5** and **+8–20% Top‑1 success** vs. votes alone.

### 9) Better reasons
- Add **margin of victory** (“XPhos 0.62 vs 0.29 RuPhos”)
- Add a **data badge** (“k=180 precedents; 24 close matches”)
- If abstaining, state what’s missing (“few close precedents for heteroaryl chloride + sulfonamide”).

---

## P2 (when ready): ML hooks that evolve with data

### 10) Core classifier head (cheap lift)
Small multiclass on (ECFP of reactants + reaction type). Blend: `score = 0.6·P_core + 0.4·vote`.

### 11) Core embeddings (generalization)
Build a mini **core graph** (edges = same family + descriptor similarity). Compute a 16–32D embedding (Node2Vec/PCA) and:
- add as ranker features,
- use for **neighbor smoothing** (weighted borrowing from similar cores).

### 12) Online learning loop
Log every recommendation + outcome; continuously refresh **priors** (family shrinkage, frequencies) and periodically refresh the **ranker**. No need to retrain the entire stack on every data update.

---

## Drop‑in snippets

**Family smoothing for votes**
```python
def smooth_votes(votes, core2fam, tau=20):
    # votes: dict core->weighted_count; core2fam: core->fam_name
    fam_means = {}
    fams = set(core2fam.get(c, '') for c in votes)
    for fam in fams:
        fam_cores = [c for c in votes if core2fam.get(c)==fam]
        if not fam_cores:
            continue
        fam_means[fam] = sum(votes[c] for c in fam_cores) / len(fam_cores)
    out = {}
    for c,v in votes.items():
        n = v  # proxy for support; or pass real neighbor count per core
        a = n/(n+tau)
        out[c] = a*v + (1-a)*fam_means.get(core2fam.get(c,''), 0.0)
    return out
```

**Abstain threshold**
```python
def should_abstain(vote_top, support_k, p_calibrated=None):
    if support_k < 50:
        return True
    if p_calibrated is not None:
        return p_calibrated < 0.45
    return vote_top < 0.35
```

**Typed guardrail example**
```python
def on_manifold(rxn_type, core, base, solvent, water_pct):
    if rxn_type=="Suzuki" and (water_pct or 0) < 5: return False
    if rxn_type=="CO_Ether" and (water_pct or 0) > 1: return False
    if rxn_type=="Borylation" and base!="KOAc": return False
    return True
```

---

## QA & monitoring

- **Holdout:** scaffold split; track Top‑1/Top‑3 core accuracy and nDCG@5.
- **Calibration:** reliability plots for confidence.
- **Drift alarms:** watch median neighbor similarity and abstain‑rate spikes.

---

## Summary

You can keep the architecture you have and get meaningful gains by:
1) **Calibrated confidence + abstain**, **family smoothing**, **typed guardrails**, and **robust T/t** (P0).  
2) Add **dual‑index retrieval** and a **tiny ranker** (P1).  
3) Later, layer a **core classifier**, **core embeddings**, and an **online loop** (P2).

This path improves accuracy, safety, and user trust without adding heavy complexity—and it positions the system to learn as your dataset grows.
