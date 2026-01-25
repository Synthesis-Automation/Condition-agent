# Shopping Recommender Workflow and HTE Lessons

## Purpose

Summarize the standard ecommerce recommendation workflow (Amazon-style) and map it to the HTE condition recommender so we can borrow proven patterns.

## Explanation: Ecommerce Workflow

1. **Data ingest**
   - User events: views, clicks, add-to-cart, purchases.
   - Item metadata: category, price, availability, brand, text/images.
   - Context: time, device, geo, referrer, session intent.

2. **Candidate generation**
   - Fast retrieval from multiple sources: co-view/co-buy, item-to-item similarity, popularity, search intent, embeddings.
   - Union of candidates from several strategies, usually thousands.

3. **Filtering / constraints**
   - Remove out-of-stock, policy violations, duplicates.
   - Business rules (promotions, category limits).

4. **Ranking**
   - Score candidates with a model using relevance + conversion + value.
   - Personalization by user and session context.

5. **Re-ranking**
   - Add diversity, novelty, or category balance.
   - Avoid near-duplicates and stale items.

6. **Feedback loop**
   - Log outcomes, run A/B tests, retrain on new data.
   - Monitor drift and coverage.

## Mapping to HTE Recommender

**User** -> target reaction (reactants, products)  
**Item** -> condition set (catalyst, ligand, base, solvent, additive)  
**Context** -> reaction type, motif specificity, substrate class, constraints

Ecommerce step -> HTE equivalent:

- Candidate generation -> retrieve conditions from similar reactant-type pairs, DRFP similarity, and reaction-type priors.
- Filtering -> remove invalid conditions, rare combinations, or disallowed reagents.
- Ranking -> score by success rate, yield, motif specificity/reactivity, and dataset support.
- Re-ranking -> diversify across catalyst/ligand families and solvent classes.
- Feedback -> evaluate on held-out reactions and update weights.

## Implementation Plan (HTE)

### 1) Candidate generation (multi-source)

Sources to union:

- **Exact reactant-type match** (current behavior).
- **Relaxed match**: drop inorganic/empty reactants, allow A/B swap, or use parent motifs.
- **DRFP similarity** from reactant + product if available.
- **Reaction-type priors** from `data/rule_db_v2/` and `chemtools/taxonomy/data/compound_logic.json`.

### 2) Candidate scoring features

Suggested features:

- **HTE stats**: success rate, median/avg yield, sample size.
- **Motif specificity**: higher specificity scores higher than generic motifs.
- **Motif reactivity**: use `reactivity_weight` for reactive motifs.
- **Reaction-type prior**: boost conditions common in the predicted reaction family.
- **Product alignment**: boost conditions seen with similar formed motifs.

### 3) Re-ranking for diversity

Greedy selection with penalties:

- Penalize repeated catalysts/ligands and solvent families.
- Keep top-k while ensuring variety in catalytic systems.

### 4) Score calibration

Normalize scores across candidate sources:

- Per-source z-score or min-max scaling.
- Blend with fixed weights (config-driven).

### 5) Evaluation and logging

Offline metrics:

- Top-k hit rate, coverage, mean yield of retrieved conditions.
- NDCG or MAP for ranking quality.

Online or batch monitoring:

- Track changes in hit rate and coverage after updates.

## Minimal Code Touch Points

Potential modules to extend:

- `chemtools/recommend/recommender.py` for multi-source retrieval and scoring.
- `chemtools/featurizers/unified.py` for product-aware features.

## Next Steps

1. Define a candidate-source list and weights in config.
2. Add a feature vector for each condition row (stats + motif + reaction type).
3. Add a diversity-aware re-ranker.
4. Add a small offline evaluation script using `data/sample500/`.
