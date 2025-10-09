# Simple Precedent Ranking - Visual Workflow

## The Three Strategies

```
┌────────────────────────────────────────────────────────────────────┐
│                    STRATEGY 1: SIMILARITY ONLY                      │
│                          (Baseline)                                 │
└────────────────────────────────────────────────────────────────────┘

Input: Reaction SMILES
   ↓
┌──────────────────────┐
│ Find k precedents    │  DRFP k-NN search
│ by similarity        │  → 30 precedents with similarity scores
└──────────────────────┘
   ↓
┌──────────────────────┐
│ Sort by similarity   │  Highest similarity = #1
└──────────────────────┘
   ↓
Output: Ranked precedents [Prec1(0.9), Prec2(0.85), Prec3(0.8), ...]


┌────────────────────────────────────────────────────────────────────┐
│                   STRATEGY 2: RULE-BASED RERANKING                  │
│                      (Chemistry-Driven)                             │
└────────────────────────────────────────────────────────────────────┘

Input: Reaction SMILES + Family
   ↓
┌──────────────────────┐
│ Find k precedents    │  DRFP k-NN search
└──────────────────────┘
   ↓
┌──────────────────────┐
│ Load rule database   │  Load family-specific rules
│ Match reaction       │  → Recommended: Cu/K3PO4/DMSO
└──────────────────────┘
   ↓
┌──────────────────────┐
│ Rerank precedents    │  For each precedent:
│ by rule alignment    │    - Catalyst match? +30% boost
│                      │    - Base match? +20% boost
│                      │    - Solvent match? +20% boost
│                      │  Score = Similarity × (1 + boost)
└──────────────────────┘
   ↓
Output: Ranked precedents [Prec2(0.85*1.7), Prec1(0.9), Prec3(0.8), ...]
        Rule-matching precedent moved to #1! ✅


┌────────────────────────────────────────────────────────────────────┐
│                 STRATEGY 3: ANALYTICS-BASED RERANKING               │
│                        (Data-Driven)                                │
└────────────────────────────────────────────────────────────────────┘

Input: Reaction SMILES + Family
   ↓
┌──────────────────────┐
│ Find k precedents    │  DRFP k-NN search
└──────────────────────┘
   ↓
┌──────────────────────┐
│ Get dataset stats    │  Top 10 catalysts: [Pd(10), PdCl2(9), ...]
│                      │  Top 10 bases: [K3PO4(10), K2CO3(9), ...]
│                      │  Top 10 solvents: [THF(10), DMF(9), ...]
└──────────────────────┘
   ↓
┌──────────────────────┐
│ Rerank precedents    │  For each precedent:
│ by popularity        │    - Popular catalyst? +0-30% boost
│                      │    - Popular base? +0-20% boost
│                      │    - Popular solvent? +0-20% boost
│                      │  Score = Similarity × (1 + boost)
└──────────────────────┘
   ↓
Output: Ranked precedents [Prec3(0.8*1.6), Prec1(0.9), Prec2(0.85), ...]
        High-yield popular conditions moved to #1! ✅
```

## Boost Calculation Examples

### Rule-Based Reranking

```
Rule match result: Cu catalyst, K3PO4 base, DMSO solvent

Precedent A:
  - Similarity: 0.50
  - Catalyst: Cu/phen ✅ matches Cu → +0.30
  - Base: K3PO4 ✅ matches → +0.20
  - Solvent: DMSO ✅ matches → +0.20
  - Total boost: 0.70
  - Final score: 0.50 × (1 + 0.70) = 0.85

Precedent B:
  - Similarity: 0.60
  - Catalyst: Pd/PPh3 ❌ no match → +0.00
  - Base: K2CO3 ❌ no match → +0.00
  - Solvent: Toluene ❌ no match → +0.00
  - Total boost: 0.00
  - Final score: 0.60 × (1 + 0.00) = 0.60

Result: Precedent A wins! (0.85 > 0.60) ✅
```

### Analytics-Based Reranking

```
Dataset statistics:
  - Pd rank: 10/10 (most popular)
  - K3PO4 rank: 10/10 (most popular)
  - THF rank: 8/10

Precedent A:
  - Similarity: 0.50
  - Catalyst: Pd (rank 10/10) → boost = (10/10) × 0.30 = 0.30
  - Base: K3PO4 (rank 10/10) → boost = (10/10) × 0.20 = 0.20
  - Solvent: THF (rank 8/10) → boost = (8/10) × 0.20 = 0.16
  - Total boost: 0.66
  - Final score: 0.50 × (1 + 0.66) = 0.83

Precedent B:
  - Similarity: 0.60
  - Catalyst: Rare catalyst (not in top 10) → boost = 0.00
  - Base: K2CO3 (rank 5/10) → boost = (5/10) × 0.20 = 0.10
  - Solvent: Rare solvent → boost = 0.00
  - Total boost: 0.10
  - Final score: 0.60 × (1 + 0.10) = 0.66

Result: Precedent A wins! (0.83 > 0.66) ✅
```

## Decision Tree: Which Strategy to Use?

```
Start: I have a reaction to plan
   ↓
   ┌─────────────────────────────────────┐
   │ Do strong chemical rules exist?     │
   │ (Ullmann→Cu, Buchwald→Pd, etc.)     │
   └─────────────────────────────────────┘
         Yes ↓                  No ↓
   ┌──────────────┐      ┌──────────────┐
   │ Use RULE     │      │ Is dataset   │
   │ reranking    │      │ large?       │
   │              │      │ (>10k rxns)  │
   │ ✅ Ensures   │      └──────────────┘
   │   correct    │         Yes ↓      No ↓
   │   chemistry  │    ┌──────────┐  ┌────────┐
   └──────────────┘    │ Use      │  │ Use    │
                       │ ANALYTICS│  │ NONE   │
                       │ reranking│  │        │
                       │          │  │ Simple │
                       │ ✅ Favors│  │ sim    │
                       │   robust │  └────────┘
                       │   popular│
                       │   reagent│
                       └──────────┘
```

## Real Example: Ullmann Reaction

```
Input: Brc1cnccn1 + Nc1ccccc1 → Product
Family: Ullmann_CN

Step 1: Find precedents
┌────────────────────────────────────────────────┐
│ Precedent 1: Cu/phen, K3PO4, DMSO     (0.50)  │
│ Precedent 2: Cu/phen, K3PO4, ethanol  (0.50)  │
│ Precedent 3: Cu/Oxine, K3PO4, DMSO    (0.50)  │
│ Precedent 4: Cu/TMEDA, Cs2CO3, DMF    (0.45)  │
└────────────────────────────────────────────────┘

Step 2: Load rules for Ullmann_CN
┌────────────────────────────────────────────────┐
│ Rule database: C_N_Coupling_Cu_db.json         │
│ Match result:                                  │
│   - Catalyst: (Cu implied by database)         │
│   - Base: K3PO4                                │
│   - Solvent: DMSO                              │
└────────────────────────────────────────────────┘

Step 3: Rerank by rule alignment
┌────────────────────────────────────────────────┐
│ Precedent 1: 0.50 × (1 + 0.3 + 0.2 + 0.2)     │
│              = 0.50 × 1.7 = 0.85  ← #1! ✅     │
│                                                │
│ Precedent 2: 0.50 × (1 + 0.3 + 0.2 + 0.0)     │
│              = 0.50 × 1.5 = 0.75               │
│                                                │
│ Precedent 3: 0.50 × (1 + 0.3 + 0.2 + 0.2)     │
│              = 0.50 × 1.7 = 0.85  ← #1! ✅     │
│                                                │
│ Precedent 4: 0.45 × (1 + 0.3 + 0.0 + 0.0)     │
│              = 0.45 × 1.3 = 0.59               │
└────────────────────────────────────────────────┘

Result:
  #1: Cu/phen, K3PO4, DMSO (perfect rule match!)
  #2: Cu/Oxine, K3PO4, DMSO (also perfect match)
  #3: Cu/phen, K3PO4, ethanol (2/3 match)
  ✅ All Cu catalysts (not Pd) - CORRECT!
```

## Code Comparison

### Old Fusion (Complex)
```python
# 100+ lines of code
evidence = collect_evidence(...)
weights = compute_adaptive_weights(evidence)  # α, β, γ, δ
candidates = generate_candidates(...)
scored = score_all_candidates(...)
final = rerank_by_alignment(...)

# Output: fusion_score=0.533, component_scores all 0.0 😕
```

### New Simple Ranker (Clear)
```python
# 10 lines of code
result = recommend_simple(
    reaction_smiles=rxn,
    family="Ullmann_CN",
    rerank_strategy='rule'
)

# Output: "Rule match: Cu/K3PO4/DMSO. Precedent #1 matches all three." ✅
```

## Summary

The simple precedent ranker achieves the goal:

✅ **Precedents are primary** - Direct ranking, not abstract fusion  
✅ **Quality control** - Rules/analytics prevent dataset errors  
✅ **Simple & transparent** - Clear workflow, no black box  
✅ **User control** - Pick strategy for use case  
✅ **Proven results** - Ullmann→Cu ✅, Suzuki→high-yield ✅
