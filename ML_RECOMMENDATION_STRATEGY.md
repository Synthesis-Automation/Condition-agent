# ML-Based Recommendation Strategy: Balancing Precedents, Analytics & Rules

## Problem Statement

Current ML recommendation has critical issues:
- ❌ **Over-reliance on random precedents**: Top-k precedents may not represent dataset well
- ❌ **Ignores dataset-level statistics**: Doesn't use common catalysts, bases, solvents
- ❌ **No rule-based integration**: Misses deterministic scheme matching
- ❌ **Weak confidence scoring**: Simple vote-share doesn't reflect data quality

**Goal**: Create a robust recommendation system that intelligently combines:
1. 🔍 **Precedent Search** (DRFP-based k-NN) - What worked for similar reactions
2. 📊 **Dataset Analytics** - What works commonly for this reaction family
3. 📋 **Rule-Based Matching** - What chemistry rules say should work
4. 🤖 **ML Yield Prediction** - What is likely to give high yields

## Core Philosophy

### The "Top-k Precedent Problem"

**Issue**: Top-k precedents are NOT a random sample of the dataset.

```python
# Current approach - PROBLEMATIC
precs = precedent.knn(family="Suzuki", k=25, reaction_smiles=rxn)
# Top 25 precedents might be:
# - All from same paper (batch effect)
# - All use expensive catalyst (researcher bias)
# - All use outdated conditions (temporal bias)
# - Not representative of successful conditions overall
```

**Example Problem**:
- Dataset has 50,215 Suzuki reactions
- 18,560 use K₂CO₃ as base (37% success rate)
- 1,129 use TEA as base (small subset)
- Top-k might return 15/25 with TEA (if query similar to TEA subset)
- **Result**: Recommends TEA when K₂CO₃ would be better!

### Solution: Multi-Source Evidence Fusion

Use **Bayesian-inspired weighting** that considers:

1. **Precedent Evidence** (PE): How similar reactions performed
   - Weight: Depends on similarity scores and diversity
   - Problem: Can be biased by query similarity

2. **Dataset Statistics** (DS): How this reaction family performs overall
   - Weight: High for well-characterized families
   - Problem: Ignores substrate specifics

3. **Rule-Based Evidence** (RE): Deterministic chemistry rules
   - Weight: High for exact SMARTS matches
   - Problem: Limited coverage

4. **ML Predictions** (ML): Yield predictions for candidate conditions
   - Weight: Depends on model confidence
   - Problem: Model may be wrong

## Proposed Architecture

### Stage 1: Evidence Collection

```python
def collect_evidence(reaction_smiles: str, family: str, k: int = 50):
    """
    Gather evidence from all sources.
    
    Returns:
        evidence = {
            'precedents': {
                'reactions': [...],  # Top-k similar reactions
                'similarities': [...],  # DRFP similarity scores
                'diversity_score': 0.65,  # How diverse are precedents?
                'coverage': 25  # Number of precedents found
            },
            'analytics': {
                'catalytic_systems': [(system, count, yield), ...],
                'bases': [(base, count, yield), ...],
                'solvents': [(solvent, count, yield), ...],
                'dataset_size': 50215,
                'high_yield_threshold': 80.0
            },
            'rules': {
                'matched_schemes': [...],  # Rule-based matches
                'confidence': 0.9,  # How confident is the match?
                'recommended': {...}  # Rule-based recommendation
            },
            'ml': {
                'model_available': True,
                'model_confidence': 0.85,
                'predictor': <DRFPYieldPredictor>
            }
        }
    """
```

### Stage 2: Candidate Generation

Generate diverse candidates from multiple sources:

```python
def generate_candidates(evidence: dict) -> List[Candidate]:
    """
    Generate candidate conditions from all evidence sources.
    
    Strategy:
        1. From precedents: Extract unique condition combinations
           - Limit to top 10 most frequent in precedents
           - Ensure diversity (different catalysts, bases, solvents)
        
        2. From analytics: Get dataset-level common conditions
           - Top 5 catalytic systems (≥80% yield)
           - Top 5 bases (by frequency)
           - Top 5 solvents (by frequency)
        
        3. From rules: Get rule-matched conditions
           - Exact SMARTS matches
           - High-confidence schemes
        
        4. Combine and deduplicate
           - Merge similar conditions
           - Limit to ~20-30 candidates
    
    Returns:
        List of Candidate objects with:
            - core, base, solvent, T_C, time_h
            - source: 'precedent', 'analytics', 'rules', or 'hybrid'
            - evidence_score: Initial score from source
    """
```

### Stage 3: Multi-Source Scoring

Score each candidate using ALL evidence:

```python
def score_candidate(candidate: Candidate, evidence: dict) -> ScoredCandidate:
    """
    Score a candidate using multi-source evidence fusion.
    
    Components:
        1. Precedent Score (PS): How well supported by similar reactions?
           PS = weighted_avg(similarity * (core_match + base_match + solvent_match))
           - Higher weight for more similar precedents
           - Penalty for low diversity (if all precedents too similar)
        
        2. Analytics Score (AS): How common in successful reactions?
           AS = (freq_score * 0.6 + yield_score * 0.4)
           - freq_score: Percentile rank in dataset frequency
           - yield_score: Average yield for this condition
        
        3. Rule Score (RS): Does it match chemistry rules?
           RS = rule_confidence * exact_match_bonus
           - 1.0 if exact SMARTS match
           - 0.5-0.8 if partial match
           - 0.0 if no match
        
        4. ML Score (MS): Predicted yield
           MS = predicted_yield / 100.0
           - Normalized to [0, 1]
           - Only if model available and confident
    
    Final Score Fusion (Bayesian-inspired):
        total_score = α*PS + β*AS + γ*RS + δ*MS
        
        Weights depend on evidence quality:
        - If k_precedents < 10: Increase β (analytics), decrease α
        - If diversity_score < 0.3: Decrease α (precedents unreliable)
        - If rule_confidence > 0.9: Increase γ
        - If model_confidence > 0.8: Increase δ
    
    Returns:
        ScoredCandidate with:
            - total_score
            - component_scores: {PS, AS, RS, MS}
            - weights: {α, β, γ, δ}
            - confidence: Overall confidence level
            - supporting_evidence: Count of supporting precedents
    """
```

### Stage 4: Adaptive Weighting

Adjust weights based on evidence quality:

```python
def compute_adaptive_weights(evidence: dict) -> dict:
    """
    Compute adaptive weights α, β, γ, δ based on evidence quality.
    
    Principles:
        1. Strong precedents (k>50, diversity>0.5, similarity>0.7):
           α=0.40, β=0.25, γ=0.15, δ=0.20
        
        2. Weak precedents (k<20, diversity<0.3):
           α=0.15, β=0.50, γ=0.20, δ=0.15
           → Trust analytics more than sparse precedents
        
        3. Rule-based match (confidence>0.9):
           α=0.25, β=0.20, γ=0.40, δ=0.15
           → Trust chemistry rules
        
        4. Large dataset (>10,000 reactions):
           α=0.30, β=0.35, γ=0.15, δ=0.20
           → Analytics very reliable
        
        5. Small dataset (<1,000 reactions):
           α=0.40, β=0.15, γ=0.25, δ=0.20
           → Precedents and rules more important
    
    Returns:
        {'α': 0.35, 'β': 0.30, 'γ': 0.20, 'δ': 0.15, 'reasoning': '...'}
    """
    
    weights = {'α': 0.35, 'β': 0.30, 'γ': 0.20, 'δ': 0.15}  # Default
    reasoning = []
    
    # Adjust based on precedent quality
    n_prec = evidence['precedents']['coverage']
    diversity = evidence['precedents']['diversity_score']
    avg_sim = mean(evidence['precedents']['similarities'])
    
    if n_prec < 10:
        weights['α'] *= 0.5  # Reduce precedent weight
        weights['β'] *= 1.5  # Increase analytics weight
        reasoning.append("Few precedents → rely more on analytics")
    
    if diversity < 0.3:
        weights['α'] *= 0.7  # Precedents not diverse
        reasoning.append("Low diversity → precedents may be biased")
    
    if avg_sim < 0.6:
        weights['α'] *= 0.8  # Precedents not very similar
        reasoning.append("Low similarity → precedents less relevant")
    
    # Adjust based on dataset size
    dataset_size = evidence['analytics']['dataset_size']
    if dataset_size > 10000:
        weights['β'] *= 1.3  # Analytics very reliable
        reasoning.append("Large dataset → trust analytics")
    elif dataset_size < 1000:
        weights['β'] *= 0.7  # Analytics less reliable
        reasoning.append("Small dataset → analytics less reliable")
    
    # Adjust based on rule match
    if evidence['rules']['confidence'] > 0.9:
        weights['γ'] *= 2.0  # Strong rule match
        reasoning.append("Strong rule match → trust chemistry")
    
    # Adjust based on ML model confidence
    if evidence['ml']['model_available']:
        ml_conf = evidence['ml'].get('model_confidence', 0.5)
        if ml_conf > 0.85:
            weights['δ'] *= 1.2
            reasoning.append("High ML confidence")
        elif ml_conf < 0.6:
            weights['δ'] *= 0.6
            reasoning.append("Low ML confidence")
    else:
        weights['δ'] = 0.0
        reasoning.append("No ML model available")
    
    # Normalize to sum to 1.0
    total = sum(weights.values())
    weights = {k: v/total for k, v in weights.items()}
    
    return {'weights': weights, 'reasoning': reasoning}
```

## Implementation Details

### Precedent Diversity Score

Measure how diverse the top-k precedents are:

```python
def compute_diversity_score(precedents: List[dict]) -> float:
    """
    Measure diversity of precedent conditions.
    
    Diversity = average pairwise distance between conditions
    
    Distance metrics:
        - Different catalyst: +1.0
        - Different base: +0.5
        - Different solvent: +0.5
        - Different temperature (>20°C apart): +0.3
    
    High diversity (>0.7): Precedents cover many conditions
    Low diversity (<0.3): Precedents are very similar (potential bias)
    
    Returns:
        Float in [0, 1] representing diversity
    """
    if len(precedents) < 2:
        return 1.0
    
    distances = []
    for i in range(len(precedents)):
        for j in range(i+1, len(precedents)):
            p1, p2 = precedents[i], precedents[j]
            
            dist = 0.0
            if p1.get('core') != p2.get('core'):
                dist += 1.0
            if p1.get('base_uid') != p2.get('base_uid'):
                dist += 0.5
            if p1.get('solvent_uid') != p2.get('solvent_uid'):
                dist += 0.5
            
            T1 = p1.get('T_C', 80)
            T2 = p2.get('T_C', 80)
            if abs(T1 - T2) > 20:
                dist += 0.3
            
            distances.append(dist)
    
    max_dist = 2.3  # Max possible distance
    avg_dist = mean(distances)
    diversity = min(1.0, avg_dist / max_dist)
    
    return diversity
```

### Analytics Integration

```python
def get_analytics_evidence(family: str, min_yield: float = 75.0) -> dict:
    """
    Get dataset analytics for condition ranking.
    
    Returns:
        {
            'catalytic_systems': [
                ('Pd(OAc)2 + XPhos', 156, 82.3),  # (system, count, avg_yield)
                ...
            ],
            'bases': [('K2CO3', 18560, 81.3), ...],
            'solvents': [('Water', 35433, 80.1), ...],
            'condition_cores': [('Pd/XPhos', 1823, 81.6), ...],
            'dataset_size': 50215,
            'high_yield_conditions': [...]  # Conditions with >min_yield
        }
    """
    from chemtools import chem
    
    stats = chem.analytics.get_stats(family)
    
    return {
        'catalytic_systems': chem.analytics.get_common_catalytic_systems(
            family, top_n=20, min_yield=min_yield
        ),
        'bases': chem.analytics.get_common_bases(
            family, top_n=15, min_yield=None
        ),
        'solvents': chem.analytics.get_common_solvents(
            family, top_n=15, min_yield=None
        ),
        'condition_cores': chem.analytics.get_condition_cores(
            family, top_n=15
        ),
        'dataset_size': stats['total_reactions'],
        'yield_stats': stats.get('yield_stats', {})
    }
```

### Candidate Scoring Example

```python
def score_with_fusion(candidate: Candidate, evidence: dict) -> float:
    """
    Example scoring for Suzuki coupling candidate:
    
    Candidate: Pd(OAc)2 + XPhos, K2CO3, Dioxane, 100°C, 12h
    
    1. Precedent Score (PS):
       - 8/25 precedents use Pd(OAc)2 + XPhos
       - 12/25 use K2CO3
       - 15/25 use Dioxane
       - Avg similarity to precedents: 0.78
       → PS = 0.78 * (8/25 + 12/25 + 15/25) / 3 = 0.78 * 0.47 = 0.36
    
    2. Analytics Score (AS):
       - "Pd(OAc)2 + XPhos": 156 reactions, 82.3% yield (rank: 15/137)
       - "K2CO3": 18,560 reactions, 81.3% yield (rank: 1/41)
       - "Dioxane": 14,870 reactions, 77.9% yield (rank: 2/137)
       → freq_score = (0.89 + 1.0 + 0.99) / 3 = 0.96
       → yield_score = (82.3 + 81.3 + 77.9) / 300 = 0.805
       → AS = 0.6 * 0.96 + 0.4 * 0.805 = 0.898
    
    3. Rule Score (RS):
       - No exact SMARTS match found
       - Partial match (aryl halide + amine → general Pd conditions)
       → RS = 0.5
    
    4. ML Score (MS):
       - Predicted yield: 79.2%
       → MS = 0.792
    
    Adaptive Weights:
       - k=25 precedents, diversity=0.55, avg_similarity=0.78
       - Dataset size=50,215 (large!)
       - No strong rule match
       → α=0.30, β=0.40, γ=0.15, δ=0.15
    
    Total Score:
       = 0.30*0.36 + 0.40*0.898 + 0.15*0.5 + 0.15*0.792
       = 0.108 + 0.359 + 0.075 + 0.119
       = 0.661
    
    Confidence: HIGH (analytics very strong, good precedent support)
    """
```

## Handling Edge Cases

### Case 1: Novel Substrate (Few Precedents)

```python
# Evidence:
n_precedents = 5
diversity_score = 0.8  # Diverse but sparse
avg_similarity = 0.45  # Not very similar

# Strategy:
weights = {
    'α': 0.15,  # Low weight on sparse precedents
    'β': 0.55,  # HIGH weight on analytics
    'γ': 0.20,  # Moderate weight on rules
    'δ': 0.10   # Low weight on ML (not enough data)
}

# Reasoning:
# "Few similar precedents found. Relying heavily on dataset-level
#  statistics for this reaction family. Use common catalytic systems,
#  bases, and solvents that work well generally."
```

### Case 2: Well-Studied Substrate (Many Precedents)

```python
# Evidence:
n_precedents = 150
diversity_score = 0.65  # Diverse and abundant
avg_similarity = 0.85  # Very similar

# Strategy:
weights = {
    'α': 0.50,  # HIGH weight on precedents
    'β': 0.20,  # Moderate weight on analytics
    'γ': 0.10,  # Low weight on rules (precedents sufficient)
    'δ': 0.20   # Moderate weight on ML
}

# Reasoning:
# "Many similar precedents found with high diversity. Strong evidence
#  for what works with this specific substrate structure."
```

### Case 3: Exact Rule Match

```python
# Evidence:
rule_confidence = 0.95  # Exact SMARTS match
matched_scheme = "Buchwald-Hartwig C-N coupling"

# Strategy:
weights = {
    'α': 0.25,  # Moderate weight on precedents
    'β': 0.20,  # Moderate weight on analytics
    'γ': 0.45,  # HIGH weight on rules
    'δ': 0.10   # Low weight on ML
}

# Reasoning:
# "Exact chemistry rule match found. Following established protocol
#  for this specific reaction type."
```

### Case 4: Analytics Dominant (Small Precedent Set)

```python
# Evidence:
n_precedents = 12
dataset_size = 41,427  # Amide formation
diversity_score = 0.25  # Low diversity (batch effect?)

# Strategy:
weights = {
    'α': 0.10,  # Very low weight (biased precedents)
    'β': 0.60,  # Very HIGH weight on analytics
    'γ': 0.15,  # Moderate weight on rules
    'δ': 0.15   # Moderate weight on ML
}

# Reasoning:
# "Precedents show low diversity and may be biased. Large dataset
#  available - trusting statistical patterns from 41,427 reactions."
```

## Output Format

```python
{
    'recommended_conditions': [
        {
            'rank': 1,
            'core': 'Pd(OAc)2 + XPhos',
            'base': 'K2CO3',
            'solvent': '1,4-Dioxane',
            'T_C': 100,
            'time_h': 12,
            'total_score': 0.661,
            'confidence': 'HIGH',
            'component_scores': {
                'precedent_score': 0.36,
                'analytics_score': 0.898,
                'rule_score': 0.5,
                'ml_score': 0.792
            },
            'weights': {'α': 0.30, 'β': 0.40, 'γ': 0.15, 'δ': 0.15},
            'evidence': {
                'supporting_precedents': 8,  # Out of 25
                'dataset_frequency': 156,  # Times seen in dataset
                'dataset_yield': 82.3,  # Average yield
                'rule_match': 'partial',
                'predicted_yield': 79.2
            },
            'reasoning': [
                "Top analytics choice: K2CO3 used in 18,560 reactions (37%)",
                "Catalytic system Pd(OAc)2+XPhos: 156 reactions, 82.3% yield",
                "8/25 precedents support this combination",
                "ML predicts 79.2% yield"
            ]
        },
        ...
    ],
    'meta': {
        'method': 'multi_source_fusion',
        'evidence_summary': {
            'precedents': '25 found, diversity=0.55, avg_similarity=0.78',
            'analytics': 'Dataset: 50,215 reactions',
            'rules': 'Partial match',
            'ml': 'Model available, confidence=0.85'
        },
        'adaptive_weights': {
            'α': 0.30,
            'β': 0.40,
            'γ': 0.15,
            'δ': 0.15,
            'reasoning': [
                "Large dataset → trust analytics (β↑)",
                "Good precedent coverage → moderate precedent weight",
                "No strong rule match → low rule weight"
            ]
        }
    }
}
```

## Implementation Roadmap

### Phase 1: Evidence Collection Infrastructure

1. ✅ **Precedent Search** (Already implemented)
   - DRFP-based k-NN
   - Binary NPZ precomputed fingerprints

2. ✅ **Dataset Analytics** (Already implemented)
   - `get_common_catalytic_systems()`
   - `get_common_bases()`, `get_common_solvents()`
   - `get_condition_cores()`

3. 🔄 **Integration Layer** (New)
   - `collect_evidence()` function
   - Unified data structures

### Phase 2: Scoring System

1. 📝 **Component Scorers** (New)
   - `score_from_precedents()`
   - `score_from_analytics()`
   - `score_from_rules()`
   - `score_from_ml()`

2. 📝 **Fusion Logic** (New)
   - `compute_adaptive_weights()`
   - `fuse_scores()`
   - `compute_confidence()`

### Phase 3: Candidate Generation

1. 📝 **Multi-Source Candidates** (New)
   - Extract from precedents (current logic)
   - Add from analytics (top systems)
   - Add from rules (matched schemes)
   - Deduplication and merging

### Phase 4: Testing & Validation

1. 📝 **Unit Tests**
   - Test each scorer independently
   - Test weight adaptation
   - Test edge cases

2. 📝 **Integration Tests**
   - End-to-end recommendations
   - Compare with current system
   - Measure improvement

3. 📝 **Benchmarking**
   - Run on test reactions
   - Measure yield prediction accuracy
   - Measure user satisfaction (if possible)

## Example: Complete Flow

```python
# Input
reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
family = "C_N_Coupling_Pd"

# Step 1: Collect evidence
evidence = collect_evidence(reaction, family, k=25)

# Step 2: Generate candidates (from all sources)
candidates = []

# From precedents
for prec in evidence['precedents']['reactions'][:10]:
    candidates.append(Candidate(
        core=prec['core'],
        base=prec['base_uid'],
        solvent=prec['solvent_uid'],
        T_C=prec.get('T_C', 100),
        time_h=prec.get('time_h', 12),
        source='precedent'
    ))

# From analytics (high-yield catalytic systems)
for system, count, avg_yield in evidence['analytics']['catalytic_systems'][:5]:
    # Parse system string "Pd(OAc)2 + XPhos" → core
    # Combine with top bases and solvents
    for base, _, _ in evidence['analytics']['bases'][:3]:
        for solvent, _, _ in evidence['analytics']['solvents'][:3]:
            candidates.append(Candidate(
                core=system,  # Complete catalytic system
                base=base,
                solvent=solvent,
                T_C=100,  # Default
                time_h=12,  # Default
                source='analytics'
            ))

# From rules (if matched)
if evidence['rules']['matched_schemes']:
    rule_rec = evidence['rules']['recommended']
    candidates.append(Candidate(
        core=rule_rec['core'],
        base=rule_rec['base'],
        solvent=rule_rec['solvent'],
        T_C=rule_rec['T_C'],
        time_h=rule_rec['time_h'],
        source='rules'
    ))

# Step 3: Deduplicate candidates
unique_candidates = deduplicate(candidates)

# Step 4: Compute adaptive weights
weights_info = compute_adaptive_weights(evidence)
weights = weights_info['weights']

# Step 5: Score each candidate
scored_candidates = []
for candidate in unique_candidates:
    ps = score_from_precedents(candidate, evidence['precedents'])
    as_ = score_from_analytics(candidate, evidence['analytics'])
    rs = score_from_rules(candidate, evidence['rules'])
    ms = score_from_ml(candidate, evidence['ml'], reaction)
    
    total_score = (
        weights['α'] * ps +
        weights['β'] * as_ +
        weights['γ'] * rs +
        weights['δ'] * ms
    )
    
    scored_candidates.append(ScoredCandidate(
        candidate=candidate,
        total_score=total_score,
        component_scores={'PS': ps, 'AS': as_, 'RS': rs, 'MS': ms},
        weights=weights
    ))

# Step 6: Rank and return top-N
ranked = sorted(scored_candidates, key=lambda x: x.total_score, reverse=True)
top_n = ranked[:5]

# Step 7: Format output with explanations
results = format_recommendations(top_n, evidence, weights_info)
```

## Key Advantages

### 1. Robust to Precedent Bias
- ✅ Doesn't blindly trust top-k precedents
- ✅ Validates against dataset-level statistics
- ✅ Detects low-diversity precedent sets

### 2. Leverages All Available Evidence
- ✅ Precedents: What worked for similar substrates
- ✅ Analytics: What works generally for this family
- ✅ Rules: What chemistry says should work
- ✅ ML: What is predicted to work

### 3. Adaptive to Data Quality
- ✅ Weak precedents → Trust analytics more
- ✅ Strong precedents → Trust precedents more
- ✅ Exact rule match → Trust chemistry rules
- ✅ Novel substrate → Rely on dataset patterns

### 4. Interpretable Results
- ✅ Shows component scores
- ✅ Explains weight choices
- ✅ Cites supporting evidence
- ✅ Quantifies confidence

### 5. Prevents Bad Recommendations
- ❌ **Prevents**: Recommending rare conditions just because they appear in top-k
- ✅ **Ensures**: Common, successful conditions are properly weighted
- ✅ **Validates**: Precedents against dataset-level success rates

## Conclusion

This multi-source fusion approach:
- 🎯 **Solves the "top-k precedent problem"** by balancing with dataset analytics
- 📊 **Uses all available evidence** intelligently with adaptive weighting
- 🔍 **Provides interpretability** through component scores and explanations
- 🚀 **Improves recommendation quality** by avoiding bias toward unrepresentative precedents

**Next Steps**:
1. Implement `collect_evidence()` function
2. Implement component scorers
3. Implement adaptive weighting logic
4. Integrate with existing recommendation pipeline
5. Test and benchmark against current system
