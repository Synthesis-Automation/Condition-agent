# Quantitative Validation Guide

## Overview

This guide explains how to assess the reliability of reaction analysis results using **quantitative validation metrics** that go beyond LLM self-reported confidence.

## Why Validation Matters

**Problem**: LLM self-reported confidence can be unreliable or overly optimistic.

**Solution**: Use multiple independent metrics that combine:
1. **Deterministic cheminformatics** (objective, tool-based)
2. **Cross-model consistency** (if multiple models agree, more likely correct)
3. **Specificity analysis** (detailed answers are more reliable than vague ones)
4. **Warning signals** (QC flags from the analysis pipeline)

---

## Validation Metrics Explained

### 1. Deterministic Quality Score (Most Important!)

**Weight in Overall Score**: 35% (highest)

**What it measures**: Quality of the cheminformatics analysis (independent of LLM)

**Components**:
- **Mapping Confidence** (50% weight): Atom mapping quality from rxnmapper
  - 0.9-1.0: Excellent
  - 0.7-0.9: Good
  - 0.5-0.7: Fair
  - <0.5: Poor (unreliable)

- **Bond Changes Detected** (25% weight): Number of bonds formed/broken/changed
  - More is better (indicates detailed analysis)

- **Reaction Center Identified** (15% weight): Core atoms involved in transformation
  - Yes/No indicator

- **Parse Quality** (10% weight): No SMILES parsing errors
  - Fewer warnings = better

**Interpretation**:
```
Deterministic Score >= 0.8  → Analysis is well-grounded in chemistry
Deterministic Score 0.6-0.8 → Reasonable but not perfect
Deterministic Score < 0.6   → Questionable chemical basis
```

**Key Insight**: This is the MOST RELIABLE metric because it's based on computational chemistry tools, not LLM opinions.

---

### 2. Specificity Score

**Weight in Overall Score**: 20%

**What it measures**: Level of detail in the LLM's analysis

**Components**:
- **Named Reactions** (30% weight): Did it identify specific reaction names (Suzuki, SNAr, etc)?
- **Events Detected** (30% weight): Number of mechanistic events described
- **Mechanism Steps** (20% weight): Number of steps in mechanism summary
- **Role Identification** (20% weight): Reagents, catalysts, electrophiles identified

**Interpretation**:
```
Specificity > 0.7  → Detailed, specific analysis (more trustworthy)
Specificity 0.4-0.7 → Generic but reasonable
Specificity < 0.4   → Vague or incomplete
```

**Rationale**: Detailed, specific answers indicate deeper understanding and are harder to fake.

---

### 3. Cross-Model Consistency (Requires 2+ Models)

**Weight in Overall Score**: 25%

**What it measures**: Agreement across different models

**Components**:
- **Class Agreement**: Do all models agree on reaction class?
- **Tag Overlap**: How many tags are shared (Jaccard similarity)?
- **Confidence Variance**: Do models have similar confidence levels?
- **Event Count Consistency**: Do models detect similar number of events?

**Interpretation**:
```
Consistency > 0.8  → Strong agreement (highly reliable)
Consistency 0.5-0.8 → Moderate agreement
Consistency < 0.5   → Models disagree (proceed with caution)
```

**Note**: This metric is only available when testing multiple models.

---

### 4. Ensemble Confidence (Requires 2+ Models)

**Weight in Overall Score**: Incorporated into consistency

**What it measures**: Weighted average confidence from multiple models

**Model Weights** (based on empirical testing):
```
o3:         1.00 (highest reliability for complex reactions)
o3-mini:    0.95
gpt-4o:     0.90
gpt-5.2:    0.85
gpt-4o-mini: 0.80
```

**Adjustment**: Penalized by disagreement factor (high variance = lower ensemble confidence)

**Interpretation**: More reliable than any single model's confidence.

---

### 5. Warning Score

**Weight in Overall Score**: 10%

**What it measures**: Presence of quality warnings

**Warnings Categorized by Severity**:
- **Critical**: mapping_failed, parse_failed (40% penalty each)
- **Moderate**: mapping_low_confidence, uncertain_bond_changes (20% penalty each)
- **Minor**: Other warnings (10% penalty each)

**Interpretation**:
```
Warning Score = 1.0  → No warnings (clean analysis)
Warning Score < 0.8  → Multiple issues detected
```

---

### 6. LLM Confidence (Least Trusted)

**Weight in Overall Score**: 10% (lowest)

**What it measures**: Self-reported confidence from the LLM

**Why Low Weight?**: LLMs can be overconfident or miscalibrated. Use other metrics to validate.

---

## Overall Validated Score

**Formula**:
```
With Cross-Model Data:
  Overall = 0.35 * Deterministic + 0.25 * Consistency + 0.20 * Specificity
          + 0.10 * Warnings + 0.10 * LLM_Confidence

Single Model:
  Overall = 0.45 * Deterministic + 0.30 * Specificity + 0.15 * Warnings
          + 0.10 * LLM_Confidence
```

**Reliability Ratings**:

| Deterministic | Overall | Rating | Recommendation |
|---------------|---------|--------|----------------|
| ≥ 0.8 | ≥ 0.7 | **HIGH** | Reliable for production use |
| ≥ 0.6 | ≥ 0.5 | **MEDIUM** | Reasonable, validate important details |
| ≥ 0.4 | ≥ 0.3 | **LOW** | Uncertain, manual review required |
| < 0.4 | < 0.3 | **VERY LOW** | Unreliable, trust with caution |

---

## Practical Examples

### Example 1: High Reliability (Rare C-N Coupling)

```
Overall Validated Score: 0.869
Reliability: MEDIUM

Individual Scores:
  deterministic_quality: 0.788  ← Good atom mapping (0.976)
  specificity:           0.733  ← Named reaction (SNAr) identified
  warning_score:         1.000  ← No warnings
  llm_confidence:        0.976  ← High (but note: lower weighted)
  consistency:           0.995  ← All 3 models agree
  ensemble_confidence:   0.955  ← Validated confidence

Key Metrics:
  Mapping Confidence: 0.976 (OK)
  Bond Changes:       1
  Reaction Center:    2 atoms
  Named Reactions:    True (SNAr)
  Events Detected:    1
  Cross-Model Consistency: 0.995

Interpretation:
✅ Strong deterministic foundation (0.788)
✅ All models agree (0.995 consistency)
✅ Specific mechanism identified (SNAr)
✅ Clean analysis (no warnings)

Recommendation: TRUST THIS RESULT
- Use for production workflows
- High confidence is validated by multiple independent metrics
```

### Example 2: Low Reliability (Complex Tandem Reaction)

```
Overall Validated Score: 0.298
Reliability: VERY LOW

Individual Scores:
  deterministic_quality: 0.412  ← Poor atom mapping (0.34)
  specificity:           0.650  ← Some detail but uncertain
  warning_score:         0.800  ← Moderate warnings
  llm_confidence:        0.300  ← Low (correctly reflects uncertainty)
  consistency:           0.450  ← Models disagree
  ensemble_confidence:   0.280  ← Low validated confidence

Key Metrics:
  Mapping Confidence: 0.340 (FAILED)  ← ⚠️ RED FLAG
  Bond Changes:       8
  Reaction Center:    12 atoms
  Named Reactions:    True (but uncertain)
  Events Detected:    2 (o3 only, others missed)
  Cross-Model Consistency: 0.450  ← ⚠️ RED FLAG

Interpretation:
❌ Weak deterministic foundation (0.412) - poor mapping
❌ Models disagree (0.450 consistency)
⚠️ Complex mechanism (8 bond changes, 2 stages)
⚠️ Some warnings present

Recommendation: DO NOT TRUST WITHOUT VALIDATION
- Manual expert review required
- Consider literature search
- May need experimental validation
- This is a genuinely difficult reaction for automated analysis
```

---

## When to Trust LLM Confidence

### ✅ Trust LLM Confidence When:

1. **Deterministic score ≥ 0.7**
   - Good atom mapping (>0.8)
   - Clear bond changes detected

2. **Cross-model consistency ≥ 0.8**
   - All models agree on class and mechanism

3. **Specificity ≥ 0.6**
   - Named reactions identified
   - Detailed events/steps provided

4. **No critical warnings**

### ❌ Distrust LLM Confidence When:

1. **Deterministic score < 0.5**
   - Poor atom mapping (<0.6)
   - Few/no bond changes detected

2. **Cross-model consistency < 0.5**
   - Models disagree on classification
   - High confidence variance

3. **Critical warnings present**
   - mapping_failed
   - parse_failed

4. **Low specificity (<0.4)**
   - Generic or vague analysis
   - No named reactions

---

## How to Use Validation in Practice

### 1. Single Model Analysis (Fast)

```bash
# Run analysis
python reaction_agent/cli.py --reaction "..." --model gpt-4o

# Check deterministic metrics only
python reaction_agent/cli.py --reaction "..." --no-llm
```

**Decision Rule**:
- Mapping confidence > 0.8 → Likely reliable
- Mapping confidence < 0.6 → Review manually

### 2. Cross-Model Validation (Recommended)

```python
from reaction_agent.scripts.quantitative_validation import validate_reaction

# Test with 3 models
validation = validate_reaction(
    rxn_smiles="...",
    models=['gpt-4o-mini', 'gpt-4o', 'o3-mini']
)

print(f"Overall Score: {validation['overall_score']:.3f}")
print(f"Reliability: {validation['reliability']}")
print(f"Recommendation: {validation['recommendation']}")
```

**Decision Rule**:
- Overall score ≥ 0.7 → Use result
- Overall score 0.5-0.7 → Validate key details
- Overall score < 0.5 → Manual review required

### 3. Batch Processing with Validation

```python
# Process reactions with automatic validation
for rxn in reactions:
    validation = validate_reaction(rxn, models=['gpt-4o', 'o3-mini'])

    if validation['reliability'] in ['HIGH', 'MEDIUM']:
        # Use result
        save_to_database(validation)
    else:
        # Flag for manual review
        flag_for_review(rxn, validation)
```

---

## Cost-Benefit Analysis

### Single Model (No Cross-Validation)

**Cost**: $0.001-0.003 per reaction
**Time**: 3-9 seconds
**Reliability**: Use deterministic score as proxy
**Best for**: High-throughput screening, well-characterized reactions

### Triple Model Validation (Cross-Validation)

**Cost**: $0.003-0.009 per reaction (3x)
**Time**: 9-27 seconds (3x)
**Reliability**: High confidence via ensemble metrics
**Best for**: Novel reactions, production decisions, literature preparation

**Recommendation**: Use single model for screening, validate hits with triple model.

---

## Key Takeaways

1. **Don't trust LLM confidence alone** - it can be miscalibrated

2. **Deterministic quality score is your most reliable indicator**
   - Based on rxnmapper and bond change analysis
   - Independent of LLM hallucinations

3. **Cross-model consistency is powerful**
   - If 3 different models agree, likely correct
   - Disagreement = proceed with caution

4. **Atom mapping quality predicts success**
   - Mapping >0.9 → Analysis confidence >0.9
   - Mapping <0.5 → Analysis unreliable

5. **Use validation for production decisions**
   - Screening: Check deterministic score (fast)
   - Important results: Run cross-model validation (thorough)

---

## Further Reading

- See `WORKFLOW_TESTING_RESULTS.md` for model comparison on complex reactions
- See `quantitative_validation.py` for implementation details
- See `test_workflows.py` for systematic testing across configurations
