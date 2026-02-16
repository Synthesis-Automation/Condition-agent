# Quick Reference: Quantitative Validation for Reaction Analysis

## TL;DR - How to Get Reliable Quality Metrics

**Question**: *"Is the LLM confidence score reliable?"*

**Answer**: **Not by itself.** Use the quantitative validation framework to get multiple independent metrics.

---

## Quick Commands

### 1. Basic Validation (3 models, comprehensive):
```bash
python reaction_agent/scripts/quantitative_validation.py
```

### 2. Custom Validation (your reaction, your models):
```python
from reaction_agent.scripts.quantitative_validation import validate_reaction

validation = validate_reaction(
    rxn_smiles="Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1",
    models=['gpt-4o-mini', 'gpt-4o', 'o3-mini'],  # 2-5 models
    max_tokens=3000
)

print(f"Overall Score: {validation['overall_score']:.3f}")
print(f"Reliability: {validation['reliability']}")
```

### 3. Fast Check (single model, deterministic only):
```bash
python reaction_agent/cli.py --reaction "..." --no-llm
```
Check the mapping confidence: >0.8 = reliable, <0.6 = questionable.

---

## Decision Matrix

| Deterministic Score | Overall Score | **Action** |
|---------------------|---------------|------------|
| **≥ 0.8** | **≥ 0.7** | ✅ **TRUST** - Use in production |
| **0.6-0.8** | **0.5-0.7** | ⚠️ **VALIDATE** - Check key details |
| **0.4-0.6** | **0.3-0.5** | ⚠️ **REVIEW** - Manual expert review |
| **< 0.4** | **< 0.3** | ❌ **CAUTION** - Don't trust, reanalyze |

---

## The 5-Second Check

```python
# Run your analysis
result = analyzer.analyze(rxn_smiles)

# Check ONE metric - mapping confidence
mapping_conf = result['tool_facts']['mapping_qc']['confidence']

if mapping_conf > 0.8:
    print("✓ Likely reliable")
elif mapping_conf > 0.6:
    print("⚠ Moderate - validate details")
else:
    print("❌ Questionable - review manually")
```

**Why this works**: Atom mapping quality is the single best predictor of overall reliability.

---

## Real Example Comparison

### Example 1: Simple SNAr Reaction
```
Reaction: Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1
LLM Confidence:         0.976  (looks good!)
Mapping Confidence:     0.976  ✓ Validates the LLM
Deterministic Score:    0.788  ✓ Strong foundation
Overall Validated:      0.795  ✓ Reliable

Conclusion: ✅ TRUST THIS RESULT
```

### Example 2: Complex Tandem Reaction
```
Reaction: c1ccc2c(c1)[nH]c1ccccc12.Fc1ccc(I)cc1.C#Cc1ccc(CC)cc1>>...
LLM Confidence:         0.200  (looks bad)
Mapping Confidence:     0.335  ❌ Poor mapping
Deterministic Score:    0.100  ❌ Weak foundation
Overall Validated:      0.352  ❌ Unreliable

Conclusion: ❌ DON'T TRUST - Manual review needed
```

**Key Insight**: In both cases, atom mapping quality predicted the outcome!

---

## What Each Metric Tells You

| Metric | What It Measures | Why It Matters |
|--------|------------------|----------------|
| **Deterministic Quality** | Cheminformatics tools (rxnmapper, RDKit) | Most objective - no LLM bias |
| **Mapping Confidence** | Atom mapping quality | Best single predictor |
| **Specificity** | Level of detail (named reactions, events) | Detailed = more trustworthy |
| **Consistency** | Agreement across models | Multiple models agree = reliable |
| **Ensemble Confidence** | Weighted average of models | Better than any single confidence |
| **LLM Confidence** | Self-reported by model | Least reliable alone |

---

## Validation Workflow

### For High-Throughput Screening (Fast)
```python
# Step 1: Run with single fast model
result = analyzer.analyze(rxn_smiles, model='gpt-4o-mini')

# Step 2: Check deterministic score
mapping = result['tool_facts']['mapping_qc']['confidence']

# Step 3: Decision
if mapping > 0.8:
    save_result(result)  # Likely reliable
else:
    flag_for_validation(result)  # Needs review
```

**Cost**: ~$0.001 per reaction
**Time**: ~5 seconds

### For Production/Important Decisions (Thorough)
```python
# Step 1: Validate with multiple models
validation = validate_reaction(
    rxn_smiles,
    models=['gpt-4o-mini', 'gpt-4o', 'o3-mini']
)

# Step 2: Check overall score
if validation['overall_score'] >= 0.7:
    save_result(validation)  # Reliable
elif validation['overall_score'] >= 0.5:
    validate_details(validation)  # Reasonable, check key points
else:
    manual_review(validation)  # Unreliable
```

**Cost**: ~$0.003-0.009 per reaction (3x models)
**Time**: ~9-27 seconds
**Benefit**: High confidence in results

---

## Common Scenarios

### Scenario 1: "LLM says 0.95 confidence, should I trust it?"

**Answer**: Check the deterministic score!

```python
validation = validate_reaction(rxn_smiles)

if validation['individual_scores']['deterministic_quality'] > 0.7:
    print("✓ Yes, the confidence is backed by good cheminformatics")
else:
    print("❌ No, LLM is overconfident - weak chemical basis")
```

### Scenario 2: "Models disagree - which is right?"

**Answer**: Look at consistency and ensemble confidence!

```python
validation = validate_reaction(rxn_smiles, models=['gpt-4o', 'o3', 'o3-mini'])

consistency = validation['individual_scores']['consistency']

if consistency > 0.8:
    print(f"✓ Models agree - trust consensus: {validation['detailed_metrics']['cross_model']['consensus_class']}")
else:
    print("❌ Models disagree - proceed with caution, manual review needed")
```

### Scenario 3: "Reaction looks complex - which model should I use?"

**Answer**: Check deterministic analysis first!

```bash
# Step 1: Check mapping quality (no LLM needed)
python reaction_agent/cli.py --reaction "..." --no-llm
```

If mapping confidence < 0.6:
- Use **o3** or **o3-mini** (4000+ tokens) for complex reactions
- Expect lower overall confidence (it's a genuinely hard problem!)
- Consider manual expert review

If mapping confidence > 0.8:
- Use **gpt-4o** or **gpt-4o-mini** (cheaper, faster)
- Expect high reliability

---

## Integration Examples

### Add to CLI:
```bash
# Add --validate flag for cross-model validation
python reaction_agent/cli.py \
  --reaction "..." \
  --validate \
  --validation-models gpt-4o-mini,gpt-4o,o3-mini
```

### Add to Batch Processing:
```python
# Process reactions with automatic quality gating
for rxn in reaction_list:
    validation = validate_reaction(rxn, models=['gpt-4o', 'o3-mini'])

    if validation['reliability'] in ['HIGH', 'MEDIUM']:
        results.append(validation)
    else:
        flagged.append((rxn, validation))

# Review flagged reactions manually
review_flagged_reactions(flagged)
```

### Add to Database:
```python
# Store validation metrics alongside results
db.insert({
    'reaction_smiles': rxn_smiles,
    'analysis': result['interpretation'],
    'llm_confidence': result['interpretation']['confidence'],
    'validated_score': validation['overall_score'],
    'reliability_rating': validation['reliability'],
    'deterministic_score': validation['individual_scores']['deterministic_quality'],
    'mapping_confidence': validation['detailed_metrics']['deterministic']['mapping_confidence'],
    'recommendation': validation['recommendation']
})
```

---

## Cost Analysis

| Validation Level | Models | Cost/Reaction | Time | Use Case |
|------------------|--------|---------------|------|----------|
| **Fast Check** | 1 (gpt-4o-mini) | ~$0.001 | 5s | High-throughput screening |
| **Standard** | 2 (gpt-4o-mini + gpt-4o) | ~$0.003 | 8s | Routine analysis |
| **Thorough** | 3 (+ o3-mini) | ~$0.006 | 15s | Important results |
| **Maximum** | 5 (+ o3, gpt-5.2) | ~$0.015 | 30s | Novel reactions, publications |

**Recommendation**: Use Fast Check for screening, validate hits with Thorough validation.

---

## Key Takeaways

1. **LLM confidence is NOT enough** - use validation framework
2. **Atom mapping quality = best single predictor** (check this first!)
3. **Deterministic score = most reliable metric** (objective, tool-based)
4. **Cross-model consistency = powerful validation** (multiple models agree?)
5. **Overall validated score > 0.7 = trustworthy** for production use

---

## Files Reference

- **Implementation**: `reaction_agent/scripts/quantitative_validation.py`
- **Detailed Guide**: `reaction_agent/docs/VALIDATION_GUIDE.md`
- **Test Script**: `reaction_agent/scripts/test_validation_comparison.py`
- **Results**: `reaction_agent/results/validation_example.json`

---

## Next Steps

1. ✅ Test on your reactions: Run `quantitative_validation.py`
2. ✅ Review metrics: See `VALIDATION_GUIDE.md` for interpretation
3. ✅ Compare reactions: Use `test_validation_comparison.py`
4. 🔧 Integrate: Add validation to your workflow (CLI, batch, API)
5. 📊 Monitor: Track validation metrics over time to improve model selection

---

**Bottom Line**: Don't trust a single LLM confidence score. Use deterministic quality + cross-model validation for reliable results.
