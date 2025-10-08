# Suzuki Reaction Test Guide

## Overview

`test_suzuki_recommendations.py` compares three approaches for recommending Suzuki coupling reaction conditions:

1. **Rule-Based**: SMARTS pattern matching via router
2. **ML (k-NN)**: Precedent-based recommendations
3. **Fusion**: Adaptive multi-evidence combination

## Quick Start

### 1. Start the API Server

```powershell
uvicorn app.main:app --reload --port 8000
```

### 2. Run the Test

```powershell
python test_suzuki_recommendations.py
```

## What It Tests

The script tests **5 Suzuki coupling reactions**:

1. **Simple Ph-Ph Coupling**
   - Phenyl bromide + phenylboronic acid → biphenyl
   - `Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1`

2. **Electron-Poor ArCl**
   - 4-chlorobenzonitrile + phenylboronic acid
   - Tests less reactive aryl chloride

3. **Heteroaryl Coupling**
   - 2-bromopyridine + phenylboronic acid
   - Tests nitrogen-containing heterocycle

4. **Ortho-Substituted**
   - 2-bromotoluene + phenylboronic acid
   - Tests steric hindrance

5. **Electron-Rich ArBr**
   - 4-bromoanisole + phenylboronic acid
   - Tests electron-rich aryl halide

## For Each Reaction, It Shows:

### [1] Rule-Based Detection
- Detected reaction family
- Confidence score
- Detection method (SMARTS match)

### [2] ML k-NN Recommendations
- Top 3 recommended conditions
- Processing time
- Dataset information

### [3] Fusion Recommendations
- **Adaptive weights** (α, β, γ, δ) for each evidence source
- **Evidence quality metrics**:
  - Precedent count
  - Diversity score (batch effect detection)
  - Dataset size
- **Weight reasoning**: Why weights were chosen
- Top 3 recommended conditions with component scores
- Processing time

### [4] Comparison
- Side-by-side comparison of ML vs Fusion top picks
- Agreement/disagreement analysis

## Expected Output

```
================================================================================
  Test 1/5: Simple Ph-Ph Coupling
================================================================================

Reaction: Phenyl bromide + phenylboronic acid → biphenyl
SMILES: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1

[1] RULE-BASED: Family Detection
--------------------------------------------------------------------------------
  Family: C_C_Coupling_Pd
  Confidence: 0.95
  Detection Method: SMARTS_match
  ✅ Correctly identified as Suzuki/C-C coupling

[2] ML-BASED: k-NN Precedent Search
--------------------------------------------------------------------------------
  Generated 5 recommendations in 0.45s

  Top 3 Recommendations:
    1. Pd + K3PO4 in toluene (confidence: high)
    2. Pd + Cs2CO3 in dioxane (confidence: medium)
    3. Pd + K2CO3 in DMF (confidence: medium)

  Metadata:
    Model: precedent_knn
    Dataset: C_C_Coupling_Pd

[3] FUSION: Adaptive Multi-Evidence
--------------------------------------------------------------------------------
  Generated 5 recommendations in 0.52s

  Adaptive Weights:
    Precedent (α): 36.4%
    Analytics (β): 50.6%
    Rules (γ):     13.0%
    ML (δ):        0.0%
    Sum: 1.000 ✅

  Evidence Quality:
    Precedents found: 15
    Diversity score: 0.142 (LOW - potential batch effect)
    Dataset size: 8234

  Weight Reasoning:
    • Low diversity (0.14) → precedents may be biased
    • Medium similarity (0.68) → precedents moderately relevant
    • No strong rule match → low rule weight
    • ML predictor not integrated → δ=0

  Top 3 Recommendations:
    1. Pd + K3PO4 in toluene (confidence: high)
       Scores: PS=0.85, AS=0.72, RS=0.45, MS=0.00
    2. Pd + Cs2CO3 in dioxane (confidence: medium)
       Scores: PS=0.78, AS=0.68, RS=0.42, MS=0.00
    3. Pd + K2CO3 in THF (confidence: medium)
       Scores: PS=0.75, AS=0.65, RS=0.40, MS=0.00

[4] COMPARISON
--------------------------------------------------------------------------------
  ML (k-NN) Top:    Pd + K3PO4 in toluene (confidence: high)
  Fusion Top:       Pd + K3PO4 in toluene (confidence: high)
  ✅ Same top recommendation

[... 4 more reactions ...]

================================================================================
  TEST SUMMARY
================================================================================

Results: 5/5 tests completed

  ✅ PASS: Simple Ph-Ph Coupling
  ✅ PASS: Electron-Poor ArCl
  ✅ PASS: Heteroaryl Coupling
  ✅ PASS: Ortho-Substituted
  ✅ PASS: Electron-Rich ArBr
```

## Method Comparison Table

The script prints a detailed comparison table:

```
┌─────────────────────┬───────────────────┬─────────────────┬──────────────────┐
│ Feature             │ Rule-Based        │ ML (k-NN)       │ Fusion           │
├─────────────────────┼───────────────────┼─────────────────┼──────────────────┤
│ Approach            │ SMARTS patterns   │ Precedent search│ Multi-evidence   │
│ Speed               │ Very Fast         │ Fast            │ Medium           │
│ Adaptability        │ Low               │ Medium          │ High             │
│ Batch Effect Aware  │ No                │ No              │ Yes ✓            │
│ Explainability      │ High (rules)      │ Medium (sim.)   │ High (reasoning) │
│ Best For            │ Known patterns    │ Common reactions│ All cases        │
└─────────────────────┴───────────────────┴─────────────────┴──────────────────┘
```

## Key Insights

### Why Fusion Is Better:

1. **Batch Effect Detection**: Low diversity score triggers weight adjustment
2. **Multi-Source Evidence**: Combines precedent, analytics, rules, and ML
3. **Adaptive Weighting**: Adjusts based on data quality
4. **Explainability**: Provides reasoning for recommendations
5. **Component Scores**: Shows contribution from each evidence source

### When to Use Each Method:

- **Rule-Based**: Quick classification, known reaction types
- **ML (k-NN)**: Fast recommendations for common reactions
- **Fusion**: Novel reactions, quality-critical applications, need explanations

## Customization

### Test Your Own Suzuki Reaction

Edit the `SUZUKI_REACTIONS` list in the script:

```python
SUZUKI_REACTIONS = [
    {
        "name": "My Custom Suzuki",
        "reaction": "YOUR_SMILES_HERE",
        "description": "Your reaction description"
    }
]
```

### Adjust Parameters

```python
# Change k value for precedent search
baseline_result = get_baseline_recommendations(reaction, k=100)

# Change max_variants for fusion
fusion_result = get_fusion_recommendations(reaction, k=50, max_variants=10)
```

## Troubleshooting

### Server Not Running
```
❌ API server not accessible: Connection refused

Please start the server:
  uvicorn app.main:app --reload --port 8000
```

**Solution**: Start the API server first.

### Timeout Errors

If reactions take too long:

```python
# Increase timeout in the script
TIMEOUT = 120  # seconds
```

### No Recommendations Generated

Check if:
1. Reaction SMILES is valid
2. Reaction family is supported
3. Dataset contains relevant precedents

## Output Files

The script only prints to console. To save results:

```powershell
# Save to file
python test_suzuki_recommendations.py > suzuki_test_results.txt

# Save and display
python test_suzuki_recommendations.py | Tee-Object -FilePath suzuki_test_results.txt
```

## Next Steps

- **Compare Results**: Analyze when fusion differs from ML k-NN
- **Tune Weights**: Understand how diversity affects weight distribution
- **Add More Reactions**: Test edge cases and novel structures
- **Integrate with Workflows**: Use insights to choose best method for your application

## Related Documentation

- `API_GUIDE.md`: All API endpoint details
- `FUSION_API_COMPLETE.md`: Deep dive into fusion system
- `test_fusion_api_simple.py`: General fusion API tests
