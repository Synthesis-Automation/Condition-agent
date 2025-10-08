# Potential Problems in C-N Coupling Test Results

**Analysis Date:** October 6, 2025  
**Model:** Buchwald DRFP v1  
**Test Set:** 62 C-N coupling reactions from sample_reactions.py

---

## 🚨 Critical Issues

### 1. **Training-Test Domain Mismatch**

**Problem:** The model was trained exclusively on **Buchwald-Hartwig reactions** (Pd-catalyzed C-N coupling), but we tested it on:

- **Buchwald-Hartwig:** 6 reactions ✅ (appropriate)
- **Generic "C-N Coupling":** 53 reactions ⚠️ (may or may not be Buchwald-Hartwig)
- **Chan-Lam Coupling:** 3 reactions ❌ (completely different mechanism)

**Why This Matters:**
- **Chan-Lam** uses Cu-mediated oxidative coupling of **boronic acids** with amines
- **Buchwald-Hartwig** uses Pd-catalyzed coupling of **aryl halides** with amines
- Different catalysts, mechanisms, and optimal conditions
- The model has NEVER seen Chan-Lam reactions during training

**Evidence:**
```
Chan-Lam reactions:
- c1ccccc1B(O)O.Nc1ccccc1.[O]>>c1ccccc1Nc1ccccc1
- Uses boronic acid, NOT halide
- Model predicted: 71.9% with Pd/SPhos (should use Cu, not Pd!)
```

**Impact:** Chan-Lam predictions are likely **invalid**. The model is hallucinating Pd-catalyzed conditions for a Cu-catalyzed reaction.

---

### 2. **Overly Optimistic Predictions**

**Problem:** ALL reactions predicted ≥65% yield, with NO failures predicted.

**Statistical Red Flags:**
- Min yield: 65.3% (4-CN-Ph-Br + aniline)
- Max yield: 90.0%
- Mean: 74.9%
- **Zero reactions <60%** ← Unrealistic!

**Reality Check:**
In real Buchwald-Hartwig chemistry:
- **Electron-poor substrates** (4-CN, 4-NO₂, 4-CHO) often give **lower yields** (30-60%)
- **Sterically hindered substrates** (ortho-substituted) can fail completely (<10%)
- **Heteroaryl chlorides** without activation are challenging (20-50%)

**Suspicious Examples:**
```
1. 4-Cyano-bromobenzene + aniline → 65.3%
   - Reality: Cyano is strongly electron-withdrawing, competes for Pd coordination
   - Expected: 30-50% (difficult substrate)

2. 4-Nitro-bromobenzene + aniline → 80.1%
   - Reality: Nitro group is VERY electron-withdrawing, often <40% yield
   - Model says: 80%! (Likely overoptimistic)

3. 2-Methylbromobenzene + aniline → 77.6%
   - Reality: Ortho-substitution causes steric hindrance, often 40-60%
   - Model says: 77.6% (Possibly overoptimistic)
```

**Hypothesis:** The model was trained on **published reactions** (Buchwald2021-2024.jsonl), which have **publication bias** toward successful reactions. Failed reactions rarely get published, so the model never learned what a "bad" reaction looks like.

---

### 3. **Missing Reaction SMILES Validation**

**Problem:** Some reaction SMILES are **incomplete or ambiguous**.

**Examples:**
```
1. "Brc1ccc(C(F)(F)F)cc1.NC1CCCCC1>>FC"
   - Product SMILES is truncated: "FC" instead of full structure
   - Should be: "FC(F)(F)c1ccc(NC2CCCCC2)cc1"

2. "Ic1ccc(C=O)cc1.Nc1ccccc1>>O=Cc1ccc"
   - Product truncated: "O=Cc1ccc" 
   - Should be: "O=Cc1ccc(Nc2ccccc2)cc1"
```

**Impact:** 
- DRFP fingerprints are generated from **incomplete products**
- Model may be making predictions based on **partial structures**
- Yields could be inaccurate due to wrong molecular representation

**Action Needed:** Validate ALL reaction SMILES before trusting predictions.

---

### 4. **Temperature and Time Fixed at Arbitrary Values**

**Problem:** All predictions use:
- **Temperature:** 100°C (fixed)
- **Time:** 12 hours (fixed)

**Reality:**
- Buchwald-Hartwig reactions range from **RT to 150°C**
- Time varies from **1-48 hours**
- Optimal conditions are substrate-dependent

**Missing Optimization:**
- Model CAN predict yields at different T/time
- Script doesn't explore this parameter space
- We're only testing 3 condition sets, missing potentially better combinations

**Example:**
```
Reaction: Ph-Br + aniline
Tested: 100°C, 12h → 74.5%
Not tested: 80°C, 24h (might be better)
Not tested: 120°C, 4h (might be faster)
```

---

### 5. **No Uncertainty Quantification**

**Problem:** Model gives **point predictions** without confidence intervals.

**Why This Matters:**
- A prediction of "74.5%" could actually be 74.5% ± 20%
- We don't know which predictions are reliable vs uncertain
- No way to prioritize experimental validation

**Test Set Performance:**
- Training: MAE 5.94%, Val MAE 9.35%
- But these are **on Buchwald reactions**
- Out-of-domain (Chan-Lam) uncertainty is UNKNOWN

---

### 6. **Condition Recommendation Logic Incomplete**

**Problem:** The script tests only 3 pre-defined condition sets:
1. Pd/XPhos, NaOtBu, Toluene
2. Pd/RuPhos, K3PO4, Toluene
3. Pd/SPhos, NaOtBu, Dioxane

**Missing:**
- No exploration of the **37 cores** in the vocabulary
- No testing of **20 bases** and **25 solvents**
- Potential optimal conditions may exist but weren't tested

**Better Approach:**
- Test all cores in vocabulary (at least top 10)
- Grid search over common bases/solvents
- Or: Use chemtools precedent-based recommendation system

---

### 7. **Substrate Scope Limitations**

**Problem:** The Buchwald dataset (2021-2024) may have **limited substrate diversity**.

**Potential Gaps:**
- **Aliphatic halides:** Model trained on aryl halides
- **F-leaving group:** Very rare in Buchwald chemistry
- **Triflates:** Less common than Br/Cl/I
- **Highly hindered substrates:** May be underrepresented

**Evidence:**
```
Test reactions with potential issues:
- tert-Butylamine (very bulky) → 70.9%
- Diethylamine (secondary) → 77.1%
- 2,6-Dimethyl substrates → 72-78%

All predicted moderately high, but this could be model extrapolation.
```

---

## ⚠️ Moderate Issues

### 8. **Duplicate Reaction Testing**

Several reactions appear multiple times with identical or very similar SMILES:

```
1. "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
   - Appears in both "Buchwald-Hartwig" and "C-N Coupling" sections
   - Tested twice with same conditions

2. Similar reactions:
   - Ph-Br + aniline (appears 3+ times)
   - 4-Br-pyridine + aniline (appears 2+ times)
```

**Impact:** Inflates success metrics, wastes computation.

---

### 9. **CAS Number Lookup Not Validated**

**Problem:** Base/solvent UIDs are assumed to be correct CAS numbers.

**Example:**
```python
base_uid: '1907-33-1'  # Assumed to be NaOtBu
solvent_uid: '108-88-3'  # Assumed to be Toluene
```

**Risk:** If CAS numbers are wrong, predictions are meaningless.

**Action:** Cross-reference with actual Buchwald training data.

---

### 10. **Missing Reaction Type Detection**

**Problem:** We manually classified reactions as "Buchwald-Hartwig" vs "C-N Coupling" based on string matching.

**Better Approach:**
- Use `chemtools.condition_core.detect_cn_coupling_type()`
- Automatically classify based on reactants (halide vs boronic acid)
- Flag reactions that don't match training distribution

---

## ✅ Things That Worked Well

1. **Model Reconstruction:** Successfully loaded and used the pickled model
2. **DRFP Fingerprinting:** All 62 reactions successfully featurized
3. **Condition Testing:** Systematic comparison of 3 condition sets
4. **Report Generation:** Comprehensive markdown report created
5. **100% Success Rate:** No prediction failures (though see Issue #2 about overoptimism)

---

## 🔧 Recommended Fixes

### Immediate Actions

1. **Remove Chan-Lam Reactions:**
   - They're out-of-domain for a Buchwald model
   - Or: Clearly label as "Exploratory/Unreliable"

2. **Validate Reaction SMILES:**
   ```python
   from chemtools.smiles import _split_reaction_smiles
   from chemtools.util.rdkit_helpers import parse_smiles
   
   # For each reaction:
   parts = _split_reaction_smiles(rxn_smiles)
   for smi in parts:
       mol = parse_smiles(smi)
       assert mol is not None, f"Invalid SMILES: {smi}"
   ```

3. **Add Warning Labels to Report:**
   ```markdown
   ## ⚠️ IMPORTANT DISCLAIMERS
   
   1. Chan-Lam predictions are UNRELIABLE (wrong catalyst type)
   2. All yields are optimistic (no failed reactions predicted)
   3. Some reaction SMILES are truncated and need validation
   4. Predictions are for 100°C, 12h only (not optimized)
   ```

4. **Flag Suspicious Predictions:**
   ```python
   # Mark as "HIGH UNCERTAINTY" if:
   - Electron-poor substrates (CN, NO2, CHO) with yield >70%
   - Ortho-substituted with yield >75%
   - Heteroaryl chlorides with yield >65%
   - Any Chan-Lam reaction
   ```

### Medium-Term Improvements

5. **Implement Uncertainty Quantification:**
   - Use LightGBM's quantile regression
   - Or: Train ensemble of models, report std dev
   - Or: Calibrate predictions with test set performance

6. **Expand Condition Search:**
   ```python
   # Test all cores in vocabulary:
   for core in predictor.core_vocab:
       for base in ['1907-33-1', '7778-53-2', '1310-73-2']:  # Top bases
           for solvent in ['108-88-3', '123-91-1', '109-99-9']:  # Top solvents
               predict_yield(rxn, core, base, solvent)
   ```

7. **Temperature/Time Optimization:**
   ```python
   # Grid search:
   for T in [80, 100, 120]:
       for time in [4, 12, 24]:
           predict_yield(rxn, core, base, solvent, T, time)
   ```

### Long-Term Enhancements

8. **Train Multi-Family Model:**
   - Combine Buchwald, Suzuki, Ullmann, Amide, Amination datasets
   - Use reaction type as a feature
   - Proper handling of different mechanisms

9. **Active Learning Loop:**
   - Run top predictions experimentally
   - Collect actual yields
   - Retrain model with new data
   - Iteratively improve

10. **Integrate with Chemtools Recommendation:**
    ```python
    from chemtools.recommend import get_recommendations
    
    # Get precedent-based recommendations
    recs = get_recommendations(reactants, products, n=5)
    
    # Compare with ML predictions
    # Use as ensemble: 50% ML + 50% precedent
    ```

---

## 📊 Summary Table

| Issue | Severity | Impact | Fix Difficulty | Priority |
|-------|----------|--------|---------------|----------|
| Chan-Lam out-of-domain | 🔴 Critical | Predictions invalid | Easy | HIGH |
| Overly optimistic yields | 🔴 Critical | Unreliable guidance | Hard | HIGH |
| Truncated SMILES | 🟠 High | Wrong predictions | Medium | HIGH |
| Fixed T/time | 🟡 Medium | Suboptimal conditions | Easy | MEDIUM |
| No uncertainty | 🟡 Medium | Can't prioritize expts | Hard | MEDIUM |
| Limited condition search | 🟡 Medium | Missing better options | Easy | MEDIUM |
| Duplicate reactions | 🟢 Low | Inflated metrics | Easy | LOW |
| Substrate scope gaps | 🟢 Low | Unknown extrapolation | Hard | LOW |

---

## 🎯 Recommended Next Steps

1. **Create Filtered Report:**
   - Remove Chan-Lam reactions
   - Add uncertainty warnings for electron-poor/hindered substrates
   - Fix truncated SMILES

2. **Experimental Validation Plan:**
   - Pick top 5 predictions (>85% yield)
   - Pick bottom 5 predictions (65-70% yield)
   - Run experiments to calibrate model

3. **Expand Testing:**
   - Test temperature range (80-120°C)
   - Test all cores in vocabulary
   - Compare with chemtools precedent system

4. **Document Limitations:**
   - Add "Model Card" describing training data
   - Specify applicability domain
   - Warn about out-of-distribution predictions

Would you like me to generate a corrected report with these fixes applied?
