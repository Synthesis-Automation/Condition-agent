# Buchwald Dataset Analysis

## Dataset Overview

**Source**: `data/reaction_dataset/Buchwald2021-2024.jsonl`  
**Total Reactions**: 1,343  
**Reaction Type**: Buchwald C-N coupling (Pd-catalyzed aryl halide + amine)

## Dataset Statistics

### Yield Distribution
- **Mean**: 75.8%
- **Median**: 77.0%
- **Std Dev**: 13.1%
- **Range**: 33.0% - 100.0%

### Reaction Conditions
- **Temperature**: 100.0°C (all reactions)
- **Time**: 15.9h ± 1.2h (mean ± std)

### Diversity
- **Unique Cores**: 39 different catalyst/ligand combinations
- **Unique Bases**: 23 different bases
- **Unique Solvents**: 26 different solvents

## Train/Val/Test Split

| Split | Count | Percentage |
|-------|-------|------------|
| Train | 940   | 70.0%      |
| Val   | 201   | 15.0%      |
| Test  | 202   | 15.0%      |

**Note**: Stratification by core was disabled due to some rare cores having <3 samples.

## Top 10 Condition Cores

| Core | Count | % |
|------|-------|---|
| Pd (no ligand specified) | 473 | 35.2% |
| Pd/Tri-tert-butylphosphinetetrafluoroborate | 166 | 12.4% |
| Pd/RuPhos | 140 | 10.4% |
| Pd/XPhos | 84 | 6.3% |
| Pd/XantPhos | 81 | 6.0% |
| Pd/PtBu3 | 56 | 4.2% |
| Pd/SPhos | 54 | 4.0% |
| Pd/triphenylphosphine | 42 | 3.1% |
| Pd/DPPP | 35 | 2.6% |
| Pd/IPr | 27 | 2.0% |

**Observations**:
- Buchwald ligands (RuPhos, XPhos, SPhos) are well-represented (278 reactions, 20.7%)
- 35% of reactions don't specify a ligand (just "Pd")
- Good diversity with 39 total cores

## Top 10 Bases (by CAS number)

| CAS | Likely Identity | Count | % |
|-----|-----------------|-------|---|
| 865-48-5 | Sodium tert-butoxide | 491 | 36.6% |
| 534-17-8 | Lithium tert-butoxide | 330 | 24.6% |
| 865-47-4 | Potassium tert-butoxide | 104 | 7.7% |
| 7778-53-2 | Potassium phosphate tribasic | 81 | 6.0% |
| 584-08-7 | Potassium carbonate | 71 | 5.3% |
| Unknown | - | 70 | 5.2% |
| 4039-32-1 | Sodium bis(trimethylsilyl)amide | 42 | 3.1% |
| 121-44-8 | Triethylamine | 36 | 2.7% |
| 1310-73-2 | Sodium hydroxide | 32 | 2.4% |
| 127-09-3 | Sodium acetate | 25 | 1.9% |

**Observations**:
- Strong bases dominate (tert-butoxides: 925 reactions, 68.9%)
- Sodium tert-butoxide is the most common (36.6%)
- Good representation of different base strengths (strong, moderate, weak)

## Top 10 Solvents (by CAS number)

| CAS | Likely Identity | Count | % |
|-----|-----------------|-------|---|
| 108-88-3 | Toluene | 804 | 59.9% |
| 123-91-1 | 1,4-Dioxane | 153 | 11.4% |
| 109-99-9 | Tetrahydrofuran (THF) | 92 | 6.9% |
| 67-68-5 | Dimethyl sulfoxide (DMSO) | 55 | 4.1% |
| 75-05-8 | Acetonitrile | 45 | 3.4% |
| 1634-04-4 | tert-Butyl methyl ether (MTBE) | 35 | 2.6% |
| 98-08-8 | Trifluorotoluene | 30 | 2.2% |
| 68-12-2 | N,N-Dimethylformamide (DMF) | 27 | 2.0% |
| 107-06-2 | 1,2-Dichloroethane | 18 | 1.3% |
| 75-09-2 | Dichloromethane | 11 | 0.8% |

**Observations**:
- Toluene dominates (59.9%) - non-polar, high boiling point
- Good representation of ethereal solvents (dioxane, THF, MTBE: 280 reactions, 20.9%)
- Some polar aprotic solvents (DMSO, DMF, acetonitrile: 127 reactions, 9.5%)

## ML Suitability

### Strengths
✅ **Good size**: 1,343 reactions is sufficient for Phase 1 (DRFP + LightGBM)  
✅ **High yields**: Mean 75.8% suggests good quality data  
✅ **Narrow variance**: Std 13.1% means model will have easier prediction task  
✅ **Diverse cores**: 39 different catalyst/ligand combinations  
✅ **Complete metadata**: All reactions have temperature and time  

### Challenges
⚠️ **Imbalanced cores**: 35% are just "Pd" without ligand specification  
⚠️ **Temperature uniformity**: All at 100°C (no variation to learn from)  
⚠️ **Rare cores**: Some cores have <3 samples (stratification disabled)  
⚠️ **CAS-based identifiers**: Need to map CAS → chemical names for interpretability  

### Recommendations
1. **Train Buchwald-specific model** first (this dataset)
2. **Compare to mixed-family model** later (all 5 datasets combined)
3. **Feature engineering**: Consider grouping similar ligands (e.g., all Buchwald ligands)
4. **Baseline**: k-NN on Buchwald should achieve ~72-76% MAE (similar to mean yield)

## Expected Performance (Phase 1)

Based on dataset characteristics:

| Metric | Expected Range | Rationale |
|--------|----------------|-----------|
| MAE | 8-12% | Narrow yield variance (13.1% std) makes prediction easier |
| RMSE | 10-15% | Same as above |
| R² | 0.70-0.80 | Good signal-to-noise ratio |
| Within ±10% | 60-70% | Most reactions cluster around 75% yield |

**Comparison**: General mixed-family models typically achieve MAE 12-15%. Buchwald-specific should be better due to:
- Single reaction mechanism (Pd-catalyzed C-N coupling)
- Consistent reaction conditions (all 100°C)
- Well-defined substrate scope (aryl halides + amines)

## Files Generated

```
data/
  buchwald_ml_train.jsonl    (940 reactions, 70%)
  buchwald_ml_val.jsonl      (201 reactions, 15%)
  buchwald_ml_test.jsonl     (202 reactions, 15%)
```

## Next Steps

1. **Install dependencies** (if not already done):
   ```bash
   pip install -r requirements-ml.txt
   ```

2. **Train Buchwald-specific model**:
   ```bash
   python scripts/ml/train_drfp_model.py \
       --train data/buchwald_ml_train.jsonl \
       --val data/buchwald_ml_val.jsonl \
       --output models/buchwald_drfp_v1.pkl
   ```

3. **Evaluate**:
   ```bash
   python scripts/ml/evaluate_models.py \
       --model models/buchwald_drfp_v1.pkl \
       --test data/buchwald_ml_test.jsonl \
       --output results/buchwald
   ```

4. **Compare to baseline**: Train k-NN baseline for comparison

5. **Later**: Train mixed-family model on all 5 datasets and compare performance

---

**Dataset**: Buchwald2021-2024.jsonl (1,343 reactions)  
**Prepared**: October 6, 2025  
**Ready for**: Phase 1 ML training (DRFP + LightGBM)
