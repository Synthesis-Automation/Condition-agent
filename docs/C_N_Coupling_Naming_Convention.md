# C-N Coupling Dataset Naming Convention

## Overview
As of October 6, 2025, we have standardized the naming convention for all C-N coupling reaction datasets to clearly indicate the metal catalyst and reaction type.

## New Naming Convention

### Dataset Files (`data/reaction_dataset/`)
| Old Name                    | New Name                         | Catalyst | Reactions |
|----------------------------|----------------------------------|----------|-----------|
| `Ullman2020-2024.jsonl`    | `C_N_coupling_Cu_Ullmann.jsonl`  | Cu       | 5,552     |
| `Buchwald2021-2024.jsonl`  | `C_N_coupling_Pd_Buchwald.jsonl` | Pd       | 1,343     |
| `Amination-Ni2014-2024.jsonl` | `C_N_coupling_Ni.jsonl`       | Ni       | 1,131     |

### Model Files (`models/`)
| Old Name                   | New Name                          | Description                 |
|---------------------------|-----------------------------------|-----------------------------|
| `ullmann_drfp_yield_v1.pkl` | `cn_coupling_cu_ullmann_v1.pkl` | Cu-catalyzed yield predictor |
| `buchwald_drfp_v1.pkl`     | `cn_coupling_pd_buchwald_v1.pkl` | Pd-catalyzed yield predictor |
| `ullmann_training_report.txt` | `cn_coupling_cu_ullmann_training_report.txt` | Cu model training report |

### Results Directories (`results/`)
| Old Name          | New Name                        | Contents                      |
|------------------|---------------------------------|-------------------------------|
| `ullmann_test/`  | `cn_coupling_cu_ullmann_test/`  | Cu model test results         |
| `cn_coupling_test/` | `cn_coupling_pd_buchwald_test/` | Pd model test results (rename pending) |

## Rationale

### Why This Convention?
1. **Clarity**: All C-N coupling reactions are grouped together with prefix `C_N_coupling_`
2. **Metal Identification**: Immediately shows which catalyst (Cu, Pd, Ni)
3. **Named Reaction Context**: Preserves historical context (Ullmann, Buchwald)
4. **Consistency**: Uniform naming across datasets, models, and results

### Benefits
- **Easier Discovery**: Searching for "C_N_coupling" finds all related files
- **Less Confusion**: Clear distinction between Cu, Pd, and Ni catalysis
- **Scalability**: Easy to add new C-N coupling methods (e.g., C_N_coupling_Fe)
- **Organization**: Groups related reactions logically in file browsers

## Updated Files

### Scripts Updated
- ✅ `scripts/train_ullmann_drfp.py` → Uses `C_N_coupling_Cu_Ullmann.jsonl` and `cn_coupling_cu_ullmann_v1.pkl`
- ✅ `scripts/test_ullmann_reactions.py` → Uses `cn_coupling_cu_ullmann_v1.pkl`, outputs to `cn_coupling_cu_ullmann_test/`
- ✅ `scripts/test_cn_coupling_reactions.py` → Uses `cn_coupling_pd_buchwald_v1.pkl`
- ✅ `scripts/verify_ml_with_rules.py` → Uses `C_N_coupling_Pd_Buchwald.jsonl` and `cn_coupling_pd_buchwald_v1.pkl`
- ✅ `scripts/ml/prepare_buchwald_dataset.py` → Uses `C_N_coupling_Pd_Buchwald.jsonl`, outputs `cn_coupling_pd_buchwald_v1.pkl`

### Files Not Updated (Legacy)
- `chemtools/featurizers/ullmann.py` - Core library, name preserved for API stability
- `tests/test_featurize_ullmann.py` - Test name preserved
- `scripts/ullmann_tester.py` - Legacy tester script
- `scripts/validate_condition_core_ullmann.py` - Validator script

## Comparison: Cu vs Pd vs Ni

| **Property**           | **Ullmann (Cu)**                 | **Buchwald (Pd)**             | **Ni Amination**              |
|------------------------|----------------------------------|-------------------------------|-------------------------------|
| **Dataset**            | C_N_coupling_Cu_Ullmann.jsonl    | C_N_coupling_Pd_Buchwald.jsonl | C_N_coupling_Ni.jsonl        |
| **Reactions**          | 5,552                            | 1,343                         | 1,131                         |
| **Catalysts**          | 27 Cu types                      | 37 Pd/ligand complexes        | TBD                           |
| **Avg Yield**          | 74.2%                            | ~73%                          | TBD                           |
| **Temp Range**         | 90-120°C                         | 80-100°C                      | TBD                           |
| **Catalyst Cost**      | Low (Cu salts)                   | High (Pd + ligands)           | Medium (Ni salts)             |
| **Ligands**            | Simple (phen, proline) or none   | Complex (XPhos, RuPhos, etc.) | TBD                           |
| **Model MAE**          | 9.61% (test)                     | 11.42% (test)                 | Not yet trained               |
| **Model Status**       | ✅ Trained & Tested              | ✅ Trained & Tested           | ⏳ Analysis pending           |

## Next Steps

### Pending Tasks
1. ⏳ **Reorganize results directories** - Rename `cn_coupling_test/` → `cn_coupling_pd_buchwald_test/`
2. ⏳ **Analyze Ni dataset** - Load `C_N_coupling_Ni.jsonl`, extract vocabulary, compare with Cu/Pd
3. 📋 **Train Ni model** - Create `cn_coupling_ni_v1.pkl`
4. 📋 **Compare all three** - Generate Cu vs Pd vs Ni comparison report

### Future Additions
- C_N_coupling_Fe (Iron-catalyzed C-N coupling, if dataset becomes available)
- C_N_coupling_Photocatalytic (Photoredox C-N coupling)
- C_N_coupling_Electrochemical (Electrochemical C-N coupling)

## Usage Examples

### Training a Model
```bash
# Train Cu model
python scripts/train_ullmann_drfp.py

# Train Pd model (requires prepare_buchwald_dataset.py first)
python scripts/ml/prepare_buchwald_dataset.py
python scripts/ml/train_drfp_model.py --train data/buchwald_ml_train.jsonl --val data/buchwald_ml_val.jsonl --output models/cn_coupling_pd_buchwald_v1.pkl
```

### Testing Models
```bash
# Test Cu model
python scripts/test_ullmann_reactions.py

# Test Pd model
python scripts/test_cn_coupling_reactions.py
```

### Verification
```bash
# Verify Pd model predictions vs precedents
python scripts/verify_ml_with_rules.py
```

## Migration Notes

### For Developers
- Update any scripts/notebooks that reference old filenames
- Check imports in `chemtools/` - core library names unchanged
- Results directories will be reorganized (see task 4)

### For Users
- Old model files are renamed, retrain or update paths
- Dataset files are renamed, update any custom scripts
- Documentation updated to reflect new convention

---

**Last Updated:** October 6, 2025  
**Status:** Phase 1 Complete (Files renamed, scripts updated)  
**Next:** Ni dataset analysis

