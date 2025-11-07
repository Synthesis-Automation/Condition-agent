# SNAr Rule Integration Issue - Root Cause Analysis

## Summary

The SNAr rule file (`SNAr_db.json`) and Reductive Amination rule file (`reductive_amination_db.json`) are **structurally valid** and can be loaded successfully, but the **agent cannot use them** because of a **feature token mismatch**.

##  Root Cause

The rule files use feature tokens like:
- `snar_applicable_electrophile_present`
- `aromatic_electrophile_present`
- `alcohol_present`
- `phenol_present`
- `thiol_present`
- `amine_nucleophile_present`
- `carbonyl_present`

But the current `chemtools.featurizers.molecular.featurize()` function **does NOT generate these tokens**. It only generates minimal features like:
- `heteroaryl`
- `electrophile_class`
- `nucleophile_class`

## Test Results

```bash
Electrophile: Clc1nc(Cl)nc(Cl)n1
Nucleophile: OC

Total features: 9

TRUE boolean features:
  heteroaryl

SNAr rule file expects these features:
  aromatic_electrophile_present                : NOT FOUND
  snar_applicable_electrophile_present         : NOT FOUND
  phenol_present                               : NOT FOUND
  alcohol_present                              : NOT FOUND ✗
  ...
```

The featurizer detects `heteroaryl` but the SNAr rule needs `snar_applicable_electrophile_present` and `alcohol_present` (for "OC").

## Solution Options

### Option 1: Update the Featurizer (Recommended)

**Extend `chemtools/featurizers/molecular.py` to generate the expected tokens.**

The functional group detection system already exists in `chemtools/util/functional_groups.py` with patterns for:
- alcohol
- phenol
- thiol
- primary_amine
- secondary_amine
- aniline
- carbonyl
- aldehyde
- ketone

These just need to be called and converted to `<name>_present` boolean tokens.

**Changes needed:**
1. Import `has_functional_group` from `chemtools.util.functional_groups`
2. In the `featurize()` function, add checks for all relevant functional groups
3. Generate `<group>_present` tokens
4. Add SNAr-specific logic for `snar_applicable_electrophile_present` and `aromatic_electrophile_present`

**Example addition to featurizer:**
```python
from chemtools.util.functional_groups import has_functional_group

def featurize(electrophile: str, nucleophile: str) -> Dict[str, Any]:
    # ... existing code ...
    
    # Add functional group detection
    nuc_groups = ["alcohol", "phenol", "thiol", "primary_amine", 
                  "secondary_amine", "aniline"]
    for group in nuc_groups:
        features[f"{group}_present"] = has_functional_group(nucleophile, group)
    
    elec_groups = ["carbonyl", "aldehyde", "ketone", "aryl_halide"]
    for group in elec_groups:
        features[f"{group}_present"] = has_functional_group(electrophile, group)
    
    # SNAr-specific detection
    features["aromatic_electrophile_present"] = is_aromatic(electrophile)
    features["snar_applicable_electrophile_present"] = (
        features.get("heteroaryl", False) or 
        has_strong_ewg(electrophile)  # Check for NO2, CN, CF3, etc.
    )
    
    return features
```

### Option 2: Update the Rule Files

**Rewrite the SNAr and Reductive Amination rule files to use only the tokens that the current featurizer generates.**

This is **not recommended** because:
1. The current featurizer is too minimal for sophisticated rule matching
2. The rule files were carefully crafted with chemistry-specific conditions
3. Other rule files (Suzuki, Buchwald-Hartwig) likely use similar tokens

## What Was Fixed

✅ Added SNAr and Reductive Amination mappings to `_FAMILY_TO_RULE_DB` in `chemtools_wrapper.py`:
```python
"snar": "SNAr_db",
"s_nar": "SNAr_db",
"aromatic_nucleophilic_substitution": "SNAr_db",
"nucleophilic_aromatic_substitution": "SNAr_db",
"reductive_amination": "reductive_amination_db",
```

✅ Added family aliases to `chemtools/analysis/reactions.py`:
```python
"S_NAr": "snar",
"SNAr": "snar",
"Aromatic_Nucleophilic_Substitution": "snar",
"Reductive_Amination": "reductive_amination",
```

✅ Verified rule file structure is valid and loads correctly

## What Still Needs Fixing

❌ **Feature detection must be expanded** to generate the tokens expected by the rule files

The minimal change needed in `chemtools/featurizers/molecular.py`:
1. Import functional group detection utilities
2. Add loops to check for all relevant functional groups
3. Generate `<group>_present` boolean flags
4. Add SNAr-specific detection logic

## Testing

Run this to verify the mappings work:
```bash
python test_snar_rule.py
```

Output shows:
- ✅ Rule files load correctly
- ✅ Family-to-database mapping works
- ❌ Feature detection doesn't generate expected tokens
- ❌ `applies_if` conditions fail

## Recommendation

**Implement Option 1 (Extend Featurizer)** because:
1. It benefits all rule-based recommendations, not just SNAr
2. The functional group library is already comprehensive
3. Other existing rule files likely depend on these same tokens
4. It aligns with the documented repository design

Would you like me to implement the featurizer extension?
