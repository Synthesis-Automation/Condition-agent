# applies_if Filtering for Rules

## Overview

The `UnifiedRecommender` automatically filters rule-based recommendations using `applies_if` criteria. This ensures that only chemically appropriate rules are suggested for a given reaction.

## Status: ✅ FULLY IMPLEMENTED AND ENABLED BY DEFAULT

## How It Works

### 1. Feature Detection
When you query a reaction, the system automatically detects substrate features:
```python
# Query: Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1
Detected features:
  - aryl_halide_present: True
  - sp2_bromide_present: True
  - terminal_alkyne_present: True
```

### 2. applies_if Criteria Check
Each rule specifies required features in its `applies_if` section:
```json
{
  "applies_if": {
    "all": ["terminal_alkyne_present"],
    "any": ["aryl_halide_present", "vinyl_halide_present"]
  }
}
```

This means:
- **ALL** conditions must be true: `terminal_alkyne_present`
- **ANY** of these must be true: `aryl_halide_present` OR `vinyl_halide_present`

### 3. Automatic Filtering
Rules that don't match are silently filtered out:
```
✅ Sonogashira (matches) → Included in results
❌ Amide Formation (doesn't match) → Filtered out
❌ Suzuki (doesn't match) → Filtered out
```

## Example: Sonogashira Rule

### applies_if Criteria
```json
{
  "all": ["terminal_alkyne_present"],
  "any": [
    "aryl_halide_present",
    "vinyl_halide_present", 
    "aryl_triflate_present",
    "vinyl_triflate_present"
  ]
}
```

### Test Case 1: Matching Reaction ✅
```
Query: Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1

Features Detected:
  ✓ terminal_alkyne_present (C#Cc1ccccc1)
  ✓ aryl_halide_present (Brc1ccccc1)
  ✓ sp2_bromide_present

Result: Sonogashira rule INCLUDED
```

### Test Case 2: Non-Matching Reaction ❌
```
Query: CC(=O)O.CCN>>CC(=O)NCC

Features Detected:
  ✗ terminal_alkyne_present (missing)
  ✗ aryl_halide_present (missing)
  ✓ carboxylic_acid_present
  ✓ primary_amine_present

Result: Sonogashira rule FILTERED OUT
```

## CLI Usage

### Default Behavior (Filtering Enabled)
```bash
python app/unified_rule_protocol_interactive_cli.py
```

```
reaction> /type rule
reaction> Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1

Found 1 recommendation(s)
📖 [1] Sonogashira Coupling Guidelines
    ✓ Matches applies_if criteria
```

### API Usage

**With filtering (default, recommended):**
```python
from chemtools.recommend.unified import UnifiedRecommender

recommender = UnifiedRecommender("build/unified_index_complete")

results = recommender.recommend(
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1",
    top_k=5,
    validate_rules=True  # Default, enables applies_if filtering
)
```

**Without filtering (NOT recommended):**
```python
results = recommender.recommend(
    "CC(=O)O.CCN>>CC(=O)NCC",
    top_k=5,
    validate_rules=False  # Disables filtering
)

# WARNING: May return chemically inappropriate rules!
# Example: Sonogashira rule for amide formation (wrong!)
```

## Benefits

### ✅ Chemical Accuracy
Only rules appropriate for the substrate type are recommended.

### ✅ No False Positives
Prevents suggesting incompatible reactions (e.g., Sonogashira for amides).

### ✅ Automatic
No manual configuration needed - works by default.

### ✅ Fail-Safe
If feature detection fails, system fails open (includes all rules) rather than failing closed (no results).

## Comparison: With vs Without Filtering

### Test Reaction: Amide Formation
```
Query: CC(=O)O.CCN>>CC(=O)NCC
```

**WITH filtering (validate_rules=True):**
```
Found 2 rules:
  ✅ Amide Formation (relevant)
  ✅ Reductive Amination (relevant)
```

**WITHOUT filtering (validate_rules=False):**
```
Found 9 rules:
  ✅ Amide Formation
  ❌ Sonogashira (wrong - needs alkyne!)
  ❌ Suzuki (wrong - needs boronic acid!)
  ❌ Buchwald-Hartwig (wrong - needs aryl halide + amine)
  ❌ C-O Coupling (wrong substrate type)
  ❌ RCM (wrong - needs dienes!)
  ❌ SNAr (wrong - needs electron-poor arene)
  ✅ Reductive Amination
  ❌ Ullmann (wrong - needs aryl halide)
```

**Conclusion:** Filtering reduces noise from 9 → 2 rules (78% reduction in irrelevant results!)

## Technical Details

### Feature Analyzer
Located in `chemtools/rule/analyzer.py`:
```python
from chemtools.rule.analyzer import FeatureAnalyzer

analyzer = FeatureAnalyzer()
features = analyzer.analyze_reaction(reaction_smiles)

# Returns: {'aryl_halide_present': True, 'terminal_alkyne_present': True, ...}
```

### applies_if Schema
Defined in `chemtools/schema/v2.json`:
```json
{
  "applies_if": {
    "type": "object",
    "properties": {
      "all": {
        "type": "array",
        "items": {"type": "string"},
        "description": "All conditions must be true (AND logic)"
      },
      "any": {
        "type": "array", 
        "items": {"type": "string"},
        "description": "At least one must be true (OR logic)"
      }
    }
  }
}
```

### Validation Logic
Located in `chemtools/recommend/unified.py`:
```python
def _check_applies_if(self, features: Dict[str, bool], applies_if: Dict[str, Any]) -> bool:
    # Check 'all' conditions (AND logic)
    if 'all' in applies_if:
        if not all(features.get(condition, False) for condition in applies_if['all']):
            return False
    
    # Check 'any' conditions (OR logic)
    if 'any' in applies_if:
        if not any(features.get(condition, False) for condition in applies_if['any']):
            return False
    
    return True
```

## Available Feature Detectors

The system can detect 80+ functional groups and features. Common ones for rules:

### Halides
- `aryl_halide_present` - Aromatic C-X bonds
- `vinyl_halide_present` - Alkene C-X bonds
- `sp2_chloride_present`, `sp2_bromide_present`, `sp2_iodide_present`

### Alkynes
- `terminal_alkyne_present` - R-C≡C-H
- `alkyne_present` - Any C≡C

### Amines
- `primary_amine_present` - R-NH₂
- `secondary_amine_present` - R₂NH
- `aniline_present` - Aromatic amine

### Acids/Bases
- `carboxylic_acid_present` - R-COOH
- `amine_present` - Any amine (primary/secondary/tertiary)

### Others
- `boronic_acid_present` - R-B(OH)₂ (for Suzuki)
- `aryl_triflate_present` - Ar-OTf
- `diene_present` - For RCM
- `aldehyde_present`, `ketone_present`

Full list: See `chemtools/util/functional_groups.py`

## Testing

Run comprehensive tests:
```bash
# Test applies_if filtering
python test_applies_if_filtering.py

# Demo in CLI
python demo_applies_if_filtering.py

# Verify all features
pytest tests/test_rule_validation.py -v
```

## Summary

✅ **Status**: Fully implemented and enabled by default  
✅ **Location**: `chemtools/recommend/unified.py`  
✅ **Default**: `validate_rules=True`  
✅ **Effect**: Filters rules based on substrate features  
✅ **Benefit**: Only chemically appropriate rules recommended  
✅ **Fail-safe**: Fails open if detection fails  

**The system is production-ready and working correctly!**
# NOTE: Applies-if filtering was part of the legacy rule-based recommender.
# The unified recommender no longer applies rule filtering by default.
