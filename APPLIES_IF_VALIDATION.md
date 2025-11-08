# Post-Similarity Validation in UnifiedRecommender

## Overview

The UnifiedRecommender now implements **two complementary validation mechanisms** to filter recommendations after DRFP similarity matching:

1. **applies_if** (Rules) - Functional group detection
2. **reaction_SMARTS** (Protocols) - Exact transformation matching

Both mechanisms act as **post-similarity filters** to ensure recommended conditions are chemically appropriate, even when DRFP similarity is high.

## How It Works

### Three-Stage Filtering Process

**Stage 1: DRFP Similarity Matching** (Discovery)
- Computes DRFP fingerprint for query reaction
- Compares with all `reference_reactions` DRFPs in the index
- Ranks by Tanimoto similarity score
- Returns top-k most similar sources

**Stage 2A: applies_if Validation for Rules** (Verification) ✨
- For each **rule** in the results, loads its full definition
- Detects chemical features from the query reaction
- Checks if detected features satisfy the rule's `applies_if` criteria
- **Filters out** rules that don't meet the criteria

**Stage 2B: reaction_SMARTS Validation for Protocols** (Verification) ✨ **NEW**
- For each **protocol** in the results, loads its full definition
- Extracts `reaction_SMARTS` patterns from protocol
- Uses RDKit to check if pattern matches query reactants
- **Filters out** protocols whose patterns don't match

## Example 1: applies_if for Rules (Amide Formation)

### Rule Definition
```json
{
  "applies_if": {
    "all": ["carboxylic_acid_present"],
    "any": ["primary_amine_present", "secondary_amine_present", "aniline_present"]
  }
}
```

### Test Case 1: Acid Chloride + Aniline
```
Reaction: ClC(=O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1
DRFP Similarity: 1.000 (perfect match - it's a reference reaction!)

Detected Features:
  ✓ aniline_present: True
  ✗ carboxylic_acid_present: False

applies_if Check:
  ALL: carboxylic_acid_present = FALSE ❌
  ANY: aniline_present = TRUE ✓
  Result: FAILS (ALL condition not met)

Without validation: amide_formation_v2 returned (rank #1)
With validation: amide_formation_v2 FILTERED OUT ✅
```

This is correct behavior! Even though the DRFP is identical (acid chloride forms amides similarly), the rule is designed for **carboxylic acids**, not acid chlorides.

## Example 2: reaction_SMARTS for Protocols (Cyanation)

### Protocol Definition
```json
{
  "name": "Copper-catalyzed Cyanation of Alkenyl Iodides",
  "reaction_SMARTS": ["IC=C.CC(O)(C#N)C>>N#CC=C"],
  "compatible_functional_groups": ["ketone", "aromatic"],
  "incompatible_functional_groups": ["enamine", "alcohol"]
}
```

### Test Case 1: Alkenyl Iodide + Acetone Cyanohydrin (Should Match)
```
Reaction: IC=CCCCC.CC(C)(O)C#N>>N#CC=CCCCC
DRFP Similarity: 1.000 (perfect - it's the reference reaction)

Pattern: IC=C.CC(O)(C#N)C>>N#CC=C
         │    │              │
         │    │              └─ Alkenyl nitrile product
         │    └─ Acetone cyanohydrin
         └─ Alkenyl iodide (I-C=C pattern)

RDKit Match:
  Reactants: IC=CCCCC ✓ matches IC=C (pattern matched on I-C=C motif)
             CC(C)(O)C#N ✓ matches CC(O)(C#N)C (acetone cyanohydrin)
  Result: MATCHES ✅

Without validation: Cyanation protocol found (rank #1, similarity: 1.000)
With validation: Cyanation protocol found (rank #1, similarity: 1.000) ✅
```

### Test Case 2: Aryl Iodide + Acetone Cyanohydrin (Should NOT Match)
```
Reaction: Ic1ccccc1.CC(C)(O)C#N>>N#Cc1ccccc1
DRFP Similarity: 0.220 (somewhat similar - both are cyanations)

Pattern: IC=C.CC(O)(C#N)C>>N#CC=C
         │
         └─ Requires alkenyl iodide (I-C=C)

RDKit Match:
  Reactants: Ic1ccccc1 ✗ does NOT match IC=C (aryl I ≠ alkenyl I)
             CC(C)(O)C#N ✓ matches CC(O)(C#N)C
  Result: NO MATCH ❌

Without validation: Cyanation protocol found (rank #15, similarity: 0.220)
With validation: Cyanation protocol FILTERED OUT ✅
```

### Test Case 3: Alkyl Iodide + Acetone Cyanohydrin (Should NOT Match)
```
Reaction: ICCCCCC.CC(C)(O)C#N>>N#CCCCCCC
DRFP Similarity: ~0.15 (low similarity)

Pattern: IC=C.CC(O)(C#N)C>>N#CC=C
         │
         └─ Requires I adjacent to C=C

RDKit Match:
  Reactants: ICCCCCC ✗ does NOT match IC=C (no C=C adjacent to I)
             CC(C)(O)C#N ✓ matches CC(O)(C#N)C
  Result: NO MATCH ❌

Without validation: Cyanation protocol not in top 20 (too low similarity)
With validation: Cyanation protocol not in results ✅
```

**Key Insight**: `reaction_SMARTS` distinguishes between:
- Alkenyl iodides (`IC=C`) ✓
- Aryl iodides (`Ic1ccccc1`) ✗
- Alkyl iodides (`ICCCC`) ✗

This level of specificity is **impossible with `applies_if`** (functional group detection alone).

### Test Case 2: Carboxylic Acid + Benzylamine
```
Reaction: O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1
DRFP Similarity: 1.000

Detected Features:
  ✓ carboxylic_acid_present: True
  ✓ primary_amine_present: True

applies_if Check:
  ALL: carboxylic_acid_present = TRUE ✓
  ANY: primary_amine_present = TRUE ✓
  Result: PASSES ✅

With validation: amide_formation_v2 returned ✅
```

### Test Case 3: Suzuki Coupling (Completely Different)
```
Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
DRFP Similarity: ~0.05 (very low - not even in top 20)

Detected Features:
  ✗ carboxylic_acid_present: False
  ✗ primary_amine_present: False
  ✗ aniline_present: False

Result: Not returned due to low DRFP similarity
(applies_if validation not even needed)
```

## Usage

### Enable Validation (Default)
```python
from chemtools.recommend import UnifiedRecommender

rec = UnifiedRecommender()
results = rec.recommend(
    "O=C(O)c1ccccc1.Nc1ccccc1>>...",
    top_k=5,
    validate_rules=True  # Default
)
```

### Disable Validation
```python
results = rec.recommend(
    "O=C(O)c1ccccc1.Nc1ccccc1>>...",
    top_k=5,
    validate_rules=False  # Get all high-similarity results
)
```

## Validation Mechanism Details

### applies_if (Rules Only)

**Syntax:**
```json
{
  "applies_if": {
    "all": ["feature1", "feature2"],  // ALL must be present (AND logic)
    "any": ["feature3", "feature4"]   // At least ONE must be present (OR logic)
  }
}
```

**Logic:**
- **ALL conditions**: Every feature must be detected (AND)
- **ANY conditions**: At least one feature must be detected (OR)
- **Both present**: ALL conditions AND ANY conditions must both pass
- **Missing applies_if**: Rule is always included (permissive)

**Common Features:**
- `carboxylic_acid_present`
- `primary_amine_present`
- `secondary_amine_present`
- `aniline_present`
- `sp2_halide_present`
- `sp2_boron_present`
- `sp2_chloride_present` (more specific)
- `alpha_chiral_acid_present`
- etc.

**Accuracy:** ⭐⭐⭐ Good - detects functional groups but can't distinguish regioisomers

### reaction_SMARTS (Protocols Only)

**Syntax:**
```json
{
  "reaction_SMARTS": [
    "IC=C.CC(O)(C#N)C>>N#CC=C",
    "IC=CC.CC(O)(C#N)C>>N#CC=CC"
  ]
}
```

**Logic:**
- Uses RDKit `ReactionFromSmarts()` and `RunReactants()` to match patterns
- Each pattern encodes: reactants → products transformation
- If ANY pattern matches, protocol is included
- Tries both forward and reverse reactant order
- **Missing reaction_SMARTS**: Protocol is always included (permissive)

**Pattern Components:**
- **Reactants**: SMARTS patterns for substrates (e.g., `IC=C` = alkenyl iodide)
- **Products**: SMARTS patterns for products (e.g., `N#CC=C` = alkenyl nitrile)
- **Separator**: `>>` indicates transformation

**Examples:**
```
"IC=C.CC(O)(C#N)C>>N#CC=C"
 │    │              │
 │    │              └─ Product: alkenyl nitrile
 │    └─ Reagent: acetone cyanohydrin
 └─ Substrate: alkenyl iodide (I adjacent to C=C)

"Ic1ccccc1.NaCN>>N#Cc1ccccc1"
 │         │       │
 │         │       └─ Product: aryl nitrile
 │         └─ Reagent: sodium cyanide
 └─ Substrate: aryl iodide (I on aromatic ring)
```

**Accuracy:** ⭐⭐⭐⭐⭐ Excellent - distinguishes exact transformations, regioisomers, stereochemistry

### Comparison

| Feature | applies_if (Rules) | reaction_SMARTS (Protocols) |
|---------|-------------------|----------------------------|
| **Accuracy** | ⭐⭐⭐ Good | ⭐⭐⭐⭐⭐ Excellent |
| **Ease of writing** | ✅ Easy | ⚠️ Moderate |
| **Distinguishes** | Functional groups | Exact transformations |
| **Regiochemistry** | ❌ No | ✅ Yes |
| **Stereochemistry** | ❌ No | ✅ Yes |
| **Substrate type** | Limited | ✅ Precise |
| **Best for** | General rules | Specific protocols |
| **Maintenance** | Low | Moderate |

## Benefits

### Overall System Benefits

1. **Chemical Correctness**: Prevents recommending inappropriate conditions despite high DRFP similarity
2. **Layered Defense**: Two complementary mechanisms catch different types of mismatches
3. **Flexibility**: Can be disabled if you want all similar reactions regardless of chemistry
4. **Fail-Open**: If validation fails, sources are still included (doesn't break on errors)
5. **Appropriate Granularity**: Rules use simpler validation, protocols use precise validation

### applies_if Benefits (Rules)

- ✅ **Easy to write**: Just list functional groups
- ✅ **Good coverage**: Catches most common mismatches (e.g., acid vs ester)
- ✅ **Low maintenance**: Robust to minor SMARTS changes
- ✅ **Appropriate for rules**: Rules cover variations, so functional group filtering is sufficient

### reaction_SMARTS Benefits (Protocols)

- ✅ **Extremely accurate**: Distinguishes regioisomers (alkenyl vs aryl iodide)
- ✅ **Stereochemistry-aware**: Can distinguish Z/E alkenes, R/S centers
- ✅ **Mechanistic specificity**: Encodes exact transformation, not just substrates
- ✅ **Worth the effort for protocols**: Protocols are curated, high-value resources
- ✅ **Safety-critical**: Prevents dangerous recommendations (e.g., wrong cyanation method)

## Implementation Details

### Feature Detection
Uses `chemtools.rule.analyzer.FeatureAnalyzer` to detect features from reaction SMILES.

### Validation Logic (in UnifiedRecommender)
```python
def _check_applies_if(features, applies_if):
    # Check ALL conditions (all must be true)
    if 'all' in applies_if:
        if not all(features.get(f, False) for f in applies_if['all']):
            return False
    
    # Check ANY conditions (at least one must be true)
    if 'any' in applies_if:
        if not any(features.get(f, False) for f in applies_if['any']):
            return False
    
    return True
```

### Source Type Validation

| Source Type | Validation Mechanism | Field Used | Accuracy |
|-------------|---------------------|------------|----------|
| **Protocols** | reaction_SMARTS matching | `reaction_SMARTS` | ⭐⭐⭐⭐⭐ Excellent |
| **Rules** | Functional group detection | `applies_if` | ⭐⭐⭐ Good |

**Why different mechanisms?**
- **Protocols** are curated, specific procedures → worth precise validation
- **Rules** cover many variations → functional group filtering is sufficient and easier to maintain

## Testing

Run the validation tests:
```bash
python test_applies_if.py
python test_acid_chloride.py
```

## Future Enhancements

Potential improvements:
- Add warning/info message when rules are filtered out
- Track which specific conditions failed
- Add confidence scoring based on feature match quality
- Support for more complex logic (NOT, nested conditions)
- User-configurable feature detection sensitivity
