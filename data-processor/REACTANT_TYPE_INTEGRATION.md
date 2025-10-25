# Reactant Type Classification Integration

## Summary
Added automatic reactant type classification to `process_reactions.py` so that generated datasets include reactant type information for use in the recommendation system.

## What Was Changed

### 1. **Import Section (Lines ~30-44)**
Added import of the `classify_reactant` module:
```python
# Import reactant classifier if available
try:
    import sys
    import os
    # Add data-processor/other_data to path for classify_reactant import
    other_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'other_data')
    if other_data_path not in sys.path:
        sys.path.insert(0, other_data_path)
    from classify_reactant import classify_reactant, load_reactant_types
    REACTANT_TYPES = load_reactant_types()
except Exception:
    classify_reactant = None
    REACTANT_TYPES = None
```

### 2. **New Function: `classify_reactants_from_smiles()` (Lines ~157-191)**
Added a helper function to classify all reactants in a dot-separated SMILES string:
```python
def classify_reactants_from_smiles(smiles_string: str) -> Tuple[List[str], List[str]]:
    """
    Classify reactants from a SMILES string (dot-separated).
    
    Args:
        smiles_string: Dot-separated SMILES of reactants (e.g., "c1ccccc1Br.CCN")
        
    Returns:
        Tuple of (reactant_types, reactant_categories)
        Each is a list matching the reactant order in SMILES
    """
```

### 3. **Classification Call in `assemble_rows()` (Lines ~1560-1563)**
After generating SMILES from MOL blocks, classify the reactants:
```python
# Classify reactant types from SMILES if available
reactant_smiles_str = '.'.join(rct_smiles_list) if rct_smiles_list else ''
reactant_types, reactant_categories = classify_reactants_from_smiles(reactant_smiles_str)
```

### 4. **Add Fields to Row Dictionary (Lines ~1775-1776)**
Added two new fields to each reaction row:
```python
# Reactant type classification (added for recommendation system)
'ReactantTypes': _json_list(reactant_types),
'ReactantCategories': _json_list(reactant_categories),
```

### 5. **Updated CSV Column Order (Lines ~1799-1801)**
Added the new columns to the output CSV:
```python
# reactant type classification (for recommendation system)
'ReactantTypes', 'ReactantCategories',
```

## New CSV Columns

### `ReactantTypes`
- **Format**: JSON array of strings
- **Example**: `["ArBr", "RNH2"]`
- **Description**: Specific reactant type IDs (member types) for each reactant in the order they appear in ReactantSMILES
- **Values**: Corresponds to `member_type` from `reactant_types.json` (e.g., "ArBr", "ArCl", "RNH2", "ArB(OH)2")

### `ReactantCategories`
- **Format**: JSON array of strings
- **Example**: `["ArX*", "Aliphatic-amine"]`
- **Description**: High-level category for each reactant in the order they appear in ReactantSMILES
- **Values**: Corresponds to category names from `reactant_types.json` (e.g., "ArX*", "ArB*", "Aliphatic-amine")

## Usage Example

When you run `process_reactions.py`, it will now automatically:
1. Extract SMILES from RDF MOL blocks (if RDKit is available)
2. Classify each reactant using SMARTS pattern matching
3. Add the classification results to the CSV output

Example output row:
```csv
ReactionID,ReactionType,ReactantSMILES,ReactantTypes,ReactantCategories,...
RXN-001,Buchwald-Hartwig,c1ccccc1Br.CCN,"[""ArBr"", ""RNH2""]","[""ArX*"", ""Aliphatic-amine""]",...
```

## Benefits

1. **Automatic Classification**: No manual labeling needed - types are inferred from structure
2. **Consistency**: Uses the same `reactant_types.json` as the recommendation system
3. **Backward Compatible**: If classify_reactant is not available, the script still works (columns will be empty)
4. **Multi-reactant Support**: Handles multiple reactants in dot-separated SMILES
5. **Ready for Recommendations**: Generated datasets can be used directly by the condition recommender

## Testing

Run the test script to verify functionality:
```bash
python data-processor/test_reactant_classification.py
```

Expected output:
```
Test: ArBr (aryl bromide)
SMILES: c1ccccc1Br
Reactant Types: ['ArBr']
Reactant Categories: ['ArX*']

Test: ArBr + RNH2 (Buchwald-Hartwig reactants)
SMILES: c1ccccc1Br.CCN
Reactant Types: ['ArBr', 'RNH2']
Reactant Categories: ['ArX*', 'Aliphatic-amine']
```

## Dependencies

- **RDKit**: Required for SMILES extraction from MOL blocks (already a dependency)
- **classify_reactant.py**: Located in `data-processor/other_data/`
- **reactant_types.json**: Located in `data-processor/other_data/`

## Integration with Recommendation System

The generated CSV with `ReactantTypes` and `ReactantCategories` can be directly used by:
1. `simple_condition_recommender.py` - for hierarchical matching (exact → category → reaction type)
2. Any future recommendation systems that need reactant type information
3. Data analysis tools to understand substrate scope

The column names follow the generic naming convention (not "Electrophile"/"Nucleophile") to work with all reaction types including those without traditional electrophiles (e.g., amide coupling, CH activation).
