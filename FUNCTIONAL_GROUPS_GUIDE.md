# Functional Groups Detection - Quick Reference

## Overview

Comprehensive functional group detection for organic molecules with **80+ functional groups** detected using SMARTS patterns (RDKit) with text-pattern fallbacks.

**Location:** `chemtools/util/functional_groups.py`

## Quick Start

### Context API (Recommended)

```python
from chemtools import chem

# Detect all functional groups
groups_dict = chem.functional_groups.detect("CC(=O)O")
# {'carboxylic_acid': True, 'carbonyl': True, 'alcohol': True, ...}

# Get list of present groups
groups_list = chem.functional_groups.get_groups("CC(=O)O")
# ['alcohol', 'carbonyl', 'carboxylic_acid']

# Check for specific group
has_ester = chem.functional_groups.has("CC(=O)OC", "ester")
# True

# Count occurrences
count = chem.functional_groups.count("O=C(O)CC(=O)O", "carboxylic_acid")
# 2

# Get categorized view
categories = chem.functional_groups.categorize("CC(=O)Oc1ccccc1")
# {'oxygen': ['ester', 'carbonyl'], 'aromatic': ['aromatic'], ...}

# Human-readable summary
summary = chem.functional_groups.summarize("CC(=O)O")
# "Oxygen: alcohol, carbonyl, carboxylic_acid"

# List all detectable groups
available = chem.functional_groups.list_available()
# ['acyl_chloride', 'alcohol', 'aldehyde', ...]
```

### Direct Import

```python
from chemtools.util.functional_groups import (
    detect_all,
    get_functional_groups,
    has_functional_group,
    count_functional_groups,
    get_group_categories,
    summarize_functional_groups,
    FUNCTIONAL_GROUP_SMARTS,
)

# Same functionality as Context API
groups = detect_all("c1ccc(Br)cc1N")
```

## Detectable Functional Groups (80+)

### Oxygen-Containing
- **Alcohols & Phenols:** `alcohol`, `phenol`, `enol`
- **Ethers:** `ether`, `epoxide`, `silyl_ether`
- **Carbonyls:** `carbonyl`, `aldehyde`, `ketone`
- **Carboxylic Derivatives:** `carboxylic_acid`, `ester`, `amide`, `anhydride`, `acyl_chloride`, `lactone`, `lactam`
- **Others:** `peroxide`, `n_oxide`, `imine_n_oxide`

### Nitrogen-Containing
- **Amines:** `amine_primary`, `amine_secondary`, `amine_tertiary`, `aniline`
- **Amides:** `amide`, `amide_primary`, `amide_secondary`, `amide_tertiary`
- **Nitriles & Imines:** `nitrile`, `imine`, `enamine`
- **Nitro & Azides:** `nitro`, `azide`, `isocyanate`
- **Others:** `hydrazine`, `hydroxylamine`, `urea`, `carbamate`, `imide`

### Sulfur-Containing
- **Thiols & Sulfides:** `thiol`, `sulfide`, `disulfide`, `thioester`
- **Sulfoxides & Sulfones:** `sulfoxide`, `sulfone`
- **Sulfonic Derivatives:** `sulfonic_acid`, `sulfonyl_chloride`, `sulfonamide`, `sulfate`

### Halides
- **Alkyl Halides:** `alkyl_fluoride`, `alkyl_chloride`, `alkyl_bromide`, `alkyl_iodide`
- **Aryl Halides:** `aryl_fluoride`, `aryl_chloride`, `aryl_bromide`, `aryl_iodide`
- **Acyl Halides:** `acyl_halide`, `acyl_chloride`

### Phosphorus-Containing
- `phosphine`, `phosphine_oxide`, `phosphate`

### Aromatic Heterocycles
- `pyridine`, `pyrrole`, `furan`, `thiophene`, `imidazole`, `thiazole`, `oxazole`

### Unsaturated
- `alkene`, `alkyne`, `aromatic`, `benzylic`, `allylic`, `propargylic`

### Special Leaving Groups
- `triflate`, `tosylate`, `mesylate`

### Protecting Groups
- `boc`, `cbz`, `fmoc`, `silyl_ether`

### Other Groups
- `carbonate`, `aziridine`, `epoxide`, `lactone`, `lactam`, `urea`, `carbamate`

## Categories

Groups are organized into chemical categories:

- **oxygen**: Alcohols, ethers, carbonyls, acids, esters, etc.
- **nitrogen**: Amines, amides, nitriles, nitro, etc.
- **sulfur**: Thiols, sulfides, sulfones, sulfonamides, etc.
- **phosphorus**: Phosphines, phosphates
- **halides**: All halogen-containing groups
- **aromatic**: Aromatic rings and heterocycles
- **unsaturated**: Alkenes, alkynes, allylic/benzylic
- **protecting_groups**: Common protecting groups
- **leaving_groups**: Triflate, tosylate, mesylate

## Examples

### Aspirin Analysis

```python
from chemtools import chem

aspirin = "CC(=O)Oc1ccccc1C(=O)O"
groups = chem.functional_groups.get_groups(aspirin)
# ['alcohol', 'ether', 'carbonyl', 'carboxylic_acid', 'ester', 'aromatic', 'benzylic']

categories = chem.functional_groups.categorize(aspirin)
# {
#   'oxygen': ['alcohol', 'ether', 'carbonyl', 'carboxylic_acid', 'ester'],
#   'aromatic': ['aromatic', 'benzylic'],
#   ...
# }
```

### Bromoaniline Analysis

```python
from chemtools import chem

bromoaniline = "c1ccc(Br)cc1N"
summary = chem.functional_groups.summarize(bromoaniline)
print(summary)
# Nitrogen: amine_primary, aniline
# Halides: aryl_bromide
# Aromatic: aromatic, benzylic
```

### Counting Groups

```python
from chemtools import chem

# Dicarboxylic acid
succinic_acid = "O=C(O)CC(=O)O"
count = chem.functional_groups.count(succinic_acid, "carboxylic_acid")
# 2
```

### Checking for Specific Groups

```python
from chemtools import chem

# Check for ester in aspirin
has_ester = chem.functional_groups.has("CC(=O)Oc1ccccc1C(=O)O", "ester")
# True

# Check for aldehyde in ketone
has_aldehyde = chem.functional_groups.has("CC(=O)C", "aldehyde")
# False
```

## Implementation Details

### SMARTS Patterns

Detection uses RDKit SMARTS patterns when available:

```python
FUNCTIONAL_GROUP_SMARTS = {
    "alcohol": "[OX2H]",
    "phenol": "c[OX2H]",
    "carboxylic_acid": "[CX3](=O)[OX2H1]",
    "ester": "[#6][CX3](=O)[OX2][#6]",
    # ... 80+ total patterns
}
```

### Fallback Mechanism

When RDKit is unavailable, falls back to text pattern matching:

```python
TEXT_PATTERNS = {
    "alcohol": ["oh", "[oh]"],
    "carboxylic_acid": ["c(=o)o", "cooh"],
    # ... simplified patterns
}
```

### Performance

- **RDKit Mode:** Fast SMARTS matching (~1-5ms per molecule)
- **Text Mode:** Very fast substring search (~0.1ms per molecule)
- **Caching:** Not cached by default (detection is fast enough)

## Integration Points

### Used By

- `chemtools/recommend/substrate_analysis.py` - Backward compatible wrapper
- Future: Constraint validation, reaction filtering, reagent compatibility

### Context API Path

```python
from chemtools import chem

# Access through unified context
chem.functional_groups.detect(smiles)
chem.functional_groups.get_groups(smiles)
chem.functional_groups.has(smiles, group_name)
chem.functional_groups.count(smiles, group_name)
chem.functional_groups.categorize(smiles)
chem.functional_groups.summarize(smiles)
chem.functional_groups.list_available()
```

## Testing

Run comprehensive test:

```bash
python test_functional_groups.py
```

Expected output: Detection of 80+ functional groups across diverse test molecules (acetic acid, aspirin, benzylamine, etc.)

## Backward Compatibility

Old `substrate_analysis.py` code still works:

```python
from chemtools.recommend.substrate_analysis import detect_functional_groups

# Returns limited set: {free_alcohol, phenol, sulfonamide, hydroxylamine}
groups = detect_functional_groups("c1ccc(O)cc1")
```

For comprehensive detection, migrate to:

```python
from chemtools import chem
groups = chem.functional_groups.detect("c1ccc(O)cc1")  # All 80+ groups
```

## Future Enhancements

- [ ] Add molecular fragment highlighting (visualize detected groups)
- [ ] Add reactivity prediction based on functional groups
- [ ] Integration with constraint system for automatic filtering
- [ ] Add functional group compatibility matrix
- [ ] Export to chemical drawing software formats

---

**Created:** October 15, 2025  
**Module:** `chemtools/util/functional_groups.py`  
**Context API:** `chem.functional_groups.*`  
**Test:** `test_functional_groups.py`
