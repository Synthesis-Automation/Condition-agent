# RXNMapper Integration for Bond Analysis

## Overview

RXNMapper has been successfully integrated into the Condition-Agent reaction analysis module, enabling automatic atom mapping and bond breaking/formation analysis for chemical reactions.

## New Features

### 1. Automatic Atom Mapping
Automatically add atom mapping to unmapped reaction SMILES using RXNMapper's attention-guided algorithm.

```python
from chemtools import add_atom_mapping

# Unmapped reaction
unmapped = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

# Add atom mapping
result = add_atom_mapping(unmapped)

if result['success']:
    print(result['mapped_smiles'])
    print(f"Confidence: {result['confidence']}")
```

### 2. Bond Breaking/Formation Analysis
Analyze which bonds break and form in a reaction, with automatic atom mapping if needed.

```python
from chemtools import analyze_bond_changes

# Works with both mapped and unmapped SMILES
result = analyze_bond_changes("Brc1ccccc1.C#C>>c1ccccc1C#C", auto_map=True)

if result['success']:
    print(f"Broken bonds: {result['broken_bonds']}")
    print(f"Formed bonds: {result['formed_bonds']}")
    print(f"Changed atoms: {result['changed_atoms']}")
    print(f"Spectator atoms: {result['spectator_atoms']}")
```

### 3. Unified API for Unmapped Reactions
One-step function to auto-map and analyze unmapped reactions.

```python
from chemtools import compare_unmapped_reaction_to_find_changes

result = compare_unmapped_reaction_to_find_changes("Br.C#C>>C#C")

if result['success']:
    print(f"Auto-mapped: {result['mapped_smiles']}")
    print(f"Confidence: {result['mapping_confidence']}")
    print(f"Broken: {result['broken_bonds']}")
    print(f"Formed: {result['formed_bonds']}")
```

### 4. Batch Processing
Efficiently process multiple reactions at once.

```python
from chemtools._atom_mapping import batch_add_atom_mapping

reactions = [
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "Ic1ccccc1.C#C>>C(#C)c1ccccc1",
    "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
]

results = batch_add_atom_mapping(reactions)

for i, result in enumerate(results, 1):
    if result['success']:
        print(f"{i}. ✅ Mapped (conf: {result['confidence']:.3f})")
```

## Installation

```bash
pip install rxnmapper
```

RXNMapper is an optional dependency. If not installed, the functions will gracefully return error messages with installation instructions.

## Public API

All functions are now exposed in the main `chemtools` namespace:

```python
from chemtools import (
    add_atom_mapping,           # Add atom mapping
    analyze_bond_changes,       # High-level bond analysis
    identify_changed_atoms_from_mapped_smiles,  # Low-level (requires mapped SMILES)
    compare_unmapped_reaction_to_find_changes,  # Auto-map then analyze
    rxnmapper_available,        # Check if RXNMapper is installed
)
```

## Integration with Existing Tools

### With Reaction Detection
```python
from chemtools import detect_reaction, analyze_bond_changes

smiles = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"

# Detect reaction type
detection = detect_reaction(smiles, use_ml=True)
print(f"Type: {detection['family']}")
print(f"Confidence: {detection['confidence']}")

# Analyze bonds
bonds = analyze_bond_changes(smiles, auto_map=True)
if bonds['success']:
    print(f"Broken: {bonds['broken_bonds']}")
    print(f"Formed: {bonds['formed_bonds']}")
```

### Interactive CLI Tool
The `app/reaction_analysis_interactive.py` tool has been updated with bond analysis:

```bash
python app/reaction_analysis_interactive.py
```

When analyzing a reaction, you'll be prompted:
```
🔬 Analyze bond breaking/formation? (y/n, default=n): y
```

This will automatically:
1. Add atom mapping using RXNMapper
2. Identify reaction center atoms
3. List broken and formed bonds
4. Provide interpretation (substitution, addition, etc.)

## Files Modified

### New Files
- `chemtools/_atom_mapping.py` - RXNMapper integration module
- `test_rxnmapper_integration.py` - Comprehensive test suite
- `RXNMAPPER_INTEGRATION.md` - This documentation

### Modified Files
- `chemtools/__init__.py` - Exposed bond analysis functions in public API
- `chemtools/util/reaction_center_detector.py` - Updated to use RXNMapper for unmapped reactions
- `app/reaction_analysis_interactive.py` - Added bond analysis option
- `requirements.txt` - Added rxnmapper as optional dependency

## Architecture

```
User Input (unmapped SMILES)
        ↓
add_atom_mapping()
        ↓
    RXNMapper
        ↓
Atom-mapped SMILES
        ↓
identify_changed_atoms_from_mapped_smiles()
        ↓
Bond Analysis Result
```

## Error Handling

All functions gracefully handle missing dependencies:

```python
result = add_atom_mapping("Br.C#C>>C#C")

if not result['success']:
    print(result['error'])  # "RXNMapper not installed..."
    print(result.get('recommendation'))  # "Install with: pip install rxnmapper"
```

## Performance Notes

- **RXNMapper initialization**: Takes a few seconds on first use (lazy loaded)
- **Batch processing**: Much faster than individual calls for multiple reactions
- **Mapping confidence**: Returned as 0-1 value when available
- **Caching**: RXNMapper instance is cached for reuse

## Examples

### Suzuki-Miyaura Coupling
```python
suzuki = "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
result = analyze_bond_changes(suzuki, auto_map=True)

# Output:
# Broken: [(6, 7)]  - Ar-Br bond
# Formed: [(6, 13)] - Ar-Ar' bond
```

### Sonogashira Coupling
```python
sono = "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1"
result = analyze_bond_changes(sono, auto_map=True)

# Output:
# Broken: [(1, 7)]  - Ar-Br bond
# Formed: [(1, 8)]  - Ar-C≡C bond
```

### Amide Coupling
```python
amide = "CC(=O)O.CCN>>CC(=O)NCC"
result = analyze_bond_changes(amide, auto_map=True)

# Output:
# Broken: [(2, 4)]  - C-OH bond
# Formed: [(2, 5)]  - C-N amide bond
```

## Testing

Run the comprehensive test suite:

```bash
python test_rxnmapper_integration.py
```

Test with interactive tool:

```bash
python app/reaction_analysis_interactive.py
```

## Future Enhancements

Potential improvements:
1. Integration with reaction mechanism prediction
2. Visualization of bond changes (highlight changed atoms)
3. Support for multi-step reactions
4. Enhanced confidence scoring
5. Alternative atom mapping tools (LocalMapper, Graphormer)

## Contributing

When adding new features to bond analysis:
1. Update `_atom_mapping.py` for new functionality
2. Expose in `chemtools/__init__.py` if public-facing
3. Add tests to `test_rxnmapper_integration.py`
4. Update this documentation

## References

- RXNMapper: https://github.com/rxn4chemistry/rxnmapper
- Paper: "Extraction of organic chemistry grammar from unsupervised learning of chemical reactions"
- Based on attention-guided atom mapping from transformer models

## License

Same as Condition-Agent project.
