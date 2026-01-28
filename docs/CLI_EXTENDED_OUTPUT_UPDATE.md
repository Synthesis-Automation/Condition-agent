# CLI Update Summary: Extended Output for Conditions Recommendation

## Overview

Updated [Cpd_rxn_featurization_cli.py](app/Cpd_rxn_featurization_cli.py) to use the new v2 featurizer output with `detailed=True` option, displaying extended analysis information useful for conditions recommendation.

## Changes Made

### 1. Enable Detailed Output
- Featurizers now called with `options={'detailed': True}` to get extended format
- Both `featurize_molecule()` and `featurize_reaction()` use detailed output

### 2. New Display Functions

#### `_print_extended_section()`
Central function to display extended analysis section with:

**For Molecules:**
- Per-motif analysis (steric/electronic properties)
- SNAr feasibility assessment
- Nearby functional groups

**For Reactions:**
- Detection matches (top candidates with confidence)
- Aggregates (reaction-wide statistics)
- Role classification (reactant/agent roles)
- Intramolecular flags

#### `_print_reaction_type_summary()`
Updated to handle v2 format where `reaction_type` is a string (not dict):
- Shows reaction type, confidence, and reaction key
- Compatible with core v2 format

### 3. Updated Display Flow

**Molecule Display:**
```
SMILES
Motifs (with ranks)
Extended Analysis:
  - Per-Motif Analysis (steric, electronic, nearby groups)
  - SNAr Feasibility
```

**Reaction Display:**
```
Reaction SMILES
Reaction Key
Reaction Type + Confidence
Extended Analysis:
  - Detection Matches
  - Aggregates (reaction statistics)
  - Role Classification (reactants + agents)
Reactants (with extended analysis)
Products (with extended analysis)
```

## Information Useful for Conditions Recommendation

The CLI now shows critical information for matching reactions to conditions:

### Molecular Features
- **Motifs with ranks**: Identify functional groups and their specificity
- **Steric properties**: Bulkiness classification (0-10 scale)
- **Electronic properties**: Electron-donating/withdrawing effects (0-10 scale)
- **Nearby groups**: Contextual functional groups affecting reactivity
- **SNAr feasibility**: Special assessment for nucleophilic aromatic substitution

### Reaction Features
- **Reaction type + confidence**: Primary classification with certainty
- **Reaction key**: Compact motif-based summary (e.g., "Ar-B(OH)2|Ar-Br -> Ar-Ar")
- **Detection matches**: Alternative reaction types considered
- **Aggregates**: 
  - Reactant count
  - Reacted/formed/spectator motifs
  - Electronic properties (avg, electron-poor flags)
  - Steric properties (max values)
- **Role classification**:
  - Reactant roles (electrophile, nucleophile, etc.)
  - Expected vs. alternative functional groups
  - Multi-functional substrate detection
- **Agent roles**:
  - Catalyst/ligand/base/solvent classification
  - Role counts and flags

## Example Output

### Molecule (Bromobenzene)
```
SMILES: c1ccccc1Br
Motifs (1):
  - Ar-Br (a=5, b=6, rank=582.00)
Extended Analysis (for Conditions Recommendation):
  Per-Motif Analysis (1):
    - Ar-Br
      steric: score=0.0, classification=
      electronic: score=5.0, description=neutral
  SNAr Feasibility:
    - Ar-Br: NO, confidence=low, score=5.0
```

### Reaction (Suzuki Coupling)
```
Reaction SMILES: c1ccccc1Br.c1cccnc1B(O)O>>c1ccccc1-c1cccnc1
Reaction Key: Ar-B(OH)2|Ar-Br -> Ar-Ar || Pyridine
Reaction Type:
  - reaction_type: Arylation_Ar_H
  - confidence: 1.0000
Extended Analysis (for Conditions Recommendation):
  Detection Matches (showing 1 of 1):
    - Arylation_Ar_H [confidence: 1.00]
  Aggregates:
    - avg_aryl_electronic: 5
    - electron_poor_aryl: False
    - formed_motifs: Ar-Ar
    - max_alkyl_steric: 0
    - max_aryl_steric: 0
    - motif_ids: Ar-B(OH)2, Ar-Br
    - reacted_motifs: Ar-B(OH)2, Ar-Br
    - spectator_motifs: Pyridine
  Role Classification:
    - reaction_type: Suzuki_miyaura
    - num_reactants: 1
    Reactants:
      - pos 0: Ar-Br (electrophile)
```

## Testing

Created test scripts to validate the new output:
- [test_cli_output.py](test_cli_output.py): Basic structure validation
- [test_cli_display.py](test_cli_display.py): Full display formatting test

All tests pass ✅

## Backward Compatibility

The CLI maintains backward compatibility:
- Still handles v2 core format (if detailed=False)
- Functions gracefully degrade if extended section not present
- All existing CLI flags work as before (--show-roles, --show-rdkit, --format)

## Usage

Run the CLI as usual:
```bash
python app/Cpd_rxn_featurization_cli.py
```

The CLI will automatically use detailed output to show extended analysis useful for conditions recommendation.
