# Automatic Reaction Type Detection with rxn-insight

## Overview

The Reaction Predictor now integrates the **rxn-insight** library for automatic reaction type detection when users select "Auto-detect" as the reaction type. This significantly improves the accuracy of reaction classification compared to the previous pattern-based detection.

## What is rxn-insight?

- **Paper**: [Rxn-INSIGHT: Fast chemical reaction analysis using bond electron-replication patterns](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980627/)
- **GitHub**: https://github.com/mrodobbe/Rxn-INSIGHT
- **Purpose**: Automated classification of chemical reactions based on reaction SMILES

## Installation

```bash
pip install rxn-insight
```

Note: You may see warnings about missing Open-Reaction-Database modules, but these don't affect functionality.

## How It Works

### 1. Auto-Detection Process

When a user selects "Auto-detect":

1. **rxn-insight Analysis**: The reaction SMILES is analyzed by rxn-insight
2. **Classification**: Returns reaction class (e.g., "Heteroatom Alkylation and Arylation") and specific name (e.g., "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)")
3. **Mapping**: The rxn-insight results are mapped to our internal reaction types
4. **Recommendations**: Appropriate ligands and solvents are suggested based on the detected type

### 2. Mapping Examples

| rxn-insight Classification | Our Internal Type |
|---------------------------|------------------|
| "Heteroatom Alkylation and Arylation" / "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)" | C-N Coupling - Buchwald-Hartwig |
| "Heteroatom Alkylation and Arylation" / "Goldberg coupling aryl amine-aryl chloride" | C-N Coupling - Ullmann |
| "C-C Coupling" / "Suzuki coupling with boronic acids" | C-C Coupling - Suzuki |
| "Acylation" / "Carboxylic acid with primary amine to amide" | Amide Formation - Acid + Amine |

### 3. Example Output

For the reaction: `Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3[nH]2)cc1`

**Auto-Detection Results:**
- **rxn-insight Class**: Heteroatom Alkylation and Arylation
- **rxn-insight Name**: Goldberg coupling aryl amine-aryl chloride  
- **Mapped Type**: C-N Coupling - Ullmann
- **Confidence**: high
- **Functional Groups**: ['Aromatic halide', 'Primary amine'] → ['Secondary amine']

## GUI Integration

### Enhanced Display

When auto-detection is used, the GUI now shows:

```
🧪 ENHANCED REACTION CONDITION RECOMMENDATIONS

Reaction: Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3[nH]2)cc1
Detected Type: Ullmann (Cu)
Auto-Detection: Goldberg coupling aryl amine-aryl chloride (confidence: high)
Classification: Heteroatom Alkylation and Arylation
Status: success
```

### Fallback Behavior

- **If rxn-insight is not installed**: Falls back to pattern-based detection
- **If rxn-insight fails**: Falls back to pattern-based detection  
- **If no mapping found**: Uses the rxn-insight result directly or falls back to pattern-based detection

## Technical Implementation

### Key Files

1. **`reaction_type_detector.py`**: New module containing rxn-insight integration
2. **`enhanced_recommendation_engine.py`**: Enhanced with auto-detection capabilities
3. **`simple_reaction_gui.py`**: Updated to display auto-detection information

### Code Structure

```python
# Import with fallback
try:
    from reaction_type_detector import detect_reaction_type
    RXN_INSIGHT_AVAILABLE = True
except ImportError:
    RXN_INSIGHT_AVAILABLE = False

# Usage in recommendation engine
if RXN_INSIGHT_AVAILABLE and detect_reaction_type:
    detection_result = detect_reaction_type(reaction_smiles)
    # Map to internal types and use for recommendations
```

## Benefits

1. **Improved Accuracy**: rxn-insight provides much more accurate reaction classification than pattern matching
2. **Detailed Information**: Users see exactly how their reaction was classified
3. **Confidence Levels**: Know how confident the system is in the detection
4. **Functional Group Analysis**: Understand what chemical transformations are occurring
5. **Fallback Safety**: System gracefully falls back if rxn-insight is unavailable

## Supported Reaction Types

The integration currently maps the following rxn-insight classifications:

### Heteroatom Alkylation and Arylation
- N-arylation (Buchwald-Hartwig/Ullmann-Goldberg) → C-N Coupling - Buchwald-Hartwig
- Goldberg coupling aryl amine-aryl chloride → C-N Coupling - Ullmann
- Chan-Lam coupling → C-N Oxidative Coupling - Chan-Lam

### C-C Coupling  
- Suzuki coupling with boronic acids/esters → C-C Coupling - Suzuki
- Stille coupling → C-C Coupling - Stille
- Heck coupling → C-C Coupling - Heck
- Sonogashira coupling → C-C Coupling - Sonogashira
- Negishi coupling → C-C Coupling - Negishi

### Acylation
- Carboxylic acid with amine to amide → Amide Formation - Acid + Amine
- Acid chloride with amine to amide → Amide Formation

### Future Extensions

The mapping can be easily extended to support additional reaction types as they become relevant for the prediction system.
