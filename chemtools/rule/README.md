# Rule-Based Recommendation System

A feature-driven condition recommendation engine for synthetic chemistry.

## Overview

This system provides automated condition recommendations based on:
- **Calculable molecular features** (v2.2 with 108 features)
- **Rule databases** in JSON format (e.g., `suzuki.json`)
- **Modifier logic** for handling symptoms and edge cases

## Quick Start

### Command Line

```bash
# Single reaction
python -m chemtools.rule.cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1"

# With verbose output
python -m chemtools.rule.cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1" --verbose

# With symptoms
python -m chemtools.rule.cli "reaction_smiles" --symptoms low_yield side_products

# Interactive mode
python -m chemtools.rule.cli --interactive

# Custom database
python -m chemtools.rule.cli "reaction" --database data/rule_db/my_reaction.json
```

### Python API

```python
from chemtools.rule import RuleEngine

# Load database and create engine
engine = RuleEngine.from_file("data/rule_db/suzuki.json")

# Get recommendation
rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")

# Print formatted output
print(rec.format_summary())

# Get JSON
rec_dict = rec.to_dict()
```

## Architecture

### Core Components

1. **`models.py`**: Data classes
   - `RuleSpec`: Rule specification with feature matching
   - `ModifierSpec`: Conditional modifier specification
   - `AppliedRule`: Applied rule with matched features
   - `ConditionRecommendation`: Complete recommendation output

2. **`analyzer.py`**: Feature detection
   - `FeatureAnalyzer`: Extracts features from reaction SMILES
   - Integrates with `chemtools.featurizers.calculable`
   - Supports multiple reactant combination methods (union, all, first)

3. **`database.py`**: Rule database management
   - `RuleDatabase`: Loads and validates JSON rule databases
   - Supports applicability checking
   - Finds matching base rules and modifiers

4. **`engine.py`**: Orchestration
   - `RuleEngine`: Main recommendation workflow
   - Combines analyzer, database, and models
   - Handles batch processing

5. **`cli.py`**: Command-line interface
   - Single reaction processing
   - Batch file processing
   - Interactive mode

## Rule Database Schema

### Basic Structure

```json
{
  "name": "Reaction Name",
  "version": "2025-11-02",
  "applies_if": {
    "all": ["required_feature_1", "required_feature_2"]
  },
  "default_rule": {
    "id": "default_id",
    "conditions": {
      "pd_source": "Pd catalyst",
      "base": "Base",
      "solvent": "Solvent"
    }
  },
  "base_rules": [...],
  "modifiers": [...]
}
```

### Base Rules

Rules are checked in order. First match wins.

```json
{
  "name": "Human-readable name",
  "id": "rule_id",
  "description": "What this rule handles",
  "reactant_features": {
    "and": ["feature1", "feature2"]  // All required
    // OR
    "any": ["feature1", "feature2"]  // At least one
    // OR
    "all": ["feature1", "feature2"]  // All required (explicit)
  },
  "conditions": {
    "pd_source": "Recommended catalyst",
    "ligand": "Recommended ligand",
    "base": "Recommended base"
  },
  "priority": 0  // Optional ordering
}
```

### Modifiers

Applied if any condition matches:

```json
{
  "when": [
    "molecular_feature",           // Feature-based
    "symptom:low_yield"           // Symptom-based
  ],
  "suggestion": "What to change",
  "rationale": "Why this helps",  // Optional
  "priority": 0                    // Optional
}
```

## Feature System Integration

The system uses **108 calculable features** from `chemtools/taxonomy/data/calculable_features.json`:

### Key Feature Categories

- **Halides**: `sp2_halide_present`, `sp2_chloride_present`, `ArBr_present`, etc.
- **Boron**: `sp2_boron_present`, `boronic_acid_present`, `boronic_ester_present`
- **Electronics**: `electron_withdrawing_present`, `strong_ewg_present`
- **Sterics**: `ortho_substitution_present`, `steric_hindrance`
- **Functional Groups**: `carbonyl_present`, `alcohol_present`, `amine_present`
- **Counts**: `halogen_count`, `sp2_halide_site_count`

### Derived Shortcuts

Some features are computed from others:
```json
{
  "name": "sp2_halide_present",
  "derive": "sp2_chloride_present OR sp2_bromide_present OR sp2_iodide_present OR sp2_fluoride_present"
}
```

## Workflow

1. **Parse reaction SMILES** → Extract reactants and products
2. **Detect features** → Use `detect_all_features()` on each reactant
3. **Check applicability** → Verify `applies_if` conditions
4. **Find matching rule** → Check `base_rules` in order, fallback to `default_rule`
5. **Apply modifiers** → Check all modifiers for feature/symptom matches
6. **Generate recommendation** → Return structured `ConditionRecommendation`

## Examples

### Example 1: Simple Suzuki Coupling

```python
engine = RuleEngine.from_file("data/rule_db/suzuki.json")
rec = engine.recommend("Brc1ccccc1.OB(O)c1ccccc1>>c1ccccc1-c1ccccc1")
```

**Output:**
- Detects: `sp2_bromide_present`, `sp2_boron_present`, `aryl_halide_present`
- Matches: Default rule (no specific base_rules match)
- Conditions: PdCl2(dtbpf), K3PO4, dioxane/H2O, 80-100°C

### Example 2: Activated Aryl Chloride

```python
rec = engine.recommend("Clc1ccc(C(=O)O)cc1.OB(O)c1ccccc1>>...")
```

**Output:**
- Detects: `sp2_chloride_present`, `strong_ewg_present`, `sp2_boron_present`
- Matches: "Aryl Chloride (Unhindered/EWG)" base rule
- Conditions: Pd(dba)2, XPhos, K3PO4, 100-120°C

### Example 3: With Symptoms

```python
rec = engine.recommend(
    "Brc1ccccc1.OB(O)c1ccccc1>>...",
    symptoms=["low_yield"]
)
```

**Output:**
- Applies symptom-based modifiers
- Suggests: Temperature increase, extended reaction time, etc.

## Testing

Run the test suite:
```bash
pytest tests/test_rule_engine.py -v
```

**Test Coverage:**
- ✅ Model serialization/deserialization
- ✅ Feature matching logic (AND/OR/ALL)
- ✅ Database loading and validation
- ✅ Reaction SMILES parsing
- ✅ Feature detection
- ✅ End-to-end recommendation workflow
- ✅ Integration with real suzuki.json

## Extending the System

### Adding a New Reaction Type

1. Create JSON database: `data/rule_db/my_reaction.json`
2. Define `applies_if` conditions (required features)
3. Add `default_rule` for fallback
4. Add `base_rules` for specific substrate classes
5. Add `modifiers` for edge cases and troubleshooting
6. Test: `python -m chemtools.rule.cli "reaction" -d my_reaction`

### Adding New Features

1. Edit `chemtools/taxonomy/data/calculable_features.json`
2. Add feature definition with SMARTS or derivation logic
3. Features are automatically available to all rule databases
4. Use in `reactant_features` or `when` conditions

## Files

```
chemtools/rule/
├── __init__.py           # Package exports
├── models.py             # Data models (RuleSpec, ModifierSpec, etc.)
├── analyzer.py           # FeatureAnalyzer (SMILES → features)
├── database.py           # RuleDatabase (JSON loader)
├── engine.py             # RuleEngine (main orchestration)
├── cli.py                # Command-line interface
└── README.md             # This file

tests/
└── test_rule_engine.py   # Comprehensive test suite (21 tests)

data/rule_db/
└── suzuki.json           # Suzuki-Miyaura coupling database
```

## Design Principles

1. **Deterministic**: No ML/randomness - pure rule-based logic
2. **Feature-driven**: All decisions based on calculable molecular features
3. **Extensible**: Easy to add new reaction types and features
4. **Transparent**: Every recommendation shows matched features and reasoning
5. **Modular**: Clean separation of concerns (analyzer, database, engine)
6. **Production-ready**: Comprehensive testing, error handling, validation

## Future Enhancements

- [ ] Confidence scoring based on literature precedent counts
- [ ] Multi-database recommendations (try multiple reaction types)
- [ ] Integration with `chemtools.precedent` for database lookups
- [ ] Web UI for interactive recommendation
- [ ] Export to reaction planning tools
- [ ] Batch optimization across reaction series

## Related Systems

- **`chemtools/featurizers/`**: Calculable features system (v2.2)
- **`chemtools/precedent/`**: Literature precedent search
- **`chemtools/protocol/`**: Protocol-based recommendations (DRFP similarity)
- **`chemtools/recommend/`**: ML-based condition recommendations

---

**Status**: ✅ Production-ready (v1.0)  
**Tests**: 21/21 passing  
**Database**: suzuki.json (Suzuki-Miyaura coupling)
