# Simplified Gradio UI - Quick Start

## Overview

`ui_simple.py` is a streamlined Gradio interface focused exclusively on **reaction condition recommendations**. Unlike the full UI (`ui_gradio.py`), this simplified version provides a clean, focused experience without unnecessary features.

## Features

### 🤖 ML-Based Recommendations
- Uses structured precedent search with DRFP similarity
- Retrieves conditions from reaction datasets
- Provides ranked recommendations with confidence scores
- Shows catalyst/ligand, base, solvent, temperature, time

### 📋 Rule-Based Recommendations  
- Uses SchemeConditionDB pattern matching
- Matches reaction SMARTS patterns to expert-curated rules
- Provides specific conditions for matched patterns
- Covers Cu, Pd, Ni C-N couplings and amide formations

### 🔄 Combined Approach
- Get both ML and rule-based recommendations simultaneously
- Compare data-driven vs. expert-rule approaches
- Comprehensive condition suggestions

## Supported Reaction Types

| Reaction Type | ML Family | Rule Database |
|--------------|-----------|---------------|
| **C-N Coupling (Cu)** | `C_N_Coupling_Cu` | `C_N_Coupling_Cu_db.json` |
| **C-N Coupling (Pd)** | `C_N_Coupling_Pd` | `C_N_Coupling_Pd_db.json` |
| **C-N Coupling (Ni)** | `C_N_Coupling_Ni` | `C_N_Coupling_Ni_db.json` |
| **Amide Formation** | `Amide_Coupling` | `Amide_formation_db.json` |
| **Suzuki Coupling** | `Suzuki_CC` | *(ML only)* |

## Usage

### Launch the UI

```bash
# From project root
python app/ui_simple.py
```

Then open your browser to: **http://127.0.0.1:7860**

### Input Format

**Reaction SMILES**: `Reactants>>Product`

Examples:
```
# Buchwald-Hartwig C-N coupling (Pd)
Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

# Ullmann C-N coupling (Cu)
Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

# Ni-catalyzed C-N coupling
Brc1ccccc1.NC1CCCC1>>c1ccc(N2CCCC2)cc1

# Amide formation
CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1
```

### Workflow

1. **Enter reaction SMILES**
2. **Select reaction type** (or use Auto-detect)
3. **Choose method**:
   - Click "🤖 ML Recommendations" for data-driven suggestions
   - Click "📋 Rule-Based Recommendations" for expert rules
   - Click "🔄 Both Methods" for comprehensive results
4. **Review results** in the summary and tables

## Key Differences from Full UI

| Feature | Simplified UI | Full UI (`ui_gradio.py`) |
|---------|---------------|--------------------------|
| **Focus** | Recommendations only | 15+ tabs with all tools |
| **Lines of Code** | ~450 | ~2700 |
| **Tabs** | 1 (Recommendations) | 15+ (Detection, Registry, etc.) |
| **Complexity** | Minimal | Comprehensive |
| **Use Case** | Quick recommendations | Full toolkit access |
| **Dependencies** | Core only | Many optional features |

## Output Format

### ML Recommendations Table

| Rank | Core | Base | Solvent | Confidence | Support |
|------|------|------|---------|------------|---------|
| 1 | Pd(OAc)2/XPhos | K3PO4 | Toluene | 89.5% | 42 |
| 2 | Pd2(dba)3/BINAP | Cs2CO3 | Dioxane | 85.2% | 38 |

### Rule-Based Conditions Table

| Component | Value | Details |
|-----------|-------|---------|
| Catalyst | Pd(OAc)2 | 5 mol% |
| Ligand | XPhos | 10 mol% |
| Base | K3PO4 | 1.5 equiv |
| Solvent | Toluene | - |
| Temperature | 110°C | - |

## Dependencies

### Required
- `gradio` - UI framework
- `chemtools` - Core recommendation engine
- Standard library modules

### Optional
- `chemtools.scdb_matcher` - For rule-based recommendations
- `rdkit` - For SMILES validation (recommended)

## Configuration

Edit these constants in `ui_simple.py` to customize:

```python
# Rule database paths
RULE_DATABASES = {
    "C-N Coupling (Cu)": "data/conditionDB/C_N_Coupling_Cu_db.json",
    # ... add more databases
}

# ML family mappings
ML_FAMILY_MAP = {
    "C-N Coupling (Cu)": "C_N_Coupling_Cu",
    # ... add more mappings
}
```

## Troubleshooting

### No ML Recommendations
- Check that reaction dataset exists in `data/reaction_dataset/`
- Verify reaction SMILES is valid
- Try Auto-detect to confirm reaction type

### No Rule-Based Recommendations
- Ensure `chemtools.scdb_matcher` is available (included in chemtools package)
- Check that rule database file exists in `data/conditionDB/`
- Verify reaction matches SMARTS patterns in database

### Auto-detect Not Working
- Ensure reaction SMILES is in correct format: `Reactants>>Product`
- Check that router model files exist
- Try manually selecting reaction type

## Performance Tips

1. **Use Auto-detect** when uncertain about reaction type
2. **Start with "Both Methods"** for comprehensive results
3. **Increase top_k slider** (1-10) for more ML suggestions
4. **Check detection confidence** - low confidence may indicate need for manual type selection

## Development

To extend the UI:

1. **Add new reaction type**:
   - Update `RULE_DATABASES` with new DB path
   - Update `ML_FAMILY_MAP` with new family name
   - Add to `THEME` choices list

2. **Modify output format**:
   - Edit `format_ml_recommendations()` for ML results
   - Edit `format_rule_recommendations()` for rule results

3. **Add new features**:
   - Create new helper functions
   - Add new Gradio components in `create_ui()`
   - Wire up event handlers

## Comparison with Full UI

**Use `ui_simple.py` when:**
- You only need condition recommendations
- You want a fast, focused interface
- You're demonstrating the recommendation system
- You're new to the Condition Agent

**Use `ui_gradio.py` when:**
- You need full access to all chemtools features
- You want to explore detection, registry, routing, etc.
- You're doing comprehensive reaction analysis
- You need advanced features like molecule visualization

## License

Same as parent Condition-agent project.

## See Also

- Full UI: `app/ui_gradio.py`
- API: `app/main.py` (FastAPI)
- CLI: `chemtools/cli/registry.py`
- Docs: `docs/` directory
