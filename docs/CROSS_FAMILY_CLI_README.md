# Cross-Family Recommendation CLI

A simple command-line tool for getting reaction condition recommendations across all reaction families without specifying catalyst or family type.

## Overview

This CLI provides a streamlined interface for cross-family recommendations:
- ✅ **Simple**: Just provide a reaction SMILES, no catalyst or family selection needed
- ✅ **Fast**: Automatically uses unified DRFP index if available
- ✅ **Readable**: Clean, emoji-enhanced output format
- ✅ **Flexible**: Supports both positional and flag-based arguments

## Quick Start

### Basic Usage

```bash
# Simplest form - just pass the reaction SMILES
python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
```

### Using --rxn Flag

```bash
python app/cross_family_recommendation_cli.py --rxn "CCBr.CCN>>CCNCC"
```

### Get More Precedents

```bash
python app/cross_family_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --k 100
```

### Disable DRFP (Feature-based Similarity Only)

```bash
python app/cross_family_recommendation_cli.py --rxn "..." --no-drfp
```

### Debug Mode

```bash
python app/cross_family_recommendation_cli.py --rxn "..." --debug
```

## Output Format

```
================================================================================
CROSS-FAMILY REACTION CONDITION RECOMMENDATION
================================================================================

🧪 Reaction: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1

✅ Using unified DRFP index for accurate cross-family search

📊 Strategy: ML
📈 Confidence: 0.85
🔍 Detected Family: C_N_Coupling (confidence: 0.92)

✨ Found 3 recommendation(s):

────────────────────────────────────────────────────────────────────────────────
Recommendation #1
────────────────────────────────────────────────────────────────────────────────
  🧬 Catalyst: CuI
  💧 Solvent: DMSO
  ⚗️  Reagents:
      base: K2CO3
      ligand: 1,10-phenanthroline
  🌡️  Temperature: 110 °C
  ⏱️  Time: 24 h
  📚 Precedent:
      Reaction ID: C_N_31-614-CAS-35265690
      Family: C_N_Coupling
      Similarity: 0.892
  ⭐ Score: 0.850

────────────────────────────────────────────────────────────────────────────────
Recommendation #2
────────────────────────────────────────────────────────────────────────────────
  ...

================================================================================
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `reaction` | Reaction SMILES (positional argument) | Required |
| `--rxn` | Reaction SMILES (alternative flag) | - |
| `--k` | Number of precedents to retrieve | 50 |
| `--no-drfp` | Disable DRFP, use feature-based similarity only | False |
| `--debug` | Print debug information | False |

## Features

### Automatic DRFP Detection

The CLI automatically detects if the unified DRFP index exists:

- ✅ **If unified index exists**: Uses DRFP for accurate similarity
- ⚠️ **If unified index missing**: Falls back to feature-based similarity

```
✅ Using unified DRFP index for accurate cross-family search
```

or

```
⚠️  Unified DRFP index not found - using feature-based similarity
   To enable DRFP: python data-processor/build_unified_recommendation_index.py
```

### Rich Output

The output includes:
- 🧪 **Reaction**: Input reaction SMILES
- 📊 **Strategy**: Recommendation strategy used
- 📈 **Confidence**: Overall confidence score
- 🔍 **Detected Family**: ML-detected reaction family
- 🧬 **Catalyst**: Recommended catalyst
- 💧 **Solvent**: Recommended solvent
- ⚗️ **Reagents**: Additional reagents (base, ligand, etc.)
- 🌡️ **Temperature**: Reaction temperature
- ⏱️ **Time**: Reaction time
- 📚 **Precedent**: Source precedent information
- ⭐ **Score**: Recommendation score

### Error Handling

Clear error messages for common issues:

```bash
# Missing reaction
python app/cross_family_recommendation_cli.py
# Output: ❌ Error: No reaction SMILES provided

# Invalid format
python app/cross_family_recommendation_cli.py "invalid"
# Output: ❌ Error: Invalid reaction SMILES format
#         Expected format: reactants>>product
```

## Comparison with local_recommendation_cli.py

| Feature | cross_family_recommendation_cli.py | local_recommendation_cli.py |
|---------|-----------------------------------|----------------------------|
| **Purpose** | Simple cross-family search | Full-featured recommendation |
| **Family selection** | ❌ No (always cross-family) | ✅ Yes (--family flag) |
| **Catalyst selection** | ❌ No | ✅ Yes (--catalyst flag) |
| **Strategy selection** | ❌ No (auto ML) | ✅ Yes (--strategy flag) |
| **Output format** | 📄 Print only | 💾 Print or save to file |
| **Simplicity** | ✅ Very simple | 🔧 Many options |
| **Use case** | Quick exploration | Detailed control |

## Examples

### C-N Coupling Reaction

```bash
python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
```

### Suzuki Coupling

```bash
python app/cross_family_recommendation_cli.py "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccccc1-c1ccccc1"
```

### Amide Formation

```bash
python app/cross_family_recommendation_cli.py "CC(=O)O.NCCc1ccccc1>>CC(=O)NCCc1ccccc1"
```

### C-O Coupling

```bash
python app/cross_family_recommendation_cli.py "Brc1ccccc1.OC>>c1ccccc1OC"
```

### With More Precedents

```bash
python app/cross_family_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --k 100
```

### Without DRFP (Faster, Less Accurate)

```bash
python app/cross_family_recommendation_cli.py --rxn "..." --no-drfp
```

## Setup

### Prerequisites

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Build unified DRFP index (recommended):**
   ```bash
   python data-processor/build_unified_recommendation_index.py
   ```
   
   Without this, the CLI will use feature-based similarity (faster but less accurate).

### Verify Installation

```bash
python app/cross_family_recommendation_cli.py --help
```

## Exit Codes

- `0`: Success (recommendations found)
- `1`: No recommendations found or error occurred

## Tips

1. **For best accuracy**: Build the unified DRFP index first
   ```bash
   python data-processor/build_unified_recommendation_index.py
   ```

2. **For faster results**: Use `--no-drfp` flag (trades accuracy for speed)

3. **For more options**: Increase `--k` value (default is 50)

4. **For troubleshooting**: Use `--debug` flag to see detailed information

## Troubleshooting

### "No recommendations found"

**Possible causes:**
- Reaction is too novel (no similar precedents)
- DRFP index not built yet
- Invalid reaction SMILES

**Solutions:**
1. Build unified index: `python data-processor/build_unified_recommendation_index.py`
2. Increase k value: `--k 100`
3. Validate reaction SMILES format

### "Unified DRFP index not found"

**Solution:**
```bash
python data-processor/build_unified_recommendation_index.py
```

Then run the recommendation again.

### Low confidence scores

**Possible causes:**
- Novel reaction without close precedents
- Limited dataset coverage for reaction type

**Solutions:**
- Try increasing `--k` to get more precedents
- Use `--debug` to see what's happening

## Integration

### Use in Scripts

```python
import subprocess
import json

# Run CLI and capture output
result = subprocess.run(
    ["python", "app/cross_family_recommendation_cli.py", 
     "--rxn", "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"],
    capture_output=True,
    text=True
)

# Check exit code
if result.returncode == 0:
    print("Success!")
    print(result.stdout)
else:
    print("No recommendations")
```

### Use Python API Directly

For programmatic access, use the Python API:

```python
from chemtools import chem

result = chem.recommend.conditions(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=50,
    search_all_families=True
)

print(result)
```

## See Also

- [`local_recommendation_cli.py`](local_recommendation_cli.py) - Full-featured CLI with all options
- [`data-processor/build_unified_recommendation_index.py`](../data-processor/build_unified_recommendation_index.py) - Build unified DRFP index
- [`docs/UNIFIED_DRFP_INDEX_GUIDE.md`](../docs/UNIFIED_DRFP_INDEX_GUIDE.md) - Unified index documentation
- [`docs/CROSS_FAMILY_SEARCH.md`](../docs/CROSS_FAMILY_SEARCH.md) - Cross-family search details

## Summary

This simple CLI provides a streamlined way to get cross-family recommendations:

```bash
# That's it! Just provide the reaction SMILES
python app/cross_family_recommendation_cli.py "reaction>>product"
```

No catalyst selection, no family specification, no output file management - just simple, fast recommendations! 🚀
# NOTE: This document refers to the legacy precedent search CLI.
# The unified recommender now uses `build_unified_recommendation_index.py`.
