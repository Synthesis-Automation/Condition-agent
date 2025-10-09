# CLI Modes Test - local_recommendation_cli.py

This document verifies that the CLI now supports both interactive and non-interactive modes.

## ✅ Non-Interactive Mode (No User Prompts)

When all required arguments are provided via command line, the script runs without prompting:

```bash
# Example 1: Full specification
python scripts/local_recommendation_cli.py \
  --rxn "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1" \
  --family Ullmann_CN \
  --strategy fusion \
  --save-dir test_output

# Example 2: Auto-detect reaction type
python scripts/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --strategy ml

# Example 3: Run only specific strategy
python scripts/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" \
  --family Buchwald_CN \
  --strategy rule
```

**Result**: ✅ No prompts shown, runs directly with provided arguments

## ✅ Interactive Mode (User Prompts)

When arguments are missing, the script prompts for input:

```bash
# No arguments - prompts for both reaction and type
python scripts/local_recommendation_cli.py

# Only strategy specified - prompts for reaction and type
python scripts/local_recommendation_cli.py --strategy fusion

# Only reaction specified - prompts for reaction type
python scripts/local_recommendation_cli.py --rxn "Brc1ccccc1>>c1ccccc1"
```

**Result**: ✅ Shows interactive prompts for missing information

## Available Command-Line Arguments

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `--rxn`, `--reaction` | str | None | Reaction SMILES (reactants>>products) |
| `--family`, `--type` | str | None | Reaction family (Ullmann_CN, Buchwald_CN, Suzuki, etc.) |
| `--k` | int | 50 | Number of precedents for ML |
| `--limit` | int | 5 | Number of ML recommendations |
| `--fusion-variants` | int | 5 | Number of fusion variants |
| `--save-dir` | str | "results" | Output directory |
| `--strategy` | str | "all" | Which strategy to run (all/rule/ml/fusion) |
| `--pretty` | flag | False | Pretty-print JSON output |

## Reaction Family Options

- `Suzuki` / `Suzuki_CC` - Suzuki-Miyaura coupling
- `Ullmann_CN` / `C_N_Coupling_Cu` - Ullmann C-N coupling (Cu)
- `Buchwald_CN` / `C_N_Coupling_Pd` - Buchwald-Hartwig C-N coupling (Pd)
- `C_N_Coupling_Ni` - Ni-catalyzed C-N coupling
- `Amide_formation` - Amide bond formation

## Examples from Testing

### ✅ Test 1: Non-interactive with full args
```bash
python scripts/local_recommendation_cli.py \
  --rxn "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1" \
  --family Ullmann_CN \
  --strategy fusion \
  --save-dir test_output
```

**Output**:
```
Local Recommendation Test
-------------------------
Reaction SMILES: Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1
Reaction type: Ullmann_CN
...
Fusion top recommendation:
  Core: Copper(I) iodide
  Base: Potassium carbonate
  Solvent: Dimethyl sulfoxide
  Confidence: LOW

Saved outputs:
  Fusion JSON: test_output\20251009_201948_ullmann_cn_fusion_local.json

Done.
```

✅ **No prompts shown** - ran directly with provided arguments
✅ **Correct Cu-based conditions** for Ullmann reaction

### ✅ Test 2: Interactive mode (no rxn argument)
```bash
python scripts/local_recommendation_cli.py --strategy fusion
```

**Output**:
```
Local Recommendation Test
-------------------------
Enter reaction SMILES: [USER PROMPTED]

Reaction Type Options:
  1) Auto-detect (server decides) (default)
  2) Suzuki Coupling
  3) Ullmann C–N (Cu)
  4) Buchwald C–N (Pd)
  5) C–N Coupling (Ni)
  6) Amide Formation
Select reaction type [1]: [USER PROMPTED]
...
```

✅ **Prompts shown** for missing arguments
✅ **Interactive selection** available

## Summary

The CLI now intelligently switches between modes:
- **Non-interactive**: All args provided → no prompts
- **Interactive**: Missing args → prompts user
- **Flexible**: Can mix (e.g., provide reaction but not type)

This makes it suitable for both:
- **Automation/scripting** (batch processing with arguments)
- **Exploratory work** (interactive testing)
