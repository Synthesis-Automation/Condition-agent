# Cross-Family Recommendation CLI - Quick Reference

## One-Line Usage

```bash
python app/cross_family_recommendation_cli.py "reaction>>product"
```

## Common Commands

```bash
# Basic usage
python app/cross_family_recommendation_cli.py "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"

# Using flag
python app/cross_family_recommendation_cli.py --rxn "CCBr.CCN>>CCNCC"

# More precedents
python app/cross_family_recommendation_cli.py --rxn "..." --k 100

# Without DRFP (faster)
python app/cross_family_recommendation_cli.py --rxn "..." --no-drfp

# Debug mode
python app/cross_family_recommendation_cli.py --rxn "..." --debug

# Help
python app/cross_family_recommendation_cli.py --help
```

## Key Features

✅ **Simple** - Just reaction SMILES, no catalyst/family selection  
✅ **Fast** - Auto-uses unified DRFP index (121,990 reactions)  
✅ **Clean** - Emoji-enhanced readable output  
✅ **Smart** - Falls back gracefully if no unified index  

## Setup (First Time)

```bash
# 1. Build unified DRFP index (optional but recommended)
python scripts/build_unified_drfp_index.py --force

# 2. Run recommendations
python app/cross_family_recommendation_cli.py "reaction>>product"
```

## Output Indicators

- ✅ **"Using unified DRFP index"** = Best accuracy (DRFP-based)
- ⚠️ **"Unified DRFP index not found"** = Using feature-based similarity

## When to Use This CLI

✅ **Use cross_family_recommendation_cli.py when:**
- Quick exploration of unknown reactions
- Don't know the reaction family
- Want simple, fast results
- Testing new reactions

🔧 **Use local_recommendation_cli.py when:**
- Need detailed control (catalyst, family, strategy)
- Want to save output to file
- Comparing different strategies
- Production workflows

## Exit Codes

- `0` = Success (recommendations found)
- `1` = No recommendations or error

## Files

- **CLI**: `app/cross_family_recommendation_cli.py`
- **Docs**: `app/CROSS_FAMILY_CLI_README.md`
- **Summary**: `docs/CROSS_FAMILY_CLI_SUMMARY.md`

## Example Output

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
```

## Quick Troubleshooting

| Problem | Solution |
|---------|----------|
| "No recommendations found" | Build unified index: `python scripts/build_unified_drfp_index.py` |
| "Unified DRFP index not found" | Same as above |
| "Invalid reaction SMILES" | Check format: `reactants>>product` |
| Slow performance | Use `--no-drfp` flag |
| Need more details | Use `--debug` flag |

## That's It!

Just run:
```bash
python app/cross_family_recommendation_cli.py "your_reaction>>product"
```

No configuration, no complexity - just fast recommendations! 🚀
