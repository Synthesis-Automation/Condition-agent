# Quick Reference: HTE Condition Summarization

## TL;DR

**Yes!** HTE tools can summarize all conditions for copper-catalyzed C-N coupling (or any reaction) for a given reactant pair.

## One-Liners

```bash
# Get all Cu-catalyzed C-N coupling conditions for a reactant pair
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" --catalyst Cu -k 20 --min-exp 1

# Export all data to CSV
python -m chemtools.HTE.analytics_cli export cu_cn.csv --catalyst Cu --reactant-a ArBr --reactant-b ArNH2

# Run comprehensive demo
python demo_hte_condition_summary.py
```

## What You Get

### 1. Condition Recommendations
- All unique condition combinations (catalyst, ligand, base, solvent)
- Statistics: avg/median yields, success rates, experiment counts
- Confidence scores for each condition
- Ranked by performance

### 2. Full Experimental Data
- Export to CSV with all 112 experiments (for ArBr + ArNH2 + Cu)
- All condition components including additives
- Individual yields and z-scores
- Reactant SMILES for each experiment

### 3. Analytics
- Component distributions (which ligands/bases are most common)
- Yield statistics across all conditions
- Success rates and confidence metrics
- Condition matrix with aggregated statistics

## Example Output

For **copper-catalyzed C-N coupling** of **aryl bromide + aryl amine**:

```
📊 Found 112 matching experiments (0.17% of database)

🏆 Top 20 conditions:

Rank  Catalyst  Ligand     Base     Solvent  Experiments  Avg_Yield_%  Success_Rate_%
1     CuI       PPBO       NaOtBu   tAmOH    1            49.8         0.0
2     CuI       DMPAO      NaOtBu   tAmOH    1            49.4         0.0
3     CuI       2-Isobu... NaOtBu   tAmOH    1            43.9         0.0
...

📊 Condition Component Analysis:
  Unique Catalysts: 2 (CuI, Cu(MeCN)4BF4)
  Unique Ligands: 25
  Unique Bases: 7 (Cs2CO3, NaOtBu, CsF, K3PO4, KOH, K2CO3, LiHMDS)
  Unique Solvents: 6 (tAmOH, DMSO, Dioxane, DMF, PhMe, MeTHF)

📈 Yield Statistics:
  Average Yield: 9.8%
  Median Yield: 6.3%
  Max Yield: 49.8%
```

## Files Generated

1. **Conditions CSV**: All experiments with full details
2. **Matrix CSV**: Aggregated statistics per condition combination

## Python API

```python
from chemtools.HTE import HTERecommender

recommender = HTERecommender()
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    catalyst_filter="Cu",
    top_k=20,
    min_experiments=1
)

for rec in result.recommendations:
    print(f"{rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent}: {rec.avg_yield:.1f}%")
```

## Key Parameters

- **`--catalyst`**: Filter by metal (Cu, Pd, Ni, copper, palladium, etc.)
- **`--reaction`**: Filter by reaction type (C_N_Coupling, Suzuki, etc.)
- **`-k`**: Number of recommendations (default: 5, use 20+ for comprehensive)
- **`--min-exp`**: Min experiments per condition (use 1 to see all)
- **`--reactant-a/b`**: Filter by reactant types (ArBr, ArNH2, etc.)

## Use Cases

✅ Find optimal conditions for a specific substrate pair  
✅ Compare different catalysts (Cu vs Pd vs Ni)  
✅ Explore ligand/base/solvent space  
✅ Export data for machine learning or DOE  
✅ Validate experimental results against HTE data  

## See Full Documentation

- **Detailed guide**: `docs/HTE_CONDITION_SUMMARY.md`
- **Comprehensive demo**: `demo_hte_condition_summary.py`
- **HTE analytics**: `docs/HTE_ANALYTICS_QUICKREF.md`
