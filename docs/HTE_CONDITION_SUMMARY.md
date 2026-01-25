# HTE Condition Summarization Guide

## Overview

Yes! Our HTE tools can comprehensively summarize all conditions for copper-catalyzed (or any metal) C-N coupling reactions for a given reactant pair. This guide shows you how.

## Quick Answer

For a specific reactant pair (e.g., aryl bromide + aryl amine with copper catalyst):

```bash
# Get all conditions with statistics
python -m chemtools.recommend.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --catalyst copper -k 20 --min-exp 1

# Export all matching experiments to CSV
python -m chemtools.recommend.analytics_cli export cu_cn_conditions.csv \
    --reaction "C_N" --catalyst Cu --reactant-a ArBr --reactant-b ArNH2
```

## Available Methods

### 1. Condition Recommendations (CLI)

Get ranked condition recommendations with statistics:

```bash
python -m chemtools.recommend.cli \
    -a "SMILES_A" \
    -b "SMILES_B" \
    --catalyst Cu \
    -k 20 \
    --min-exp 1
```

**Output includes:**
- All matching condition combinations
- Success rate (% experiments with yield > 50%)
- Average/median yields
- Number of experiments per condition
- Confidence scores

**Key parameters:**
- `-k`: Number of recommendations (default: 5)
- `--min-exp`: Minimum experiments per condition (use 1 to see all)
- `--catalyst`: Filter by metal (Cu, Pd, Ni, copper, palladium, etc.)
- `--reaction`: Filter by reaction type (C_N_Coupling, Suzuki, etc.)

### 2. Export Filtered Experiments (CLI)

Export all matching experiments to CSV for detailed analysis:

```bash
python -m chemtools.recommend.analytics_cli export output.csv \
    --catalyst Cu \
    --reactant-a ArBr \
    --reactant-b ArNH2 \
    --reaction "C_N"
```

**Output:** CSV file with all experimental data including:
- Reactant SMILES
- All condition components (catalyst, ligand, base, solvent, etc.)
- Yields and z-scores
- Reaction types

**Key parameters:**
- `--catalyst`: Filter by metal
- `--reactant-a`, `--reactant-b`: Filter by reactant types
- `--reaction`: Filter by reaction type
- `--min-yield`: Minimum yield threshold

### 3. Programmatic Access (Python API)

#### Get Condition Recommendations

```python
from chemtools.recommend import HTERecommender

recommender = HTERecommender()
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    catalyst_filter="Cu",
    top_k=20,
    min_experiments=1
)

# Access results
print(f"Found {len(result.recommendations)} conditions")
for rec in result.recommendations:
    print(f"{rec.catalyst}, {rec.ligand}, {rec.base}, {rec.solvent}")
    print(f"  Avg Yield: {rec.avg_yield:.1f}%")
    print(f"  Experiments: {rec.num_experiments}")
```

#### Export and Analyze Data

```python
from chemtools.recommend.analytics import HTEAnalytics
import pandas as pd

# Export filtered data
analytics = HTEAnalytics()
analytics.export_subset(
    output_path="cu_cn_conditions.csv",
    catalyst_filter="Cu",
    reactant_a_type="ArBr",
    reactant_b_type="ArNH2"
)

# Load and analyze
df = pd.read_csv("cu_cn_conditions.csv")

# Group by conditions
conditions = df.groupby(['Catalyst', 'Ligand', 'Base', 'Solvent']).agg({
    'AREA_TOTAL_REDUCED': ['count', 'mean', 'median', 'max']
}).reset_index()

# Analyze distributions
print(f"Unique catalysts: {df['Catalyst'].nunique()}")
print(f"Unique ligands: {df['Ligand'].nunique()}")
print(f"Top ligands:\n{df['Ligand'].value_counts().head()}")
```

## Complete Example

The `demo_hte_condition_summary.py` script demonstrates all capabilities:

```bash
python demo_hte_condition_summary.py
```

This comprehensive demo:
1. Gets ranked condition recommendations
2. Exports all matching experiments
3. Analyzes condition component distributions
4. Creates a condition matrix with statistics

**Generated files:**
- `demo_cu_cn_conditions.csv` - All experiments
- `demo_cu_cn_conditions_matrix.csv` - Condition combinations with aggregated statistics

## Use Cases

### 1. Finding Optimal Conditions

```bash
# Get top conditions ranked by success rate
python -m chemtools.recommend.cli -a "ArBr_SMILES" -b "ArNH2_SMILES" \
    --catalyst Cu -k 10 --min-exp 2
```

### 2. Exploring Condition Space

```bash
# Export all conditions for manual analysis
python -m chemtools.recommend.analytics_cli export all_cu_cn.csv \
    --catalyst Cu --reaction "C_N"
```

### 3. Component Analysis

```python
from chemtools.recommend.analytics import HTEAnalytics

analytics = HTEAnalytics()

# Get catalyst statistics
catalyst_stats = analytics.get_catalyst_stats(
    reaction_type="C_N",
    reactant_a_type="ArBr",
    reactant_b_type="ArNH2"
)

print(catalyst_stats)
```

### 4. Comparing Catalysts

```bash
# Export copper conditions
python -m chemtools.recommend.analytics_cli export cu_cn.csv \
    --catalyst Cu --reaction "C_N"

# Export palladium conditions  
python -m chemtools.recommend.analytics_cli export pd_cn.csv \
    --catalyst Pd --reaction "C_N"

# Compare in Python/Excel
```

## Output Format

### Recommendation Output

```
================================================================================
Recommendation #1
================================================================================
Confidence Score: 10.3/100
Success Rate: 0.0% (1 experiments)
Avg Yield: 49.8% | Median: 49.8%

ðŸ§ª CONDITIONS:
  Catalyst: CuI
  Ligand: PPBO
  Base: NaOtBu
  Solvent: tAmOH

ðŸ“Š STATISTICS:
  Reaction Type: C_N_Coupling
  Z-Score Range: 2.61 to 2.61
  Reactant Types: ArBr + ArNH2
```

### CSV Export Columns

- `Reaction_Type_Standardized`
- `AREA_TOTAL_REDUCED` (yield)
- `z-Score`
- `Reactant_A`, `Reactant_B` (SMILES)
- `Reactant_A_Type`, `Reactant_B_Type`
- `Catalyst`, `Ligand`, `Base`, `Solvent`
- `Secondary Solvent`, `Additive`, `Coupling Reagent`

## Key Features

âœ?**Supports all metal catalysts**: Cu, Pd, Ni, Ir, Rh, Ru, Pt, Au, Ag, Fe, Co, Zn  
âœ?**Reactant type detection**: Automatic classification (ArBr, ArNH2, etc.)  
âœ?**Statistical ranking**: Success rates, confidence scores, yield statistics  
âœ?**Flexible filtering**: By reaction type, catalyst, reactant types, yield  
âœ?**Multiple output formats**: Terminal display, JSON, CSV export  
âœ?**Batch processing**: Process multiple reactant pairs from file  

## Tips

1. **Use `--min-exp 1`** to see all conditions, even rare ones
2. **Export to CSV** for custom analysis in Excel/Python/R
3. **Filter early** with `--catalyst` and `--reaction` to reduce noise
4. **Check reactant types** - the system automatically detects ArBr, ArNH2, etc.
5. **Look at z-scores** - high positive z-scores indicate better-than-average conditions

## Limitations

- Single experiment conditions have low confidence scores
- Success rate threshold is fixed at 50% yield
- Reactant type detection requires RDKit
- Database coverage depends on available HTE data

## Related Tools

- `chemtools.recommend.cli` - Condition recommendations
- `chemtools.recommend.analytics_cli` - Database analytics
- `demo_hte_condition_summary.py` - Comprehensive demo
- `chemtools.recommend.analytics.HTEAnalytics` - Python API for analytics

## References

- HTE database: `data/HTE_db/HTE_canonical.csv` (65,468 experiments)
- Reactant detection: `chemtools.featurizers.analysis.reactants.classify_reactant_smiles`
- Implementation: `chemtools/recommend/recommender.py`, `chemtools/recommend/analytics.py`
