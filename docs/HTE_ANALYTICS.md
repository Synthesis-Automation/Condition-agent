# HTE Analytics Tools

Comprehensive analytics interface for exploring and analyzing the HTE database.

## Overview

The HTE Analytics module provides tools to:

- **List reactant pairs** by reaction type and/or catalyst
- **Analyze catalyst statistics** for specific reactions or substrates
- **Summarize reaction types** in the database
- **Analyze metal usage patterns** across reactions
- **Find similar reactant pairs** based on reaction type or catalyst
- **Export filtered datasets** for further analysis

## Installation

The analytics tools are part of the `chemtools.recommend` module:

```python
from chemtools.recommend import HTEAnalytics

analytics = HTEAnalytics()
```

## Quick Start

### Python API

```python
from chemtools.recommend import HTEAnalytics

# Initialize analytics
analytics = HTEAnalytics()

# List all Suzuki reactant pairs with Pd catalysts
pairs = analytics.list_reactant_pairs(
    reaction_type="Suzuki",
    catalyst_filter="Pd",
    min_experiments=50
)

# Analyze Cu catalysts for C-N coupling
catalysts = analytics.get_catalyst_stats(
    reaction_type="C-N",
    catalyst_filter="Cu"
)

# Get reaction type summary
reactions = analytics.get_reaction_type_summary()

# Analyze metal usage
metals = analytics.analyze_metal_usage()
```

### Command Line Interface

```bash
# List Suzuki reactant pairs with Pd catalysts
python -m chemtools.recommend.analytics_cli pairs --reaction Suzuki --catalyst Pd --top 10

# Analyze Cu catalysts
python -m chemtools.recommend.analytics_cli catalysts --reaction "C-N" --catalyst Cu --compact

# View reaction type summary
python -m chemtools.recommend.analytics_cli reactions --top 20

# Analyze metal usage
python -m chemtools.recommend.analytics_cli metals --detailed

# Export filtered dataset
python -m chemtools.recommend.analytics_cli export suzuki_pd.csv --reaction Suzuki --catalyst Pd --min-yield 50
```

## API Reference

### HTEAnalytics Class

#### `__init__(hte_db_path="data/HTE_db/HTE_canonical.csv")`

Initialize analytics with HTE database.

**Parameters:**

- `hte_db_path` (str): Path to HTE database JSONL file

**Example:**

```python
analytics = HTEAnalytics("data/HTE_db/HTE_canonical.csv")
```

---

#### `list_reactant_pairs(reaction_type=None, catalyst_filter=None, min_experiments=1, sort_by="count")`

List all reactant type pairs in the database with optional filters.

**Parameters:**

- `reaction_type` (str, optional): Filter by reaction type (e.g., "Suzuki", "C-N Coupling")
- `catalyst_filter` (str, optional): Filter by catalyst metal (e.g., "Pd", "Cu", "palladium")
- `min_experiments` (int): Minimum number of experiments for inclusion (default: 1)
- `sort_by` (str): Sort by "count" or "success_rate" (default: "count")

**Returns:**

- `pd.DataFrame` with columns:
  - `Reactant_A_Type`: Type of reactant A
  - `Reactant_B_Type`: Type of reactant B
  - `Reaction_Type`: Reaction type
  - `Num_Experiments`: Number of experiments
  - `Avg_Yield`: Average yield (%)
  - `Median_Yield`: Median yield (%)
  - `Success_Rate`: Percentage of experiments with yield > 50%
  - `Top_Catalyst`: Most frequently used catalyst

**Example:**

```python
# Find common Suzuki pairs with Pd catalysts
df = analytics.list_reactant_pairs(
    reaction_type="Suzuki",
    catalyst_filter="Pd",
    min_experiments=100,
    sort_by="count"
)

print(f"Found {len(df)} pairs")
for _, row in df.head(5).iterrows():
    print(f"{row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
    print(f"  Experiments: {row['Num_Experiments']}")
    print(f"  Success: {row['Success_Rate']:.1f}%")
```

**Output:**

```
Found 14 pairs
ArCl + ArB(OR)2
  Experiments: 2528
  Success: 25.9%
ArBr + ArB(OH)2
  Experiments: 1908
  Success: 32.6%
...
```

---

#### `get_catalyst_stats(reaction_type=None, reactant_a_type=None, reactant_b_type=None)`

Get catalyst usage statistics with optional filters.

**Parameters:**

- `reaction_type` (str, optional): Filter by reaction type
- `reactant_a_type` (str, optional): Filter by reactant A type
- `reactant_b_type` (str, optional): Filter by reactant B type

**Returns:**

- `pd.DataFrame` with columns:
  - `Catalyst`: Catalyst name
  - `Metal`: Extracted metal symbol (Pd, Cu, Ni, etc.)
  - `Num_Experiments`: Number of experiments
  - `Avg_Yield`: Average yield (%)
  - `Success_Rate`: Percentage of experiments with yield > 50%
  - `Reaction_Types`: Comma-separated list of reaction types

**Example:**

```python
# Find top Pd catalysts for Suzuki
df = analytics.get_catalyst_stats(reaction_type="Suzuki")
df_pd = df[df['Metal'] == 'Pd']

for _, row in df_pd.head(10).iterrows():
    print(f"{row['Catalyst']}: {row['Num_Experiments']} exp, {row['Avg_Yield']:.1f}% avg")
```

**Output:**

```
dtbpfPdCl2: 1635 exp, 33.4% avg
P(tBu)3 Pd(crotyl)Cl: 999 exp, 26.2% avg
SPhos Pd(crotyl)Cl: 815 exp, 35.2% avg
...
```

---

#### `get_reaction_type_summary()`

Get summary statistics for all reaction types in the database.

**Returns:**

- `pd.DataFrame` with columns:
  - `Reaction_Type`: Reaction type name
  - `Num_Experiments`: Total experiments
  - `Num_Reactant_Pairs`: Number of unique reactant pairs
  - `Num_Catalysts`: Number of unique catalysts
  - `Avg_Yield`: Average yield (%)
  - `Success_Rate`: Percentage of experiments with yield > 50%
  - `Top_Catalyst`: Most frequently used catalyst
  - `Top_Reactant_Pair`: Most common reactant pair

**Example:**

```python
df = analytics.get_reaction_type_summary()

for _, row in df.head(10).iterrows():
    print(f"{row['Reaction_Type']}: {row['Num_Experiments']:,} exp")
    print(f"  Pairs: {row['Num_Reactant_Pairs']}, Catalysts: {row['Num_Catalysts']}")
    print(f"  Success: {row['Success_Rate']:.1f}%, Top: {row['Top_Catalyst']}")
```

**Output:**

```
C_N_Coupling: 24,012 exp
  Pairs: 38, Catalysts: 111
  Success: 16.7%, Top: CuI
Suzuki: 11,588 exp
  Pairs: 22, Catalysts: 74
  Success: 27.4%, Top: dtbpfPdCl2
...
```

---

#### `analyze_metal_usage()`

Analyze metal usage across the entire database.

**Returns:**

- `dict` with:
  - `metal_distribution`: DataFrame with metal counts and percentages
  - `by_reaction_type`: Dict mapping metals to reaction type counts
  - `total_experiments`: Total number of experiments

**Example:**

```python
result = analytics.analyze_metal_usage()

print(f"Total: {result['total_experiments']:,} experiments\n")

for _, row in result['metal_distribution'].iterrows():
    print(f"{row['Metal']}: {row['Num_Experiments']:,} ({row['Percentage']:.1f}%)")

# Top reactions for Pd
pd_reactions = result['by_reaction_type']['Pd']
for rxn, count in sorted(pd_reactions.items(), key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {rxn}: {count:,}")
```

**Output:**

```
Total: 66,308 experiments

Pd: 45,205 (68.2%)
  C_N_Coupling: 20,265
  Suzuki: 11,588
  Arylation-acidic-C-H: 4,152
Cu: 5,478 (8.3%)
  C_N_Coupling: 3,747
  CO-Coupling: 1,239
...
```

---

#### `find_similar_pairs(reactant_a_type, reactant_b_type, similarity_criteria="reaction_type")`

Find reactant pairs similar to the given pair.

**Parameters:**

- `reactant_a_type` (str): Reactant A type to search for
- `reactant_b_type` (str): Reactant B type to search for
- `similarity_criteria` (str): How to define similarity:
  - `"reaction_type"`: Same reaction type (default)
  - `"catalyst"`: Same catalyst metal
  - `"both"`: Same reaction type AND catalyst metal

**Returns:**

- `pd.DataFrame` of similar reactant pairs with statistics

**Example:**

```python
# Find pairs similar to ArCl + ArB(OH)2 by reaction type
df = analytics.find_similar_pairs(
    reactant_a_type="ArCl",
    reactant_b_type="ArB(OH)2",
    similarity_criteria="reaction_type"
)

print(f"Found {len(df)} similar pairs\n")
for _, row in df.head(5).iterrows():
    print(f"{row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
    print(f"  {row['Reaction_Type']}: {row['Num_Experiments']} exp")
```

---

#### `export_subset(output_path, reaction_type=None, catalyst_filter=None, reactant_a_type=None, reactant_b_type=None, min_yield=None)`

Export a filtered subset of the database to CSV.

**Parameters:**

- `output_path` (str): Path to save the CSV file
- `reaction_type` (str, optional): Filter by reaction type
- `catalyst_filter` (str, optional): Filter by catalyst metal
- `reactant_a_type` (str, optional): Filter by reactant A type
- `reactant_b_type` (str, optional): Filter by reactant B type
- `min_yield` (float, optional): Minimum yield threshold

**Returns:**

- `int`: Number of experiments exported

**Example:**

```python
# Export high-performing Suzuki data
count = analytics.export_subset(
    output_path="suzuki_high_yield.csv",
    reaction_type="Suzuki",
    catalyst_filter="Pd",
    min_yield=80.0
)
print(f"Exported {count:,} experiments")
```

---

## CLI Reference

### Commands

#### `pairs` - List Reactant Pairs

```bash
python -m chemtools.recommend.analytics_cli pairs [OPTIONS]
```

**Options:**

- `--reaction TEXT`: Filter by reaction type
- `--catalyst TEXT`: Filter by catalyst metal (e.g., Pd, Cu)
- `--min-experiments INT`: Minimum number of experiments (default: 1)
- `--sort {count,success_rate}`: Sort by count or success_rate (default: count)
- `--top INT`: Number of results to show (default: 20)
- `--compact`: Use compact output format
- `-o, --output PATH`: Save results to CSV

**Examples:**

```bash
# Top 10 Suzuki pairs with Pd catalysts
python -m chemtools.recommend.analytics_cli pairs --reaction Suzuki --catalyst Pd --top 10 --compact

# All C-N coupling pairs with Cu, sorted by success rate
python -m chemtools.recommend.analytics_cli pairs --reaction "C-N" --catalyst Cu --sort success_rate

# Export to CSV
python -m chemtools.recommend.analytics_cli pairs --reaction Suzuki --catalyst Pd -o suzuki_pairs.csv
```

---

#### `catalysts` - Analyze Catalysts

```bash
python -m chemtools.recommend.analytics_cli catalysts [OPTIONS]
```

**Options:**

- `--reaction TEXT`: Filter by reaction type
- `--reactant-a TEXT`: Filter by reactant A type
- `--reactant-b TEXT`: Filter by reactant B type
- `--top INT`: Number of results to show (default: 20)
- `--compact`: Use compact output format
- `-o, --output PATH`: Save results to CSV

**Examples:**

```bash
# Top Pd catalysts for Suzuki
python -m chemtools.recommend.analytics_cli catalysts --reaction Suzuki --compact

# Cu catalysts for ArBr + RNH2
python -m chemtools.recommend.analytics_cli catalysts --reactant-a ArBr --reactant-b RNH2 --compact
```

---

#### `reactions` - Analyze Reaction Types

```bash
python -m chemtools.recommend.analytics_cli reactions [OPTIONS]
```

**Options:**

- `--top INT`: Number of results to show (default: 20)
- `--compact`: Use compact output format
- `-o, --output PATH`: Save results to CSV

**Examples:**

```bash
# View top 15 reaction types
python -m chemtools.recommend.analytics_cli reactions --top 15 --compact

# Export full summary
python -m chemtools.recommend.analytics_cli reactions -o reaction_summary.csv
```

---

#### `metals` - Analyze Metal Usage

```bash
python -m chemtools.recommend.analytics_cli metals [OPTIONS]
```

**Options:**

- `--detailed`: Show detailed breakdown by reaction type
- `-o, --output PATH`: Save results to CSV

**Examples:**

```bash
# Basic metal distribution
python -m chemtools.recommend.analytics_cli metals

# Detailed breakdown by reaction
python -m chemtools.recommend.analytics_cli metals --detailed
```

---

#### `export` - Export Filtered Dataset

```bash
python -m chemtools.recommend.analytics_cli export OUTPUT_PATH [OPTIONS]
```

**Options:**

- `--reaction TEXT`: Filter by reaction type
- `--catalyst TEXT`: Filter by catalyst metal
- `--reactant-a TEXT`: Filter by reactant A type
- `--reactant-b TEXT`: Filter by reactant B type
- `--min-yield FLOAT`: Minimum yield threshold

**Examples:**

```bash
# Export high-yield Suzuki data
python -m chemtools.recommend.analytics_cli export suzuki_high.csv --reaction Suzuki --min-yield 80

# Export Cu-catalyzed C-N coupling
python -m chemtools.recommend.analytics_cli export cn_copper.csv --reaction "C-N" --catalyst Cu

# Export specific substrate combination
python -m chemtools.recommend.analytics_cli export arbr_arb.csv --reactant-a ArBr --reactant-b "ArB(OH)2"
```

---

## Use Cases

### 1. Catalyst Selection

**Goal:** Find the most successful Pd catalyst for ArCl + ArB(OH)2 Suzuki coupling.

```python
from chemtools.recommend import HTEAnalytics

analytics = HTEAnalytics()

# Get catalyst stats for specific substrates
df = analytics.get_catalyst_stats(
    reaction_type="Suzuki",
    reactant_a_type="ArCl",
    reactant_b_type="ArB(OH)2"
)

# Filter for Pd and sort by success rate
df_pd = df[df['Metal'] == 'Pd'].sort_values('Success_Rate', ascending=False)

print("Top 5 Pd catalysts:")
for _, row in df_pd.head(5).iterrows():
    print(f"{row['Catalyst']}: {row['Success_Rate']:.1f}% success")
```

### 2. Reaction Scope Analysis

**Goal:** What reactant pairs work well with copper catalysts?

```python
# List all Cu-catalyzed pairs with good success rates
df = analytics.list_reactant_pairs(
    catalyst_filter="Cu",
    min_experiments=50,
    sort_by="success_rate"
)

print(f"Found {len(df)} Cu-catalyzed reactant pairs\n")
for _, row in df.head(10).iterrows():
    print(f"{row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
    print(f"  {row['Reaction_Type']}: {row['Success_Rate']:.1f}% success")
    print(f"  {row['Num_Experiments']} experiments")
```

### 3. Database Exploration

**Goal:** What are the most common reactions in the database?

```bash
# CLI: View top reactions
python -m chemtools.recommend.analytics_cli reactions --top 15 --compact
```

Output:

```
1. C_N_Coupling
   Experiments: 24,012, Pairs: 38, Catalysts: 111
   Performance: 16.7% success, 20.0% avg
   Top: CuI
   Most Common: ArBr + R2NH-a-branch

2. Suzuki
   Experiments: 11,588, Pairs: 22, Catalysts: 74
   Performance: 27.4% success, 30.4% avg
   Top: dtbpfPdCl2
   Most Common: ArCl + ArB(OR)2
...
```

### 4. Literature Search Preparation

**Goal:** Export all successful Pd-catalyzed Buchwald-Hartwig data for analysis.

```bash
# Export filtered dataset
python -m chemtools.recommend.analytics_cli export buchwald_pd_success.csv \
    --reaction "C-N" \
    --catalyst Pd \
    --min-yield 70
```

### 5. Metal Comparison

**Goal:** Compare Pd vs Cu usage in C-N coupling.

```python
# Analyze metal usage for C-N coupling
df = analytics.get_catalyst_stats(reaction_type="C-N")

pd_count = df[df['Metal'] == 'Pd']['Num_Experiments'].sum()
cu_count = df[df['Metal'] == 'Cu']['Num_Experiments'].sum()

print(f"C-N Coupling Metal Usage:")
print(f"  Palladium: {pd_count:,} experiments")
print(f"  Copper: {cu_count:,} experiments")
print(f"  Ratio: {pd_count/cu_count:.2f}:1 (Pd:Cu)")
```

---

## Database Schema

The HTE database (`HTE_canonical.csv`) contains the following fields:

- `reaction_type`: Standardized reaction type name
- `reactant_types`: Classified reactant types
- `catalyst_type`: Catalyst type tags (Pd, Cu, Ni, organocatalyst, etc.)
- `conditions`: Condition components (catalyst, ligand, base, solvent, secondary_solvent, additive, coupling_reagent)
- `metrics`: Yield/response metrics (area_total_reduced, z_score)

**Total:** 66,308 experiments across 41 reaction types

---

## Tips and Best Practices

### Filtering

1. **Reaction types** are case-insensitive substring matches:
   - `"Suzuki"` matches "Suzuki", "Suzuki-Miyaura", etc.
   - `"C-N"` matches "C_N_Coupling", "C-N-Coupling", etc.

2. **Catalyst filters** support both symbols and full names:
   - `"Pd"` or `"palladium"`
   - `"Cu"` or `"copper"`
   - `"Ni"` or `"nickel"`

3. **Minimum experiments** helps filter noise:
   - Use `min_experiments=10` for exploration
   - Use `min_experiments=50` for statistical significance
   - Use `min_experiments=100` for high-confidence trends

### Performance

1. **Cache the analytics object** if making multiple queries:

   ```python
   analytics = HTEAnalytics()  # Load once
   
   # Multiple queries
   pairs1 = analytics.list_reactant_pairs(reaction_type="Suzuki")
   pairs2 = analytics.list_reactant_pairs(reaction_type="C-N")
   ```

2. **Use CLI for quick exploration**, Python API for automation:

   ```bash
   # Quick check
   python -m chemtools.recommend.analytics_cli pairs --reaction Suzuki --top 5
   
   # Automated analysis
   python my_analysis_script.py
   ```

3. **Export filtered datasets** for external analysis (Excel, R, etc.):

   ```bash
   python -m chemtools.recommend.analytics_cli export my_data.csv --reaction Suzuki --catalyst Pd
   ```

---

## Troubleshooting

### No results found

**Problem:** Query returns 0 results

**Solutions:**

1. Check reaction type spelling (use `reactions` command to see available types)
2. Verify catalyst filter (use `metals` command to see available metals)
3. Lower `min_experiments` threshold
4. Check reactant type names match database (use `pairs` without filters to see all)

### Slow performance

**Problem:** Queries take a long time

**Solutions:**

1. Add more specific filters (reaction type, catalyst, etc.)
2. Use `min_experiments` to reduce result size
3. Limit results with `--top` in CLI
4. For large exports, consider filtering by yield threshold

### Unexpected catalyst counts

**Problem:** Catalyst appears in multiple metal categories

**Explanation:** Some catalysts contain multiple metals (e.g., "CuI, PPh3)2PdCl2")

- Metal extraction finds the first match
- Use `get_catalyst_stats()` to see full catalyst names
- Filter by exact catalyst name if needed

---

## Future Enhancements

Planned features:

- [ ] Ligand usage analysis
- [ ] Base/solvent combination analysis
- [ ] Time-series analysis (if timestamp data added)
- [ ] Automated statistical significance testing
- [ ] Visualization integration (plots, charts)
- [ ] Similarity search by chemical structure
- [ ] Machine learning-based predictions

---

## Related Documentation

- [HTE Recommender](./HTE_RECOMMENDER.md) - Condition recommendation system
- [Catalyst Filtering](../CATALYST_FILTER_FEATURE.md) - Catalyst filtering feature
- [Repository Guidelines](../AGENTS.md) - Project structure and conventions
