# HTE Analytics - Quick Reference Card

## Installation

```python
from chemtools.HTE import HTEAnalytics
analytics = HTEAnalytics()
```

## Common Operations

### 1. List Reactant Pairs

```python
# Find Suzuki pairs with Pd catalysts (≥50 experiments)
pairs = analytics.list_reactant_pairs(
    reaction_type="Suzuki",
    catalyst_filter="Pd",
    min_experiments=50
)
```

**CLI:**
```bash
python -m chemtools.HTE.analytics_cli pairs --reaction Suzuki --catalyst Pd --min-experiments 50 --compact
```

---

### 2. Analyze Catalysts

```python
# Get Cu catalyst statistics for C-N coupling
catalysts = analytics.get_catalyst_stats(
    reaction_type="C-N",
    catalyst_filter="Cu"
)
```

**CLI:**
```bash
python -m chemtools.HTE.analytics_cli catalysts --reaction "C-N" --catalyst Cu --compact
```

---

### 3. Reaction Type Summary

```python
# Get summary of all reaction types
reactions = analytics.get_reaction_type_summary()
print(f"Found {len(reactions)} reaction types")
```

**CLI:**
```bash
python -m chemtools.HTE.analytics_cli reactions --top 20 --compact
```

---

### 4. Metal Usage Analysis

```python
# Analyze metal distribution
result = analytics.analyze_metal_usage()
print(result['metal_distribution'])
```

**CLI:**
```bash
python -m chemtools.HTE.analytics_cli metals --detailed
```

---

### 5. Export Filtered Data

```python
# Export high-performing Suzuki data
count = analytics.export_subset(
    output_path="suzuki_high.csv",
    reaction_type="Suzuki",
    catalyst_filter="Pd",
    min_yield=80
)
```

**CLI:**
```bash
python -m chemtools.HTE.analytics_cli export suzuki_high.csv --reaction Suzuki --catalyst Pd --min-yield 80
```

---

## Filtering Options

All functions support:

| Filter | Type | Example | Description |
|--------|------|---------|-------------|
| `reaction_type` | str | `"Suzuki"`, `"C-N"` | Case-insensitive substring match |
| `catalyst_filter` | str | `"Pd"`, `"copper"` | Metal symbol or name |
| `reactant_a_type` | str | `"ArBr"`, `"ArCl"` | Exact type match |
| `reactant_b_type` | str | `"RNH2"`, `"ArB(OH)2"` | Exact type match |
| `min_experiments` | int | `50`, `100` | Minimum sample size |
| `min_yield` | float | `50.0`, `80.0` | Minimum yield threshold |

---

## Output Columns

### Reactant Pairs

- `Reactant_A_Type`, `Reactant_B_Type` - Substrate types
- `Reaction_Type` - Reaction name
- `Num_Experiments` - Sample size
- `Avg_Yield`, `Median_Yield` - Performance metrics
- `Success_Rate` - % experiments with yield > 50%
- `Top_Catalyst` - Most frequently used catalyst

### Catalysts

- `Catalyst` - Full catalyst name
- `Metal` - Extracted metal symbol (Pd, Cu, Ni, etc.)
- `Num_Experiments` - Number of uses
- `Avg_Yield` - Average performance
- `Success_Rate` - % yield > 50%
- `Reaction_Types` - Comma-separated list

### Reactions

- `Reaction_Type` - Reaction name
- `Num_Experiments` - Total experiments
- `Num_Reactant_Pairs` - Unique substrate combinations
- `Num_Catalysts` - Unique catalysts used
- `Avg_Yield`, `Success_Rate` - Performance
- `Top_Catalyst` - Most common catalyst
- `Top_Reactant_Pair` - Most common substrates

---

## Database Statistics

- **Total Experiments**: 66,308
- **Reaction Types**: 41
- **Unique Catalysts**: 229
- **Reactant Type Combinations**: 71

### Top Metals
- Pd: 45,205 experiments (68.2%)
- Cu: 5,478 experiments (8.3%)
- Ni: 328 experiments (0.5%)

### Top Reactions
1. C_N_Coupling: 24,012 experiments
2. Suzuki: 11,588 experiments
3. Arylation-acidic-C-H: 4,152 experiments

---

## Quick Examples

### Find Best Catalyst for Specific Substrates

```python
df = analytics.get_catalyst_stats(
    reaction_type="Suzuki",
    reactant_a_type="ArCl",
    reactant_b_type="ArB(OH)2"
)
df = df.sort_values('Success_Rate', ascending=False)
print(f"Best catalyst: {df.iloc[0]['Catalyst']}")
```

### Compare Pd vs Cu for C-N Coupling

```python
df = analytics.get_catalyst_stats(reaction_type="C-N")
pd_count = df[df['Metal'] == 'Pd']['Num_Experiments'].sum()
cu_count = df[df['Metal'] == 'Cu']['Num_Experiments'].sum()
print(f"Pd: {pd_count:,}, Cu: {cu_count:,}")
```

### Export High-Success Experiments

```bash
python -m chemtools.HTE.analytics_cli export high_success.csv --min-yield 80
```

---

## Help Commands

```bash
# Main help
python -m chemtools.HTE.analytics_cli --help

# Command-specific help
python -m chemtools.HTE.analytics_cli pairs --help
python -m chemtools.HTE.analytics_cli catalysts --help
python -m chemtools.HTE.analytics_cli reactions --help
python -m chemtools.HTE.analytics_cli metals --help
python -m chemtools.HTE.analytics_cli export --help
```

---

## Documentation

- **Full API Reference**: `docs/HTE_ANALYTICS.md`
- **Demo Script**: `python demo_hte_analytics.py`
- **Main README**: `chemtools/HTE/README.md`
- **Implementation Summary**: `HTE_ANALYTICS_IMPLEMENTATION.md`

---

## Cheat Sheet

| Task | Command |
|------|---------|
| List all Suzuki pairs | `python -m chemtools.HTE.analytics_cli pairs --reaction Suzuki` |
| Top Pd catalysts | `python -m chemtools.HTE.analytics_cli catalysts --catalyst Pd --top 10` |
| Reaction summary | `python -m chemtools.HTE.analytics_cli reactions --compact` |
| Metal distribution | `python -m chemtools.HTE.analytics_cli metals` |
| Export Pd Suzuki | `python -m chemtools.HTE.analytics_cli export out.csv --reaction Suzuki --catalyst Pd` |

---

**Last updated**: November 15, 2025
