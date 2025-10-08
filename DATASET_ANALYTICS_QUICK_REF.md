# Dataset Analytics - Quick Reference

## Installation

No installation needed - part of ChemTools core.

```python
from chemtools import chem
```

## Quick Start

```python
# List available families
families = chem.analytics.get_all_families()

# Get basic stats
stats = chem.analytics.get_stats("Suzuki")
print(f"Total: {stats['total_reactions']:,} reactions")
print(f"Avg yield: {stats['yield_stats']['mean']:.1f}%")

# Get top catalysts
catalysts = chem.analytics.get_common_catalysts("C_N_Coupling_Pd", top_n=10)
for name, count, avg_yield in catalysts:
    print(f"{name}: {count} reactions, {avg_yield:.1f}% yield")

# Design 96-well plate
plate = chem.analytics.get_plate_recommendations(
    "C_N_Coupling_Pd", 
    n_conditions=96, 
    min_yield=60.0
)
```

## Common Use Cases

### 1. Get Top Reagents

```python
# Top 10 bases for Suzuki
bases = chem.analytics.get_common_bases("Suzuki", top_n=10)

# Top 10 solvents
solvents = chem.analytics.get_common_solvents("Suzuki", top_n=10)

# Top 10 ligands
ligands = chem.analytics.get_common_ligands("C_N_Coupling_Pd", top_n=10)
```

### 2. Filter by Yield

```python
# Only high-yield conditions (≥80%)
high_yield_cats = chem.analytics.get_common_catalysts(
    "C_N_Coupling_Pd",
    top_n=20,
    min_yield=80.0
)
```

### 3. Plate Design

```python
# 96-well plate optimized for diversity
plate = chem.analytics.get_plate_recommendations(
    family="Suzuki",
    n_conditions=96,
    min_yield=70.0,
    optimize_for='diversity'  # or 'yield' or 'frequency'
)

# Export to CSV
import csv
with open('plate.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=plate[0].keys())
    writer.writeheader()
    writer.writerows(plate)
```

### 4. Print Summary Report

```python
# Comprehensive analytics summary
chem.analytics.print_summary("Suzuki", top_n=10)
```

## API Functions

| Function | Purpose | Returns |
|----------|---------|---------|
| `get_all_families()` | List reaction families | `List[str]` |
| `get_stats(family)` | Dataset statistics | `Dict` |
| `get_common_catalysts(family, top_n, min_yield)` | Ranked catalysts | `List[Tuple[str, int, float]]` |
| `get_common_ligands(family, top_n, min_yield)` | Ranked ligands | `List[Tuple[str, int, float]]` |
| `get_common_bases(family, top_n, min_yield)` | Ranked bases | `List[Tuple[str, int, float]]` |
| `get_common_solvents(family, top_n, min_yield)` | Ranked solvents | `List[Tuple[str, int, float]]` |
| `get_common_reagents(family, role, top_n, min_yield)` | Ranked reagents | `List[Tuple[str, str, int, float]]` |
| `get_condition_cores(family, top_n, min_yield)` | Ranked condition cores | `List[Tuple[str, int, float]]` |
| `get_plate_recommendations(family, n, min_yield, optimize_for)` | HTE plate design | `List[Dict]` |
| `print_summary(family, top_n)` | Print report | `None` |

## Output Formats

### Reagent Rankings
Tuple: `(name, count, avg_yield)`
```python
('Potassium carbonate', 18560, 81.3)
```

### Reagent Rankings (with role)
Tuple: `(name, role, count, avg_yield)`
```python
('EDC·HCl', 'COUPLING_REAGENT', 2345, 89.2)
```

### Plate Recommendations
Dict with keys:
- `condition_id`: Unique ID
- `catalyst`, `ligand`, `base`, `solvent`: Reagent names (or None)
- `temperature_c`, `time_h`: Conditions (or None)
- `avg_yield`: Average yield %
- `frequency`: Number of precedents
- `score`: Recommendation score

## Available Datasets

- `Suzuki` (50,215 reactions)
- `Amide_formation` (41,427 reactions)
- `C_N_Coupling_Pd` (1,343 reactions)
- `C_N_Coupling_Cu` (5,552 reactions)
- `C_N_Coupling_Ni` (large dataset)

## Demo & Tests

```bash
# Run comprehensive demo
python demo_dataset_analytics.py

# Run quick test
python test_analytics_quick.py

# Run full test suite
pytest tests/test_dataset_analytics.py -v
```

## Tips

1. **Start with `print_summary()`** - Gets you an overview quickly
2. **Use `min_yield` filtering** - Focus on successful conditions
3. **Try different `optimize_for` strategies** - Diversity, yield, or frequency
4. **Export plates to CSV** - Easy integration with plate readers
5. **Compare datasets** - Use `get_stats()` across multiple families
