# Reagent Database Analytics

## Overview

The reagent analytics module provides statistical analysis and insights about your reagent database. Query counts, families, missing data, and more.

## Quick Start

### Python API

```python
from chemtools.reagent import (
    get_database_statistics,
    count_reagents_by_type,
    get_family_statistics,
    find_reagents_by_family,
    print_database_summary,
)

# How many ligands?
count = count_reagents_by_type('ligand')
print(f"Ligands: {count}")

# How many bases?
count = count_reagents_by_type('base')
print(f"Bases: {count}")

# Full statistics
stats = get_database_statistics()
print(f"Total: {stats['total_reagents']}")
print(f"By type: {stats['by_type']}")

# Pretty summary
print_database_summary()
```

### Command Line

```bash
# Database summary
python -m chemtools.reagent.analytics summary

# Specific role
python -m chemtools.reagent.analytics role ligand
python -m chemtools.reagent.analytics role base

# Missing data report
python -m chemtools.reagent.analytics missing
```

### Demo Script

```bash
python scripts/demo_reagent_analytics.py
```

## Available Functions

### Basic Counts

**`get_all_reagent_types() -> List[str]`**

Get list of all reagent types in database.

```python
from chemtools.reagent import get_all_reagent_types

types = get_all_reagent_types()
# Returns: ['acid', 'additive', 'base', 'ligand', 'solvent', ...]
```

**`count_reagents_by_type(reagent_type: str) -> int`**

Count reagents of a specific type.

```python
from chemtools.reagent import count_reagents_by_type

ligand_count = count_reagents_by_type('ligand')
base_count = count_reagents_by_type('base')
solvent_count = count_reagents_by_type('solvent')

print(f"Ligands: {ligand_count}")
print(f"Bases: {base_count}")
print(f"Solvents: {solvent_count}")
```

**`get_all_reagents_by_type(reagent_type: str) -> List[Dict]`**

Get all reagent entries of a specific type.

```python
from chemtools.reagent import get_all_reagents_by_type

all_ligands = get_all_reagents_by_type('ligand')
for ligand in all_ligands[:5]:
    print(f"  {ligand['name']} (CAS: {ligand.get('cas', 'N/A')})")
```

### Database Statistics

**`get_database_statistics(registry_dir=None) -> Dict`**

Get comprehensive statistics about the entire database.

```python
from chemtools.reagent import get_database_statistics

stats = get_database_statistics()

print(f"Total reagents: {stats['total_reagents']}")
print(f"With CAS: {stats['total_with_cas']} ({stats['percent_with_cas']:.1f}%)")
print(f"With InChIKey: {stats['total_with_inchikey']} ({stats['percent_with_inchikey']:.1f}%)")

# By type
for reagent_type, count in stats['by_type'].items():
    print(f"  {reagent_type}: {count}")

# Top families
for family, count in stats['top_families'][:10]:
    print(f"  {family}: {count}")
```

Returns:
```python
{
    "total_reagents": 449,
    "by_type": {"ligand": 155, "base": 47, ...},
    "total_with_cas": 400,
    "total_with_inchikey": 1,
    "total_with_smiles": 0,
    "percent_with_cas": 89.1,
    "families_by_role": {...},
    "top_families": [("trialkyl_triaryl_phosphines", 29), ...],
    "id_format_stats": {"inchikey": 4, "cas": 45, "other": 400},
    "multi_role_reagents": 0,
    ...
}
```

**`print_database_summary(registry_dir=None) -> None`**

Print formatted database summary.

```python
from chemtools.reagent import print_database_summary

print_database_summary()
```

Output:
```
======================================================================
REAGENT DATABASE SUMMARY
======================================================================

Registry: C:\...\data\reagents
Total reagents: 449

======================================================================
BY TYPE
======================================================================
  ligand              :  155 reagents
  solvent             :   67 reagents
  metal_precursor     :   50 reagents
  base                :   47 reagents
  ...
```

### Family Analysis

**`get_family_statistics(role: str, registry_dir=None) -> Dict`**

Get statistics about families for a specific role.

```python
from chemtools.reagent import get_family_statistics

ligand_stats = get_family_statistics('ligand')

print(f"Total ligands: {ligand_stats['total_reagents']}")
print(f"Total families: {ligand_stats['total_families']}")

for family_data in ligand_stats['families']:
    print(f"  {family_data['name']}: {family_data['count']}")
```

Returns:
```python
{
    "role": "ligand",
    "total_reagents": 155,
    "total_families": 32,
    "families": [
        {"name": "trialkyl_triaryl_phosphines", "count": 29},
        {"name": "dialkylbiaryl_phosphines", "count": 29},
        ...
    ]
}
```

**`find_reagents_by_family(role: str, family: str) -> List[Dict]`**

Find all reagents in a specific family.

```python
from chemtools.reagent import find_reagents_by_family

# Find all phosphine ligands
phosphines = find_reagents_by_family('ligand', 'trialkyl_triaryl_phosphines')

print(f"Found {len(phosphines)} phosphine ligands:")
for ligand in phosphines:
    name = ligand['name']
    cas = ligand.get('cas', 'N/A')
    print(f"  - {name} (CAS: {cas})")
```

**`print_role_summary(role: str, registry_dir=None) -> None`**

Print formatted summary for a specific role.

```python
from chemtools.reagent import print_role_summary

print_role_summary('ligand')
print_role_summary('base')
```

### Missing Data Analysis

**`get_missing_data_report(registry_dir=None) -> Dict`**

Generate report of missing data.

```python
from chemtools.reagent import get_missing_data_report

report = get_missing_data_report()

print(f"Total reagents: {report['total_reagents']}")

# Missing by field
for field, data in report['by_field'].items():
    print(f"{field}: {data['missing']} missing ({data['percent']:.1f}%)")

# Missing by type
for reagent_type, data in report['by_type'].items():
    print(f"{reagent_type}:")
    print(f"  Missing CAS: {data['missing_cas']}")
    print(f"  Missing InChIKey: {data['missing_inchikey']}")
```

**`print_missing_data_report(registry_dir=None) -> None`**

Print formatted missing data report.

```python
from chemtools.reagent import print_missing_data_report

print_missing_data_report()
```

## Common Use Cases

### 1. Count Reagents by Type

```python
from chemtools.reagent import count_reagents_by_type, get_all_reagent_types

types = get_all_reagent_types()

print("Reagent inventory:")
for reagent_type in types:
    count = count_reagents_by_type(reagent_type)
    print(f"  {reagent_type:20s}: {count:4d}")
```

### 2. Find All Ligands in a Family

```python
from chemtools.reagent import find_reagents_by_family

# Find all NHC ligands
nhc_ligands = find_reagents_by_family('ligand', 'nhc_imidazolylidenes')

print(f"Found {len(nhc_ligands)} NHC ligands:")
for ligand in nhc_ligands:
    print(f"  - {ligand['name']}")
    
    # Get role-specific data
    role_data = ligand['roles']['ligand']
    donors = role_data.get('donors', [])
    denticity = role_data.get('denticity')
    print(f"    Donors: {donors}, Denticity: {denticity}")
```

### 3. Analyze Base Database

```python
from chemtools.reagent import get_family_statistics

base_stats = get_family_statistics('base')

print(f"Base Database:")
print(f"  Total: {base_stats['total_reagents']}")
print(f"  Families: {base_stats['total_families']}")

print("\nTop 5 base families:")
for family_data in base_stats['families'][:5]:
    name = family_data['name']
    count = family_data['count']
    percent = count / base_stats['total_reagents'] * 100
    print(f"  {name:30s}: {count:2d} ({percent:4.1f}%)")
```

### 4. Data Quality Check

```python
from chemtools.reagent import get_database_statistics, get_missing_data_report

stats = get_database_statistics()
report = get_missing_data_report()

print("Data Quality Report:")
print(f"  Total reagents: {stats['total_reagents']}")
print(f"  With CAS: {stats['percent_with_cas']:.1f}%")
print(f"  With InChIKey: {stats['percent_with_inchikey']:.1f}%")
print(f"  With SMILES: {stats['percent_with_smiles']:.1f}%")

# Find types with most missing data
print("\nMissing InChIKeys by type:")
for reagent_type, data in report['by_type'].items():
    if data['missing_inchikey'] > 0:
        print(f"  {reagent_type}: {data['missing_inchikey']} missing")
```

### 5. Export Family List for a Role

```python
from chemtools.reagent import get_family_statistics
import json

ligand_stats = get_family_statistics('ligand')

# Export to JSON
output = {
    "role": "ligand",
    "total": ligand_stats['total_reagents'],
    "families": ligand_stats['families']
}

with open('ligand_families.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"Exported {len(ligand_stats['families'])} families to ligand_families.json")
```

## Current Database Status

Based on latest analysis:

- **Total reagents**: 449
- **Ligands**: 155 (32 families)
- **Bases**: 47 (16 families)
- **Solvents**: 67 (23 families)
- **Metal precursors**: 50 (23 families)
- **Additives**: 39 (12 families)

**Data completeness**:
- CAS numbers: 89.1%
- InChIKeys: 0.2% (needs improvement)
- SMILES: 0.0% (needs improvement)
- Abbreviations: 93.8%

**Top ligand families**:
1. Trialkyl/triaryl phosphines: 29
2. Dialkylbiaryl phosphines: 29
3. Bidentate diphosphines: 17
4. Diamines: 8
5. NHC imidazolylidenes: 8

## Files

- **Implementation**: `chemtools/reagent/analytics.py`
- **Demo script**: `scripts/demo_reagent_analytics.py`
- **This guide**: `chemtools/reagent/ANALYTICS_README.md`

## See Also

- `VALIDATOR_README.md` - Database validation
- `lookup.py` - Runtime reagent lookup
- `taxonomy_store.py` - Reagent taxonomy management
