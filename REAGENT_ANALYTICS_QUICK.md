# Quick Answer: Reagent Analytics

## Your Question: "How many ligands, bases, etc.?"

### Quick Answer

```python
from chemtools.reagent import count_reagents_by_type

# How many ligands?
ligands = count_reagents_by_type('ligand')
print(f"Ligands: {ligands}")
# Output: Ligands: 155

# How many bases?
bases = count_reagents_by_type('base')
print(f"Bases: {bases}")
# Output: Bases: 47
```

### All Counts at Once

```python
from chemtools.reagent import get_all_reagent_types, count_reagents_by_type

types = get_all_reagent_types()

for reagent_type in types:
    count = count_reagents_by_type(reagent_type)
    print(f"{reagent_type:25s}: {count:3d}")
```

**Current Database:**
```
acid                     :   2
additive                 :  39
base                     :  47
condensation_agent       :  26
enzyme                   :   1
ligand                   : 155  ← Answer!
metal_precursor          :  50
organo_catalyst          :   2
other_reagent            :   3
oxidant                  :  14
preformed_metal_catalyst :  34
reductant                :   9
solvent                  :  67

Total                    : 449
```

## Yes, We Have Analytics Tools! ✅

### What's Available

**Basic Counts:**
- ✅ `get_all_reagent_types()` - List all types
- ✅ `count_reagents_by_type(type)` - Count by type
- ✅ `get_all_reagents_by_type(type)` - Get all entries

**Advanced Analytics:**
- ✅ `get_database_statistics()` - Full database stats
- ✅ `get_family_statistics(role)` - Family breakdown
- ✅ `find_reagents_by_family(role, family)` - Search by family
- ✅ `get_missing_data_report()` - Data quality report

**Pretty Print:**
- ✅ `print_database_summary()` - Formatted summary
- ✅ `print_role_summary(role)` - Role-specific summary
- ✅ `print_missing_data_report()` - Missing data report

### Command Line

```bash
# Quick summary
python -m chemtools.reagent.analytics summary

# Specific role (e.g., ligands)
python -m chemtools.reagent.analytics role ligand

# Missing data
python -m chemtools.reagent.analytics missing
```

### Demo Script

```bash
python scripts/demo_reagent_analytics.py
```

## More Examples

### Find All Phosphine Ligands

```python
from chemtools.reagent import find_reagents_by_family

phosphines = find_reagents_by_family('ligand', 'trialkyl_triaryl_phosphines')
print(f"Found {len(phosphines)} phosphine ligands")
# Output: Found 29 phosphine ligands
```

### Ligand Family Distribution

```python
from chemtools.reagent import get_family_statistics

stats = get_family_statistics('ligand')
print(f"Total ligands: {stats['total_reagents']}")
print(f"Total families: {stats['total_families']}")

# Top families
for family in stats['families'][:5]:
    print(f"  {family['name']}: {family['count']}")
```

Output:
```
Total ligands: 155
Total families: 32
  trialkyl_triaryl_phosphines: 29
  dialkylbiaryl_phosphines: 29
  bidentate_diphosphines: 17
  diamines: 8
  nhc_imidazolylidenes: 8
```

### Full Database Summary

```python
from chemtools.reagent import print_database_summary

print_database_summary()
```

## Documentation

- **Full API**: `chemtools/reagent/ANALYTICS_README.md`
- **Demo script**: `scripts/demo_reagent_analytics.py`
- **Implementation**: `ANALYTICS_IMPLEMENTATION.md`
- **This guide**: `REAGENT_ANALYTICS_QUICK.md`

## Summary

✅ **Yes, we have comprehensive analytics tools!**
- Count reagents by type
- Analyze families
- Find reagents by family
- Get database statistics
- Data quality reports
- Command line interface
- Pretty-printed summaries

**Most common use:**
```python
from chemtools.reagent import count_reagents_by_type

count = count_reagents_by_type('ligand')  # 155
```
