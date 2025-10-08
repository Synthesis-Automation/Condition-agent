# Get All Reagents by Type - Feature Summary

**Date**: January 2025  
**Status**: ✅ **IMPLEMENTED**

## Overview

Added new functions to `chemtools/reagent_lookup.py` to retrieve **all reagents of a specific type**, enabling users to explore the full reagent database, filter by properties, and build custom reagent lists.

---

## New Functions

### 1. `get_all_reagents_by_type(reagent_type: str) -> List[Dict[str, Any]]`

**Purpose**: Get all reagents from a specific database.

**Parameters**:
- `reagent_type` (str): Type of reagent database
  - Available types: `'ligand'`, `'base'`, `'solvent'`, `'metal_precursor'`, `'additive'`, `'acid'`, `'oxidant'`, `'reductant'`, etc.

**Returns**: List of reagent dictionaries with full details

**Example**:
```python
from chemtools.reagent_lookup import get_all_reagents_by_type

# Get all solvents
solvents = get_all_reagents_by_type('solvent')
print(f"Found {len(solvents)} solvents")

for s in solvents[:5]:
    print(f"  - {s['name']} (CAS: {s['cas']})")
```

**Output**:
```
Found 66 solvents
  - 1-Butanol (CAS: 71-36-3)
  - 1-Propanol (CAS: 71-23-8)
  - Ethanol (CAS: 64-17-5)
  - Isobutanol (CAS: 78-83-1)
  - Methanol (CAS: 67-56-1)
```

---

### 2. `count_reagents_by_type(reagent_type: str) -> int`

**Purpose**: Quick count of reagents in a database.

**Parameters**:
- `reagent_type` (str): Type of reagent database

**Returns**: Integer count

**Example**:
```python
from chemtools.reagent_lookup import count_reagents_by_type

count = count_reagents_by_type('ligand')
print(f"Available ligands: {count}")
# Output: Available ligands: 152
```

---

## Database Statistics

Based on current data (January 2025):

| Database Type | Count | Examples |
|--------------|-------|----------|
| **ligand** | 152 | BINAP, XPhos, PPh3, DavePhos |
| **solvent** | 66 | DMF, DMSO, Toluene, THF |
| **metal_precursor** | 49 | Pd(OAc)2, Pd2(dba)3, CuI |
| **base** | 47 | K3PO4, Cs2CO3, NaOtBu, DIPEA |
| **additive** | 39 | CuI, CuBr, NaI |
| **preformed_metal_catalyst** | 34 | Pd(PPh3)4, etc. |
| **condensation_agent** | 26 | DCC, EDC, etc. |
| **oxidant** | 14 | Oxone, H2O2, etc. |
| **reductant** | 9 | NaBH4, LiAlH4, etc. |
| **acid** | 2 | Acetic acid, etc. |
| **organo_catalyst** | 2 | Proline, etc. |
| **other_reagent** | 2 | Miscellaneous |
| **enzyme** | 1 | Enzymatic catalysts |

**Total**: 443+ reagents across 13 databases

---

## Usage Examples

### Example 1: List All Ligands

```python
from chemtools.reagent_lookup import get_all_reagents_by_type

ligands = get_all_reagents_by_type('ligand')

print(f"Total ligands: {len(ligands)}")
for lig in ligands[:10]:
    name = lig.get('name')
    abbr = lig.get('abbreviation', [])
    abbr_str = f" ({', '.join(abbr[:2])})" if abbr else ""
    print(f"  • {name}{abbr_str}")
```

**Output**:
```
Total ligands: 152
  • 1,10-Phenanthroline (phen)
  • BINAP (BINAP)
  • XPhos (XPhos)
  • DavePhos (DavePhos)
  • ...
```

---

### Example 2: Filter Solvents by Properties

```python
from chemtools.reagent_lookup import get_all_reagents_by_type

solvents = get_all_reagents_by_type('solvent')

# Filter by boiling point
high_bp = [s for s in solvents 
           if s.get('properties', {}).get('boiling_point', 0) > 150]

print(f"Solvents with bp > 150°C: {len(high_bp)}")
for s in high_bp[:5]:
    bp = s['properties']['boiling_point']
    print(f"  • {s['name']}: {bp}°C")
```

---

### Example 3: Get All Bases with pKa Values

```python
from chemtools.reagent_lookup import get_all_reagents_by_type

bases = get_all_reagents_by_type('base')

# Filter bases with pKa data
bases_with_pka = [b for b in bases if b.get('pka') is not None]

# Sort by pKa (strongest to weakest)
sorted_bases = sorted(bases_with_pka, key=lambda x: x.get('pka', 0), reverse=True)

print("Strongest bases (by pKa):")
for b in sorted_bases[:10]:
    print(f"  • {b['name']}: pKa {b['pka']}")
```

---

### Example 4: Build Custom Reagent Library

```python
from chemtools.reagent_lookup import get_all_reagents_by_type

# Build a custom ligand library for Buchwald reactions
all_ligands = get_all_reagents_by_type('ligand')

# Filter for phosphine ligands (example)
buchwald_ligands = [
    lig for lig in all_ligands 
    if any(name in lig.get('name', '').lower() 
           for name in ['phos', 'binap', 'xantphos'])
]

print(f"Buchwald-type ligands: {len(buchwald_ligands)}")
for lig in buchwald_ligands[:10]:
    print(f"  • {lig['name']}")
```

---

### Example 5: Export Reagent List to CSV

```python
import csv
from chemtools.reagent_lookup import get_all_reagents_by_type

solvents = get_all_reagents_by_type('solvent')

# Export to CSV
with open('solvents.csv', 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['Name', 'CAS', 'Abbreviation', 'SMILES', 'Boiling Point'])
    
    for s in solvents:
        abbr = ', '.join(s.get('abbreviation', []))
        bp = s.get('properties', {}).get('boiling_point', '')
        writer.writerow([
            s.get('name', ''),
            s.get('cas', ''),
            abbr,
            s.get('smiles', ''),
            bp
        ])

print(f"Exported {len(solvents)} solvents to solvents.csv")
```

---

### Example 6: Compare All Available Options

```python
from chemtools.reagent_lookup import get_all_reagent_types, count_reagents_by_type

print("Reagent Database Summary:")
print("-" * 50)

types = get_all_reagent_types()
for reagent_type in types:
    count = count_reagents_by_type(reagent_type)
    print(f"{reagent_type:30s} {count:3d} reagents")
```

**Output**:
```
Reagent Database Summary:
--------------------------------------------------
acid                           2 reagents
additive                       39 reagents
base                           47 reagents
condensation_agent             26 reagents
enzyme                         1 reagents
ligand                         152 reagents
metal_precursor                49 reagents
organo_catalyst                2 reagents
other_reagent                  2 reagents
oxidant                        14 reagents
preformed_metal_catalyst       34 reagents
reductant                      9 reagents
solvent                        66 reagents
```

---

## Reagent Dictionary Structure

Each reagent dictionary contains:

```python
{
    "name": str,              # Full chemical name
    "cas": str,               # CAS registry number
    "abbreviation": [str],    # List of abbreviations
    "smiles": str,            # SMILES notation (if available)
    "inchi_key": str,         # InChI key (if available)
    "aliases": [str],         # Alternative names
    "roles": dict,            # Reagent roles/applications
    "properties": {           # Physical/chemical properties
        "boiling_point": float,
        "melting_point": float,
        "density": float,
        "protic": bool,
        # ... other properties
    },
    "pka": float,             # pKa value (for bases)
    # ... other fields
}
```

---

## Integration with Existing API

The new functions complement the existing API:

```python
from chemtools.reagent_lookup import (
    find_reagent,              # Find ONE reagent by name
    get_all_reagents_by_type,  # Get ALL reagents of a type (NEW)
    count_reagents_by_type,    # Count reagents (NEW)
    get_all_reagent_types      # List available types
)

# OLD: Find specific reagent
dmf = find_reagent('DMF', 'solvent')

# NEW: Get all solvents
all_solvents = get_all_reagents_by_type('solvent')

# NEW: Count ligands
num_ligands = count_reagents_by_type('ligand')

# List all database types
types = get_all_reagent_types()
```

---

## Performance Notes

- **Caching**: Uses `@lru_cache` to avoid reloading databases
- **Speed**: First call loads from disk (~10-50ms), subsequent calls are instant (<1ms)
- **Memory**: All databases together use ~2-5MB RAM (negligible)

**Recommendation**: Safe to call repeatedly without performance concerns.

---

## Advanced Use Cases

### 1. Build Reaction Condition Matrix

```python
from chemtools.reagent_lookup import get_all_reagents_by_type

ligands = get_all_reagents_by_type('ligand')
bases = get_all_reagents_by_type('base')
solvents = get_all_reagents_by_type('solvent')

# Generate all combinations for screening
combinations = []
for lig in ligands[:10]:
    for base in bases[:5]:
        for solv in solvents[:5]:
            combinations.append({
                'ligand': lig['name'],
                'base': base['name'],
                'solvent': solv['name']
            })

print(f"Generated {len(combinations)} screening conditions")
```

---

### 2. Search Across All Databases

```python
from chemtools.reagent_lookup import get_all_reagent_types, get_all_reagents_by_type

def search_all_databases(search_term: str):
    """Search for a reagent across all databases."""
    results = []
    
    for db_type in get_all_reagent_types():
        reagents = get_all_reagents_by_type(db_type)
        
        for r in reagents:
            if search_term.lower() in r.get('name', '').lower():
                results.append({
                    'name': r['name'],
                    'type': db_type,
                    'cas': r.get('cas')
                })
    
    return results

# Example: Find all reagents with "methyl" in name
results = search_all_databases('methyl')
print(f"Found {len(results)} reagents containing 'methyl'")
```

---

### 3. Generate Reagent Catalog

```python
from chemtools.reagent_lookup import get_all_reagent_types, get_all_reagents_by_type

def generate_catalog():
    """Generate a complete reagent catalog."""
    catalog = {}
    
    for db_type in get_all_reagent_types():
        reagents = get_all_reagents_by_type(db_type)
        
        catalog[db_type] = {
            'count': len(reagents),
            'names': [r['name'] for r in reagents],
            'with_smiles': sum(1 for r in reagents if r.get('smiles')),
            'with_cas': sum(1 for r in reagents if r.get('cas'))
        }
    
    return catalog

catalog = generate_catalog()
for db_type, info in catalog.items():
    print(f"{db_type}: {info['count']} total, {info['with_smiles']} with SMILES")
```

---

## Testing

### Test Script: `test_get_all_reagents.py`

Run the provided test script to see all functionality:

```bash
python test_get_all_reagents.py
```

This demonstrates:
1. ✅ Listing all available reagent types
2. ✅ Getting all solvents with details
3. ✅ Getting all ligands
4. ✅ Getting all bases with pKa values
5. ✅ Filtering reagents by properties

---

## API Reference

### Function Signatures

```python
def get_all_reagents_by_type(reagent_type: str) -> List[Dict[str, Any]]:
    """
    Get all reagents of a specific type.
    
    Args:
        reagent_type: Type of reagent (e.g., 'ligand', 'base', 'solvent')
        
    Returns:
        List of all reagent dictionaries in that database
    """

def count_reagents_by_type(reagent_type: str) -> int:
    """
    Count reagents in a specific database.
    
    Args:
        reagent_type: Type of reagent database
        
    Returns:
        Number of reagents in that database
    """

def get_all_reagent_types() -> List[str]:
    """
    Get list of all available reagent database types.
    
    Returns:
        Sorted list of database type names
    """
```

---

## Comparison: Old vs New API

### OLD API (find specific reagent)
```python
from chemtools.reagent_lookup import find_reagent

# Find ONE specific reagent
dmf = find_reagent('DMF', 'solvent')
if dmf:
    print(dmf['name'])
```

**Limitation**: Only works if you know the exact reagent name.

---

### NEW API (explore all reagents)
```python
from chemtools.reagent_lookup import get_all_reagents_by_type

# Get ALL reagents and explore
solvents = get_all_reagents_by_type('solvent')

# Filter, sort, analyze
polar_solvents = [s for s in solvents if 'polar' in s.get('properties', {})]
by_name = sorted(solvents, key=lambda x: x['name'])
```

**Benefit**: Browse entire database, filter by any criteria, build custom lists.

---

## Conclusion

The new `get_all_reagents_by_type()` function enables:

✅ **Exploration**: Browse all available reagents by type  
✅ **Filtering**: Find reagents matching specific criteria  
✅ **Analysis**: Count, sort, and compare reagents  
✅ **Export**: Generate custom reagent lists/catalogs  
✅ **Screening**: Build condition matrices for experimentation  

**Total Coverage**: 443+ reagents across 13 databases, all accessible programmatically!

---

## Related Functions

- `find_reagent(name, type)` - Find specific reagent by name
- `enrich_reagent_info(name, type)` - Get enriched reagent details
- `load_reagent_database(type)` - Load database (internal/cached)
- `get_all_reagent_types()` - List available database types
- `get_all_reagents_by_type(type)` - **NEW**: Get all reagents of type
- `count_reagents_by_type(type)` - **NEW**: Count reagents

---

**Status**: ✅ **READY FOR USE**

See `test_get_all_reagents.py` for working examples!
