# Catalytic System Analysis Feature

## Overview

Added `get_common_catalytic_systems()` function to analyze complete catalytic systems (catalyst + ligand combinations) rather than individual components. This preserves the important relationship between catalysts and ligands used together.

## Problem Statement

Previously, the analytics module only extracted individual catalyst names from the `catalytic_system` array:

```python
# Old approach - loses pairing information
for cat in rxn.get('catalytic_system', []):
    name = cat.get('name')
    if name:
        catalysts.add(name)
```

**Problem**: We couldn't see which catalyst-ligand pairs are used together and how successful those combinations are.

## Solution

New `get_common_catalytic_systems()` function analyzes complete systems:

```python
# New approach - preserves complete systems
cat_sys = rxn.get('catalytic_system', [])
if cat_sys:
    component_names = [comp.get('name', '') for comp in cat_sys if comp.get('name')]
    if component_names:
        system_str = ' + '.join(component_names)
        # Track frequency and yields for complete system
```

## API

### Function Signature

```python
def get_common_catalytic_systems(
    family: str,
    top_n: int = 10,
    min_yield: Optional[float] = None
) -> List[Tuple[str, int, float]]
```

### Parameters

- **family**: Reaction family name (e.g., "C_N_Coupling_Pd", "Suzuki")
- **top_n**: Number of top catalytic systems to return (default: 10)
- **min_yield**: Optional minimum yield filter (0-100)

### Returns

List of tuples: `(system_string, count, avg_yield)`
- **system_string**: Components joined with " + " (e.g., "Pd(OAc)2 + RuPhos")
- **count**: Number of reactions using this system
- **avg_yield**: Average yield (%) or None

Sorted by count (descending)

## Usage Examples

### Basic Usage

```python
from chemtools import chem

# Get top 10 catalytic systems for C-N coupling
systems = chem.analytics.get_common_catalytic_systems("C_N_Coupling_Pd", top_n=10)

for system_str, count, avg_yield in systems:
    print(f"{system_str}: {count} reactions, {avg_yield:.1f}% avg yield")
```

**Output**:
```
Tris(dibenzylideneacetone)dipalladium(0) + RuPhos: 71 reactions, 73.4% avg yield
Tris(dibenzylideneacetone)dipalladium(0) + Tri-tert-butylphosphine tetrafluoroborate: 62 reactions, 78.2% avg yield
Palladium(II) acetate: 61 reactions, 71.9% avg yield
...
```

### High-Yield Systems

```python
# Find catalytic systems with ≥85% average yield
high_yield = chem.analytics.get_common_catalytic_systems(
    "C_N_Coupling_Pd",
    top_n=20,
    min_yield=85.0
)

for system_str, count, avg_yield in high_yield:
    print(f"{system_str}: {count} reactions, {avg_yield:.1f}%")
```

**Output**:
```
Palladium(II) acetate + 76189-56-5: 41 reactions, 90.9%
Palladium(II) acetate + 210169-40-7: 39 reactions, 91.4%
Tris(dibenzylideneacetone)dipalladium(0) + SPhos: 26 reactions, 86.3%
...
```

### Comparison: Systems vs Individual Components

```python
# Individual components (old way)
catalysts = chem.analytics.get_common_catalysts("C_N_Coupling_Pd", top_n=5)
# → Loses pairing: Pd(OAc)2 (474 reactions), RuPhos (127 reactions)

# Complete systems (new way)
systems = chem.analytics.get_common_catalytic_systems("C_N_Coupling_Pd", top_n=5)
# → Preserves pairing: Pd(OAc)2 + RuPhos (71 reactions, 73.4% yield)
```

## Key Insights from C_N_Coupling_Pd Dataset

### Most Common Systems
1. **Tris(dibenzylideneacetone)dipalladium(0) + RuPhos** (71 reactions, 73.4% yield)
2. **Tris(dibenzylideneacetone)dipalladium(0) + Tri-tert-butylphosphine tetrafluoroborate** (62 reactions, 78.2% yield)
3. **Palladium(II) acetate** (61 reactions, 71.9% yield) - single component

### Highest-Yield Systems
1. **Palladium(II) chloride + 2241598-31-0** (21 reactions, 94.4% yield)
2. **Palladium(II) acetate + 210169-40-7** (39 reactions, 91.4% yield)
3. **Palladium(II) acetate + 76189-56-5** (41 reactions, 90.9% yield)

### Multi-Component Systems
- **14 out of 15** top systems contain multiple components (catalyst + ligand)
- Shows importance of analyzing complete systems rather than individual parts

## Implementation Details

### Dataset Structure

Catalytic systems are stored as arrays in the JSONL dataset:

```json
{
  "catalytic_system": [
    {"name": "Palladium(II) acetate", "cas": "3375-31-3"},
    {"name": "RuPhos", "cas": "787618-22-8"}
  ]
}
```

### System String Format

Components are joined with ` + `:
- Single component: `"Palladium(II) acetate"`
- Two components: `"Palladium(II) acetate + RuPhos"`
- Three components: `"Pd(OAc)2 + Ligand1 + Ligand2"`

### Performance

- **C_N_Coupling_Pd** (1,343 reactions): ~0.05 seconds
- **Suzuki** (50,215 reactions): ~2.2 seconds
- Comparable to other analytics functions

## Testing

### Test Coverage

Added **Test 8** to `tests/test_analytics_module.py`:

```python
def test_8_catalytic_systems():
    """Test 8: Get common catalytic systems (catalyst + ligand combinations)."""
    systems = chem.analytics.get_common_catalytic_systems("C_N_Coupling_Pd", top_n=15)
    
    # Validation
    assert isinstance(systems, list)
    assert len(systems) > 0
    
    # Check multi-component systems
    multi_component = [s for s in systems if " + " in s[0]]
    assert len(multi_component) > 0  # Should find paired systems
```

### Test Results

```
Test 8: Common Catalytic Systems
  C_N_Coupling_Pd - Top 15 Catalytic Systems
  ⏱️  Time: 0.0471 seconds
  📊 Found 15 catalytic systems
  ℹ️  Found 14 multi-component systems
  ✅ Test passed: Found and validated 15 catalytic systems
```

## Demo Script

Run `demo_catalytic_systems.py` to see:

1. **Individual Component Analysis** (old way) - shows limitations
2. **Complete System Analysis** (new way) - shows advantages
3. **High-Yield Systems** - filtered by ≥85% average yield
4. **Suzuki Comparison** - larger dataset example

```bash
python demo_catalytic_systems.py
```

## Benefits

### For Chemists
- ✅ **Preserves catalyst-ligand relationships**
- ✅ **Shows complete systems as used in practice**
- ✅ **Identifies successful combinations**
- ✅ **Supports yield filtering for high-performance systems**

### For Applications
- 🧪 **HTE Plate Design**: Use proven catalyst/ligand pairs
- 📊 **Condition Recommendations**: Suggest complete systems
- 🔍 **Gap Analysis**: Identify untested combinations
- 📚 **Learning from Precedents**: See what works in literature

## Integration with ChemTools API

Accessible via `chem.analytics.get_common_catalytic_systems()`:

```python
from chemtools import chem

# Use alongside other analytics functions
stats = chem.analytics.get_stats("Suzuki")
systems = chem.analytics.get_common_catalytic_systems("Suzuki", top_n=10)
catalysts = chem.analytics.get_common_catalysts("Suzuki", top_n=10)
bases = chem.analytics.get_common_bases("Suzuki", top_n=10)
```

## Files Modified

### Core Implementation
- `chemtools/dataset_analytics.py`: Added `get_common_catalytic_systems()` function (48 lines)

### API Integration
- `chemtools/context.py`: Added `get_common_catalytic_systems()` method to `DatasetAnalyticsNamespace`

### Testing
- `tests/test_analytics_module.py`: Added Test 8 for catalytic systems (55 lines)

### Documentation & Demos
- `demo_catalytic_systems.py`: Interactive demo script (150 lines)
- `CATALYTIC_SYSTEM_ANALYSIS.md`: This documentation

## Next Steps

### Potential Enhancements
1. **System Similarity**: Find similar catalytic systems (e.g., same catalyst, different ligands)
2. **Component Analysis**: Analyze which ligands pair best with which catalysts
3. **Trend Analysis**: Track catalytic system popularity over time
4. **Recommendation Integration**: Suggest complete systems in condition recommendations

### Usage in Recommendations
```python
# Future: Integrate with recommendation engine
systems = chem.analytics.get_common_catalytic_systems(
    family="C_N_Coupling_Pd",
    top_n=5,
    min_yield=80.0
)

# Use these proven systems in plate design
for system_str, count, avg_yield in systems:
    # Parse system_str to extract components
    # Add to recommendation candidates
    pass
```

## Summary

The `get_common_catalytic_systems()` function addresses a critical gap in the analytics module by preserving and analyzing complete catalytic systems as they're used in practice. This enables better condition recommendations, more effective HTE plate design, and deeper insights into successful catalyst-ligand combinations.

**Key Achievement**: We now analyze chemistry as it's actually performed—with complete catalytic systems working together—rather than isolated components.
