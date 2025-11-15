# Catalyst Filtering Feature - Summary

## ✅ Feature Added: Catalyst Metal Type Filtering

The HTE recommendation system now supports **filtering by catalyst metal type** (e.g., palladium, copper, nickel).

## What Was Added

### 1. New Parameter: `catalyst_filter`

**Python API**:
```python
from chemtools.HTE import HTERecommender

recommender = HTERecommender()
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="CCN",
    catalyst_filter="Pd"  # Filter for palladium catalysts
)
```

**CLI**:
```bash
# Palladium catalysts only
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" --catalyst Pd

# Copper catalysts only
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "CCN" --catalyst copper

# Combine with other filters
python -m chemtools.HTE.cli -a "c1ccc(Cl)cc1" -b "c1ccc(B(O)O)cc1" \
    --reaction Suzuki --catalyst palladium -k 5
```

### 2. Supported Metals

The filter accepts both **metal symbols** and **full names**:

| Metal     | Symbol | Experiments | % of DB |
|-----------|--------|-------------|---------|
| Palladium | Pd     | 45,205      | 68.2%   |
| Copper    | Cu     | 5,907       | 8.9%    |
| Nickel    | Ni     | 619         | 0.9%    |
| Iridium   | Ir     | 120         | 0.2%    |
| Ruthenium | Ru     | -           | -       |
| Platinum  | Pt     | -           | -       |
| Gold      | Au     | -           | -       |
| Others    | Fe, Co, Zn, Ag, Rh | - | - |

## Implementation Details

### Code Changes

**File: `chemtools/HTE/recommender.py`**

1. Added `catalyst_filter` parameter to `recommend()` method
2. Added `_filter_by_catalyst()` helper method
3. Integrated filtering into recommendation pipeline

```python
def _filter_by_catalyst(self, df: pd.DataFrame, catalyst_filter: str) -> pd.DataFrame:
    """
    Filter dataframe by catalyst metal type.
    
    Maps common names to symbols (e.g., 'copper' -> 'Cu')
    Performs case-insensitive substring search in Catalyst column
    """
    # ... implementation
```

**File: `chemtools/HTE/cli.py`**

1. Added `--catalyst` CLI argument
2. Updated help examples
3. Passed filter to `recommend()` calls in both single and batch modes

### Testing

Created `demo_catalyst_filtering.py` demonstrating:
- ✅ Filtering C-N coupling by Pd vs Cu
- ✅ Filtering Suzuki reactions by Pd
- ✅ CLI usage examples

## Usage Examples

### Example 1: Suzuki with Palladium

```python
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(B(O)O)cc1",
    reaction_type_filter="Suzuki",
    catalyst_filter="Pd",
    top_k=5
)

# Result: 904 matching experiments (all Pd-based)
# Top: Pd-PEPPSI-IPent Cl, 97.5% avg yield
```

### Example 2: C-N Coupling Comparison

```python
# Palladium conditions
result_pd = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    catalyst_filter="palladium",
    top_k=3
)
# 1,997 matching experiments

# Copper conditions
result_cu = recommender.recommend(
    reactant_a_smiles="c1ccc(Cl)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    catalyst_filter="copper",
    top_k=3
)
# 0 matching experiments (Cu not commonly used for ArCl + ArNH2)
```

### Example 3: CLI Usage

```bash
# Suzuki with palladium
python -m chemtools.HTE.cli \
    -a "c1ccc(Cl)cc1" \
    -b "c1ccc(B(O)O)cc1" \
    --reaction Suzuki \
    --catalyst Pd \
    -k 3 \
    --compact

# Output:
# Reactant A: ArCl (ArX*)
# Reactant B: ArB(OH)2 (ArB*)
# Predicted: Suzuki (100%)
# Matches: 904 experiments
# 
# TOP RECOMMENDATION (Score: 70.1/100)
#   Catalyst: Pd-PEPPSI-IPent Cl o-picoline
#   Ligand: IPENT Cl
#   Base: Cs2CO3
#   Solvent: Brij 35
#   Success: 100.0% (2 exp, avg 97.5%)
```

## Database Coverage by Metal

Analysis of HTE_0.csv (66,308 experiments):

### Palladium (68.2% of database)
- **Most common**: Suzuki, C-N Coupling, Heck, Negishi
- **Top catalysts**: tBuBrettPhos Pd(allyl)OTf (2,395), dtbpfPdCl2 (2,158)
- **191 unique Pd catalyst formulations**

### Copper (8.9% of database)
- **Most common**: C-N Coupling (3,747), CO-Coupling (1,239), Sonogashira (269)
- **Top catalysts**: CuI (3,520), Cu(MeCN)4BF4 (771)
- **66 unique Cu catalyst formulations**
- **Often used with**: ArI, ArBr + various nucleophiles

### Nickel (0.9% of database)
- **619 experiments** with Ni catalysts
- Less common but present for specific transformations

### Other Metals (rare)
- Iridium: 120 experiments
- Rhodium: 0 experiments
- Ruthenium, Platinum, Gold: Present but minimal

## Important Notes

### 1. Reactant Type Mismatch

⚠️ **Database uses different type names than our classifier**

Example:
- Our classifier: `terminal-alkyne`
- Database: `alkyne`
- Our classifier: `indole_present`
- Database: `arom-NH`

**Impact**: Some combinations may not find matches even though data exists.

**Workaround**: The system is based on exact type matching. For best results:
- Use common reactant types (ArBr, ArCl, ArI, RNH2, ArNH2, ArB(OH)2)
- These have good coverage and consistent naming

### 2. Copper Catalyst Usage

Copper catalysts are prevalent in:
- **C-N Coupling**: ArBr + RNH2-a-branch (672 exp), ArBr + arom-NH (671 exp)
- **CO-Coupling**: ArI + ArOH, ArBr + ROH-a-branch
- **Sonogashira**: ArBr + alkyne (230 exp)

But not common for:
- **Suzuki** (primarily Pd)
- **Heck** (primarily Pd)
- **ArCl substrates** (Cu less reactive with ArCl)

### 3. Filter Behavior

- Case-insensitive substring match
- Accepts both symbols ("Pd", "Cu") and names ("palladium", "copper")
- Can be combined with `reaction_type_filter`
- Returns empty result if no matches found

## Files Modified

1. `chemtools/HTE/recommender.py` - Added filtering logic
2. `chemtools/HTE/cli.py` - Added CLI support
3. `demo_catalyst_filtering.py` - Demo script (new)
4. `check_catalyst_types.py` - Analysis script (new)

## Documentation Updates Needed

- [ ] Update `docs/HTE_RECOMMENDER.md` with catalyst filtering
- [ ] Update `chemtools/HTE/README.md` with examples
- [ ] Add to API reference section

## Testing Checklist

- [x] Python API with Pd filter
- [x] Python API with Cu filter
- [x] CLI with --catalyst Pd
- [x] CLI with --catalyst copper
- [x] Combination with --reaction filter
- [x] Edge case: metal with no matches
- [x] Demo script runs successfully

## Summary

✅ **Catalyst filtering is fully functional**

Users can now filter recommendations by:
- Metal type (Pd, Cu, Ni, etc.)
- Full metal name (palladium, copper, nickel, etc.)
- Case-insensitive matching
- Works in both Python API and CLI
- Combines with existing reaction type filters

**Most useful for**:
- Comparing Pd vs Cu conditions
- Finding metal-specific protocols
- Screening within specific catalyst families
- Optimizing based on metal availability/cost
