# C_N_Coupling Dataset - Complete Catalyst Analysis

## Executive Summary

**Total Reactions:** 15,967

The "missing" reactions you noticed are classified as:
1. **Organocatalysts** (33.03%) - Base-promoted or acid-catalyzed reactions
2. **Other/No catalyst** (10.54%) - Reactions with no identifiable catalyst

## Complete Catalyst Distribution

| Catalyst Class | Count | Percentage |
|----------------|-------|------------|
| **Organocatalyst** | 5,274 | 33.03% |
| **Cu (Copper)** | 4,887 | 30.61% |
| **Pd (Palladium)** | 2,669 | 16.72% |
| **Other/No catalyst** | 1,683 | 10.54% |
| **Ni (Nickel)** | 1,365 | 8.55% |
| **Ir (Iridium)** | 35 | 0.22% |
| **Ru (Ruthenium)** | 22 | 0.14% |
| **Zn (Zinc)** | 20 | 0.13% |
| **Fe (Iron)** | 11 | 0.07% |
| **Rh (Rhodium)** | 1 | 0.01% |
| **TOTAL** | **15,967** | **100.00%** ✓ |

---

## Grouped by Major Categories

| Category | Count | Percentage |
|----------|-------|------------|
| **Palladium (Pd)** | 2,669 | 16.72% |
| **Copper (Cu)** | 4,887 | 30.61% |
| **Nickel (Ni)** | 1,365 | 8.55% |
| **Other Metals** (Ir, Ru, Zn, Fe, Rh) | 89 | 0.56% |
| **Non-metal** (Organocatalyst + Other) | 6,957 | 43.57% |

---

## Why Was I Wrong Before?

My initial numbers were **incorrect** because I was counting only the reactions with **metal symbols in the core** (like "Pd", "Pd/XantPhos"). I missed:

1. **Organocatalyzed reactions** (5,274) - These have cores like:
   - `Base: DIPEA` (1,155 reactions)
   - `Base: K2CO3` (962 reactions)
   - `Base: Et3N` (628 reactions)
   - `Acid: HCl` (367 reactions)
   - `Base: Cs2CO3` (259 reactions)

2. **No catalyst reactions** (1,683) - These have empty `condition_core` fields

3. **Rare metals** (89 total):
   - Iridium (Ir): 35 reactions
   - Ruthenium (Ru): 22 reactions
   - Zinc (Zn): 20 reactions
   - Iron (Fe): 11 reactions
   - Rhodium (Rh): 1 reaction

---

## Detailed Breakdown by Catalyst

### 1. Organocatalyst (5,274 reactions - 33.03%)

**What it is:** Base-promoted or acid-catalyzed C-N coupling reactions without transition metals.

**Top condition cores:**
- `Base: DIPEA` - 1,155 reactions (N,N-Diisopropylethylamine)
- `Base: K2CO3` - 962 reactions (Potassium carbonate)
- `Base: Et3N` - 628 reactions (Triethylamine)
- `Acid: HCl` - 367 reactions (Hydrochloric acid)
- `Base: Cs2CO3` - 259 reactions (Cesium carbonate)

**Chemistry:** These are typically:
- SNAr reactions (nucleophilic aromatic substitution)
- Base-promoted amine alkylations
- Acid-catalyzed condensations

### 2. Copper (4,887 reactions - 30.61%)

**What it is:** Cu-catalyzed C-N coupling (Ullmann-Goldberg reaction family).

**Top condition cores:**
- `Cu` - 3,716 reactions (generic copper)
- `Cu/N,N'-Dimethylethylenediamine` - 305 reactions
- `Cu/phen` - 284 reactions (1,10-phenanthroline)
- `Cu/L-Proline` - 175 reactions
- `Cu/1,2-Diaminocyclohexane,trans` - 80 reactions

**Chemistry:** Classic Ullmann coupling and modern variations.

### 3. Palladium (2,669 reactions - 16.72%)

**What it is:** Pd-catalyzed C-N coupling (Buchwald-Hartwig reaction family).

**Top condition cores:**
- `Pd` - 555 reactions (generic palladium)
- `Pd/XantPhos` - 350 reactions
- `Pd/rac-BINAP` - 296 reactions
- `Pd/XPhos` - 243 reactions
- `Pd/Tri-tert-butylphosphinetetrafluoroborate` - 225 reactions

**Chemistry:** Modern cross-coupling with bulky phosphine ligands.

### 4. No Catalyst / Other (1,683 reactions - 10.54%)

**What it is:** Reactions with no identifiable catalyst core (empty `condition_core`).

**Possible reasons:**
- Direct nucleophilic substitution (no catalyst needed)
- Missing data in original dataset
- Novel/unusual reaction conditions

### 5. Nickel (1,365 reactions - 8.55%)

**What it is:** Ni-catalyzed C-N coupling (emerging alternative to Pd).

**Top condition cores:**
- `Ni` - 800 reactions (generic nickel)
- `Ni/dtbbpy` - 250 reactions (4,4'-di-tert-butyl-2,2'-bipyridine)
- `Ni/13-propanediyl` - 84 reactions
- `Ni/4,4'-Dimethyl-2,2'-bipyridine` - 84 reactions
- `Ni/24-difluorophenyl` - 41 reactions

**Chemistry:** Cost-effective alternative to Pd for specific substrates.

### 6. Rare Metals (89 reactions - 0.56%)

**Iridium (Ir)** - 35 reactions (0.22%)
- Used in photoredox C-N coupling
- Often combined with photoredox catalysis

**Ruthenium (Ru)** - 22 reactions (0.14%)
- Used in specialized C-N coupling conditions

**Zinc (Zn)** - 20 reactions (0.13%)
- Used as reducing agent or co-catalyst

**Iron (Fe)** - 11 reactions (0.07%)
- Emerging earth-abundant metal catalyst

**Rhodium (Rh)** - 1 reaction (0.01%)
- Very rare for C-N coupling

---

## Implications for Catalyst Filtering

### 1. All Catalyst Classes Are Supported

You can filter to any of these catalyst types:

```python
# Metal catalysts
relax={"catalyst_class": "Pd"}   # 2,669 reactions
relax={"catalyst_class": "Cu"}   # 4,887 reactions
relax={"catalyst_class": "Ni"}   # 1,365 reactions
relax={"catalyst_class": "Ir"}   # 35 reactions
relax={"catalyst_class": "Ru"}   # 22 reactions
relax={"catalyst_class": "Zn"}   # 20 reactions
relax={"catalyst_class": "Fe"}   # 11 reactions
relax={"catalyst_class": "Rh"}   # 1 reaction

# Non-metal catalysts
relax={"catalyst_class": "organo_catalyst"}  # 5,274 reactions
relax={"catalyst_class": "other"}            # 1,683 reactions
```

### 2. Practical Use Cases

**Laboratory with limited catalyst inventory:**
```python
# Only have Pd and Cu available
for metal in ["Pd", "Cu"]:
    result = recommend_from_reaction(rxn, relax={"catalyst_class": metal})
```

**Green chemistry preference (non-metal):**
```python
# Prefer base-promoted reactions (no transition metal)
result = recommend_from_reaction(rxn, relax={"catalyst_class": "organo_catalyst"})
```

**Cost optimization:**
```python
# Compare Cu (cheap) vs Pd (expensive) vs Ni (moderate)
costs = {"Cu": "cheap", "Ni": "moderate", "Pd": "expensive"}
for metal in ["Cu", "Ni", "Pd"]:
    result = recommend_from_reaction(rxn, relax={"catalyst_class": metal})
    print(f"{metal} ({costs[metal]}): {result['recommendation']['core']}")
```

**Emerging methods:**
```python
# Try iron catalysis (earth-abundant, sustainable)
result = recommend_from_reaction(rxn, relax={"catalyst_class": "Fe"})
```

---

## Updated Statistics Summary

### Corrected Numbers

| Category | Count | Percentage | Notes |
|----------|-------|------------|-------|
| Cu-catalyzed | 4,887 | 30.61% | **Largest group** |
| Organocatalyst | 5,274 | 33.03% | **Base/acid promoted** |
| Pd-catalyzed | 2,669 | 16.72% | Was 1,664 in my initial count |
| No catalyst | 1,683 | 10.54% | Empty condition_core |
| Ni-catalyzed | 1,365 | 8.55% | Was 800 in my initial count |
| Other metals | 89 | 0.56% | Ir, Ru, Zn, Fe, Rh |
| **TOTAL** | **15,967** | **100.00%** | ✓ Complete |

### Why the Discrepancy?

My initial analysis counted only reactions where the `condition_core` **started with** the metal symbol:
- `Pd`, `Pd/XantPhos` → Counted as Pd
- `Cu`, `Cu/phen` → Counted as Cu

But I **missed** reactions where the classifier had to parse the full catalyst system:
- `Base: DIPEA` → Classified as `organo_catalyst`
- Empty core → Classified as `other`
- Complex catalyst names → Parsed and classified

The **correct** approach uses `_row_catalyst_class()` which:
1. Parses `condition_core` field
2. Scans `catalytic_system` list for metal names
3. Falls back to enzyme/organo_catalyst/other classification

---

## Key Insights

1. **Organocatalysis is huge** (33%!) - Don't overlook base-promoted reactions
2. **Cu is the most common** (31%) - Not Pd! Ullmann coupling is very popular
3. **Pd is only 17%** - Important but not dominant
4. **Rare metals exist** - Fe, Ru, Ir for specialized applications
5. **No catalyst needed** (11%) - SNAr and direct substitution are viable

---

## Updated Documentation

I'll update the documentation to reflect these corrected numbers:
- ✅ `CATALYST_FILTERING_GUIDE.md` - Need to update statistics
- ✅ `CATALYST_FILTERING_QUICK_REF.txt` - Need to update percentages
- ✅ `scripts/analyze_catalyst_distribution.py` - Already correct!

---

## Testing

Run the analysis yourself:

```bash
python scripts/analyze_catalyst_distribution.py
```

This will show you the **complete, accurate** catalyst distribution with all 10 catalyst classes accounting for 100% of reactions.
