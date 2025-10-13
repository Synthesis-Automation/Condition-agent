# Reagent Database Validation Errors Report

**Date**: October 13, 2025  
**Total Entries**: 407  
**Valid Entries**: 403 (99.0%)  
**Invalid Entries**: 4 (1.0%)  
**Critical Errors**: 4

---

## Summary of Errors (🔴 Critical Issues Only)

### 1. **ACID** Role - 2 Errors

Both acid entries are missing the required `acidity` field.

| Index | ID | Compound Name | Error |
|-------|--------|---------------|-------|
| 0 | cas-109-63-7 | **BF3.Et2O** (Boron trifluoride etherate) | Missing `roles.acid.acidity` field |
| 1 | cas-7647-01-0 | **HCl** (Hydrochloric Acid) | Missing `roles.acid.acidity` field |

**Fix Required**: Add `"acidity": "<value>"` to the `roles.acid` section for both compounds.

---

### 2. **LIGAND** Role - 1 Error

One ligand entry is missing the required `smiles` field.

| Index | ID | Compound Name | Error |
|-------|--------|---------------|-------|
| 16 | ZEMZPXWZVTUONV-UHFFFAOYSA-N | **DavePhos** (2-(Dicyclohexylphosphino)-2'-(dimethylamino)biphenyl) | Missing `smiles` field |

**Fix Required**: Add SMILES structure to this ligand entry.

---

### 3. **SOLVENT** Role - 1 Error

One solvent entry is missing the required `polarity` field.

| Index | ID | Compound Name | Error |
|-------|--------|---------------|-------|
| 18 | cas-1076-43-3 | **C6D6** (Benzene-d6) | Missing `roles.solvent.polarity` field |

**Fix Required**: Add `"polarity": "<value>"` to the `roles.solvent` section.

---

## How to Fix These Errors

### Fix ACID entries (2 compounds)

Edit `data/reagents/taxonomy_acid.json` (or similar file):

```json
{
  "id": "cas-109-63-7",
  "name": "1,4-Dimethoxyethane",
  "roles": {
    "acid": {
      "families": ["..."],
      "acidity": "weak"  // ← ADD THIS FIELD
    }
  }
}
```

```json
{
  "id": "cas-7647-01-0",
  "name": "Hydrochloric acid",
  "roles": {
    "acid": {
      "families": ["..."],
      "acidity": "strong"  // ← ADD THIS FIELD
    }
  }
}
```

### Fix LIGAND entry (1 compound)

Edit `data/reagents/taxonomy_ligand.json`:

1. Find entry at index 16
2. Add the `smiles` field:

```json
{
  "id": "...",
  "name": "...",
  "smiles": "C1=CC=CC=C1P..."  // ← ADD THIS FIELD
}
```

### Fix SOLVENT entry (1 compound)

Edit `data/reagents/taxonomy_solvent.json`:

```json
{
  "id": "cas-1076-43-3",
  "name": "Cyclohexyl methyl carbonate",
  "roles": {
    "solvent": {
      "families": ["..."],
      "polarity": "polar_aprotic"  // ← ADD THIS FIELD (or appropriate value)
    }
  }
}
```

---

## Recommended Polarity Values for Solvents

- `"polar_protic"` - e.g., water, methanol, ethanol
- `"polar_aprotic"` - e.g., DMF, DMSO, acetonitrile
- `"nonpolar"` - e.g., hexane, toluene, benzene

---

## Recommended Acidity Values for Acids

- `"superacid"` - e.g., triflic acid
- `"strong"` - e.g., HCl, H2SO4
- `"moderate"` - e.g., acetic acid
- `"weak"` - e.g., phenol

---

## Quick Fix Commands

After editing the JSON files, re-run validation:

```bash
# Validate entire database
python scripts/validate_reagent_db.py --verbose

# Validate specific role
python scripts/validate_reagent_db.py --role acid --verbose
python scripts/validate_reagent_db.py --role ligand --verbose
python scripts/validate_reagent_db.py --role solvent --verbose
```

---

## Files to Edit

1. `data/reagents/taxonomy_acid.json` (or `acid.json`) - Fix 2 entries
2. `data/reagents/taxonomy_ligand.json` (or `ligand.json`) - Fix 1 entry
3. `data/reagents/taxonomy_solvent.json` (or `solvent.json`) - Fix 1 entry

---

## Note on Warnings (⚠️)

The validation also shows **230 warnings**, mostly about:
- Non-standard ID formats (e.g., `cas-123-45-6` instead of InChIKey)
- These are **not critical** and won't break functionality
- They're just recommendations for better data consistency

**Focus on the 4 critical errors first** (🔴 shown above).
