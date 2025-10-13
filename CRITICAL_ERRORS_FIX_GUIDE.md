# Critical Errors Summary - Reagent Database

**Generated**: October 13, 2025  
**Database**: `data/reagents/`  
**Total Errors**: 4 compounds

---

## 🔴 Compounds with Critical Errors

### 1. **BF3.Et2O** (Boron trifluoride etherate)
- **Role**: ACID
- **ID**: cas-109-63-7
- **Error**: Missing `roles.acid.acidity` field
- **Fix**: Add `"acidity": "lewis_acid"` to the acid role section
- **File**: `data/reagents/acid.json` (index 0)

### 2. **HCl** (Hydrochloric Acid)
- **Role**: ACID
- **ID**: cas-7647-01-0
- **Error**: Missing `roles.acid.acidity` field
- **Fix**: Add `"acidity": "strong"` to the acid role section
- **File**: `data/reagents/acid.json` (index 1)

### 3. **DavePhos**
- **Full Name**: 2-(Dicyclohexylphosphino)-2'-(dimethylamino)biphenyl
- **Role**: LIGAND
- **ID**: ZEMZPXWZVTUONV-UHFFFAOYSA-N
- **Error**: Missing `smiles` field
- **Fix**: Add SMILES structure for this phosphine ligand
- **File**: `data/reagents/ligand.json` (index 16)
- **Suggested SMILES**: `CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4`

### 4. **C6D6** (Benzene-d6)
- **Role**: SOLVENT
- **ID**: cas-1076-43-3
- **Error**: Missing `roles.solvent.polarity` field
- **Fix**: Add `"polarity": "nonpolar"` to the solvent role section
- **File**: `data/reagents/solvent.json` (index 18)

---

## 📝 How to Fix

### For ACID entries (BF3.Et2O and HCl)

Edit `data/reagents/acid.json`:

```json
{
  "id": "cas-109-63-7",
  "name": "Boron trifluoride etherate",
  "roles": {
    "acid": {
      "families": ["lewis_acids"],
      "acidity": "lewis_acid"  // ← ADD THIS
    }
  }
}
```

```json
{
  "id": "cas-7647-01-0",
  "name": "Hydrochloric Acid",
  "roles": {
    "acid": {
      "families": ["mineral_acids"],
      "acidity": "strong"  // ← ADD THIS
    }
  }
}
```

### For LIGAND entry (DavePhos)

Edit `data/reagents/ligand.json` at index 16:

```json
{
  "id": "ZEMZPXWZVTUONV-UHFFFAOYSA-N",
  "name": "2-(Dicyclohexylphosphino)-2'-(dimethylamino)biphenyl",
  "abbreviation": ["DavePhos"],
  "smiles": "CN(C)C1=CC=CC=C1C2=CC=CC=C2P(C3CCCCC3)C4CCCCC4",  // ← ADD THIS
  "roles": {
    "ligand": {
      // ... existing fields
    }
  }
}
```

### For SOLVENT entry (C6D6)

Edit `data/reagents/solvent.json` at index 18:

```json
{
  "id": "cas-1076-43-3",
  "name": "Benzene-d6",
  "abbreviation": ["C6D6"],
  "roles": {
    "solvent": {
      "families": ["aromatic_solvents"],
      "polarity": "nonpolar"  // ← ADD THIS
    }
  }
}
```

---

## ✅ Verify Fixes

After making the changes, re-run validation:

```bash
python scripts/validate_reagent_db.py
```

Expected output:
```
✅ NO CRITICAL ERRORS - All entries are valid!
```

---

## 📊 Current Database Status

- **Total Entries**: 407
- **Valid Entries**: 403 (99.0%)
- **Invalid Entries**: 4 (1.0%)
- **Warnings**: 230 (mostly non-standard ID formats - not critical)

After fixing these 4 errors, the database will be 100% valid! 🎉
