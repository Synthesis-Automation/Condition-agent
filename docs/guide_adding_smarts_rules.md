# Guide: Adding SMARTS-Based Rules to Suzuki Database

**Date:** October 5, 2025  
**Purpose:** Demonstrate how to easily expand the Suzuki condition database with new SMARTS-based rules

---

## Overview of the Hybrid Approach

Each rule combines:
1. **Simple SMARTS** - Easy to understand reaction pattern
2. **Feature requirements** - Precise selectivity filters  
3. **Priority** - Control which rule fires first

---

## New Rules Added (Examples)

### 1. Meta-Disubstituted Aryl Iodides → tBuXPhos

**Use Case:** 1,3,5-Trisubstituted aromatics, meta-hindered substrates

```json
{
  "id": "SCDB-SUZ-ARI-META-TBUXPHOS",
  "rxn_smiles_min": "[c:1]-[I:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 55,
  "feature_requirements": {
    "electrophile.lg_class": "I",
    "electrophile.meta_sub_count": ">= 2"
  }
}
```

**SMARTS Breakdown:**
- `[c:1]-[I:2]` → Aromatic C-I bond (any aryl iodide)
- `[B:3](-[O])-[c:4]` → Boron with oxygen (boronic acid or Bpin)

**Feature Logic:**
- Requires ArI (not ArBr or ArCl)
- Requires ≥2 meta substituents (detected by feature code)

**Why This Works:**
- ✅ SMARTS keeps it simple (just "aryl iodide + boron")
- ✅ Features handle the complex part (counting meta positions)
- ✅ Priority 55 ensures it beats defaults but loses to higher-priority specialty rules

---

### 2. 2-Heteroaryl Halides → SPhos

**Use Case:** 2-Bromopyridine, 2-chlorothiazole, electron-deficient heteroaryls

```json
{
  "id": "SCDB-SUZ-2HETARYL-SPHOS",
  "rxn_smiles_min": "[n,s:1][c:2]-[Br,Cl:3].[B:4](-[O])-[c:5]>>[n,s:1][c:2]-[c:5]",
  "priority": 65,
  "feature_requirements": {
    "electrophile.is_2_hetaryl": true,
    "electrophile.lg_class": ["Br", "Cl"]
  }
}
```

**SMARTS Breakdown:**
- `[n,s:1]` → Nitrogen or sulfur heteroatom
- `[c:2]-[Br,Cl:3]` → Adjacent aromatic carbon with Br/Cl
- Pattern matches 2-halopyridines, 2-halothiazoles, etc.

**Feature Logic:**
- `is_2_hetaryl` flag confirms C-X is alpha to heteroatom
- Accepts both Br and Cl

**Why This Works:**
- ✅ SMARTS explicitly matches heteroaryl structure
- ✅ Feature confirms it's the 2-position (not 3- or 4-)
- ✅ Priority 65 beats most rules (heteroaryls need special care)

---

### 3. Bulky Nucleophiles → XPhos

**Use Case:** 2,6-Dimethylphenylboronic acid, tetra-ortho biaryls

```json
{
  "id": "SCDB-SUZ-BULKY-NUC-XPHOS",
  "rxn_smiles_min": "[c:1]-[Br,I:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 48,
  "feature_requirements": {
    "electrophile.lg_class": ["Br", "I"],
    "nucleophile.ortho_sub_count": ">= 2"
  }
}
```

**SMARTS Breakdown:**
- Generic aryl halide + boron pattern
- No structural constraints in SMARTS!

**Feature Logic:**
- Filters for Br/I only (Cl too unreactive for this)
- Requires nucleophile has ≥2 ortho substituents

**Why This Works:**
- ✅ SMARTS is maximally simple
- ✅ ALL selectivity comes from features
- ✅ Priority 48 is lower than electrophile-ortho (50) but higher than defaults

---

## Template for Adding New Rules

### Step 1: Identify the Chemical Pattern

Ask yourself:
- What makes this substrate special? (steric, electronic, heteroatom?)
- Can SMARTS capture the core structure easily?
- What needs feature detection? (counting, positions, properties?)

### Step 2: Write Simple SMARTS

**Examples:**
```smarts
# Generic aryl halide
[c:1]-[Br:2]

# Specific heteroaryl
[n:1]cccc[c:2]-[Br:3]

# Vinyl halide
C=[C:1]-[Br:2]

# Triflate
[c:1]-[O:2]S(=O)(=O)C(F)(F)F
```

**Tips:**
- Keep it readable!
- Use atom mapping `:1`, `:2` for anchors
- Don't try to encode everything in SMARTS

### Step 3: Define Feature Requirements

**Common features available:**
```json
{
  "electrophile.lg_class": ["I", "Br", "Cl", "OTf"],
  "electrophile.ortho_sub_count": ">= 1",
  "electrophile.meta_sub_count": ">= 2",
  "electrophile.ring_hetero_count": ">= 1",
  "electrophile.electronics": "electron_poor",
  
  "nucleophile.ortho_sub_count": ">= 2",
  "nucleophile.is_hindered": true,
  
  "boron.partner_type": "Bpin",
  "boron.prone_to_protodeboronation": true
}
```

### Step 4: Set Priority

**Priority guidelines:**
- 0-10: Default conditions (catch-all)
- 20-40: General optimized conditions
- 45-55: Substrate-specific (steric, positional)
- 60-70: Special cases (heteroaryls, sensitive groups)
- 75-85: Very specific (rare substrates)

**Rule of thumb:** More specific = higher priority

### Step 5: Write the Entry

```json
{
  "id": "SCDB-SUZ-YOUR-RULE-NAME",
  "reaction_type": "Suzuki_Miyaura",
  "name": "Descriptive name for chemists",
  "rxn_smiles_min": "[your SMARTS here]",
  "priority": 50,
  "token_signature": [
    "Suzuki",
    "substrate_class",
    "Pd_source",
    "ligand",
    "base",
    "solvent"
  ],
  "conditions": {
    "pd_source": ["Pd2(dba)3"],
    "ligands": ["XPhos"],
    "boron_partner": ["aryl-Bpin 1.3 eq"],
    "base": ["K3PO4 2.0 eq"],
    "solvent": ["toluene/H2O (3:1)"],
    "temperature_C": [80, 95],
    "time_h": [6, 12],
    "loadings": {
      "Pd_mol%": [1.0, 2.5],
      "ligand_mol%": [2.0, 5.0],
      "base_eq": [2.0, 2.5]
    }
  },
  "env": {
    "anchors": {
      "electrophile_ring_c_mapno": 1,
      "leaving_group_mapno": 2
    },
    "features_from_smiles": {
      "electrophile.lg_class": "Br/I",
      "your.custom.feature": "description"
    },
    "feature_requirements": {
      "electrophile.lg_class": ["Br", "I"],
      "your.custom.feature": ">= 1"
    }
  },
  "evidence": {
    "provenance": "expert_default",
    "last_updated": "2025-10-05T03:00:00Z"
  },
  "notes": [
    "Explain why these conditions work",
    "Mention key mechanistic considerations"
  ]
}
```

---

## More Rule Ideas to Try

### 4. Aryl Chlorides with Strong EDG → Ligand 95
```json
{
  "id": "SCDB-SUZ-ARCL-STRONG-EDG-L95",
  "rxn_smiles_min": "[c:1]-[Cl:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 52,
  "feature_requirements": {
    "electrophile.lg_class": "Cl",
    "electrophile.electronics": "very_electron_rich"
  }
}
```

### 5. Indole/Pyrrole 2-Halides → RuPhos
```json
{
  "id": "SCDB-SUZ-INDOLE-2-HAL-RUPHOS",
  "rxn_smiles_min": "[nH,n]1[c:2][c:3]-[Br:4]>>[nH,n]1[c:2][c:3]-[c:5]",
  "priority": 68,
  "feature_requirements": {
    "electrophile.is_indole_2_position": true,
    "electrophile.lg_class": "Br"
  }
}
```

### 6. Tetra-Substituted Aryl Bromides → BrettPhos
```json
{
  "id": "SCDB-SUZ-TETRASUB-BRETTPHOS",
  "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 70,
  "feature_requirements": {
    "electrophile.lg_class": "Br",
    "electrophile.total_sub_count": ">= 4"
  }
}
```

### 7. Alpha-Halopyrimidines → XPhos Low Water
```json
{
  "id": "SCDB-SUZ-PYRIMIDINE-XPHOS",
  "rxn_smiles_min": "[n:1][c:2][n:3]~[c:4]-[Cl,Br:5]>>[n:1][c:2][n:3]~[c:4]-[c:6]",
  "priority": 72,
  "feature_requirements": {
    "electrophile.is_pyrimidine": true,
    "electrophile.lg_class": ["Cl", "Br"]
  }
}
```

---

## Best Practices

### ✅ DO:
- Keep SMARTS simple and readable
- Use features for complex logic (counting, properties)
- Set appropriate priority based on specificity
- Add clear notes explaining the chemistry
- Test with representative substrates

### ❌ DON'T:
- Overload SMARTS with too many constraints
- Forget to set priority (will default to 0)
- Duplicate existing rules unnecessarily
- Use features that aren't implemented yet
- Set priority too high (causes false matches)

---

## Testing Your New Rule

Create a simple test script:

```python
from chemtools.scdb_matcher import loader, match

db = loader.load_db("data/conditionDB/suzuki_db.json")

# Test reaction
rxn = "Brc1c(C)cccc1C.c1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccccc1"

result = match(db, rxn)

print(f"Matched Entry: {result.entry_id}")
print(f"Priority: {result.priority}")
print(f"Match Type: {result.match_type}")
```

---

## Summary

**The hybrid SMARTS + features approach gives you:**

1. **Simplicity** - SMARTS patterns are intuitive for chemists
2. **Power** - Features handle complex positional/electronic logic
3. **Flexibility** - Easy to add new rules without breaking existing ones
4. **Control** - Priority system ensures correct rule precedence
5. **Maintainability** - Clear separation between structure and properties

**Adding a new rule takes ~10 minutes:**
1. Write simple SMARTS (2 min)
2. Define feature requirements (3 min)
3. Set conditions and priority (3 min)
4. Add notes and metadata (2 min)

**You're ready to expand the database!** 🚀

---

**Files modified:**
- `data/conditionDB/suzuki_db.json` - Added 3 new example rules
- `docs/smarts_ortho_patterns.md` - Theory and limitations
- `docs/guide_adding_smarts_rules.md` - This practical guide

**Next steps:**
- Test the new rules with real substrates
- Add more rules for your specific use cases
- Consider implementing new features if needed
