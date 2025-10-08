# Ortho Substitution Detection in SCDB Matcher

**Date:** October 5, 2025  
**Entry:** SCDB-SUZ-ARBR-ORTHO-XPhos

---

## Overview

Yes, the SCDB matcher **does have methods to detect ortho substitution** in starting materials. The detection happens at two levels:

1. **SMARTS Pattern Matching** (in `rxn_smiles_min`)
2. **Feature Detection** (in `env.features_from_smiles` and matcher code)

---

## Detection Methods

### 1. SMARTS-Based Detection (SCDB-SUZ-ARBR-ORTHO-XPhos)

**Pattern:** `[c:1]([cH0]):[cH0]-[Br:2]`

**What it matches:**
- `[c:1]` - aromatic carbon (ipso position)
- `([cH0])` - first ortho position must have 0 hydrogens (substituted)
- `:[cH0]` - second ortho position must have 0 hydrogens (substituted)
- `-[Br:2]` - bromide leaving group

**Requirement:** BOTH ortho positions must be substituted (no H atoms)

**Example matches:**
```
✓ Brc1cc(C)ccc1C     (2-Bromo-4,5-dimethylbenzene)
✗ Brc1ccccc1C        (2-Bromotoluene - only ONE ortho substituted)
✗ Brc1c(C)cccc1C     (2,6-Dimethylbromobenzene - pattern doesn't match this geometry)
✗ Brc1ccccc1         (Bromobenzene - no ortho substituents)
```

**Limitation:** The current SMARTS pattern is very strict and only matches substrates with both ortho positions fully substituted. Many ortho-hindered substrates are not matched.

---

### 2. Feature-Based Detection (Matcher Code)

**Location:** `scdb_matcher/matcher.py` lines 268-296

**Algorithm:**
```python
def _compile_feature_sets(normalization):
    ortho_sub_count = 0.0
    
    # For each halogen atom (Cl, Br, I):
    for halogen in halogen_atoms:
        for neighbor in halogen.GetNeighbors():
            # Find aromatic carbon attached to halogen
            if neighbor.GetIsAromatic():
                # Find the 6-membered ring
                ring = get_ring_containing(neighbor)
                
                # Check ortho positions
                ortho = 0
                for nbr in neighbor.GetNeighbors():
                    if nbr != halogen and nbr.GetIsAromatic():
                        # Count if ortho carbon has no H (substituted)
                        if nbr.GetTotalNumHs() == 0:
                            ortho += 1
                
                ortho_sub_count = max(ortho_sub_count, float(ortho))
```

**Output Feature:** `electrophile.ortho_sub_count`

**Values:**
- `0.0` = No ortho substituents
- `1.0` = One ortho position is substituted
- `2.0` = Both ortho positions are substituted

---

### 3. Alternative Detection (Ullmann Featurizer)

**Location:** `chemtools/featurizers/ullmann.py` lines 95-136

**Algorithm:**
```python
def detect_ortho_substitution(mol, ipso_atom):
    ortho_count = 0
    
    # Find 6-membered aromatic ring containing ipso carbon
    rings = [r for r in ri.AtomRings() if atom_idx in r and len(r) == 6]
    
    if rings:
        ring = rings[0]
        pos = ring.index(atom_idx)
        
        # Ortho positions are pos-1 and pos+1 (modulo 6)
        ortho_atoms = [
            mol.GetAtomWithIdx(ring[(pos - 1) % 6]), 
            mol.GetAtomWithIdx(ring[(pos + 1) % 6])
        ]
        
        # Check if ortho position has substituents
        def is_substituted(ar_atom):
            # Count non-ring heavy neighbors
            for nb in ar_atom.GetNeighbors():
                if nb.GetIdx() not in ring and nb.GetAtomicNum() > 1:
                    return True
            return False
        
        ortho_count = sum(1 for a in ortho_atoms if is_substituted(a))
    
    return ortho_count  # 0, 1, or 2
```

**Output Feature:** `ortho_count` (string: "0", "1", or "2+")

---

## Entry Configuration (SCDB-SUZ-ARBR-ORTHO-XPhos)

```json
{
  "id": "SCDB-SUZ-ARBR-ORTHO-XPhos",
  "rxn_smiles_min": "[c:1]([cH0]):[cH0]-[Br:2].[B:3](Oc)(Oc)-[c:4]>>[c:1]-[c:4]",
  "env": {
    "features_from_smiles": {
      "electrophile.ortho_sub_count": ">=1 indicates steric congestion",
      "electrophile.lg_class": "Br"
    },
    "feature_requirements": {
      "electrophile.lg_class": "Br",
      "electrophile.ortho_sub_count": ">= 1"
    }
  }
}
```

**Matching Logic:**
1. SMARTS pattern must match (very strict - requires both ortho positions substituted)
2. OR feature requirements must be satisfied:
   - `electrophile.lg_class` = "Br"
   - `electrophile.ortho_sub_count` >= 1

---

## Why Ortho-Substituted Reactions Match Default Conditions

**Problem:** In the test run, all ortho-substituted substrates matched `SCDB-SUZ-DEFAULT-ArI_ArBr` instead of `SCDB-SUZ-ARBR-ORTHO-XPhos`.

**Reasons:**

1. **SMARTS Pattern Too Strict:**
   - Pattern `[c:1]([cH0]):[cH0]-[Br:2]` requires BOTH ortho positions to be substituted with specific geometry
   - Most ortho-hindered substrates have only ONE ortho substituent
   - Pattern fails on 2,6-disubstituted benzenes due to geometry

2. **Feature Requirements Not Enforced:**
   - The matcher compiles the feature `electrophile.ortho_sub_count`
   - But the entry requires SMARTS match OR feature match
   - Since SMARTS fails, features should trigger, but default entries have priority 0 and match first

3. **Priority System:**
   - SCDB-SUZ-ARBR-ORTHO-XPhos has no explicit priority (defaults to 0?)
   - SCDB-SUZ-DEFAULT-ArI_ArBr also has priority 0
   - Default entries may be checked first

---

## Recommendations

### Option 1: Fix the SMARTS Pattern

Make it less strict to match substrates with at least one ortho substituent:

**Current (too strict):**
```
[c:1]([cH0]):[cH0]-[Br:2]
```

**Proposed (matches ≥1 ortho):**
```
[c:1]([cH0,cH1])-[Br:2]
```
This matches any aromatic carbon with at least one ortho position having reduced H count.

### Option 2: Rely Only on Feature Detection

Remove the strict SMARTS pattern and rely on `feature_requirements`:

```json
{
  "rxn_smiles_min": "[c:1]-[Br:2].[B:3]>>[c:1]-[B:3]",
  "feature_requirements": {
    "electrophile.lg_class": "Br",
    "electrophile.ortho_sub_count": ">= 1"
  }
}
```

### Option 3: Increase Entry Priority

Set a higher priority for the ortho-hindered entry:

```json
{
  "id": "SCDB-SUZ-ARBR-ORTHO-XPhos",
  "priority": 50,
  ...
}
```

This ensures it's checked before default conditions.

---

## Test Results

### SMARTS Pattern Matching

| Substrate | SMARTS Matches | Actual Ortho Count |
|-----------|----------------|-------------------|
| Brc1ccccc1 (bromobenzene) | ✗ | 0 |
| Brc1ccccc1C (2-bromotoluene) | ✗ | 1 |
| Brc1cc(C)ccc1C (2-Br-4,5-dimethyl) | ✓ | 2 |
| Brc1c(C)cccc1C (2,6-dimethyl) | ✗ | 2 |

**Conclusion:** SMARTS pattern only matches 1 out of 3 ortho-substituted substrates.

---

## Implementation Details

### Files Involved

1. **scdb_matcher/matcher.py**
   - Lines 268-309: `_compile_feature_sets()` - detects ortho substitution
   - Feature: `electrophile.ortho_sub_count` (float: 0.0, 1.0, 2.0)

2. **chemtools/featurizers/ullmann.py**
   - Lines 95-136: Ortho detection for Ullmann C-N coupling
   - Feature: `ortho_count` (string: "0", "1", "2+")

3. **data/conditionDB/suzuki_db.json**
   - Line 176+: SCDB-SUZ-ARBR-ORTHO-XPhos entry definition

---

## Conclusion

✅ **YES**, the SCDB matcher has methods to detect ortho substitution:

1. **Feature detection** in `scdb_matcher/matcher.py` calculates `electrophile.ortho_sub_count`
2. **Ullmann featurizer** in `chemtools/featurizers/ullmann.py` provides `ortho_count`
3. **SMARTS matching** in database entries can explicitly match ortho-substituted patterns

⚠️ **HOWEVER**, the current SCDB-SUZ-ARBR-ORTHO-XPhos entry has a **very strict SMARTS pattern** that fails to match many ortho-hindered substrates.

### To fix:
1. Relax the SMARTS pattern, OR
2. Remove SMARTS and rely on feature requirements, OR
3. Increase entry priority to ensure it's checked before defaults

The feature detection code is robust and correctly identifies ortho substituents - the issue is in how the database entry is configured to use this information.

---

**Test Script:** `scripts/test_ortho_detection.py`  
**Documentation Generated:** October 5, 2025
