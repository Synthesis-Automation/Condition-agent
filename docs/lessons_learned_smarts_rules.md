# Lessons Learned: Adding SMARTS Rules

**Date:** October 5, 2025  
**Context:** Real debugging session while adding new rules

---

## Lesson 1: SMARTS Alone Isn't Enough - Feature Requirements Matter!

### What Happened

We added a new rule for 2-heteroaryl halides:

```json
{
  "id": "SCDB-SUZ-2HETARYL-SPHOS",
  "rxn_smiles_min": "[n,s:1][c:2]-[Br,Cl:3].[B:4](-[O])-[c:5]>>[n,s:1][c:2]-[c:5]",
  "priority": 65,
  "feature_requirements": {
    "electrophile.is_2_hetaryl": true,  // ← Feature doesn't exist yet!
    "electrophile.lg_class": ["Br", "Cl"]
  }
}
```

### What We Expected

2-Bromopyridine should match this rule because:
- SMARTS `[n:1][c:2]-[Br:3]` matches perfectly ✓
- Priority 65 > 50 (higher than ortho rule) ✓

### What Actually Happened

**Entry was FILTERED OUT** before it could match! 

**Reason:** The feature `electrophile.is_2_hetaryl` doesn't exist in the feature detector code, so the feature requirement check fails and the entire entry is rejected.

### The Fix

**Option 1: Remove the non-existent feature requirement**

```json
{
  "feature_requirements": {
    "electrophile.lg_class": ["Br", "Cl"]
    // Removed "electrophile.is_2_hetaryl": true
  }
}
```

This works because SMARTS already filters for the correct structure!

**Option 2: Implement the feature**

Add detection code in `scdb_matcher/matcher.py` to calculate `is_2_hetaryl`.

**Option 3: Use existing features**

```json
{
  "feature_requirements": {
    "electrophile.lg_class": ["Br", "Cl"],
    "electrophile.ring_hetero_count": ">= 1"
  }
}
```

---

## Lesson 2: Why Did 2-Bromopyridine Match the Ortho Rule?

### The Ortho Rule

```json
{
  "id": "SCDB-SUZ-ARBR-ORTHO-XPhos",
  "feature_requirements": {
    "electrophile.lg_class": "Br",
    "electrophile.ortho_sub_count": ">= 1"
  }
}
```

### Why It Matched

In 2-bromopyridine (`Brc1ccccn1`):
- The nitrogen atom is in the ortho position relative to Br
- The feature detector counts **any non-hydrogen atom** at ortho positions
- Ring heteroatoms count as "ortho substituents"!

**This is actually chemically reasonable** - the nitrogen IS an ortho substituent that creates steric/electronic effects.

### The Key Insight

**Priority determines which rule wins:**
- If hetaryl rule had priority 65 and worked, it would match first
- But since it was filtered (missing feature), it fell through to ortho rule
- Ortho rule (priority 50) then matched correctly

---

## Lesson 3: Feature Requirements Are Gates

### How Matching Works

1. **Parse reaction** → Extract SMARTS fragments
2. **Calculate features** → Run feature detection code
3. **For each entry:**
   - Check `feature_requirements` first (GATE)
   - If any required feature is missing or doesn't match → **FILTER OUT**
   - If all features pass → Try SMARTS matching
   - If SMARTS match → **SUCCESS**

### The Gate Analogy

```
[Reaction] → [Feature Detection] → [Feature Gate] → [SMARTS Matching] → [Result]
                                         ↓
                                    Missing feature?
                                         ↓
                                      REJECT!
```

**Critical Rule:** Never reference a feature that doesn't exist!

---

## Lesson 4: SMARTS vs Features - When to Use What

### Use SMARTS Alone (No Features)

**When:** The structure is specific enough

```json
{
  "rxn_smiles_min": "[c]-OS(=O)(=O)C(F)(F)F",
  "feature_requirements": {}  // No features needed!
}
```

Triflate pattern is so specific, SMARTS alone is sufficient.

### Use SMARTS + Existing Features

**When:** You need simple filters (LG type, electronics)

```json
{
  "rxn_smiles_min": "[c:1]-[Cl:2]",
  "feature_requirements": {
    "electrophile.lg_class": "Cl",  // Already implemented
    "electrophile.electronics": "electron_poor"  // Already implemented
  }
}
```

### Use SMARTS + New Features

**When:** You need complex detection not expressible in SMARTS

**Step 1:** Implement the feature detector first!

```python
# In scdb_matcher/matcher.py
def _compile_feature_sets(normalization):
    # ... existing code ...
    
    # Add new feature
    for mol in normalization.masked_mols:
        is_indole = detect_indole_scaffold(mol)
        set_features["electrophile.is_indole"] = {"true"} if is_indole else set()
```

**Step 2:** Then use it in database

```json
{
  "feature_requirements": {
    "electrophile.is_indole": true  // Now it exists!
  }
}
```

---

## Lesson 5: Priority Is Your Friend

### Priority Ranges (Reminder)

```
0-10   → Defaults (catch-all, low specificity)
20-40  → General optimized conditions
45-55  → Substrate-specific (steric, positional)
60-70  → Special cases (heteroaryls, sensitive groups)
75-85  → Very specific (rare substrates)
```

### Example Conflict Resolution

**Scenario:** Both rules could match

```
Rule A: Generic ArBr + ortho sub → Priority 50
Rule B: Heteroaryl Br → Priority 65
```

**Result:** Rule B wins (higher priority)

### Best Practice

When adding a new rule, ask:
- "Could this overlap with existing rules?"
- "Is my substrate more specific or less specific?"
- "What priority ensures correct precedence?"

---

## Lesson 6: Testing Is Essential

### Quick Test Template

```python
from chemtools.scdb_matcher import loader, match

db = loader.load_db("data/conditionDB/suzuki_db.json")

# Test your reaction
result = match(db, "Your.SMILES.Here")

# Check what matched
print(f"Entry: {result.entry_id}")
print(f"Priority: {result.priority}")

# Inspect the trace
for eval in result.trace['evaluations']:
    if eval['id'] == 'YOUR-RULE-ID':
        print(f"Status: {eval.get('matched', False)}")
        print(f"Filtered: {eval.get('filtered', False)}")
        print(f"Matches: {eval.get('matches', [])}")
```

### What to Check

- ✅ Did my rule match?
- ✅ Did it match the right reactions?
- ✅ Did it correctly reject wrong reactions?
- ✅ Is the priority correct vs other rules?

---

## Summary: The Right Way to Add Rules

### ✅ Correct Process

1. **Write SMARTS** - Keep it simple, match the core structure
2. **Check existing features** - Use `grep "electrophile\." scdb_matcher/matcher.py`
3. **Only use features that exist** - Don't reference unimplemented features!
4. **Set appropriate priority** - Based on specificity
5. **Test immediately** - Run a quick test before committing
6. **Debug with trace** - Use `result.trace` to see what happened

### ❌ Common Mistakes

1. **Referencing non-existent features** → Entry gets filtered
2. **Priority too low** → Gets overridden by defaults
3. **Priority too high** → Overrides more specific rules
4. **SMARTS too specific** → Misses valid substrates
5. **No testing** → Silent failures in production

### 🎯 Golden Rule

**SMARTS matches structure, Features filter properties, Priority resolves conflicts.**

All three must work together!

---

## Fixed Version of 2-Hetaryl Rule

**Option A: Simple (SMARTS only)**

```json
{
  "id": "SCDB-SUZ-2HETARYL-SPHOS",
  "rxn_smiles_min": "[n,s:1][c:2]-[Br,Cl:3].[B:4](-[O])-[c:5]>>[n,s:1][c:2]-[c:5]",
  "priority": 65,
  "feature_requirements": {
    "electrophile.lg_class": ["Br", "Cl"]
  }
}
```

**Option B: With existing features**

```json
{
  "id": "SCDB-SUZ-2HETARYL-SPHOS",
  "rxn_smiles_min": "[n,s:1][c:2]-[Br,Cl:3].[B:4](-[O])-[c:5]>>[n,s:1][c:2]-[c:5]",
  "priority": 65,
  "feature_requirements": {
    "electrophile.lg_class": ["Br", "Cl"],
    "electrophile.ring_hetero_count": ">= 1"
  }
}
```

Both work! Option A is simpler since SMARTS already filters for heteroaryls.

---

**Key Takeaway:** Start simple, test early, iterate based on real behavior. The hybrid approach gives you flexibility, but you must understand how the pieces fit together!
