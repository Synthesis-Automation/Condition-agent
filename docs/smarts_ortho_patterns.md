# SMARTS Patterns for Ortho-Substituted Aryl Bromides

## The Challenge

Detecting ortho substitution purely with SMARTS is tricky because we need to:
1. Identify the ipso carbon (attached to Br)
2. Check that **neighboring** aromatic carbons have substituents (not just any carbon in the ring)

## Recommended SMARTS-Based Approach

### Option 1: Explicit Pattern (Most Reliable)
```smarts
[c:1]([c,C,N,O,S,F,Cl,Br,I,P])-[Br:2]
```

**What it matches:**
- Aromatic carbon `[c:1]` bonded to Br
- That carbon has a direct neighbor that is:
  - Aromatic carbon `c` (part of ring)
  - OR any carbon/heteroatom substituent `C,N,O,S,F,Cl,Br,I,P`

**Problem:** This matches the ortho carbons in benzene ring too!

### Option 2: Use Recursive SMARTS (Complex but Precise)
```smarts
[c:1]1[c;$(c([!H])[!H])][c][c][c][c]1-[Br:2]
```

**What it matches:**
- 6-membered aromatic ring with Br
- Position 2 (ortho) has a carbon with non-H substituents

**Problem:** Only checks one ortho position, misses 2,6-disubstituted cases

### Option 3: Combine with Feature Requirements (CURRENT SOLUTION)

**Best Practice:**
```json
{
  "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "feature_requirements": {
    "electrophile.lg_class": "Br",
    "electrophile.ortho_sub_count": ">= 1"
  }
}
```

**Advantages:**
- Simple, readable SMARTS for the core reaction
- Feature detection handles ortho counting correctly
- Easy to extend to other leaving groups

## Why SMARTS Alone is Limited Here

SMARTS patterns match **local structure** but have difficulty with:
- **Positional relationships** in rings (ortho/meta/para)
- **Counting** specific positions
- **Excluding** ring atoms while including substituents

For ortho detection specifically:
- Hard to distinguish "this carbon has an ortho substituent" from "this carbon IS an ortho substituent"
- Pattern `[c]([C])-Br` matches both the ipso carbon AND the ortho carbon

## When to Use SMARTS vs Features

### Use SMARTS when:
✅ Matching specific functional groups (triflates, vinyl halides, heteroaryls)
✅ Detecting core transformations (Suzuki, Buchwald, etc.)
✅ Identifying reactant classes (ArI, ArBr, ArCl)
✅ Matching specific ring systems (pyridines, indoles)

### Use Feature Detection when:
✅ Counting positions (ortho/meta/para substituents)
✅ Electronic effects (EDG/EWG counting)
✅ Steric congestion assessment
✅ Complex property calculations

## Recommended Hybrid Approach

For the **SCDB-SUZ-ARBR-ORTHO-XPhos** entry:

1. **SMARTS** for core reaction matching:
   - Aryl bromide: `[c:1]-[Br:2]`
   - Boron partner: `[B:3](-[O])-[c:4]`

2. **Features** for selectivity:
   - `electrophile.lg_class`: "Br"
   - `electrophile.ortho_sub_count`: ">= 1"

3. **Priority** for precedence:
   - `"priority": 50` (higher than defaults)

This gives you:
- ✅ Clear, readable SMARTS patterns
- ✅ Accurate ortho detection via features
- ✅ Easy to extend to new cases (just add more feature requirements)
- ✅ Better than pure SMARTS for positional logic

## Expanding to More Rules

### Example 1: Meta-Substituted Aryl Chlorides
```json
{
  "id": "SCDB-SUZ-ARCL-META-RuPhos",
  "rxn_smiles_min": "[c:1]-[Cl:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 55,
  "feature_requirements": {
    "electrophile.lg_class": "Cl",
    "electrophile.meta_sub_count": ">= 2"
  }
}
```

### Example 2: Electron-Poor Heterocycles
```json
{
  "id": "SCDB-SUZ-HETBR-EPoor-SPhos",
  "rxn_smiles_min": "[n,s,o:1]-[c:2]-[Br:3].[B:4](-[O])-[c:5]>>[n,s,o:1]-[c:2]-[c:5]",
  "priority": 60,
  "feature_requirements": {
    "electrophile.ring_hetero_count": ">= 1",
    "electrophile.electronics": "electron_poor"
  }
}
```

### Example 3: Bulky Nucleophiles
```json
{
  "id": "SCDB-SUZ-BULKY-NUC-XPhos",
  "rxn_smiles_min": "[c:1]-[Br:2].[B:3](-[O])-[c:4]>>[c:1]-[c:4]",
  "priority": 45,
  "feature_requirements": {
    "electrophile.lg_class": "Br",
    "nucleophile.ortho_sub_count": ">= 2",
    "nucleophile.is_hindered": true
  }
}
```

## Summary

**SMARTS are great for general matching**, but **features are better for nuanced selectivity**. 

The **hybrid approach** (simple SMARTS + feature requirements + priority) gives you:
1. Easy-to-read patterns
2. Accurate positional/electronic detection
3. Clear rule precedence
4. Maximum flexibility for expansion

You're right that SMARTS is more intuitive - use it for what it's good at (structure matching), and let features handle the complex logic (counting, positions, properties)!
