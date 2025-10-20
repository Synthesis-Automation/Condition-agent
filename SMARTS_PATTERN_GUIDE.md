# SMARTS Pattern Matching Guide for Protocol Validation

## The Problem

When validating protocols, you may encounter errors like:
```
❌ Suzuki_protocols.json[1]
   ERROR: Reaction SMILES does not match any of the 1 SMARTS pattern(s)
   Patterns:
     - [c,n,o,s][#35,#53].[c,n,o,s]B(O[H])O[H]>>[c,n,o,s][c,n,o,s]
```

This happens when there's a mismatch between the **SMARTS pattern** and the **reaction SMILES**.

## Common Causes & Solutions

### 1. **Implicit vs Explicit Hydrogens**

**Problem**: Boronic acid representation mismatch

- **Reaction SMILES**: `B(O)O` (implicit hydrogens)
- **SMARTS Pattern**: `B(O[H])O[H]` (explicit hydrogens)

**Solution A - Match the SMILES (Recommended)**:
```json
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
]
```

**Solution B - Make SMARTS Flexible**:
```json
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B([OH,O])([OH,O])>>[c,n,o,s][c,n,o,s]"
]
```

**Solution C - Use Wildcard for H**:
```json
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B(O[H,*])O[H,*]>>[c,n,o,s][c,n,o,s]"
]
```

### 2. **Too Specific SMARTS Patterns**

**Problem**: Pattern too restrictive

- **Reaction**: `FC1=CC=C(C=C1)B(O)O` (fluorinated boronic acid)
- **SMARTS**: `c1ccccc1B(O)O` (only matches plain phenylboronic acid)

**Solution - Use Generic Aromatic**:
```json
"reaction_SMARTS": [
  "[c,n,o,s]B(O)O"  // Matches any aromatic boronic acid
]
```

### 3. **Missing Functional Group Tolerance**

**Problem**: Pattern doesn't account for functional groups

- **Reaction**: `BrC1=CC=C(C=C1)C(OC)=O` (ester group present)
- **SMARTS**: `[c]Br` (doesn't specify ester tolerance)

**Solution - Add Wildcards**:
```json
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
]
```
This allows ANY aromatic/heteroatom with the halide/boronic acid.

## SMARTS Pattern Best Practices

### ✅ Good SMARTS Patterns

1. **Generic Heteroatom Class**: `[c,n,o,s]` - Matches carbon, nitrogen, oxygen, sulfur in aromatic
2. **Halogen by Atomic Number**: `[#35,#53]` - Matches Br (35) or I (53)
3. **Flexible Hydrogen**: `O` instead of `O[H]` - Matches both implicit and explicit H
4. **Functional Group Neutral**: Don't overconstrain unless necessary

### ❌ Problematic SMARTS Patterns

1. **Explicit Hydrogens**: `B(O[H])O[H]` - Won't match `B(O)O`
2. **Too Specific**: `c1ccccc1` - Only matches benzene, not substituted aromatics
3. **Too Restrictive**: `[C][Br]` - Only matches aliphatic C-Br, not aromatic

## Quick Fix Script

If you have SMARTS validation errors, use the provided fix script:

```powershell
python fix_suzuki_smarts.py
```

This automatically:
1. Detects `B(O[H])O[H]` patterns
2. Replaces with `B(O)O` to match SMILES
3. Updates the JSON file
4. Reports changes made

## Manual Fix Process

### Step 1: Identify the Mismatch
```powershell
python -m chemtools.protocol.validate_protocols --verbose
```

Look for errors showing:
- Reaction SMILES
- SMARTS patterns that failed to match

### Step 2: Compare Reaction SMILES vs SMARTS

**Example Error**:
```
Reaction: BrC1=CC=C(C=C1)C(OC)=O.FC2=CC=C(C=C2)B(O)O>>...
SMARTS:   [c,n,o,s][#35,#53].[c,n,o,s]B(O[H])O[H]>>[c,n,o,s][c,n,o,s]
                                         ↑ ↑  ↑ ↑
                            Explicit H doesn't match implicit H in SMILES
```

### Step 3: Fix the Pattern

Edit `data/protocol_db/Suzuki_protocols.json`:

```json
{
  "reaction": {
    "reaction_smiles": "BrC1=CC=C(C=C1)C(OC)=O.FC2=CC=C(C=C2)B(O)O>>...",
    "reaction_SMARTS": [
      "[c,n,o,s][#35,#53].[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
    ]
  }
}
```

### Step 4: Rebuild Index and Validate

```powershell
# Rebuild protocol index
python rebuild_protocol_index.py

# Validate all protocols
python -m chemtools.protocol.validate_protocols
```

## Common SMARTS Pattern Library

### Suzuki-Miyaura Coupling

```json
// Aryl halide + Boronic acid → Biaryl
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
]

// Aryl halide + Boronic ester → Biaryl
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[c,n,o,s]B(O[C])O[C]>>[c,n,o,s][c,n,o,s]"
]

// Aryl triflate + Boronic acid → Biaryl
"reaction_SMARTS": [
  "[c,n,o,s]OS(=O)(C(F)(F)F)=O.[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
]
```

### Sonogashira Coupling

```json
// Aryl halide + Terminal alkyne → Aryl alkyne
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[C]#C>>[c,n,o,s]C#C"
]
```

### Buchwald-Hartwig Amination

```json
// Aryl halide + Amine → Aryl amine
"reaction_SMARTS": [
  "[c,n,o,s][#35,#53].[N,NH,NH2]>>[c,n,o,s][N]"
]
```

### C-N Coupling (Generic)

```json
// Heteroatom-halide + Nitrogen nucleophile → C-N bond
"reaction_SMARTS": [
  "[C,c][#35,#53].[N;H2,H1,H0]>>[C,c][N]"
]
```

## Testing SMARTS Patterns

You can test if a SMARTS pattern matches a reaction using RDKit:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# Test reaction
rxn_smiles = "BrC1=CC=C(C=C1)C(OC)=O.FC2=CC=C(C=C2)B(O)O>>FC3=CC=C(C4=CC=C(C=C4)C(OC)=O)C=C3"
rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)

# Test SMARTS pattern
pattern = "[c,n,o,s][#35,#53].[c,n,o,s]B(O)O>>[c,n,o,s][c,n,o,s]"
pattern_rxn = AllChem.ReactionFromSmarts(pattern)

# Check if reaction matches pattern
rxn_reactants = rxn.GetReactants()
pattern_reactants = pattern_rxn.GetReactants()

# Each pattern reactant should match at least one actual reactant
matches = all(
    any(rxn_r.HasSubstructMatch(pat_r) for rxn_r in rxn_reactants)
    for pat_r in pattern_reactants
)

print(f"Pattern matches: {matches}")
```

## Summary

**The Fix Applied**:
- Changed `B(O[H])O[H]` → `B(O)O` in 3 Suzuki protocols
- Result: 100% validation success (20/20 protocols valid)

**Key Takeaway**: 
Always ensure SMARTS patterns match the hydrogen representation (implicit vs explicit) used in your reaction SMILES. When in doubt, use implicit hydrogens (`B(O)O`) as they're more common in SMILES notation.

## Additional Resources

- **RDKit SMARTS Tutorial**: https://www.rdkit.org/docs/GettingStartedInPython.html#smarts-support
- **SMARTS Theory**: http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
- **Protocol Validation**: `python -m chemtools.protocol.validate_protocols --help`
