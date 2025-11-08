# reaction_SMARTS vs applies_if: Validation Strategy Analysis

## Current State

### Two Validation Mechanisms Exist

1. **applies_if** (Rules only) - Functional Group Detection
   - Location: `data/rule_db_v2/*.json`
   - Checks: Presence of functional groups in substrates
   - Example: `{"all": ["carboxylic_acid_present"], "any": ["primary_amine_present", "aniline_present"]}`
   - **Status**: ✅ **Currently implemented** in `UnifiedRecommender._validate_rule_applicability()`
   - Pros: Easy to write, human-readable, covers common cases
   - Cons: Limited specificity, can't distinguish regiochemistry or stereochemistry

2. **reaction_SMARTS** (Protocols only) - Reaction Transformation Matching
   - Location: `data/protocol_db_v2/*.json`
   - Checks: **Exact reaction transformation pattern** matches query
   - Example: `["IC=C.CC(O)(C#N)C>>N#CC=C"]` (alkenyl iodide + acetone cyanohydrin → nitrile)
   - **Status**: ❌ **Currently NOT used for validation** (only schema-validated)
   - Pros: **Much more accurate** - encodes exact transformation, stereochemistry, regiochemistry
   - Cons: Harder to build, requires RDKit reaction matching

### Example: Cyanation Protocol

```json
{
  "reaction_SMARTS": ["IC=C.CC(O)(C#N)C>>N#CC=C"],
  "family": "Cyanation_Cu_Alkenyl_Iodide",
  "compatible_functional_groups": ["ketone", "aromatic"],
  "incompatible_functional_groups": ["enamine", "alcohol", "carboxylic acid"]
}
```

This pattern encodes:
- **Substrate**: Alkenyl iodide (`IC=C`) specifically at sp² carbon
- **Cyanide source**: Acetone cyanohydrin (`CC(O)(C#N)C`)
- **Product**: Alkenyl nitrile (`N#CC=C`)
- **Transformation**: I → CN replacement preserving alkene geometry

## Accuracy Comparison

### applies_if (Functional Group Detection)
```python
# Can detect:
✓ "sp2_halide_present" → Matches both IC=C (alkenyl) AND IBr (aryl)
✓ "ketone_present" → Matches any C=O
✗ Cannot distinguish: alkenyl vs aryl halides
✗ Cannot distinguish: reaction mechanism (substitution vs coupling)
✗ Cannot encode: stereochemistry, regiochemistry
```

### reaction_SMARTS (Transformation Matching)
```python
# Can detect:
✓ Exact transformation: IC=C >> N#CC=C (alkenyl iodide cyanation)
✓ Distinguishes: IC=C (alkenyl) from Ic1ccccc1 (aryl)
✓ Preserves: Stereochemistry (Z vs E alkenes)
✓ Encodes: Specific reagent requirements (acetone cyanohydrin)
✓ Guards: Can include "guards_forbid" to exclude conflicting groups
✗ More complex: Requires RDKit reaction matching
```

## Real-World Example: Cyanation Specificity

### Query 1: Alkenyl Iodide Cyanation (Should Match)
```
Reaction: IC=CCCC.CC(C)(O)C#N >> N#CC=CCCC
```

**applies_if approach** (if we had it for protocols):
- Would need: `{"all": ["sp2_halide_present"], "any": ["cyanide_source"]}`
- Problem: Also matches aryl iodides, alkyl iodides with C=C elsewhere

**reaction_SMARTS approach**:
```python
Pattern: "IC=C.CC(O)(C#N)C>>N#CC=C"
Match: ✅ YES - exact transformation matches
```

### Query 2: Aryl Iodide Cyanation (Should NOT Match)
```
Reaction: Ic1ccccc1.NaCN >> N#Cc1ccccc1
```

**applies_if approach**:
- Detects: `sp2_halide_present` = TRUE
- Result: ❌ FALSE POSITIVE - would recommend alkenyl iodide protocol

**reaction_SMARTS approach**:
```python
Pattern: "IC=C.CC(O)(C#N)C>>N#CC=C"
Match: ✅ NO - aryl (Ic1ccc) doesn't match alkenyl (IC=C)
```

### Query 3: Alkyl Iodide with Remote C=C (Should NOT Match)
```
Reaction: ICCCC=C.CC(C)(O)C#N >> N#CCCCC=C
```

**applies_if approach**:
- Detects: `sp2_halide_present` = FALSE, `halide_present` = TRUE, has alkene elsewhere
- Result: ✅ Correctly excludes (but for wrong reasons)

**reaction_SMARTS approach**:
```python
Pattern: "IC=C.CC(O)(C#N)C>>N#CC=C"
Match: ✅ NO - C-I not adjacent to C=C
```

## Agreement

**YES, I completely agree** that `reaction_SMARTS` is:

1. ✅ **More accurate** than `applies_if`
   - Encodes exact transformation, not just substrate features
   - Distinguishes regioisomers (alkenyl vs aryl halides)
   - Preserves stereochemistry (Z/E, R/S)

2. ✅ **Harder to build**
   - Requires understanding reaction mechanism
   - Needs RDKit expertise to write patterns
   - Must test against diverse substrates
   - More prone to overfitting (too specific)

3. ✅ **Should be used for protocol validation**
   - Protocols are already curated, high-value resources
   - Worth the effort to add accurate matching
   - Can prevent false positives from DRFP similarity alone

## Proposed Implementation Strategy

### Three-Tier Validation System

```
Query Reaction
  ↓
[Stage 1: DRFP Similarity Search]
  - Compute query DRFP
  - Find top-k similar protocols/rules
  - Fast, broad search
  ↓
[Stage 2A: reaction_SMARTS Validation] ← NEW for protocols
  - For each protocol with reaction_SMARTS
  - Use RDKit reaction matching: rxn.RunReactants(query_reactants)
  - Check if products match or pattern matches reactants
  - FILTER OUT if no match
  ↓
[Stage 2B: applies_if Validation] ← EXISTING for rules
  - For each rule with applies_if
  - Detect features from query
  - Check all/any conditions
  - FILTER OUT if criteria not met
  ↓
[Stage 3: Condition Selection]
  - For remaining sources, select best base_rule/conditions
  - Use feature specificity scoring
  ↓
Display: Validated, substrate-optimized recommendations
```

### Benefits of Three-Tier System

1. **Stage 1 (DRFP)**: Fast discovery of potentially relevant chemistry
2. **Stage 2A (reaction_SMARTS)**: High-accuracy protocol filtering
3. **Stage 2B (applies_if)**: Functional-group based rule filtering
4. **Stage 3**: Within-source optimization

### Implementation Priorities

1. ✅ **DONE**: applies_if validation for rules
2. ⏳ **NEXT**: reaction_SMARTS validation for protocols
3. ⏳ **FUTURE**: Add applies_if to protocols as fallback (easier to write than SMARTS)
4. ⏳ **FUTURE**: Add reaction_SMARTS to rules (for high-value, well-characterized rules)

## Code Design

### Proposed Method Signature

```python
def _validate_protocol_smarts(
    self,
    protocol: Dict[str, Any],
    query_smiles: str
) -> bool:
    """
    Validate if a protocol's reaction_SMARTS pattern matches the query reaction.
    
    Args:
        protocol: Protocol dictionary with 'reaction' containing 'reaction_SMARTS'
        query_smiles: Query reaction SMILES (reactants>>products)
    
    Returns:
        True if pattern matches OR no reaction_SMARTS present (permissive)
        False if pattern exists but doesn't match (strict filtering)
    """
```

### Matching Strategy

Two approaches (both should be supported):

#### Approach 1: Reactant Substructure Matching
```python
# Check if query reactants match pattern reactants
rxn_pattern = AllChem.ReactionFromSmarts(smarts)
query_reactants = [Chem.MolFromSmiles(r) for r in reactants]

# Try to match reactants
matches = rxn_pattern.RunReactants(tuple(query_reactants))
if len(matches) > 0:
    return True  # Pattern can apply to these reactants
```

#### Approach 2: Full Reaction Matching
```python
# Check if entire reaction (reactants>>products) matches pattern
# More stringent - requires products to also match
# Useful for validating actual transformations occurred
```

### Error Handling

```python
# Fail-open behavior (like applies_if)
try:
    if 'reaction_SMARTS' not in protocol.get('reaction', {}):
        return True  # No SMARTS → always include
    
    # Try to match
    matches = self._check_reaction_smarts(smarts, query_smiles)
    return matches
except Exception as e:
    # If matching fails (RDKit error, invalid SMARTS), include the protocol
    # Better to over-recommend than miss valid protocols
    return True
```

## Gradual Rollout Plan

### Phase 1: Implementation
1. Add `_validate_protocol_smarts()` method to `UnifiedRecommender`
2. Implement RDKit reaction matching logic
3. Add `validate_smarts` parameter (default=True) to `recommend()`
4. Call validation in recommendation pipeline after DRFP search

### Phase 2: Testing
1. Test with existing protocols that have `reaction_SMARTS`
2. Verify filtering reduces false positives
3. Check that true positives are retained
4. Compare with/without validation

### Phase 3: Documentation
1. Document validation behavior in protocol schema
2. Add examples of good `reaction_SMARTS` patterns
3. Create guide for protocol authors
4. Update `APPLIES_IF_VALIDATION.md` to cover both mechanisms

### Phase 4: Expansion
1. Add `reaction_SMARTS` to more protocols
2. Consider adding `applies_if` to protocols as simpler fallback
3. Provide tooling to generate SMARTS from reaction SMILES
4. Consider adding `reaction_SMARTS` to high-value rules

## Comparison Table

| Feature | applies_if | reaction_SMARTS |
|---------|-----------|-----------------|
| **Accuracy** | Moderate | **Very High** |
| **Ease of writing** | **Easy** | Hard |
| **Specificity** | Functional groups only | **Exact transformation** |
| **Current status** | ✅ Implemented | ❌ Not used |
| **Best for** | Rules (many variants) | Protocols (curated) |
| **Can distinguish** | Presence/absence | **Regio/stereochemistry** |
| **Maintenance** | **Low** | High |
| **False positives** | Higher | **Lower** |
| **False negatives** | Lower | Higher (if too specific) |

## Recommendation

**Implement both validation mechanisms** with appropriate usage:

1. **Protocols**: Use `reaction_SMARTS` (more accurate, worth the effort)
2. **Rules**: Use `applies_if` (easier to maintain, good enough for variations)
3. **Fallback**: Always fail-open if validation fails or field missing
4. **Future**: Consider adding both fields to all sources for defense-in-depth

This hybrid approach maximizes accuracy while maintaining practicality.
