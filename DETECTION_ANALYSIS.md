# Reaction Type Detection Analysis & Recommendations

## Current Detection Method

**Hybrid approach: Rule-based (primary) + rxn_insight ML (optional fallback)**

### 1. Rule-Based Detection (`chemtools/router.py`)

- **Location**: `detect_family()` function (lines 175-215)
- **Method**: SMARTS pattern matching on reactants
- **Coverage**: **Only 7 reaction patterns hardcoded**
  1. C-N coupling (aryl/vinyl halide + amine)
  2. C-O coupling (aryl/vinyl halide + alcohol)
  3. C-S coupling (aryl/vinyl halide + thiol)
  4. Suzuki (aryl halide + boron)
  5. Sonogashira (aryl/vinyl halide + terminal alkyne)
  6. Amide coupling (carboxylic acid + amine)
  7. **Default: Unknown** for everything else

### 2. rxn_insight ML Detection (optional)

- **Status**: ✓ Installed, available
- **Accuracy on Phase 2 types**: ⚠️ **Poor** (see test results below)
- **Method**: Trained ML classifier from rxn_insight library
- **Used as**: Fallback/augmentation to rule-based

### 3. Taxonomy JSON Files (metadata only)

- **Files**: `reactant_types.json` (49 types), `reaction_types.json` (72 types)
- **Purpose**: Define reactant categories and reaction metadata
- **Usage**: Classification of individual reactants, NOT reaction family detection
- **Current state**: ✓ Phase 1+2 additions validated (49 reactant types, 72 reaction types)

---

## Phase 2 Detection Test Results

| Reaction Type | Rule-Based | rxn_insight | Final Result | Status |
|--------------|------------|-------------|--------------|--------|
| Esterification (Fischer) | Unknown | Acylation→amide (❌ wrong) | **Unknown** | ❌ FAIL |
| Esterification (benzoic) | Unknown | Not detected | **Unknown** | ❌ FAIL |
| Grignard addition | co_coupling (❌ wrong) | Not detected | **co_coupling** | ❌ WRONG |
| Hydroboration | Unknown | Not detected | **Unknown** | ❌ FAIL |
| Nitrile formation (SN2) | Unknown | Not detected | **Unknown** | ❌ FAIL |
| Finkelstein (Br→I) | Unknown | Not detected | **Unknown** | ❌ FAIL |
| Williamson ether | Unknown | ullmann_cn (❌ wrong) | **ullmann_cn** | ❌ WRONG |
| Claisen condensation | Unknown | C-C Coupling→None | **Unknown** | ❌ FAIL |
| Michael addition | Unknown | C-C Coupling→None | **Unknown** | ❌ FAIL |

**Result: 0/9 Phase 2 reactions correctly detected** (rxn_insight provides incorrect or no mappings)

---

## Evaluation of 3 Options

### **Option A: Add Detection Rules to `router.py`** ⭐⭐⭐⭐

**Approach**: Manually add SMARTS patterns and logic for each Phase 2 reaction type

**Implementation Steps**:

1. Add SMARTS patterns to `_compile_smarts()`:
   - Ketone/aldehyde: `[CX3](=O)[C,H]`
   - Ester: `[CX3](=O)[OX2][C,H]`
   - Grignard: `[Mg][Br,Cl,I]`
   - Borane: `[BH3,B2H6,BX3]`
   - Cyanide: `[C-]#N`
   - Iodide: `[I-]`
   - Alkoxide: `[O-][C,H]`
   - etc.

2. Add detection logic to `detect_family()`:

   ```python
   # Esterification
   if h.get("acid") and h.get("alcohol"):
       fam, conf = "esterification", 0.85
   
   # Grignard addition
   if h.get("carbonyl") and h.get("grignard"):
       fam, conf = "grignard_addition", 0.90
   
   # Nitrile formation
   if h.get("alkyl_halide") and h.get("cyanide"):
       fam, conf = "nitrile_formation", 0.85
   ```

**Pros**:

- ✅ Deterministic, reproducible
- ✅ Fast (no ML overhead)
- ✅ No external dependencies
- ✅ Full control over logic
- ✅ Easy to debug and validate
- ✅ Aligns with existing codebase architecture

**Cons**:

- ⚠️ Manual work (2-4 hours for all Phase 2 types)
- ⚠️ Requires chemistry knowledge for SMARTS
- ⚠️ Each new type needs manual addition

**Estimated effort**: 3-5 hours for Phase 2 completion

---

### **Option B: Refactor to Use Taxonomy JSON** ⭐⭐⭐

**Approach**: Rewrite `detect_family()` to iterate through `reaction_types.json`

**Implementation Steps**:

1. Load taxonomy JSON files at startup
2. For each reaction type in JSON:
   - Check if all required reactant types are present
   - Match SMARTS patterns from `reactant_types.json`
3. Return best match based on confidence scoring

**Pros**:

- ✅ Systematic, scalable
- ✅ Taxonomy-driven (add types via JSON, not code)
- ✅ Easier to maintain long-term
- ✅ Supports dynamic expansion

**Cons**:

- ❌ Major refactoring (200+ lines of new code)
- ❌ Need to design matching algorithm
- ❌ Performance overhead (JSON parsing, iteration)
- ❌ Risk of breaking existing detection
- ❌ Requires extensive testing
- ⚠️ SMARTS in JSON may need refinement

**Estimated effort**: 8-12 hours for full implementation + testing

---

### **Option C: Rely on rxn_insight ML** ⭐

**Approach**: Use existing rxn_insight, add mapping rules for new types

**Implementation Steps**:

1. Add mapping rules in `reaction_type_detector.py` `_map_to_family()`
2. Trust rxn_insight classifications
3. Add manual overrides for edge cases

**Pros**:

- ✅ Minimal code changes
- ✅ ML may generalize to unseen reactions

**Cons**:

- ❌ **Test results show poor accuracy** (0/9 correct)
- ❌ External dependency (black box)
- ❌ No control over ML model
- ❌ Misclassifications hard to fix
- ❌ May not cover all chemistry types
- ❌ Version updates may break mappings

**Estimated effort**: 2-3 hours, but **LOW SUCCESS RATE**

---

## **Recommendation: Option A (Add Detection Rules to router.py)** ⭐⭐⭐⭐⭐

### Why Option A is Best

1. **Proven Architecture**: Current system already uses rule-based detection successfully for 7 reaction types
2. **High Accuracy**: Deterministic rules can achieve >95% accuracy when properly designed
3. **Fast Implementation**: 3-5 hours vs 8-12 hours for Option B
4. **Low Risk**: Incremental additions, easy to validate each rule
5. **Debuggable**: Simple logic, easy to trace failures
6. **No Dependencies**: Works offline, no ML black boxes
7. **Aligns with Codebase**: Minimal disruption to existing architecture

### Implementation Plan (3-5 hours)

**Phase 2A: Add SMARTS patterns** (1 hour)

- Carbonyl, ester, Grignard, borane, cyanide, iodide, alkoxide, alkene

**Phase 2B: Add detection logic** (2-3 hours)

- Esterification, Grignard addition, hydroboration
- Nitrile formation, Finkelstein, Williamson ether
- Claisen condensation, Michael addition, organolithium addition

**Phase 2C: Validation** (1 hour)

- Re-run test suite
- Measure UNKNOWN reduction (expect 107 → 40-60)

### Expected Outcome

- Coverage improvement: **74.9% → 85-90%** (320 → ~360-380 classified)
- UNKNOWN reduction: **107 → ~40-60** reactions

---

## Why NOT Option B or C

**Option B (JSON refactoring)**:

- Too complex for current needs
- Better suited for future v3.0 rewrite
- Save for when you have >200 reaction types

**Option C (rxn_insight)**:

- Test results prove **0% accuracy** on Phase 2 types
- Unreliable for production use
- Keep as optional augmentation only

---

## Next Steps (if approved)

1. Implement Option A detection rules in `router.py`
2. Validate with test suite
3. Document new patterns in code comments
4. Generate final coverage report
