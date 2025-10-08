# Recommend/ML Refactoring - Work In Progress

**Date**: October 7, 2025  
**Status**: ⏸️ **PAUSED** - Awaiting decision on approach  
**Progress**: 20% Complete

---

## ✅ Completed Work

### 1. Directory Structure Created
```
chemtools/
├── recommend/              ✅ Created
│   ├── families/           ✅ Created
│   ├── utils.py            ✅ Created (130 lines)
│   ├── substrate_analysis.py  ✅ Created (165 lines)
│   └── plate_design.py     ✅ Created (190 lines)
```

### 2. Modules Completed

#### ✅ `recommend/utils.py` (130 lines)
**Extracted from**: recommend.py lines 19-40, 817-856  
**Functions**:
- `canonical_family()` - Family name normalization
- `pick_electrophile_nucleophile()` - Reactant classification
- `median()` - Median calculation
- `pick_with_constraints()` - Constraint application

#### ✅ `recommend/substrate_analysis.py` (165 lines)
**Extracted from**: recommend.py lines 462-506  
**Functions**:
- `has_free_alcohol()` - Free alcohol detection
- `has_phenol()` - Phenol detection
- `has_sulfonamide()` - Sulfonamide detection
- `has_hydroxylamine()` - Hydroxylamine detection
- `detect_functional_groups()` - Combined functional group analysis

#### ✅ `recommend/plate_design.py` (190 lines)
**Extracted from**: recommend.py lines 1339-1454  
**Functions**:
- `design_plate_from_reaction()` - Main plate design function
- `_well_ids()` - Well ID generation

**Note**: Has import dependency on `core.py` (not yet created)

---

## ⏸️ Remaining Work

### Large Extraction Tasks

#### 1. **families/base.py** (~100 lines)
**Purpose**: Base interface for family-specific analyzers  
**Status**: ❌ Not started

**Needs**:
- Abstract `FamilyAnalyzer` class
- `analyze_substrates()` method
- `build_features()` method
- `get_default_conditions()` method

---

#### 2. **families/amide.py** (~550 lines) ⚠️ **LARGEST MODULE**
**Extract from**: recommend.py lines 297-836  
**Status**: ❌ Not started

**Code sections to extract**:
- Lines 297-411: `_analyze_carboxylic_acid()` - Acid substrate analysis (115 lines)
- Lines 412-461: `_analyze_amine()` - Amine substrate analysis (50 lines)
- Lines 507-543: `_acid_classification()` - Acid classification (37 lines)
- Lines 544-595: `_amine_classification()` - Amine classification (52 lines)
- Lines 596-601: `_amide_substrate_class()` - Substrate class mapping (6 lines)
- Lines 602-626: `_infer_amide_category()` - Category inference (25 lines)
- Lines 627-637: `_water_management_for_category()` - Water tolerance (11 lines)
- Lines 638-650: `_looks_like_aromatic_aniline()` - Aniline detection (13 lines)
- Lines 729-794: `_amide_rule_feature_builder()` - Feature builder (66 lines)

**Plus**: Create `AmideAnalyzer` class integrating all above

---

#### 3. **families/cn_coupling.py** (~200 lines)
**Extract from**: recommend.py lines 651-728, scattered C-N logic  
**Status**: ❌ Not started

**Code sections to extract**:
- Lines 178-210: `_electrophile_class_text()` - Electrophile description
- Lines 211-229: `_nucleophile_class_text()` - Nucleophile description
- Lines 651-728: `_default_rule_feature_builder()` - C-N feature builder
- Lines 817-835: `_pick_electrophile_nucleophile()` - Already in utils.py ✅

**Plus**: Create `CNCouplingAnalyzer` class

---

#### 4. **families/cc_coupling.py** (~150 lines)
**Extract from**: Scattered Suzuki/Sonogashira logic in recommend.py  
**Status**: ❌ Not started

**Needs**:
- `SuzukiAnalyzer` class
- `SonogashiraAnalyzer` class
- Boron/halide matching logic

---

#### 5. **families/__init__.py** (~50 lines)
**Status**: ❌ Not started

**Needs**:
- `FAMILY_REGISTRY` dict
- `get_family_analyzer()` factory function
- Exports

---

#### 6. **core.py** (~400-500 lines) ⚠️ **CRITICAL - MOST COMPLEX**
**Extract from**: recommend.py lines 858-1338  
**Status**: ❌ Not started

**Main functions**:
- Lines 858-1263: `recommend_from_reaction()` - Main recommendation (405 lines)
  - Reaction normalization
  - Family detection
  - Role featurization integration
  - Precedent search
  - Output formatting
  
- Lines 1264-1338: `recommend_conditions_structured()` - Structured output (74 lines)
  - Multiple variants
  - Reagent enrichment
  - Confidence scoring

**Supporting code**:
- Lines 42-89: Optional role-aware featurization
- Lines 90-177: Role featurization logic
- Lines 230-240: `_RuleFeatureContext` dataclass
- Lines 241-296: Role pack utilities
- Lines 795-816: `_compose_rule_features()` - Feature composition

---

#### 7. **recommend/__init__.py** (~30 lines)
**Status**: ❌ Not started

**Needs**:
- Public API exports
- `recommend_from_reaction`
- `recommend_conditions_structured`
- `design_plate_from_reaction`

---

### ML Package Enhancement

#### 8. **ml/predictors/** subpackage
**Status**: ❌ Not started

**Tasks**:
- Create `ml/predictors/__init__.py`
- Move `ml/drfp_predictor.py` → `ml/predictors/drfp_predictor.py`
- Update imports

---

#### 9. **ml/hybrid.py** (~250 lines)
**Extract from**: recommend_ml.py (entire file)  
**Status**: ❌ Not started

**Functions**:
- `hybrid_recommend()` - ML+k-NN hybrid (150 lines)
- `recommend_with_yield_prediction()` - Simplified API (100 lines)

---

#### 10. **ml/__init__.py** update
**Status**: ❌ Not started

**Needs**:
- Update exports to include `hybrid`
- Import from `predictors/`

---

### Integration & Testing

#### 11. **context.py** updates
**Status**: ❌ Not started

**Changes needed**:
- Update `RecommendNamespace` imports
- Change `from . import recommend` → `from .recommend import core`
- Change `from . import recommend_ml` → `from .ml import hybrid`

---

#### 12. **Delete old files**
**Status**: ❌ Not started

**Files to remove**:
- `chemtools/recommend.py` (1,454 lines)
- `chemtools/recommend_ml.py` (226 lines)

---

#### 13. **Testing**
**Status**: ❌ Not started

**Needs**:
- Create `test_recommend_refactoring.py`
- Test all public APIs
- Test family analyzers
- Test ML hybrid

---

#### 14. **Documentation**
**Status**: ❌ Not started

**Updates needed**:
- `CHEMTOOLS_FILE_SUMMARY.md` - Update directory structure
- `RECOMMEND_ML_REFACTORING_COMPLETE.md` - Final summary

---

## 📊 Progress Summary

| Task | Lines | Status | Complexity |
|------|-------|--------|------------|
| ✅ utils.py | 130 | **DONE** | Low |
| ✅ substrate_analysis.py | 165 | **DONE** | Low |
| ✅ plate_design.py | 190 | **DONE** | Medium |
| ❌ families/base.py | 100 | Not started | Medium |
| ❌ families/amide.py | ~550 | Not started | **High** |
| ❌ families/cn_coupling.py | ~200 | Not started | Medium |
| ❌ families/cc_coupling.py | ~150 | Not started | Medium |
| ❌ families/__init__.py | 50 | Not started | Low |
| ❌ core.py | ~450 | Not started | **Very High** |
| ❌ __init__.py | 30 | Not started | Low |
| ❌ ml/predictors/ | - | Not started | Low |
| ❌ ml/hybrid.py | 250 | Not started | Medium |
| ❌ ml/__init__.py | 40 | Not started | Low |
| ❌ context.py updates | - | Not started | Medium |
| ❌ Delete old files | - | Not started | Low |
| ❌ Testing | - | Not started | High |
| ❌ Documentation | - | Not started | Low |

**Overall Progress**: 485 / ~2,700 lines = **~18% complete**

---

## 🤔 Decision Points

### Option 1: Continue Full Automated Refactoring
**Pros**:
- Complete implementation
- All code properly extracted
- Ready to use immediately

**Cons**:
- Will take significant time (~15-20 more file operations)
- Risk of errors in complex extractions (especially `core.py` and `amide.py`)
- Harder to review all changes

**Estimate**: 2-3 hours of work, 40+ more tool calls

---

### Option 2: Create Detailed Extraction Guide
**Pros**:
- Clear documentation of what goes where
- You can review before implementing
- Easier to spot issues

**Cons**:
- Requires manual implementation
- Takes longer overall

**Estimate**: 30 minutes to document, your time to implement

---

### Option 3: Hybrid Approach - Critical Modules Only
**Pros**:
- Focus on most important parts (core.py, families/)
- Less risky
- Can iterate

**Cons**:
- Incomplete refactoring
- Still need to finish later

**Estimate**: 1-2 hours for critical modules

---

### Option 4: Python Script-Assisted Extraction
**Pros**:
- Automated but reviewable
- Can run multiple times
- Less error-prone

**Cons**:
- Need to write the script first
- Still need to review output

**Estimate**: 1 hour for script + execution

---

## 💡 Recommendation

Given the complexity, I recommend **Option 3: Hybrid Approach**:

1. ✅ **Complete** the three simple modules (utils, substrate_analysis, plate_design) - **DONE**
2. ⏸️ **Create** a Python extraction script for the complex modules
3. ⏸️ **Focus** on core.py and families/amide.py (the two largest)
4. ⏸️ **Test** incrementally
5. ⏸️ **Finish** remaining modules if tests pass

This balances automation with safety and allows us to validate the approach before completing everything.

---

## ❓ Next Steps - Your Decision

**Please choose**:
- **A**: Continue with full automated refactoring (Option 1)
- **B**: Create detailed extraction guide (Option 2)  
- **C**: Hybrid approach - critical modules only (Option 3) ⭐ **RECOMMENDED**
- **D**: Create Python extraction script (Option 4)
- **E**: Pause and review what's done so far

**Current state**: 3 modules created, directory structure ready, ~485 lines extracted from 1,454-line monolith.
