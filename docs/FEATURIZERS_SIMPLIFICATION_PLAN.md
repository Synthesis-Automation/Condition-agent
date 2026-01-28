# Featurizers Module Simplification Plan

## Executive Summary

The `/chemtools/featurizers` directory has grown to **~6,000 lines** across 21 files with overlapping responsibilities, legacy APIs, and unclear module boundaries. This plan outlines a phased approach to simplify, consolidate, and modernize the featurizers architecture.

---

## Current State Analysis

### Module Inventory (by size)

#### Core Modules (5,115 lines)
1. **reaction_detection.py** (820 lines) - Reaction type matching
2. **unified.py** (735 lines) - Main public API for molecule/reaction features
3. **molecule.py** (592 lines) - Molecule-level featurization
4. **reaction_pair.py** (528 lines) - Electrophile/nucleophile pair features
5. **motif_detect.py** (478 lines) - SMARTS pattern matching
6. **motif_registry.py** (467 lines) - Build compound registry from taxonomy
7. **calculable.py** (410 lines) - Generate calculable features from spec
8. **aryl_electronics.py** (344 lines) - Aryl electronic effects
9. **alkyl_steric.py** (234 lines) - Alkyl steric analysis
10. **aryl_steric.py** (218 lines) - Aryl steric analysis

#### Support Modules (266 lines)
11. **spectator_rank.py** (114 lines) - Rank spectator groups by relevance
12. **reaction.py** (87 lines) - Legacy reaction featurization
13. **nearby_groups.py** (65 lines) - Detect nearby functional groups

#### Shims (10 lines)
14. **structural.py** (10 lines) - Compatibility shim (LEGACY)

#### Analysis Submodule (1,861 lines)
15. **analysis/reactants.py** (803 lines) - Reactant classification
16. **analysis/reaction_context.py** (340 lines) - Context-aware reactant roles
17. **analysis/reactions.py** (219 lines) - Reaction family normalization
18. **analysis/smiles.py** (175 lines) - SMILES normalization
19. **analysis/feasibility.py** (125 lines) - SNAr feasibility analysis
20. **analysis/__init__.py** (157 lines) - Analysis module orchestration
21. **analysis/_registry.py** (42 lines) - Private registry helpers

---

## Problems Identified

### 1. **API Fragmentation**
- **3 overlapping entry points**: `molecule.py`, `reaction.py`, `unified.py`
- **Usage**: Most code now imports from `unified`, but legacy imports remain
- **Impact**: Confusion about which API to use

### 2. **Duplicate Featurization Logic**
- `reaction.py` (87 lines) is a thin wrapper around `unified.py` functionality
- `molecule.py` and `unified.py` have overlapping molecule featurization
- `structural.py` is a 10-line compatibility shim (dead weight)

### 3. **Oversized Modules**
- **reaction_detection.py** (820 lines): Monolithic reaction type matcher
- **unified.py** (735 lines): Does too much (orchestration + formatting + aggregation)
- **molecule.py** (592 lines): Complex feature extraction
- **analysis/reactants.py** (803 lines): Reactant classification is too large

### 4. **Unclear Responsibilities**
- **molecule.py** vs **unified.py**: Both featurize molecules, unclear split
- **reaction.py** vs **unified.py**: Duplicate reaction featurization
- **motif_detect.py** vs **motif_registry.py**: Registry building vs detection split is unclear

### 5. **Legacy Compatibility Overhead**
- `analysis/reactants.py` has extensive legacy ID mapping (lines 149-365)
- Legacy APIs maintained for backward compatibility
- No clear deprecation path

### 6. **Analysis Submodule Coupling**
- `analysis/` is tightly coupled to taxonomy data files
- Hard to test in isolation
- Mixed concerns: SMILES normalization, reactant classification, reaction analysis

---

## Simplification Strategy

### **Phase 1: API Consolidation** (1 week)

#### Goal: Single, clear entry point for featurization

**Actions:**
1. **Deprecate legacy APIs**:
   - Mark `reaction.py::featurize_reaction` as deprecated (redirect to `unified`)
   - Mark `structural.py` as deprecated (delete after migration window)
   - Add deprecation warnings with migration guide

2. **Consolidate molecule featurization**:
   - Keep `unified.py::featurize_molecule` as **primary API**
   - Keep `molecule.py::featurize_molecule` as **implementation**
   - Move `molecule.py::analyze_smiles` to `unified.py` (user-facing)

3. **Update imports**:
   ```python
   # Before (multiple sources)
   from chemtools.featurizers.molecule import featurize_molecule
   from chemtools.featurizers.reaction import featurize_reaction
   from chemtools.featurizers.structural import featurize_molecule
   
   # After (single source)
   from chemtools.featurizers.unified import featurize_molecule, featurize_reaction
   ```

**Expected Impact**:
- Remove 100 lines (delete `structural.py`, simplify `reaction.py`)
- Single import pattern across codebase

---

### **Phase 2: Module Refactoring** (2 weeks)

#### Goal: Split oversized modules, clarify responsibilities

**2.1. Split `reaction_detection.py` (820 lines)**

Current structure is monolithic. Split into:

```
featurizers/
  detection/
    __init__.py          # Public API
    core.py              # Main detect_reaction_types function (150 lines)
    matchers.py          # Reaction type matching logic (300 lines)
    motif_analysis.py    # Motif set analysis (200 lines)
    scoring.py           # Confidence scoring (150 lines)
    utils.py             # Helper functions (100 lines)
```

**2.2. Split `unified.py` (735 lines)**

Separate orchestration from formatting:

```
featurizers/
  unified.py           # Public API + orchestration (300 lines)
  formatters/
    __init__.py
    molecule.py        # Molecule output formatting (150 lines)
    reaction.py        # Reaction output formatting (200 lines)
    aggregation.py     # Aggregate feature extraction (100 lines)
```

**2.3. Split `analysis/reactants.py` (803 lines)**

Too much responsibility. Split into:

```
featurizers/analysis/
  reactants/
    __init__.py         # Public API
    classification.py   # Core classification (250 lines)
    legacy_mapping.py   # Legacy ID compatibility (200 lines)
    taxonomy.py         # Taxonomy integration (200 lines)
    formatters.py       # Output formatting (150 lines)
```

**Expected Impact**:
- Reduce max file size from 820 → 300 lines
- Clearer module responsibilities
- Easier testing and maintenance

---

### **Phase 3: Registry Consolidation** (1 week)

#### Goal: Simplify motif detection pipeline

**3.1. Merge `motif_detect.py` + `motif_registry.py`**

These are tightly coupled (registry builds, detect uses). Merge into:

```
featurizers/motifs/
  __init__.py          # Public API
  registry.py          # Build registry from taxonomy (400 lines)
  detection.py         # SMARTS pattern matching (400 lines)
  patterns.py          # Pattern definitions (100 lines)
```

**3.2. Move taxonomy interaction to dedicated module**

Currently scattered across multiple files:

```
featurizers/motifs/
  taxonomy_loader.py   # Single source for loading taxonomy data
```

**Expected Impact**:
- Remove 100 lines of duplicate loading logic
- Clearer separation: registry building vs detection

---

### **Phase 4: Legacy Cleanup** (1 week)

#### Goal: Remove legacy compatibility layer

**4.1. Remove legacy ID mapping**

After Phase 1 deprecation period (3 months):
- Delete `analysis/reactants.py::_pick_legacy_alias` (50 lines)
- Delete legacy ID fields from output
- Update documentation

**4.2. Remove deprecated files**

Delete after migration:
- `structural.py` (10 lines)
- `reaction.py` (87 lines) - merge into `unified.py`

**4.3. Simplify analysis module**

```
featurizers/analysis/
  __init__.py          # Public API (50 lines, down from 157)
  reactants/           # New structure (from Phase 2)
  context.py           # Rename reaction_context.py for clarity
  families.py          # Rename reactions.py for clarity
  feasibility.py       # Keep as-is (125 lines)
  smiles.py            # Keep as-is (175 lines)
```

**Expected Impact**:
- Remove 300+ lines of legacy code
- Clearer module naming

---

### **Phase 5: Documentation & Testing** (1 week)

#### Goal: Comprehensive migration guide and tests

**5.1. Documentation**

Create:
- `docs/FEATURIZERS_API.md` - Complete API reference
- `docs/MIGRATION_GUIDE_V3.md` - Migration from old APIs
- Inline docstrings for all public functions

**5.2. Testing**

Add:
- Integration tests for new structure
- Deprecation warning tests
- Performance benchmarks (ensure no regression)

---

## Success Metrics

### Before (Current State)
- **21 files**, 6,000+ lines
- **3 entry points** (molecule, reaction, unified)
- **Max file size**: 820 lines
- **Legacy overhead**: ~500 lines
- **Import confusion**: Multiple valid import paths

### After (Target State)
- **~25 files** (more files, but smaller/focused), 5,500 lines
- **1 entry point** (unified)
- **Max file size**: 300 lines
- **Legacy overhead**: 0 lines (after deprecation period)
- **Import clarity**: Single canonical path

### Quality Improvements
- ✅ No file exceeds 400 lines
- ✅ Clear module responsibilities
- ✅ Easier to test (smaller units)
- ✅ Better performance (remove legacy overhead)
- ✅ Modern Python structure (subpackages)

---

## Risk Mitigation

### Breaking Changes
- **Phase 1**: Add deprecation warnings, maintain old APIs
- **Grace period**: 3 months before removal
- **Documentation**: Clear migration examples

### Performance
- Benchmark before/after each phase
- Target: No regression, potential 10-20% improvement from legacy removal

### Testing
- Full test coverage before refactoring
- Run existing tests after each change
- Add new tests for refactored modules

---

## Implementation Timeline

| Phase | Duration | Effort | Dependencies |
|-------|----------|--------|--------------|
| Phase 1: API Consolidation | 1 week | 2 days | None |
| Phase 2: Module Refactoring | 2 weeks | 5 days | Phase 1 |
| Phase 3: Registry Consolidation | 1 week | 3 days | Phase 2 |
| Phase 4: Legacy Cleanup | 1 week | 2 days | Phase 1 (after grace period) |
| Phase 5: Documentation & Testing | 1 week | 3 days | All phases |
| **Total** | **6 weeks** | **15 days** | - |

*Note: Duration includes review, testing, and buffer time*

---

## Next Steps

1. **Review this plan** with team/stakeholders
2. **Create issues** for each phase in version control
3. **Start Phase 1** (low risk, high value)
4. **Monitor usage** during deprecation period
5. **Iterate** based on feedback

---

## Questions for Decision

1. **Deprecation period**: 3 months enough, or should we extend?
2. **Breaking changes**: Accept for v3.0, or maintain full backward compatibility?
3. **Analysis submodule**: Keep separate, or merge into main featurizers?
4. **Priority**: Urgent cleanup, or can we stage over 2-3 releases?

---

**Document Version**: 1.0  
**Date**: 2026-01-28  
**Author**: AI Assistant  
**Status**: Draft for Review
