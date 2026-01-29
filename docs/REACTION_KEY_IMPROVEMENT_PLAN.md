# Reaction Key Generation Improvement Plan

## Overview

The reaction key system generates a standardized representation of chemical transformations:

```
[Reacted Motifs] -> [Formed Motifs] || [Spectator Groups]
```

This document outlines the current architecture, known issues, and a phased improvement plan.

---

## Current Architecture

### Pipeline Flow

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────────┐
│  Parse SMILES   │────▶│  Detect Motifs   │────▶│  Count Comparison   │
│  (RDKit)        │     │  (SMARTS)        │     │  (reactant vs prod) │
└─────────────────┘     └──────────────────┘     └──────────┬──────────┘
                                                            │
┌─────────────────┐     ┌──────────────────┐     ┌──────────▼──────────┐
│  Reaction Key   │◀────│  Pattern Filter  │◀────│  Classify: Reacted  │
│  Formatting     │     │  (NEW)           │     │  Formed / Spectator │
└─────────────────┘     └──────────────────┘     └─────────────────────┘
```

### Key Files

| File | Purpose |
|------|---------|
| `chemtools/featurizers/motifs/detection.py` | SMARTS-based motif detection |
| `chemtools/featurizers/formatters/aggregation.py` | Motif change analysis & filtering |
| `chemtools/featurizers/formatters/reaction.py` | Reaction featurization & key generation |
| `chemtools/taxonomy/data/organic_compounds.v1.3.json` | Motif definitions (scaffold + substituent) |
| `chemtools/taxonomy/data/transformation_patterns.json` | Pattern-based filtering rules (NEW) |

---

## Current Issues

### Issue 1: Over-Detection of Reacted Motifs

**Problem**: Multiple motifs detected on the same carbon center all appear as "reacted" even when only one bond transforms.

**Example**:

```
Reaction: IC=CCCCCCC.CC(C)(O)C#N >> N#CC=CCCCCCC

Raw Detection:
  Reacted: [Alkenyl-I, R_acidic-CN, R_acidic-OH]  ❌ R_acidic-OH shouldn't be here
  
Expected:
  Reacted: [Alkenyl-I, R_acidic-CN]  ✓
```

**Root Cause**: Counter-based comparison doesn't track which specific bond transformed.

---

### Issue 2: No Bond-Level Tracking

**Problem**: System counts motif IDs globally, not per-bond transformation.

**Impact**: When a molecule has multiple instances of the same motif, the system can't distinguish which specific bond reacted.

---

### Issue 3: Scaffold Ambiguity

**Problem**: Same substituent with different scaffolds can cause confusion.

**Example**: `Ar-Br` vs `HeteroAr-Br` vs `Alkenyl-Br` all have `-Br` but represent different reactive sites.

---

### Issue 4: Context-Dependent Motifs

**Problem**: Some motifs are only reactive under specific conditions.

**Example**: `R_acidic-H` (acidic C-H) is only relevant when using strong base conditions.

---

## Improvement Plan

### Phase 1: Pattern-Based Filtering ✅ COMPLETED

**Status**: Implemented (January 2026)

**What was done**:

- Created `transformation_patterns.json` with reaction-specific filtering rules
- Added `filter_reacted_by_pattern()` to aggregation logic
- Each reaction type defines expected consumed substituents

**Impact**: Medium - Filters out irrelevant motifs based on reaction type

**Files Changed**:

- `chemtools/taxonomy/data/transformation_patterns.json` (NEW)
- `chemtools/featurizers/formatters/aggregation.py`
- `chemtools/featurizers/formatters/reaction.py`

---

### Phase 2: Primary Reactive Site Selection

**Status**: Planned

**Goal**: For each reactant molecule, select only the most reactive motif.

**Approach**:

```python
def select_primary_motif(motifs: List[Dict]) -> Optional[Dict]:
    """Select the most reactive motif per molecule based on reactivity_weight."""
    if not motifs:
        return None
    return max(motifs, key=lambda m: m.get("reactivity_weight", 0))
```

**Changes Required**:

1. Add `reactivity_weight` to motif taxonomy entries
2. Update aggregation to use primary motif selection
3. Fallback to highest priority if weights are equal

**Effort**: Low  
**Impact**: Medium

---

### Phase 3: Bond-Level Motif Tracking

**Status**: ✅ COMPLETED (January 2026)

**Goal**: Track motifs by their bond indices, not just IDs.

**What was done**:

1. Added `select_primary_motifs_by_atom()` - selects one motif per attachment atom
2. Added `extract_motif_with_bond_info()` - preserves full motif dict with bond info
3. Added `analyze_motif_changes_with_bonds()` - bond-level change comparison
4. Updated `aggregate_reaction_features()` to preserve `reactant_motifs_full` with bond info
5. Added `primary_motif_ids` to output - shows motifs after per-atom selection

**Data Structure** (now preserved through aggregation):

```python
reactant_motifs_full = [
    {"id": "Alkenyl-I", "bond": (0, 1), "a_idx": 0, "b_idx": 1, "reactivity_weight": 0.9},
    {"id": "R_acidic-CN", "bond": (5, 6), "a_idx": 5, "b_idx": 6, "reactivity_weight": 0.8},
    # R_acidic-OH filtered out - same a_idx=5, lower reactivity_weight
]
```

**Logic**:

- Group motifs by `a_idx` (attachment carbon)
- If multiple motifs share same `a_idx`, select highest `reactivity_weight`
- Compare product bonds to determine which specific bond transformed

**Files Changed**:

- `chemtools/featurizers/formatters/aggregation.py` (added functions, updated aggregation)

**Impact**: High - Correctly handles multiple motifs on same carbon

---

### Phase 4: Substituent-Centric Comparison

**Status**: ✅ COMPLETED (January 2026)

**Goal**: Focus on substituent changes rather than full motif IDs.

**What was done**:

- Added `get_scaffold()` helper function
- Created `analyze_substituent_changes()` for cross-scaffold comparison
- Builds substituent → motifs mapping for both reactants and products
- Identifies truly consumed substituents (gone from ALL scaffolds)
- Refines reacted/formed classification based on substituent presence

**Key Functions**:

```python
def get_substituent(motif_id: str) -> str:
    """'Ar-Br' -> '-Br', 'Alkenyl-I' -> '-I'"""
    return "-" + motif_id.split("-", 1)[1] if "-" in motif_id else ""

def analyze_substituent_changes(reactant_motif_ids, product_motif_ids):
    """Cross-scaffold substituent comparison."""
    # Build substituent -> motifs mapping
    reactant_subs = {sub: [motifs...] for each substituent}
    product_subs = {sub: [motifs...] for each substituent}
    
    # Identify truly consumed/appeared substituents
    consumed_subs = set(reactant_subs) - set(product_subs)
    appeared_subs = set(product_subs) - set(reactant_subs)
    
    # Refine classification
    return (reacted_set, formed_set, spectator_set)
```

**Benefits**:

- Works across scaffold types (Ar-Br, HeteroAr-Br, Alkenyl-Br all have -Br)
- Better spectator detection (if substituent exists on other scaffold, might be spectator)
- More chemically intuitive

**Impact**: High

**Files Changed**:

- `chemtools/featurizers/formatters/aggregation.py`

---

### Phase 5: Atom Mapping Integration

**Status**: Future

**Goal**: Use atom-to-atom mapping (AAM) for precise bond change tracking.

**Approach**:

1. Use existing `_atom_mapping.py` module (RXNMapper integration)
2. Map atoms from reactants to products
3. Identify which bonds changed partners

**Example**:

```
Mapped: [C:1]-[Br:2] + [B:3]-[C:4] >> [C:1]-[C:4]
Analysis: 
  - C:1 lost Br:2 partner
  - C:4 lost B:3 partner  
  - C:1 gained C:4 partner
  => Reacted: C-Br, C-B; Formed: C-C
```

**Prerequisites**:

- RXNMapper availability
- Fallback for unmappable reactions

**Effort**: High  
**Impact**: Very High

---

### Phase 6: Taxonomy-Driven Transformation Rules

**Status**: Partially Implemented

**Goal**: Define expected transformations in reaction taxonomy itself.

**Current** (in `transformation_patterns.json`):

```json
{
  "Suzuki_miyaura": {
    "pattern": "coupling",
    "electrophile_consumed": ["-I", "-Br", "-Cl", "-OTf"],
    "nucleophile_consumed": ["-B(OH)2", "-Bpin", "-BF3K"]
  }
}
```

**Future Enhancement**:

```json
{
  "Suzuki_miyaura": {
    "transformation_logic": {
      "bond_formed": ["C(sp2)-C(sp2)"],
      "bond_broken": ["C-X", "C-B"],
      "atoms_lost": ["X", "B(OR)2"]
    }
  }
}
```

**Effort**: Medium  
**Impact**: Medium

---

## Priority Matrix

| Phase | Effort | Impact | Priority | Status |
|-------|--------|--------|----------|--------|
| 1. Pattern Filtering | Low | Medium | P0 | ✅ Done |
| 2. Primary Site Selection | Low | Medium | P1 | Planned |
| 3. Bond-Level Tracking | Medium | High | P1 | ✅ Done |
| 4. Substituent-Centric | Medium | High | P0 | ✅ Done |
| 5. Atom Mapping | High | Very High | P3 | Future |
| 6. Taxonomy Rules | Medium | Medium | P2 | Partial |

---

## Testing Strategy

### Unit Tests

```python
# tests/test_reaction_key_generation.py

def test_coupling_reacted_motifs():
    """Verify only coupling-relevant motifs marked as reacted."""
    rxn = "Brc1ccc(C(=O)O)cc1.Nc1ccccc1>>O=C(O)c1ccc(Nc2ccccc2)cc1"
    result = featurize_reaction(rxn)
    
    reacted = result["aggregates"]["reacted_motifs"]
    assert "Ar-Br" in reacted
    assert "Ar-NH2" in reacted
    assert "Ar-CO2H" not in reacted  # Should be spectator

def test_same_carbon_multiple_motifs():
    """Only primary motif from same carbon should be reacted."""
    rxn = "IC=CCCCCCC.CC(C)(O)C#N>>N#CC=CCCCCCC"
    result = featurize_reaction(rxn)
    
    reacted = result["aggregates"]["reacted_motifs"]
    # Should not include both R_acidic-CN and R_acidic-OH from same carbon
    assert len([m for m in reacted if "R_acidic" in m]) <= 1
```

### Integration Tests

```python
def test_reaction_key_format():
    """Verify reaction key follows expected format."""
    rxn = "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    result = featurize_reaction(rxn)
    
    key = result["reaction_key"]
    assert "->" in key
    assert "||" in key
    assert "Ar-Br" in key
    assert "Ar-B(OH)2" in key
```

---

## Success Metrics

1. **Accuracy**: % of reactions with correct reacted/formed classification
2. **Precision**: False positive rate for reacted motifs
3. **Coverage**: % of reaction types with pattern definitions
4. **Performance**: Processing time per reaction (target: <100ms)

---

## Timeline

| Milestone | Target Date | Deliverable |
|-----------|-------------|-------------|
| Phase 1 Complete | ✅ Jan 2026 | Pattern filtering implemented |
| Phase 4 Complete | ✅ Jan 2026 | Substituent-centric logic |
| Phase 2 Complete | Feb 2026 | Primary site selection |
| Phase 3 Complete | Mar 2026 | Bond-level tracking |
| Phase 5 Evaluation | Q3 2026 | Atom mapping feasibility |

---

## References

- `chemtools/taxonomy/data/transformation_patterns.json` - Pattern definitions (v1.3)
- `chemtools/featurizers/formatters/aggregation.py` - Substituent-centric analysis
- `chemtools/taxonomy/data/organic_compounds.v1.3.json` - Motif taxonomy
- `chemtools/taxonomy/data/reaction_types.v4.0.json` - Reaction type definitions
- `chemtools/_atom_mapping.py` - Atom mapping utilities (for Phase 5)
