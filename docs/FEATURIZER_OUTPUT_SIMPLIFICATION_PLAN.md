# Featurizer Output Simplification Plan

## Executive Summary

The current compound and reaction featurizers produce overly complex nested output structures with significant redundancy, making them difficult to consume. This plan proposes a **two-tier architecture** with simplified core outputs and optional extended analysis.

---

## Current Output Complexity

### **Molecule Output (13 top-level fields)**

```python
{
  "schema_version": "v1",
  "kind": "molecule",
  "molecule": {
    "schema_version": "v2",
    "smiles": "...",
    "motifs": [...],              # List of detected motifs with metadata
    "context_motifs": [...],       # Background/scaffold motifs
    "ranked_motifs": [...],        # Sorted by importance
    "steric": {                    # Steric analysis per motif
      "aryl": [...],
      "alkyl": [...]
    },
    "electronics": {               # Electronic analysis per motif
      "aryl": [...]
    },
    "nearby": [...],               # Nearby functional groups
    "aryl_analysis": [...],        # Duplicate of electronics.aryl?
    "analyses": [...],             # Unified per-motif analysis
    "meta": {...},                 # Error/warning metadata
    "rdkit_properties": {...},     # 10+ RDKit descriptors
    "snar_feasibility": [...]      # SNAr reaction check
  },
  "meta": {...}                    # Top-level errors
}
```

**Problems:**
1. **Double nesting**: `payload.molecule.motifs` instead of `payload.motifs`
2. **Redundancy**: `steric`/`electronics`/`aryl_analysis` vs `analyses` (same data, different groupings)
3. **Excessive metadata**: Two `meta` sections, two `schema_version` fields
4. **Mixed concerns**: Detection + analysis + feasibility + properties in one bundle
5. **Unclear hierarchy**: 13+ fields with overlapping purposes

---

### **Reaction Output (16+ top-level fields)**

```python
{
  "schema_version": "v1",
  "kind": "reaction",
  "reaction": {
    "reaction_smiles": "...",
    "reaction_key": "...",         # CRK-v1 summary
    "normalized": {...},           # Normalization metadata
    "reaction_type": {...},        # Best match + alternatives
    "detection": {...},            # All detection matches + evidence
    "aggregates": {                # Reaction-wide statistics
      "reactant_count": 2,
      "motif_ids": [...],
      "reacted_motifs": [...],
      "formed_motifs": [...],
      "spectator_motifs": [...],
      "spectator_groups_combined": [...],
      "spectator_groups_ranked": [...],
      "max_aryl_steric": 0.5,
      "max_alkyl_steric": 2.0,
      "avg_aryl_electronic": 5.2,
      "electron_poor_aryl": false,
      "has_multi_functional_substrates": false
    },
    "roles": {...},                # Reactant role classification
    "agent_roles": {...},          # Reagent role classification
    "intramolecular": {...},       # Intramolecular check
    "feasibility": {...},          # Feasibility analysis
    "reactants": [                 # Full molecule bundles
      {
        "smiles": "...",
        "motifs": [...],
        "analyses": [...],
        "rdkit_properties": {...}
      },
      ...
    ],
    "products": [...],             # Full molecule bundles
    "meta": {...}
  },
  "meta": {...}
}
```

**Problems:**
1. **Redundant summaries**: `aggregates` duplicates info from `reactants[].motifs` + `products[].motifs`
2. **Multiple detection sections**: `reaction_type` + `detection` + `roles` all describe the same classification
3. **Excessive nesting**: `payload.reaction.reactants[0].motifs` is 4 levels deep
4. **Feature bloat**: 16+ fields when most consumers only need 5-7
5. **Unclear priorities**: No indication of what fields are essential vs optional

---

## Proposed Simplification Strategy

### **Phase 1: Two-Tier Architecture** (Recommended)

#### **Core Output** (Default, ~6-8 fields)
Essential information for 80% of use cases:

```python
# Molecule Core
{
  "smiles": "...",
  "motifs": [                      # Primary functional groups
    {
      "id": "Ar-Br",
      "rank": 450,
      "a_atom_idx": 5,
      "b_atom_idx": 6
    }
  ],
  "properties": {                   # Key steric/electronic summary
    "max_steric": 2.0,
    "avg_electronic": 5.2
  },
  "rdkit": {                        # Standard descriptors
    "mw": 157.0,
    "logP": 2.4,
    "tpsa": 29.1
  }
}

# Reaction Core
{
  "reaction_smiles": "...",
  "reaction_type": "Suzuki-Miyaura",   # Best match only
  "confidence": 0.92,
  "reaction_key": "|Ar-Br|Ar-B(OR)2 -> Ar-Ar | bond_formed: C(ar)-C(ar) | bond_broken: Br-C(ar); B-C(ar) | spectators: []",
  "reactants": [                       # Simplified molecules
    {
      "smiles": "...",
      "motifs": [{"id": "Ar-Br", "rank": 450}]
    }
  ],
  "products": [                        # Simplified molecules
    {
      "smiles": "...",
      "motifs": [{"id": "Ar-Ar", "rank": 380}]
    }
  ],
  "feasibility": "high"                # Simple enum: high/medium/low
}
```

#### **Extended Output** (On-demand, opt-in via `detailed=True`)
Advanced analysis for specialized needs:

```python
{
  # Core fields (always included)
  ...
  
  # Extended analysis (opt-in)
  "extended": {
    "detection": {                     # All matches, not just top
      "matches": [...]
    },
    "per_motif_analysis": [            # Detailed steric/electronic breakdown
      {
        "motif_id": "Ar-Br",
        "steric": {
          "score": 2.0,
          "classification": "secondary",
          "description": "moderately hindered"
        },
        "electronic": {
          "score": 4.8,
          "description": "electron-neutral"
        },
        "nearby_groups": ["CN", "OCH3"]
      }
    ],
    "aggregates": {                    # Reaction-wide statistics
      "spectator_motifs": [...],
      "spectator_groups_ranked": [...]
    },
    "role_classification": {...},     # Advanced reactant/agent roles
    "snar_feasibility": [...],         # Detailed feasibility check
    "normalization_log": {...}         # Input normalization details
  }
}
```

**Benefits:**
- **80/20 rule**: Core covers most use cases, extended available when needed
- **Reduced complexity**: Core has 6-8 fields vs current 13-16
- **Clear separation**: Detection/analysis separated from core structure
- **Backward compatible**: Extended matches current output structure

---

### **Phase 2: Field Consolidation** (Aggressive simplification)

#### Merge Redundant Sections

**Molecule:**
```python
# BEFORE: 3 separate sections
"steric": {"aryl": [...], "alkyl": [...]},
"electronics": {"aryl": [...]},
"analyses": [...]                      # Contains duplicates

# AFTER: Single unified section
"motif_analysis": [                    # One analysis per motif
  {
    "motif_id": "Ar-Br",
    "steric": {...},                   # Inline
    "electronic": {...},               # Inline
    "nearby_groups": [...]             # Inline
  }
]
```

**Reaction:**
```python
# BEFORE: 3 separate detection sections
"reaction_type": {...},                # Best match
"detection": {"matches": [...]},       # All matches
"roles": {...}                         # Role classification

# AFTER: Single detection section
"detection": {
  "reaction_type": "Suzuki-Miyaura",
  "confidence": 0.92,
  "alternatives": [                    # Top 3 alternatives
    {"type": "Negishi", "confidence": 0.68}
  ],
  "slot_evidence": {                   # What motifs matched
    "electrophile": ["Ar-Br"],
    "nucleophile": ["Ar-B(OR)2"]
  }
}
```

#### Remove Double Nesting

```python
# BEFORE: Unnecessary wrapper
{
  "kind": "molecule",
  "molecule": {                        # Extra level
    "smiles": "...",
    "motifs": [...]
  }
}

# AFTER: Flat structure
{
  "kind": "molecule",
  "smiles": "...",
  "motifs": [...]
}
```

#### Simplify Metadata

```python
# BEFORE: Two meta sections
{
  "meta": {...},                       # Top level
  "molecule": {
    "meta": {...}                      # Nested
  }
}

# AFTER: Single meta section
{
  "meta": {
    "errors": [],
    "warnings": [],
    "rdkit_available": true
  }
}
```

---

## Implementation Plan

### **Week 1: Core/Extended Split**

**Goal:** Introduce two-tier architecture without breaking existing code

1. **Create new formatters** (`formatters/simplified.py`):
   ```python
   def featurize_molecule_core(smiles, options=None) -> dict:
       """Return core molecule features (6 fields)."""
       
   def featurize_molecule_extended(smiles, options=None) -> dict:
       """Return core + extended analysis (13 fields)."""
   ```

2. **Update `unified.py`** to support `detailed` flag:
   ```python
   def featurize_molecule(smiles, options=None):
       options = options or {}
       if options.get('detailed', False):
           return featurize_molecule_extended(smiles, options)
       return featurize_molecule_core(smiles, options)
   ```

3. **Add backward compatibility wrapper**:
   ```python
   def featurize_molecule_legacy(smiles, options=None):
       """Legacy output format (deprecated)."""
       warnings.warn("Legacy format deprecated, use detailed=True", DeprecationWarning)
       return featurize_molecule_extended(smiles, options)
   ```

**Testing:**
- Add 20+ test cases comparing core/extended/legacy outputs
- Validate schema with JSON Schema validators
- Benchmark performance (core should be ~30% faster)

---

### **Week 2: Field Consolidation**

**Goal:** Merge redundant sections in extended output

1. **Merge `analyses` sections**:
   - Combine `steric`/`electronics`/`nearby` into single `motif_analysis` list
   - Keep backward compatibility by exporting both formats temporarily

2. **Consolidate detection**:
   - Merge `reaction_type` + `detection` into single `detection` section
   - Move `roles` into `detection.role_classification`

3. **Remove double nesting**:
   - Flatten `molecule` wrapper in core output
   - Keep nested format in legacy mode

**Testing:**
- Update all test fixtures to use new format
- Run integration tests against recommender system
- Validate CLI output formatting

---

### **Week 3: Migration & Cleanup**

**Goal:** Update consumers and deprecate legacy format

1. **Update consumers**:
   - `app/Cpd_rxn_featurization_cli.py` → use core output by default
   - `chem_assistant/chemtools_wrapper.py` → add `detailed` parameter to tools
   - `chemtools/recommend/recommender.py` → use core motifs/properties
   - `app/A_convert_to_hte_format.py` → use core output

2. **Documentation**:
   - Add migration guide (`docs/FEATURIZER_MIGRATION.md`)
   - Update API examples in `docs/HTE_RECOMMENDER_API.md`
   - Create comparison table (current vs simplified)

3. **Deprecation timeline**:
   - **Month 1**: Core/extended available, legacy warns
   - **Month 2**: Legacy requires explicit flag `legacy=True`
   - **Month 3**: Legacy removed entirely

---

## Expected Impact

### **Code Metrics**

| Metric | Current | After Phase 1 | After Phase 2 |
|--------|---------|---------------|---------------|
| **Molecule fields** | 13 | 6 (core) / 13 (ext) | 8 (unified) |
| **Reaction fields** | 16 | 7 (core) / 16 (ext) | 10 (unified) |
| **Nesting depth** | 4 levels | 3 levels | 2 levels |
| **Typical JSON size** | 5-10 KB | 1-2 KB (core) | 2-4 KB |
| **Consumer LOC** | 50-100 lines | 20-30 lines | 15-25 lines |

### **User Experience**

**Before:**
```python
# Extract motifs requires 4 levels of nesting
motifs = payload["reaction"]["reactants"][0]["molecule"]["motifs"]

# Find steric score requires searching 3 sections
steric = None
for entry in payload["reaction"]["reactants"][0]["molecule"]["analyses"]:
    if entry.get("steric"):
        steric = entry["steric"].get("score_0_10")
```

**After (Phase 1 Core):**
```python
# Direct access to essentials
motifs = payload["reactants"][0]["motifs"]
steric = payload["reactants"][0]["properties"]["max_steric"]
```

**After (Phase 2 Unified):**
```python
# Even simpler with consolidated structure
motifs = payload["reactants"][0]["motifs"]
analysis = payload["reactants"][0]["motif_analysis"][0]
steric = analysis["steric"]["score"]
```

---

## Risk Mitigation

### **Breaking Changes**

**Risk:** Existing consumers break when output structure changes

**Mitigation:**
- Maintain legacy format for 3 months with deprecation warnings
- Provide automated migration script to update consumer code
- Version output with `schema_version` field
- Add validation tests for all output formats

### **Performance Regression**

**Risk:** Simplified output generation slower than current

**Mitigation:**
- Benchmark core output (should be faster, not slower)
- Extended output can be slightly slower (more processing)
- Cache expensive computations (already in place)
- Profile before/after with realistic data

### **Lost Information**

**Risk:** Simplified output removes useful data

**Mitigation:**
- Extended output preserves all current fields
- Survey consumers to identify essential vs optional fields
- Add `include_*` flags for granular control
- Document what information moved to extended format

---

## Success Criteria

1. **Adoption**: 80% of consumers use core output within 1 month
2. **Performance**: Core output generation ≤ 60% of current time
3. **Clarity**: New contributor can extract motifs in < 5 minutes
4. **Compatibility**: Zero breakage with legacy flag enabled
5. **Size**: Core JSON payload ≤ 30% of current size

---

## Alternatives Considered

### **A. Do Nothing**
- **Pros**: No migration cost
- **Cons**: Technical debt compounds, new consumers struggle

### **B. Complete Redesign**
- **Pros**: Perfect structure from scratch
- **Cons**: High breakage risk, 2-3 month timeline

### **C. Incremental Cleanup Only**
- **Pros**: Lower risk than full redesign
- **Cons**: Doesn't address fundamental complexity

**Selected: Phase 1 (Two-Tier) + Phase 2 (Consolidation)** balances simplification with compatibility.

---

## Open Questions

1. **Field naming**: Should simplified fields use `motifs` or `functional_groups`?
2. **Feasibility format**: Simple enum ("high/medium/low") vs score (0-10) vs dict?
3. **RDKit properties**: Include in core or extended? (Currently core)
4. **Schema versioning**: Bump to v2 or keep v1 with sub-versions?
5. **Migration tooling**: Auto-convert consumer code or manual updates?

---

## Next Steps

1. **Review with team** - Get feedback on two-tier architecture
2. **Create prototype** - Implement core output for molecule featurizer
3. **User testing** - Have 2-3 consumers try simplified output
4. **Iterate** - Adjust based on feedback before full rollout
5. **Execute** - Follow 3-week implementation plan

---

**Document Status:** Draft  
**Author:** GitHub Copilot  
**Date:** 2026-01-28  
**Review:** Pending team feedback
