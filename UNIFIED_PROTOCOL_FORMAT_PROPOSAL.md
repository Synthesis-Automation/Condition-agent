# UNIFIED PROPOSAL: Adopt Protocol Format for Rules

## Key Insight 💡

**Protocols already have the perfect structure for automated execution!**

Instead of creating a new format, we should:
1. **Add `reaction_setup` to rule outputs** (borrowing protocol structure)
2. **Keep rule authoring simple** (current condition dictionaries)
3. **Generate protocol-like structure dynamically** when needed

## Current Protocol Structure (Already Perfect)

```json
{
  "schema_version": "2.0",
  "source_type": "protocol",
  "reaction_setup": [
    {
      "chemicals": [
        {
          "name": "2-Iodobenzamide",
          "abbreviation": null,
          "cas": "3930-83-4",
          "smiles": "NC(=O)c1ccccc1I",
          "amount": {
            "weight_g": 19.8,
            "mmol": 80.0,
            "volume_ml": null,
            "equivalents": 1.0
          },
          "role": "starting_material"
        },
        {
          "name": "Triphenylphosphine",
          "role": "ligand",
          "amount": {"equivalents": 0.02}
        },
        // ... more chemicals in addition order
      ],
      "conditions": [
        {
          "temperature_C": 60,
          "time_h": 16,
          "atmosphere": "N2"
        }
      ]
    }
  ]
}
```

**✅ This structure has everything automation needs:**
- Explicit addition order (array position)
- Role-based organization
- Precise amounts with units
- Reaction conditions grouped together

## Current Rule Structure (Simple & Flexible)

```json
{
  "schema_version": "2.0",
  "source_type": "rule",
  "default_rule": {
    "conditions": {
      "catalyst": "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM",
      "catalyst_loading_molpct": "0.5-2.0",
      "base": "Et3N or DIPEA",
      "base_equiv": "2.0-3.0",
      "solvent": "THF, toluene, or DMF",
      "temperature_C": "40-80",
      "time_h": "1-8",
      "atmosphere": "N2 or Ar",
      "additives": ["CuI (2-5 mol%)"]
    }
  }
}
```

**✅ Advantages of keeping this:**
- Human-readable
- Easy to author
- Flexible (ranges, options)
- No redundancy

## Proposed Solution: Dynamic Format Conversion

### Strategy

**At recommendation time**, convert rule conditions → protocol-like `reaction_setup`:

```python
# chemtools/formatters/rule_to_protocol.py

def rule_conditions_to_reaction_setup(
    conditions: Dict[str, Any],
    user_substrates: Optional[List[Dict[str, Any]]] = None,
    scale_mmol: float = 1.0
) -> Dict[str, Any]:
    """
    Convert rule conditions to protocol-like reaction_setup.
    
    Borrows protocol structure but generates from rule data.
    """
    chemicals = []
    
    # Standard addition order (like protocols)
    # 1. Solvent
    if "solvent" in conditions:
        chemicals.append({
            "name": _pick_first_option(conditions["solvent"]),
            "abbreviation": _get_abbreviation(conditions["solvent"]),
            "role": "solvent",
            "amount": {
                "volume_ml": None,  # To be calculated
                "note": "Volume based on target concentration"
            }
        })
    
    # 2. Base
    if "base" in conditions:
        base_equiv = _parse_range_midpoint(conditions.get("base_equiv", "2.0"))
        chemicals.append({
            "name": _pick_first_option(conditions["base"]),
            "abbreviation": _get_abbreviation(conditions["base"]),
            "role": "base",
            "amount": {
                "mmol": base_equiv * scale_mmol,
                "equivalents": base_equiv
            }
        })
    
    # 3. Ligand (if separate from catalyst)
    if "ligand" in conditions and "built-in" not in conditions["ligand"].lower():
        ligand_loading = _parse_range_midpoint(conditions.get("ligand_loading_molpct", "2.0"))
        chemicals.append({
            "name": _pick_first_option(conditions["ligand"]),
            "role": "ligand",
            "amount": {
                "mmol": (ligand_loading / 100.0) * scale_mmol,
                "equivalents": ligand_loading / 100.0
            }
        })
    
    # 4. Metal catalyst
    if "catalyst" in conditions:
        cat_loading = _parse_range_midpoint(conditions.get("catalyst_loading_molpct", "1.0"))
        chemicals.append({
            "name": _pick_first_option(conditions["catalyst"]),
            "role": "metal_catalyst",
            "amount": {
                "mmol": (cat_loading / 100.0) * scale_mmol,
                "equivalents": cat_loading / 100.0
            }
        })
    
    # 5. Additives
    if "additives" in conditions:
        for additive in conditions["additives"]:
            chemicals.append({
                "name": additive,
                "role": "additive",
                "amount": {
                    "note": "As specified in rule"
                }
            })
    
    # 6. User substrates (if provided)
    if user_substrates:
        for substrate in user_substrates:
            chemicals.append({
                "name": substrate.get("name", "Substrate"),
                "smiles": substrate.get("smiles"),
                "role": substrate.get("role", "starting_material"),
                "amount": {
                    "mmol": substrate.get("mmol", scale_mmol),
                    "equivalents": substrate.get("equivalents", 1.0)
                }
            })
    else:
        # Placeholder
        chemicals.append({
            "name": "Substrate (user-provided)",
            "role": "starting_material",
            "amount": {
                "mmol": scale_mmol,
                "equivalents": 1.0
            }
        })
    
    # Reaction conditions (like protocol format)
    reaction_conditions = {
        "temperature_C": _parse_range_midpoint(conditions.get("temperature_C", "25")),
        "time_h": _parse_range_midpoint(conditions.get("time_h", "4")),
        "atmosphere": conditions.get("atmosphere", "N2 or Ar")
    }
    
    return {
        "reaction_setup": [
            {
                "chemicals": chemicals,
                "conditions": [reaction_conditions]
            }
        ],
        "metadata": {
            "generated_from": "rule",
            "format": "protocol-compatible",
            "scale_mmol": scale_mmol
        }
    }
```

## Standardization Plan - REVISED

### Phase 1: Quick Fixes (30 minutes) ✅

Fix the one file with mixed types:

```bash
python fix_suzuki_types.py
```

**Result:** All 9 rule files now have consistent string types.

### Phase 2: Adopt Protocol Field Naming (1 hour)

Create mapping between rule fields and protocol structure:

```python
# chemtools/formatters/field_mapping.py

PROTOCOL_ROLE_MAPPING = {
    # Protocol roles (from reaction_setup)
    "starting_material": ["substrate", "electrophile", "nucleophile"],
    "metal_catalyst": ["catalyst", "pd_precatalyst", "pd_source", "ru_catalyst", "cu_source"],
    "ligand": ["ligand"],
    "base": ["base"],
    "solvent": ["solvent"],
    "additive": ["additives", "additive"],
    "reagent": ["coupling_system", "reducing_agent"],
}

PROTOCOL_AMOUNT_FIELDS = {
    # Protocol amount structure
    "weight_g": None,
    "mmol": None,
    "volume_ml": None,
    "equivalents": None
}

# Mapping rule fields to protocol amount
RULE_TO_PROTOCOL_AMOUNT = {
    "catalyst_loading_molpct": "equivalents",  # Convert mol% to equivalents
    "ligand_loading_molpct": "equivalents",
    "base_equiv": "equivalents",
    # ... etc
}
```

### Phase 3: Create Documentation (1 hour)

#### File: `data/rule_db_v2/AUTHORING_GUIDE.md`

```markdown
# Rule Authoring Guide

## Overview

Rules use a simplified format for easy authoring. When output to users or 
automation systems, they are converted to protocol-like `reaction_setup` format.

## Rule Format (What You Write)

```json
{
  "conditions": {
    "catalyst": "PdCl2(PPh3)2",
    "catalyst_loading_molpct": "1.0-3.0",
    "base": "Et3N",
    "base_equiv": "2.0",
    "solvent": "THF",
    "temperature_C": "60-80"
  }
}
```

**Field Naming Standards:**
- Use protocol-compatible names where possible
- Add units in suffix: `_molpct`, `_equiv`, `_C`, `_h`
- All numeric fields as strings (even single values)
- Ranges: "min-max" format
- Options: " or " separator

## Output Format (What Users See)

When `format_for_automation=True`, rules are converted to:

```json
{
  "reaction_setup": [
    {
      "chemicals": [
        {
          "name": "THF",
          "role": "solvent",
          "amount": {"volume_ml": null}
        },
        {
          "name": "Et3N",
          "role": "base",
          "amount": {"equivalents": 2.0, "mmol": 2.0}
        },
        {
          "name": "PdCl2(PPh3)2",
          "role": "metal_catalyst",
          "amount": {"equivalents": 0.02, "mmol": 0.02}
        }
      ],
      "conditions": [
        {"temperature_C": 70, "time_h": 4}
      ]
    }
  ]
}
```

Same structure as protocols - ready for automation!

## Role Mapping

| Rule Field | Protocol Role |
|------------|---------------|
| `catalyst`, `pd_precatalyst`, etc. | `metal_catalyst` |
| `ligand` | `ligand` |
| `base` | `base` |
| `solvent` | `solvent` |
| `additives` | `additive` |
| User substrates | `starting_material` |

## Standard Fields

### Core Fields (Most Rules)
- `solvent`: String with options
- `temperature_C`: String, range or single
- `time_h`: String, range or single
- `atmosphere`: String (e.g., "N2 or Ar")

### Catalyst Fields
- `catalyst` / `pd_precatalyst` / `ru_catalyst` / `cu_source`
- `catalyst_loading_molpct`: String
- `ligand`: String (if not built-in)
- `ligand_loading_molpct`: String

### Other Common Fields
- `base`: String
- `base_equiv`: String
- `additives`: Array of strings
- `notes`: String
```

### Phase 4: Implement Converter (2 hours)

Create `chemtools/formatters/rule_to_protocol.py` with full implementation.

### Phase 5: Integration (1 hour)

Update `UnifiedRecommender.recommend()`:

```python
def recommend(
    self,
    reaction_smiles: str,
    top_k: int = 10,
    format_for_automation: bool = False,
    user_substrates: Optional[List[Dict]] = None,
    scale_mmol: float = 1.0,
    ...
) -> List[Dict[str, Any]]:
    """..."""
    
    results = self._search(...)
    
    if format_for_automation:
        from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup
        
        for result in results:
            if result["source_type"] == "rule":
                # Load full rule data
                rule_data = self._load_source(result["source_file"])
                conditions = self._get_conditions(rule_data)
                
                # Convert to protocol format
                result["reaction_setup"] = rule_conditions_to_reaction_setup(
                    conditions=conditions,
                    user_substrates=user_substrates,
                    scale_mmol=scale_mmol
                )
                result["format"] = "protocol-compatible"
            elif result["source_type"] == "protocol":
                # Already has reaction_setup
                result["format"] = "native-protocol"
    
    return results
```

## Benefits of This Approach

### 1. Reuse Proven Format ✅
- Protocols already validated for automation
- No need to invent new structure
- Users see consistent format

### 2. Simple Rule Authoring ✅
- Keep current simple dictionary format
- Easy to read and edit
- No redundancy (don't repeat role info)

### 3. Unified Output ✅
Both protocols and rules produce:
```json
{
  "reaction_setup": [...],  // Standard structure
  "source_type": "protocol" | "rule",
  "format": "native-protocol" | "protocol-compatible"
}
```

### 4. Minimal Changes ✅
- Fix 1 file (Suzuki) ← 10 minutes
- Add converter ← 2 hours
- Update recommender ← 1 hour
- Documentation ← 1 hour
**Total: 4-5 hours**

### 5. Future-Proof ✅
- Protocol format is already stable
- Any improvements to protocol format automatically benefit rules
- Clear migration path if we want to make rules more protocol-like

## Comparison with Original Proposals

### Original Proposal
- Create new `addition_sequence` format
- Custom field ordering logic
- Separate from protocols

### New Proposal (Adopt Protocol Format)
- ✅ Reuse existing `reaction_setup` format
- ✅ Leverage protocol structure
- ✅ Single format for all outputs
- ✅ Less code to maintain
- ✅ Proven for automation

## Implementation Timeline

### Week 1: Foundation
- [x] Analysis complete (done ✅)
- [ ] Fix Suzuki_db.json types (30 min)
- [ ] Create field mapping (1 hour)
- [ ] Write AUTHORING_GUIDE.md (1 hour)

### Week 2: Implementation
- [ ] Implement `rule_to_protocol.py` converter (2 hours)
- [ ] Update UnifiedRecommender (1 hour)
- [ ] Write tests (2 hours)

### Week 3: Integration
- [ ] Update chem_assistant wrapper (30 min)
- [ ] Add examples to README (1 hour)
- [ ] Validate on all 9 rule files (1 hour)

**Total Time: ~10 hours** (vs 15+ hours for custom format)

## Next Steps

1. **Approve approach**: Using protocol `reaction_setup` format?
2. **Fix Suzuki file**: Run `fix_suzuki_types.py`
3. **Implement converter**: Create `rule_to_protocol.py`
4. **Test**: Verify output matches protocol structure
5. **Ship**: Update recommender and tools

## Questions?

1. Should generated `reaction_setup` include CAS numbers?
   - **Recommendation:** No, keep it simple initially. Can add lookup later.

2. Should we calculate exact volumes/weights?
   - **Recommendation:** Provide equivalents/mmol, let automation system calculate based on its constraints.

3. Handle multi-stage reactions (like protocols can)?
   - **Recommendation:** Start with single stage, add later if needed.

4. Support workup/purification sections?
   - **Recommendation:** No - that's protocol-specific. Rules just provide reaction conditions.

---

**This approach gives you automation-ready output while keeping rule authoring simple!** 🎯
