##LLM Synthesis Standard Format for Robotic Execution

## Overview

The LLM-enhanced recommendation system now outputs in **two formats**:

1. **Analysis Format** (`*_llm_analysis.json`) - Rich LLM synthesis with consensus analysis, warnings, and reasoning
2. **Standard Format** (`*_llm_local.json`) - Standardized structure for robotic execution systems

This ensures LLM recommendations are **compatible with the same robotic workflows** as ML, Rule-based, and Protocol recommendations.

## Problem Solved

**Before**: LLM synthesis returned only analysis-focused JSON:
```json
{
  "status": "success",
  "synthesis": {
    "consensus_analysis": {...},
    "recommended_condition": {...},
    "backup_conditions": [...]
  }
}
```

**Issue**: Robotic systems couldn't parse this format - they expect `meta`, `input`, `detection`, `recommended_conditions` structure.

**After**: LLM synthesis now returns **both** formats:
- Analysis format for human review
- Standard format for robotic execution

## Standard Format Structure

The standard format matches `output_formatter.py` specifications:

```json
{
  "meta": {
    "generated_at": "2025-10-12T14:23:17.107278Z",
    "model": "LLM-synthesis",
    "model_version": "1.0.0",
    "status": "success",
    "processing_time_ms": 45823.0
  },
  "input": {
    "reaction_smiles": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "requested_reaction_type": "Suzuki_Miyaura"
  },
  "detection": {
    "family": "Suzuki_Miyaura",
    "detected_reaction_type": "Suzuki_Miyaura",
    "method": "llm-multi-source",
    "confidence": 0.95
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "chemicals": [
        {
          "name": null,
          "abbreviation": null,
          "cas": null,
          "smiles": "Brc1ccccc1",
          "equivalents": null,
          "role": "electrophile"
        },
        {
          "name": null,
          "abbreviation": null,
          "cas": null,
          "smiles": "c1ccc(B(O)O)cc1",
          "equivalents": null,
          "role": "nucleophile"
        },
        {
          "name": "Tetrakis(triphenylphosphine)palladium(0)",
          "abbreviation": "Pd(PPh3)4",
          "cas": "14221-01-3",
          "smiles": "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.[Pd]",
          "equivalents": null,
          "role": "metal_precursor"
        },
        {
          "name": "Triphenylphosphine",
          "abbreviation": "PPh3",
          "cas": "603-35-0",
          "smiles": "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3",
          "equivalents": null,
          "role": "ligand"
        },
        {
          "name": "Potassium carbonate",
          "abbreviation": "K2CO3",
          "cas": "584-08-7",
          "smiles": "[K+].[K+].[O-]C(=O)[O-]",
          "equivalents": 2.0,
          "role": "base"
        },
        {
          "name": "Tetrahydrofuran",
          "abbreviation": "THF",
          "cas": "109-99-9",
          "smiles": "C1CCOC1",
          "equivalents": null,
          "role": "solvent"
        }
      ],
      "conditions": {
        "temperature": [80.0, 80.0],
        "time": [6.0, 24.0],
        "atmosphere": null
      },
      "confidence": 0.95,
      "summary": {
        "rationale": "High consensus across all sources for standard Suzuki coupling with electron-neutral aryl bromide.",
        "warnings": [
          "Monitor for proto-debromination if reaction stalls",
          "Degassed solvent recommended for Pd(0) catalyst",
          "Boronic acid may oxidize if exposed to air for prolonged periods"
        ]
      },
      "source": {
        "method": "llm-multi-source-synthesis",
        "confidence_level": "high"
      }
    },
    {
      "rank": 2,
      "chemicals": [...],
      "conditions": {...},
      "confidence": 0.85,
      "summary": {
        "rationale": "If Pd(PPh3)4 gives <30% conversion after 12h at 80°C",
        "warnings": [...]
      },
      "source": {
        "method": "llm-multi-source-synthesis",
        "confidence_level": "high"
      }
    },
    {
      "rank": 3,
      "chemicals": [...],
      "conditions": {...},
      "confidence": 0.75,
      "summary": {
        "rationale": "If backup 1 fails (<20% after 12h) - more active catalyst system",
        "warnings": [...]
      },
      "source": {
        "method": "llm-multi-source-synthesis",
        "confidence_level": "high"
      }
    }
  ],
  "extras": {
    "llm_synthesis": {
      "consensus_analysis": {...},
      "confidence_level": "high",
      "confidence_reasoning": "...",
      "recommended_condition": {...},
      "backup_conditions": [...],
      "warnings": [...],
      "source_comparison": {...}
    },
    "sources_used": {
      "ml_precedents": 3,
      "rule_matches": 1,
      "protocol_procedures": 2
    },
    "llm_metadata": {
      "model": "deepseek-v3.2-exp",
      "tokens": 1480,
      "latency_ms": 43600,
      "processing_time_ms": 45823
    }
  }
}
```

## Key Features

### 1. Reagent Enrichment

All chemicals are enriched with database information:

```python
{
  "name": "Tetrakis(triphenylphosphine)palladium(0)",
  "abbreviation": "Pd(PPh3)4",
  "cas": "14221-01-3",
  "smiles": "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.[Pd]",
  "equivalents": 0.05,  # Auto-calculated for catalysts
  "role": "metal_precursor"
}
```

The `reagent_lookup` module provides:
- Full IUPAC names
- CAS numbers
- SMILES strings
- Typical equivalents based on role

### 2. Multiple Recommendations

LLM synthesis includes:
- **Primary recommendation** (rank 1): Best consensus from all sources
- **Backup condition 1** (rank 2): Alternative if primary fails
- **Backup condition 2** (rank 3): Different catalyst family if backups fail

Each recommendation has:
- Complete chemical list with enrichment
- Conditions (temperature, time, atmosphere)
- Confidence score
- Rationale (when to use)
- Warnings

### 3. Preserved Analysis Data

Original LLM synthesis stored in `extras`:

```json
{
  "extras": {
    "llm_synthesis": {
      "consensus_analysis": {...},
      "confidence_reasoning": "All sources agree on Pd(PPh3)4...",
      "warnings": [...],
      "source_comparison": {
        "ml_contribution": "High - similarity 0.95...",
        "rule_contribution": "Primary - high-priority SCDB match...",
        "protocol_contribution": "Supporting - confirms 80°C..."
      }
    },
    "sources_used": {...},
    "llm_metadata": {...}
  }
}
```

This allows:
- Human review of reasoning
- Debugging recommendation logic
- Understanding source contributions
- Audit trail for decision-making

## Implementation

### New Function: `convert_llm_synthesis_to_standard_format()`

**Location**: `llmtools/recommendation_llm.py`

**Signature**:
```python
def convert_llm_synthesis_to_standard_format(
    reaction_smiles: str,
    synthesis_result: Dict[str, Any],
    requested_type: Optional[str] = None,
    processing_time_ms: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Convert LLM synthesis result to standard output format for robotic execution.
    
    Transforms the analysis-focused LLM output into the same structure as
    ML/Rule/Protocol outputs, making it compatible with robotic systems.
    """
```

**Process**:

1. **Extract LLM Analysis**
   ```python
   synthesis = synthesis_result.get("synthesis", {})
   recommended_condition = synthesis.get("recommended_condition", {})
   backup_conditions = synthesis.get("backup_conditions", [])
   ```

2. **Build Chemical Entries**
   ```python
   # Start with reactants
   chemicals = _starting_material_entries(reaction_smiles)
   
   # Add catalyst system with enrichment
   catalyst_name = recommended_condition.get("catalyst")
   chemicals.append(_normalize_rule_string_value(catalyst_name, "metal_precursor"))
   ```

3. **Parse Conditions**
   ```python
   # Handle formats like "80°C", "80-90°C", "80-100 °C"
   temp_match = re.search(r'(\d+)(?:\s*-\s*(\d+))?\s*°?C?', temperature)
   temp_value = (temp_low + temp_high) / 2
   temp_range = [temp_low, temp_high]
   ```

4. **Create Multiple Recommendations**
   ```python
   recommendations = []
   
   # Primary (rank 1)
   primary_rec = _build_recommendation_from_llm(
       rank=1,
       condition_dict=recommended_condition,
       confidence=0.95,
       reasoning=synthesis.get("confidence_reasoning"),
   )
   
   # Backups (rank 2, 3)
   for i, backup in enumerate(backup_conditions[:2], start=2):
       backup_rec = _build_recommendation_from_llm(rank=i, ...)
   ```

5. **Preserve Original Data**
   ```python
   output["extras"] = {
       "llm_synthesis": synthesis,
       "sources_used": sources_used,
       "llm_metadata": llm_metadata,
   }
   ```

### Helper Function: `_build_recommendation_from_llm()`

**Responsibilities**:
- Extract chemicals from LLM condition dict
- Enrich with database info (CAS, SMILES, etc.)
- Parse temperature ranges
- Build conditions block
- Add rationale and warnings

**Example**:
```python
recommendation = _build_recommendation_from_llm(
    rank=1,
    condition_dict={
        "catalyst": "Pd(PPh3)4",
        "ligand": "PPh3 (pre-complexed)",
        "solvent": "THF",
        "temperature": "80°C",
        "base": "K2CO3 (2 equiv)",
        "rationale": "High consensus..."
    },
    reaction_smiles="Brc1ccccc1...",
    confidence=0.95,
    warnings=["Monitor for proto-debromination..."]
)
```

**Output**:
```python
{
    "rank": 1,
    "chemicals": [
        # Reactants
        {"smiles": "Brc1ccccc1", "role": "electrophile"},
        {"smiles": "c1ccc(B(O)O)cc1", "role": "nucleophile"},
        # Enriched catalyst system
        {"name": "Tetrakis(triphenylphosphine)palladium(0)", "abbreviation": "Pd(PPh3)4", "cas": "14221-01-3", ...},
        {"name": "Triphenylphosphine", "abbreviation": "PPh3", "cas": "603-35-0", ...},
        {"name": "Potassium carbonate", "abbreviation": "K2CO3", "equivalents": 2.0, ...},
        {"name": "Tetrahydrofuran", "abbreviation": "THF", ...},
    ],
    "conditions": {
        "temperature": [80.0, 80.0],
        "time": [6.0, 24.0],
        "atmosphere": None
    },
    "confidence": 0.95,
    "summary": {
        "rationale": "High consensus...",
        "warnings": ["Monitor for proto-debromination..."]
    }
}
```

## CLI Tool Updates

### Modified: `local_llm_synthesis()`

**Old Return**:
```python
return result  # Single dict
```

**New Return**:
```python
return analysis_result, standard_result  # Tuple of two dicts
```

**Usage**:
```python
llm_analysis, llm_standard = local_llm_synthesis(
    reaction=reaction,
    ml_result=ml_result,
    rule_result=rule_result,
    protocol_result=protocol_result,
    constraints=constraints,
    requested_type=args.rxn_type,  # NEW: For standard format
)
```

### File Outputs

**Two Files Saved**:

1. `{timestamp}_{label}_llm_analysis.json`
   - Rich LLM synthesis
   - Consensus analysis
   - Source comparison
   - Warnings and reasoning
   - **For human review**

2. `{timestamp}_{label}_llm_local.json`
   - Standard format
   - Enriched chemicals
   - Multiple recommendations
   - Robotic-compatible structure
   - **For robotic execution**

**Console Output**:
```
Saved outputs:
  Rule JSON:            results/20251012_143821_suzuki_rule_local.json
  ML JSON:              results/20251012_143821_suzuki_ml_local.json
  Protocol JSON:        results/20251012_143821_suzuki_protocol_local.json
  LLM Analysis JSON:    results/20251012_143821_suzuki_llm_analysis.json
  LLM Standard JSON:    results/20251012_143821_suzuki_llm_local.json  (for robotic execution)
```

## Testing

### Test Suite: `tests/test_llm_standard_format.py`

**Coverage**:
1. ✅ Structure validation (meta, input, detection, recommended_conditions)
2. ✅ Chemical enrichment (CAS, SMILES, names)
3. ✅ Multiple recommendations (primary + backups)
4. ✅ Conditions parsing (temperature ranges, time)
5. ✅ Confidence scores
6. ✅ Rationale and warnings
7. ✅ Compatibility with `ensure_standard_output()`
8. ✅ Robotic execution readiness

**Run**:
```powershell
python tests/test_llm_standard_format.py
```

**Expected Output**:
```
======================================================================
Testing LLM Synthesis Standard Format Conversion
======================================================================

1. Converting LLM synthesis to standard format...
✅ Conversion successful

2. Verifying standard structure...
✅ All required keys present: {'meta', 'recommended_conditions', 'detection', 'input'}
✅ Meta section valid
✅ Input section valid
✅ Detection section valid (confidence: 0.95)
✅ Recommendations present: 3 recommendation(s)
✅ Recommendation 1 has required keys: {'conditions', 'chemicals', 'rank'}
✅ Chemicals array present: 6 chemical(s)
   - Catalyst: Pd(PPh3)4
✅ Conditions present: {'temperature': [80.0, 80.0], 'time': [6.0, 24.0], 'atmosphere': None}
✅ Original synthesis data preserved in extras

3. Checking robotic execution compatibility...
✅ Output passes ensure_standard_output validation

4. Robotic execution readiness:
   Reaction: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1
   Family: Suzuki_Miyaura
   Confidence: 95.00%
   Recommendations: 3

   Recommendation 1:
     Rank: 1
     Chemicals: 6
       - Metal_precursor: Pd(PPh3)4
       - Ligand: Triphenylphosphine
       - Base: Potassium carbonate
       - Solvent: Tetrahydrofuran
       - Temperature: [80.0, 80.0]
       - Rationale: High consensus across all sources...

======================================================================
✅ All tests passed! LLM synthesis standard format ready for robots.
======================================================================
```

## Usage Examples

### Example 1: Basic LLM Synthesis

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --rxn-type suzuki
```

**Output**:
- `*_llm_analysis.json` - LLM reasoning and warnings
- `*_llm_local.json` - **Standard format for robot** ✅

### Example 2: With Constraints

```powershell
python scripts/local_recommendation_cli.py `
  --strategy llm `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --rxn-type suzuki `
  --constraints '{\"cost\": \"low\", \"scale\": \"multigram\"}'
```

**Result**: LLM considers constraints and outputs cost-effective conditions in standard format.

### Example 3: Compare All Methods

```powershell
python scripts/local_recommendation_cli.py `
  --strategy all `
  --rxn "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1" `
  --rxn-type suzuki
```

**Output Files**:
- `*_rule_local.json` - Rule-based (standard format) ✅
- `*_ml_local.json` - ML precedent (standard format) ✅
- `*_fusion_local.json` - Fusion (standard format) ✅
- `*_protocol_local.json` - Protocol (standard format) ✅
- `*_llm_analysis.json` - LLM analysis
- `*_llm_local.json` - **LLM synthesis (standard format)** ✅

**All 5 methods now output in the same format for robotic systems!**

## Robotic System Integration

### Workflow

1. **Generate Recommendations**
   ```python
   # Run LLM synthesis
   analysis, standard = local_llm_synthesis(...)
   
   # Save standard format for robot
   with open("robot_conditions.json", "w") as f:
       json.dump(standard, f)
   ```

2. **Robot Reads Standard Format**
   ```python
   # Robotic controller
   with open("robot_conditions.json", "r") as f:
       conditions = json.load(f)
   
   # Extract primary recommendation
   primary = conditions["recommended_conditions"][0]
   ```

3. **Execute Reaction**
   ```python
   # Prepare reagents
   for chemical in primary["chemicals"]:
       role = chemical["role"]
       cas = chemical["cas"]
       equivalents = chemical.get("equivalents", 1.0)
       
       if role == "electrophile":
           robot.load_reagent(position=1, cas=cas, amount=1.0)
       elif role == "metal_precursor":
           robot.load_reagent(position=3, cas=cas, amount=equivalents)
       # ... etc
   
   # Set conditions
   temp_range = primary["conditions"]["temperature"]
   robot.set_temperature(temp_range[0])
   
   # Run reaction
   robot.start_reaction()
   ```

4. **Handle Failures with Backups**
   ```python
   # Check yield
   yield_pct = robot.measure_yield()
   
   if yield_pct < 30:
       # Try backup condition 1
       backup = conditions["recommended_conditions"][1]
       robot.setup_new_reaction(backup)
   ```

## Benefits

### 1. **Universal Compatibility**
- LLM output now matches ML/Rule/Protocol format
- Single parser for all recommendation methods
- Robotic systems don't need LLM-specific code

### 2. **Reagent Enrichment**
- Automatic CAS lookup for reagent ordering
- SMILES for molecular weight calculations
- Abbreviations for laboratory labels

### 3. **Multiple Backup Strategies**
- Primary recommendation (best consensus)
- Backup 1 (alternative conditions)
- Backup 2 (different catalyst family)
- Clear "when to use" guidance

### 4. **Audit Trail**
- Original LLM synthesis preserved in `extras`
- Source contributions documented
- Confidence reasoning explained
- Warnings captured

### 5. **Decision Tree Support**
- Recommendations ranked by confidence
- Rationale explains when to escalate
- Warnings flag potential issues
- Robotic systems can auto-escalate

## Migration Guide

### For Existing Code

**Before** (old format):
```python
result = local_llm_synthesis(reaction, ...)
recommended = result["synthesis"]["recommended_condition"]
catalyst = recommended["catalyst"]  # Just string
```

**After** (dual format):
```python
analysis, standard = local_llm_synthesis(reaction, ...)

# Use standard format for robot
primary = standard["recommended_conditions"][0]
for chem in primary["chemicals"]:
    if chem["role"] == "metal_precursor":
        catalyst_cas = chem["cas"]  # Full CAS number
        catalyst_smiles = chem["smiles"]  # SMILES string
        break

# Use analysis for debugging
print(analysis["synthesis"]["confidence_reasoning"])
```

### For Robotic Systems

**No Changes Required!**

If your system already parses ML/Rule/Protocol outputs:
```python
def parse_recommendation(json_file):
    data = json.load(open(json_file))
    
    # Works for ALL methods now (ML, Rule, Protocol, LLM)
    primary = data["recommended_conditions"][0]
    chemicals = primary["chemicals"]
    conditions = primary["conditions"]
    
    return chemicals, conditions
```

LLM synthesis outputs are now drop-in compatible!

## Future Enhancements

### Planned Features

1. **Time Estimation**
   - Add typical reaction times to conditions
   - Based on precedent database
   - Include "check after X hours" guidance

2. **Yield Prediction**
   - Add `expected_outcome` section
   - Typical yield range from precedents
   - Based on similarity and consensus

3. **Scale Factors**
   - Add scaling guidance
   - Adjust catalyst loading for scale
   - Heat/mass transfer considerations

4. **Safety Data**
   - Add chemical hazard info
   - Equipment requirements
   - Protective measures needed

5. **Cost Estimation**
   - Add reagent costs
   - Total reaction cost estimate
   - Cost-effective alternatives

## Conclusion

The LLM synthesis system now provides **dual outputs**:

1. **Analysis Format** - Rich LLM reasoning for human review
2. **Standard Format** - Robotic-compatible structure for execution

This makes LLM-enhanced recommendations **fully compatible with existing robotic systems** while preserving the valuable LLM analysis and reasoning.

**Key Achievement**: All recommendation methods (ML, Rule, Protocol, LLM) now output in the same standard format, enabling universal robotic integration.

---

**Status**: ✅ **COMPLETE AND TESTED**

**Files Modified**:
- `llmtools/recommendation_llm.py` (+200 lines)
- `scripts/local_recommendation_cli.py` (+30 lines)
- `tests/test_llm_standard_format.py` (+250 lines, NEW)

**Documentation**:
- [LLM Standard Format Guide](./LLM_STANDARD_FORMAT.md) (this file)
- [CLI Usage Guide](./CLI_LLM_SYNTHESIS_USAGE.md)
- [Step 4 Completion](./STEP4_PROMPT_REFINEMENT_COMPLETE.md)

**Date**: October 12, 2025
