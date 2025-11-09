# Automation-Ready Condition Format

## Overview

The condition recommendation system now supports **automation-ready output format** with clear, ordered addition sequences suitable for automated execution workflows.

## Key Features

### 1. Protocol-Compatible Structure

Both rules and protocols now produce the same standardized format:

```python
{
  "reaction_setup": [
    {
      "chemicals": [
        {"name": "THF", "role": "solvent", "amount": {...}},
        {"name": "Et3N", "role": "base", "amount": {...}},
        {"name": "Pd(PPh3)4", "role": "metal_catalyst", "amount": {...}},
        ...
      ],
      "conditions": [
        {"temperature_C": 80, "time_h": 4, "atmosphere": "N2"}
      ]
    }
  ],
  "metadata": {
    "generated_from": "rule",  # or "protocol"
    "format": "protocol-compatible",
    "scale_mmol": 1.0
  }
}
```

### 2. Standard Addition Order

Chemicals are ordered by standard laboratory addition sequence:

1. **Solvent** - Added first to establish reaction medium
2. **Base** - Added to prepare reaction conditions
3. **Ligand** - Coordinates with catalyst
4. **Metal Catalyst** - Active catalytic species
5. **Catalyst** - Other catalysts (if not metal-based)
6. **Additive** - Additional reagents
7. **Starting Material** - Substrates (user-provided or placeholder)
8. **Reagent** - Other reagents

### 3. Intelligent Value Parsing

The converter handles common rule formats automatically:

- **Ranges**: `"60-80"` → `70.0` (midpoint)
- **Options**: `"THF or toluene"` → `"THF"` (first option)
- **Equivalents**: Calculated based on scale_mmol
- **Mol%**: Converted to equivalents (e.g., `"1.5 mol%"` → `0.015 equiv`)

## Usage

### In Python

```python
from chemtools.recommend.unified import UnifiedRecommender

# Initialize recommender
recommender = UnifiedRecommender()

# Get recommendations with automation format
results = recommender.recommend(
    reaction_smiles="CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    top_k=3,
    format_for_automation=True,
    scale_mmol=1.0  # Reaction scale
)

# Access ordered chemicals
for result in results:
    if result.full_data:
        setup = result.full_data["reaction_setup"][0]
        
        print(f"Addition sequence for {result.name}:")
        for i, chem in enumerate(setup["chemicals"], 1):
            print(f"  {i}. {chem['name']} ({chem['role']})")
            if "equivalents" in chem["amount"]:
                print(f"     Amount: {chem['amount']['equivalents']} equiv")
```

### With LangChain Agent

```python
from chem_assistant.chemtools_wrapper import unified_recommender_tool

# Use in agent workflow
result = unified_recommender_tool(
    reaction_smiles="CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    top_k=3,
    format_for_automation=True,
    scale_mmol=1.0
)

# Result includes ordered reaction_setup for each recommendation
for rec in result["recommendations"]:
    if "reaction_setup" in rec:
        chemicals = rec["reaction_setup"][0]["chemicals"]
        print(f"Addition order: {[c['name'] for c in chemicals]}")
```

## Direct Converter Usage

For custom workflows, use the converter directly:

```python
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

# Load rule conditions
conditions = {
    "catalyst": "Pd(PPh3)4",
    "catalyst_loading_molpct": "1.0",
    "base": "K2CO3",
    "base_equiv": "2.0",
    "solvent": "DMF",
    "temperature_C": "80",
    "time_h": "4"
}

# Convert to automation format
result = rule_conditions_to_reaction_setup(
    conditions=conditions,
    scale_mmol=1.0,
    reaction_family="Suzuki_Miyaura"
)

# Use result for automation
setup = result["reaction_setup"][0]
for chem in setup["chemicals"]:
    print(f"Add {chem['name']}: {chem['amount']}")
```

## Benefits

### 1. Consistency

- Same format whether source is rule or protocol
- Predictable structure for automation systems
- Standard field names across all families

### 2. Automation-Ready

- Ordered addition sequences
- Calculated amounts and equivalents
- Units clearly specified
- Compatible with robotic systems

### 3. Flexible

- Handles multiple catalyst types (Pd, Cu, Ru, etc.)
- Supports ranges and alternatives
- Scales to any reaction size
- Extensible to new reaction families

### 4. Simple Authoring

- Rule authors continue using simple key-value format
- No need to specify addition order manually
- Conversion happens automatically on output
- Existing rules work without modification

## Technical Details

### Converter: `chemtools/formatters/rule_to_protocol.py`

**Main function:**
```python
def rule_conditions_to_reaction_setup(
    conditions: Dict[str, Any],
    user_substrates: Optional[List[Dict]] = None,
    scale_mmol: float = 1.0,
    reaction_family: str = "Unknown"
) -> Dict[str, Any]:
    """Convert rule conditions to protocol-compatible reaction_setup format."""
```

**Key mappings:**

```python
# Addition order priorities
ADDITION_ORDER_PRIORITY = {
    "solvent": 1,
    "base": 2,
    "ligand": 3,
    "metal_catalyst": 4,
    "catalyst": 5,
    "additive": 6,
    "starting_material": 7,
    "reagent": 8
}

# Field-to-role mapping
RULE_FIELD_TO_ROLE = {
    "catalyst": "metal_catalyst",
    "pd_precatalyst": "metal_catalyst",
    "pd_source": "metal_catalyst",
    "ru_catalyst": "metal_catalyst",
    "cu_source": "metal_catalyst",
    "base": "base",
    "ligand": "ligand",
    "solvent": "solvent",
    "additives": "additive"
}
```

### Integration Points

1. **UnifiedRecommender** (`chemtools/recommend/unified.py`):
   - New `format_for_automation` parameter
   - Automatically loads and converts rule conditions
   - Protocols pass through with existing structure

2. **LangChain Wrapper** (`chem_assistant/chemtools_wrapper.py`):
   - `unified_recommender_tool` supports automation format
   - Includes `reaction_setup` in response when enabled

## Example Output

### Suzuki-Miyaura (Rule → Protocol)

**Input (Rule):**
```json
{
  "pd_source": "PdCl2(dtbpf)",
  "ligand": "dtbpf",
  "catalyst_loading_molpct": "1.5",
  "base": "K3PO4 (aq, 3.25 M)",
  "base_equiv": "2.0",
  "solvent": "1,4-dioxane/H2O 4:1",
  "temperature_C": "80–100",
  "time_h": "1–6"
}
```

**Output (Protocol-Compatible):**
```json
{
  "reaction_setup": [{
    "chemicals": [
      {"name": "1,4-dioxane/H2O 4:1", "role": "solvent", "amount": {...}},
      {"name": "K3PO4 (aq, 3.25 M)", "role": "base", "amount": {"equivalents": 2.0, "mmol": 2.0}},
      {"name": "dtbpf", "role": "ligand", "amount": {"equivalents": 0.02, "mmol": 0.02}},
      {"name": "PdCl2(dtbpf)", "role": "metal_catalyst", "amount": {"equivalents": 0.015, "mmol": 0.015}},
      {"name": "Substrate (user-provided)", "role": "starting_material", "amount": {"mmol": 1.0}}
    ],
    "conditions": [
      {"temperature_C": 90.0, "time_h": 3.5, "atmosphere": "N2 or Ar"}
    ]
  }],
  "metadata": {
    "generated_from": "rule",
    "format": "protocol-compatible",
    "scale_mmol": 1.0
  }
}
```

## Testing

Run tests to verify functionality:

```bash
# Test converter with all rule files
python test_rule_to_protocol.py

# Test end-to-end conversion
python test_simple_conversion.py

# Test UnifiedRecommender integration (requires DRFP + built index)
python test_unified_automation.py
```

All tests should pass with ✅.

## Future Extensions

1. **User Substrates**: Specify actual substrates to replace placeholders
2. **Volume Calculations**: Auto-calculate solvent volumes from concentration
3. **Stock Solutions**: Handle pre-made reagent stock solutions
4. **Multi-Step**: Support sequential addition steps
5. **Equipment**: Specify vessels, stirring, etc.

## References

- Protocol format: `data/protocol_db_v2/*.json`
- Rule format: `data/rule_db_v2/*.json`
- Converter implementation: `chemtools/formatters/rule_to_protocol.py`
- Integration: `chemtools/recommend/unified.py`
