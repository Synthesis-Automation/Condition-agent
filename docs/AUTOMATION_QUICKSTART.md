# Quick Start: Automation-Ready Condition Format

## What's New?

Condition recommendations now include **ordered addition sequences** for automated execution!

## Simple Example

```python
from chemtools.recommend.unified import UnifiedRecommender

# Get recommendations with automation format
recommender = UnifiedRecommender()
results = recommender.recommend(
    reaction_smiles="CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    format_for_automation=True,  # ← Enable automation format
    scale_mmol=1.0
)

# Print addition sequence
for result in results:
    if result.full_data:
        setup = result.full_data["reaction_setup"][0]
        print(f"\n{result.name} - Addition sequence:")
        for i, chem in enumerate(setup["chemicals"], 1):
            print(f"  {i}. Add {chem['name']} ({chem['role']})")
```

## Output Example

```
Suzuki-Miyaura Coupling - Addition sequence:
  1. Add 1,4-dioxane/H2O 4:1 (solvent)
  2. Add K3PO4 (base)
  3. Add dtbpf (ligand)
  4. Add PdCl2(dtbpf) (metal_catalyst)
  5. Add Substrate (starting_material)
```

## Key Features

✅ **Standard Addition Order**: Solvent → Base → Ligand → Catalyst → Substrate  
✅ **Automatic Conversion**: Rules converted to protocol format  
✅ **Calculated Amounts**: Equivalents computed from scale  
✅ **Consistent Structure**: Same format for rules and protocols  

## Use Cases

### 1. Robot Integration

```python
# Get conditions for automated synthesis
results = recommender.recommend(
    reaction_smiles=your_reaction,
    format_for_automation=True,
    scale_mmol=5.0  # 5 mmol scale
)

# Send to robot controller
for chem in results[0].full_data["reaction_setup"][0]["chemicals"]:
    robot.add_chemical(
        name=chem["name"],
        amount_mmol=chem["amount"]["mmol"]
    )
```

### 2. LLM Agent

```python
from chem_assistant.chemtools_wrapper import unified_recommender_tool

# Agent can request automation format
result = unified_recommender_tool(
    reaction_smiles="CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    format_for_automation=True
)

# Result includes ordered sequence
print(result["recommendations"][0]["reaction_setup"])
```

### 3. Direct Converter

```python
from chemtools.formatters.rule_to_protocol import rule_conditions_to_reaction_setup

# Convert any rule to automation format
conditions = {
    "catalyst": "Pd(PPh3)4",
    "base": "K2CO3",
    "solvent": "DMF"
}

result = rule_conditions_to_reaction_setup(conditions, scale_mmol=1.0)
# → Protocol-compatible format with addition order
```

## When to Use

**Enable `format_for_automation=True` when:**
- Integrating with robotic systems
- Need explicit addition order
- Generating lab protocols
- Automating reaction setup

**Keep default (`False`) when:**
- Just browsing recommendations
- Only need similarity scores
- Want faster responses (no extra data loading)

## Learn More

- Full documentation: `docs/AUTOMATION_FORMAT.md`
- Test examples: `test_rule_to_protocol.py`, `test_simple_conversion.py`
- Rule format: `data/rule_db_v2/*.json`
- Protocol format: `data/protocol_db_v2/*.json`

## Quick Test

```bash
# Test the converter
python test_simple_conversion.py

# Should output:
# ✅ Has reaction_setup
# ✅ Setup has chemicals (5) and conditions
# ✅ Addition order: solvent → base → ligand → metal_catalyst → starting_material
# SUCCESS! Rule → Protocol conversion working perfectly! 🎉
```

---

**Ready to automate your chemistry!** 🤖⚗️
