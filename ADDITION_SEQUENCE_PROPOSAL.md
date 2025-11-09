# Addition Sequence for Automated Execution - Proposal

## Problem Statement

For automated robotic execution of recommended reaction conditions, we need a **well-defined addition sequence** with clear ordering and timing information. Currently:

- ✅ **Protocols** have full `reaction_setup` with ordered chemicals
- ❌ **Rules** only have condition dictionaries (unordered key-value pairs)

## Current State Analysis

### Protocol Format (Has Addition Sequence)
```json
{
  "reaction_setup": [
    {
      "chemicals": [
        {"name": "2-Iodobenzamide", "role": "starting_material", "amount": {...}},
        {"name": "Triphenylphosphine", "role": "ligand", "amount": {...}},
        {"name": "Copper(I) iodide", "role": "metal_catalyst", "amount": {...}},
        {"name": "Palladium(II) acetate", "role": "metal_catalyst", "amount": {...}},
        {"name": "N,N-Dimethylformamide", "role": "solvent", "amount": {...}},
        {"name": "Triethylamine", "role": "base", "amount": {...}},
        {"name": "1-Hexyne", "role": "starting_material", "amount": {...}}
      ],
      "conditions": [
        {"temperature_C": 60, "time_h": 16, "atmosphere": "N2"}
      ]
    }
  ]
}
```

**Advantages:**
- Complete addition order preserved from literature
- Exact amounts with units
- Role-based grouping
- Multiple setup stages supported
- Real experimental procedure captured

### Rule Format (No Addition Sequence)
```json
{
  "conditions": {
    "catalyst": "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM",
    "catalyst_loading_molpct": "0.5-2.0",
    "ligand": "Built-in (PPh3 or dppf)",
    "base": "Et3N or DIPEA",
    "base_equiv": "2.0-3.0",
    "solvent": "THF, toluene, or DMF",
    "temperature_C": "40-80",
    "time_h": "1-8",
    "atmosphere": "N2 or Ar; thoroughly degassed",
    "additives": ["CuI (2-5 mol%)"]
  }
}
```

**Challenges:**
- Unordered dictionary (JSON objects have no inherent order)
- Ranges instead of exact values ("0.5-2.0", "40-80")
- Multiple options ("PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM")
- No explicit addition order
- Ambiguous amounts ("2-5 mol%")

## Solution Options

### Option 1: Generate Synthetic `reaction_setup` from Rules ⭐ **RECOMMENDED**

**Concept:** Create a standardized addition sequence when formatting rule-based recommendations for output.

**Implementation:**
1. Add a post-processing step in `UnifiedRecommender.recommend()`
2. Convert rule conditions → protocol-like `reaction_setup` structure
3. Apply **role-based addition order heuristics**
4. Pick **midpoint values** from ranges
5. Select **first option** from choices

**Advantages:**
- ✅ Minimal schema changes (rules stay simple)
- ✅ Backward compatible (rules remain human-readable)
- ✅ Format conversion at output time only
- ✅ Can be improved iteratively
- ✅ Keeps rule authoring simple

**Role-Based Addition Order Template:**
```python
STANDARD_ADDITION_ORDER = [
    "solvent",           # 1. Add solvent first (establish medium)
    "base",              # 2. Add base (if needed for catalyst activation)
    "ligand",            # 3. Add ligand (before metal catalyst)
    "metal_catalyst",    # 4. Add metal catalyst (forms active complex)
    "catalyst",          # 5. Add other catalysts
    "additive",          # 6. Add additives (promoters, cocatalysts)
    "starting_material", # 7. Add starting materials (main substrates)
    "reagent",           # 8. Add reagents (coupling partners, etc.)
]
```

**Example Transformation:**

Input (Rule):
```json
{
  "catalyst": "PdCl2(PPh3)2",
  "catalyst_loading_molpct": "0.5-2.0",
  "base": "Et3N",
  "base_equiv": "2.0-3.0",
  "solvent": "THF",
  "temperature_C": "40-80"
}
```

Output (Generated `reaction_setup`):
```json
{
  "reaction_setup": [
    {
      "chemicals": [
        {
          "name": "THF",
          "role": "solvent",
          "amount": {"volume_ml": "to_be_scaled"},
          "addition_order": 1
        },
        {
          "name": "Et3N",
          "role": "base",
          "amount": {"equivalents": 2.5},  # midpoint of 2.0-3.0
          "addition_order": 2
        },
        {
          "name": "PdCl2(PPh3)2",
          "role": "metal_catalyst",
          "amount": {"mol_percent": 1.25},  # midpoint of 0.5-2.0
          "addition_order": 3
        },
        {
          "name": "Substrate 1",
          "role": "starting_material",
          "amount": {"equivalents": 1.0},
          "addition_order": 4,
          "note": "User-provided substrate"
        }
      ],
      "conditions": [
        {
          "temperature_C": 60,  # midpoint of 40-80
          "atmosphere": "N2 or Ar"
        }
      ]
    }
  ],
  "metadata": {
    "generated_from": "rule",
    "original_rule_id": "BR_ArI_ArBr_standard",
    "note": "Addition sequence generated from rule conditions using standard heuristics"
  }
}
```

---

### Option 2: Add `addition_sequence` Field to Rule Schema

**Concept:** Extend rule schema with optional explicit addition order.

```json
{
  "conditions": {
    "catalyst": "PdCl2(PPh3)2",
    "base": "Et3N",
    "solvent": "THF"
  },
  "addition_sequence": [
    {"role": "solvent", "timing": "first"},
    {"role": "base", "timing": "after_solvent"},
    {"role": "catalyst", "timing": "after_base"},
    {"role": "starting_material", "timing": "after_catalyst"}
  ]
}
```

**Disadvantages:**
- ❌ Requires schema migration
- ❌ Need to update all 9 rule files manually
- ❌ More verbose (duplicates role information)
- ❌ Still need to handle ranges and choices

---

### Option 3: Hybrid Approach

**Concept:** 
1. Use **Option 1** as default (automatic generation)
2. Allow rules to optionally specify `addition_sequence_override` for special cases
3. Best of both worlds: simple rules + flexibility when needed

```json
{
  "conditions": {...},
  "addition_sequence_override": [
    {"component": "base", "order": 1, "note": "Must be added before catalyst"},
    {"component": "catalyst", "order": 2, "note": "Add after base activation"}
  ]
}
```

---

## Recommended Implementation Plan

### Phase 1: Add Output Formatter (1-2 hours)

**File:** `chemtools/formatters/addition_sequence.py` (NEW)

```python
"""
Generate standardized addition sequences for automated execution.

Converts rule-based conditions into protocol-like reaction_setup format
with explicit addition ordering suitable for robotic systems.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


@dataclass
class ChemicalAddition:
    """Represents a single chemical addition step."""
    name: str
    role: str
    amount: Dict[str, Any]
    order: int
    timing: Optional[str] = None
    note: Optional[str] = None


# Standard addition order based on chemistry best practices
ROLE_ORDER = {
    "solvent": 1,
    "base": 2,
    "ligand": 3,
    "metal_catalyst": 4,
    "catalyst": 5,
    "additive": 6,
    "starting_material": 7,
    "reagent": 8,
}


def generate_addition_sequence(
    conditions: Dict[str, Any],
    substrate_smiles: Optional[List[str]] = None,
    scale_mmol: float = 1.0
) -> Dict[str, Any]:
    """
    Generate protocol-like reaction_setup from rule conditions.
    
    Args:
        conditions: Rule conditions dict (from rule JSON)
        substrate_smiles: Optional list of substrate SMILES for auto-detection
        scale_mmol: Reaction scale in mmol (default 1.0)
    
    Returns:
        Dict with reaction_setup structure suitable for automation
    """
    chemicals = []
    
    # Map conditions to chemicals
    if "solvent" in conditions:
        chemicals.append(
            _make_chemical(
                name=_pick_first_option(conditions["solvent"]),
                role="solvent",
                amount={"volume_ml": "to_be_scaled", "note": "Typical: 0.1-0.5 M concentration"},
                order=ROLE_ORDER["solvent"]
            )
        )
    
    if "base" in conditions:
        base_equiv = _pick_midpoint(conditions.get("base_equiv", "2.0"))
        chemicals.append(
            _make_chemical(
                name=_pick_first_option(conditions["base"]),
                role="base",
                amount={"equivalents": base_equiv},
                order=ROLE_ORDER["base"]
            )
        )
    
    if "ligand" in conditions and "Built-in" not in conditions["ligand"]:
        ligand_loading = _pick_midpoint(conditions.get("ligand_loading_molpct", "2.0"))
        chemicals.append(
            _make_chemical(
                name=_pick_first_option(conditions["ligand"]),
                role="ligand",
                amount={"mol_percent": ligand_loading},
                order=ROLE_ORDER["ligand"]
            )
        )
    
    if "catalyst" in conditions:
        cat_loading = _pick_midpoint(conditions.get("catalyst_loading_molpct", "1.0"))
        chemicals.append(
            _make_chemical(
                name=_pick_first_option(conditions["catalyst"]),
                role="metal_catalyst" if "Pd" in conditions["catalyst"] or "Ni" in conditions["catalyst"] else "catalyst",
                amount={"mol_percent": cat_loading},
                order=ROLE_ORDER.get("metal_catalyst", ROLE_ORDER["catalyst"])
            )
        )
    
    if "additives" in conditions and conditions["additives"]:
        for idx, additive in enumerate(conditions["additives"]):
            chemicals.append(
                _make_chemical(
                    name=additive,
                    role="additive",
                    amount={"note": "As specified in conditions"},
                    order=ROLE_ORDER["additive"] + idx * 0.1
                )
            )
    
    # Add placeholder for substrates
    chemicals.append(
        _make_chemical(
            name="Substrate (user-provided)",
            role="starting_material",
            amount={"equivalents": 1.0, "mmol": scale_mmol},
            order=ROLE_ORDER["starting_material"],
            note="Replace with actual substrate"
        )
    )
    
    # Sort by order
    chemicals.sort(key=lambda x: x["addition_order"])
    
    # Format conditions
    reaction_conditions = {
        "temperature_C": _pick_midpoint(conditions.get("temperature_C", "25")),
        "time_h": _pick_midpoint(conditions.get("time_h", "4")),
        "atmosphere": conditions.get("atmosphere", "N2 or Ar"),
    }
    
    return {
        "reaction_setup": [
            {
                "chemicals": chemicals,
                "conditions": [reaction_conditions]
            }
        ],
        "metadata": {
            "generated_from": "rule_conditions",
            "scale_mmol": scale_mmol,
            "note": "Addition sequence generated using standard heuristics. Verify before execution."
        }
    }


def _make_chemical(name: str, role: str, amount: Dict[str, Any], order: float, note: Optional[str] = None) -> Dict[str, Any]:
    """Helper to create chemical entry."""
    chem = {
        "name": name,
        "role": role,
        "amount": amount,
        "addition_order": order
    }
    if note:
        chem["note"] = note
    return chem


def _pick_first_option(value: str) -> str:
    """Pick first option from 'A or B or C' string."""
    if " or " in value:
        return value.split(" or ")[0].strip()
    return value.strip()


def _pick_midpoint(value: str | float) -> float:
    """Pick midpoint from range like '0.5-2.0' or '40-80'."""
    if isinstance(value, (int, float)):
        return float(value)
    
    value_str = str(value).strip()
    
    if "-" in value_str:
        try:
            parts = value_str.split("-")
            low = float(parts[0].strip())
            high = float(parts[1].strip())
            return (low + high) / 2.0
        except (ValueError, IndexError):
            return 1.0  # Fallback
    
    try:
        return float(value_str)
    except ValueError:
        return 1.0  # Fallback
```

### Phase 2: Integrate into UnifiedRecommender (30 minutes)

Update `chemtools/recommend/unified.py`:

```python
from chemtools.formatters.addition_sequence import generate_addition_sequence

class UnifiedRecommender:
    def recommend(
        self,
        reaction_smiles: str,
        top_k: int = 10,
        min_similarity: float = 0.0,
        source_types: Optional[List[str]] = None,
        validate_rules: bool = True,
        format_for_automation: bool = False,  # NEW parameter
    ) -> List[Dict[str, Any]]:
        """
        ...
        Args:
            format_for_automation: If True, convert rule conditions to
                protocol-like reaction_setup format for robotic execution
        """
        results = self._search(...)
        
        # Existing formatting...
        
        # NEW: Add automation-friendly format
        if format_for_automation:
            for result in results:
                if result["source_type"] == "rule":
                    # Get conditions from rule
                    rule_data = self._load_source(result["source_file"])
                    conditions = rule_data.get("default_rule", {}).get("conditions", {})
                    
                    # Generate addition sequence
                    result["reaction_setup"] = generate_addition_sequence(
                        conditions=conditions,
                        substrate_smiles=self._extract_substrates(reaction_smiles)
                    )
                    result["formatted_for_automation"] = True
        
        return results
```

### Phase 3: Update LangChain Tool (15 minutes)

Update `chem_assistant/chemtools_wrapper.py`:

```python
class UnifiedRecommenderInput(BaseModel):
    # ... existing fields ...
    format_for_automation: bool = Field(
        default=False,
        description="Generate protocol-like addition sequences for robotic execution"
    )

@tool("unified_recommender", args_schema=UnifiedRecommenderInput, return_direct=False)
def unified_recommender_tool(..., format_for_automation: bool = False) -> Dict[str, Any]:
    """
    ... (updated docstring to mention automation format)
    """
    results = recommender.recommend(
        reaction_smiles=reaction_smiles,
        top_k=top_k,
        min_similarity=min_similarity,
        source_types=source_types,
        validate_rules=validate_rules,
        format_for_automation=format_for_automation,  # NEW
    )
```

### Phase 4: Testing & Validation (1 hour)

Create `tests/test_addition_sequence.py`:

```python
def test_generate_addition_sequence_from_sonogashira_rule():
    """Test addition sequence generation for Sonogashira rule."""
    conditions = {
        "catalyst": "PdCl2(PPh3)2 or Pd(dppf)Cl2·DCM",
        "catalyst_loading_molpct": "0.5-2.0",
        "base": "Et3N or DIPEA",
        "base_equiv": "2.0-3.0",
        "solvent": "THF, toluene, or DMF",
        "temperature_C": "40-80",
        "time_h": "1-8",
        "additives": ["CuI (2-5 mol%)"]
    }
    
    result = generate_addition_sequence(conditions)
    
    # Check structure
    assert "reaction_setup" in result
    assert len(result["reaction_setup"]) > 0
    
    chemicals = result["reaction_setup"][0]["chemicals"]
    
    # Check addition order
    orders = [c["addition_order"] for c in chemicals]
    assert orders == sorted(orders), "Chemicals should be sorted by addition order"
    
    # Check first is solvent
    assert chemicals[0]["role"] == "solvent"
    
    # Check midpoints picked
    base_chem = next(c for c in chemicals if c["role"] == "base")
    assert base_chem["amount"]["equivalents"] == 2.5  # midpoint of 2.0-3.0
```

## Benefits

1. **Automation-Ready**: Output can be directly fed to robotic systems
2. **Backward Compatible**: No changes to existing rule files needed
3. **Flexible**: Can override with custom sequences if needed
4. **Validated**: Follows chemistry best practices for addition order
5. **Standardized**: Both protocols and rules output the same format
6. **Extensible**: Easy to add more sophistication later

## Timeline

- **Phase 1:** 2 hours (create formatter module)
- **Phase 2:** 30 minutes (integrate into UnifiedRecommender)
- **Phase 3:** 15 minutes (update LangChain tool)
- **Phase 4:** 1 hour (testing)

**Total:** ~4 hours of development time

## Next Steps

1. Get approval on this approach
2. Create `chemtools/formatters/addition_sequence.py`
3. Add `format_for_automation` parameter to UnifiedRecommender
4. Update chem_assistant tool wrapper
5. Write comprehensive tests
6. Document in AGENTS.md and README.md
