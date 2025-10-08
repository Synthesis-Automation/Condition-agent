# Reaction Type Parameters - Quick Answer

## ✅ YES! Both Systems Accept Reaction Type Parameters

### Quick Overview

Both ML and rule-based recommendation functions accept a `reaction_type` parameter in addition to SMILES:

```python
# ML Recommendations
get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>product",
    reaction_type="C-N Coupling (Pd)",  # ← Reaction type parameter
    top_k=3
)

# Rule-Based Recommendations
get_rule_recommendations(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>product",
    reaction_type="C-N Coupling (Pd)"  # ← Reaction type parameter
)
```

## Two Modes of Operation

### 1. **Auto-Detection Mode**

```python
# Let system detect reaction type automatically
get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>product",
    reaction_type="Auto-detect",  # ← System will detect as Suzuki
    top_k=3
)
```

**How it works:**
- Uses rxn-insight ML classifier (primary)
- Falls back to rule-based functional group detection
- Provides confidence scores
- Maps detected type to both ML and rule systems

### 2. **Manual Specification Mode**

```python
# Explicitly specify reaction type
get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>product",
    reaction_type="C-N Coupling (Pd)",  # ← Manually specified
    top_k=3
)
```

**When to use:**
- You know the exact reaction type
- Auto-detection confidence is low
- Want reproducible results
- Testing specific models

## Supported Reaction Types

### UI Names (User Input)
```
"Auto-detect"
"C-N Coupling (Pd)"
"C-N Coupling (Cu)"
"C-N Coupling (Ni)"
"Amide Formation"
"Suzuki Coupling"
```

### ML API Names (Internal)
```python
None  # Auto-detect
"C_N_Coupling_Pd"
"C_N_Coupling_Cu"
"C_N_Coupling_Ni"
"Amide_Coupling"
"Suzuki_CC"
```

## Full Function Signatures

### ML Recommendations

```python
def get_ml_recommendations(
    reaction_smiles: str,
    reaction_type: str,      # ← Accepts UI names or "Auto-detect"
    top_k: int = 3,
) -> Tuple[str, List[List[Any]]]:
    """Get ML-based recommendations with enhanced formatting."""
```

### Rule-Based Recommendations

```python
def get_rule_recommendations(
    reaction_smiles: str,
    reaction_type: str,      # ← Accepts UI names or "Auto-detect"
) -> Tuple[str, List[List[Any]]]:
    """Get rule-based (SchemeConditionDB) recommendations."""
```

### Core ML API

```python
def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,  # ← Accepts ML family names or None
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Structured condition recommendations for API/UI consumers."""
```

## Parameter Flow

```
User Input: reaction_type="Auto-detect" or "C-N Coupling (Pd)"
    ↓
detect_and_map_reaction_type()
    ↓
Maps to appropriate naming convention:
    ├─ ML family: "C_N_Coupling_Pd"
    └─ Rule DB: "C-N Coupling (Pd)"
    ↓
Passed to respective systems
```

## Usage Examples

### Example 1: Auto-Detect (Recommended)

```python
from app.ui_simple import get_ml_recommendations

summary, table = get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    reaction_type="Auto-detect",  # ← System detects Suzuki
    top_k=3
)

# Output shows detection:
# **Auto-detected (rxn-insight):** Suzuki_CC
#   Class: C-C Coupling
#   Name: Suzuki coupling with boronic acids
#   Confidence: 80.00%
```

### Example 2: Manual Specification

```python
from app.ui_simple import get_rule_recommendations

summary, table = get_rule_recommendations(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_type="C-N Coupling (Pd)",  # ← Manually specified
)

# Uses specified type, skips detection
```

### Example 3: Direct API Call

```python
from chemtools import recommend

data = recommend.recommend_conditions_structured(
    reaction="Brc1ccccc1.Nc1ccccc1>>product",
    reaction_type="C_N_Coupling_Pd",  # ← ML family name
    k=10,
    limit=3
)

recommendations = data["recommendations"]
```

### Example 4: Check Confidence First

```python
from app.ui_simple import detect_and_map_reaction_type

# Check detection confidence
result = detect_and_map_reaction_type(
    "Brc1ccccc1.Nc1ccccc1>>product",
    "Auto-detect"
)

if result['confidence'] > 0.6:
    # High confidence, proceed with auto-detect
    summary, table = get_ml_recommendations(
        "Brc1ccccc1.Nc1ccccc1>>product",
        "Auto-detect"
    )
else:
    # Low confidence, use manual selection
    print(f"⚠ Low confidence: {result['confidence']:.1%}")
    summary, table = get_ml_recommendations(
        "Brc1ccccc1.Nc1ccccc1>>product",
        "C-N Coupling (Pd)"  # Manual fallback
    )
```

## Key Features

✅ **Flexible Input**
- Auto-detection OR manual specification
- Supports both modes seamlessly

✅ **Intelligent Mapping**
- Automatically converts between naming conventions
- Maps to both ML and rule-based systems

✅ **Catalyst-Aware**
- Detects Pd vs Cu vs Ni catalysts
- Extracts from SMILES agents section

✅ **Confidence Scoring**
- Provides detection confidence (0-100%)
- Warns when confidence < 50%

✅ **Error Handling**
- Helpful messages for unsupported types
- Always suggests alternatives
- Graceful fallbacks

## When to Use Each Mode

| Scenario | Recommendation |
|----------|---------------|
| **Standard reactions** | Auto-detect |
| **Know exact type** | Manual specification |
| **Low confidence (<60%)** | Manual specification |
| **Testing/development** | Manual specification |
| **Production/user-facing** | Auto-detect with confidence check |
| **Unusual reactions** | Manual specification |

## Complete Documentation

For detailed information, see:
- **`docs/REACTION_TYPE_PARAMETER_GUIDE.md`** - Complete guide with all details
- **`docs/AUTO_DETECTION_GUIDE.md`** - Auto-detection system documentation
- **`docs/ML_ERROR_HANDLING.md`** - Error handling when type detection fails

## Summary

**Yes, both ML and rule-based systems accept reaction type parameters!**

The system provides:
- ✅ **Two modes**: Auto-detect or manual specification
- ✅ **Flexible API**: Accepts both UI names and internal names
- ✅ **Intelligent mapping**: Converts between naming conventions
- ✅ **High accuracy**: ML-based detection with rule-based fallback
- ✅ **User-friendly**: Clear error messages and alternatives

You can use either mode depending on your needs, and the system handles all the complexity of mapping between different naming conventions and routing to the appropriate recommendation engine.
