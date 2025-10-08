# Reaction Type Parameter Flow - Complete Guide

## Overview

Yes! Both ML and rule-based recommendation systems accept reaction type parameters in addition to SMILES. The system supports:

1. **Manual specification** - User explicitly selects reaction type
2. **Auto-detection** - System automatically detects from SMILES
3. **Hybrid approach** - Auto-detect with confidence thresholds and fallback

## API Signatures

### 1. UI-Level Functions (app/ui_simple.py)

#### ML Recommendations

```python
def get_ml_recommendations(
    reaction_smiles: str,
    reaction_type: str,
    top_k: int = 3,
) -> Tuple[str, List[List[Any]]]:
    """
    Get ML-based recommendations with enhanced formatting.
    
    Args:
        reaction_smiles: Reaction SMILES string (e.g., "reactants>>products")
        reaction_type: Reaction type name OR "Auto-detect"
                      Valid types: "C-N Coupling (Pd)", "Suzuki Coupling", etc.
        top_k: Number of recommendations to return (default: 3)
    
    Returns:
        (summary_markdown, table_data)
    """
```

#### Rule-Based Recommendations

```python
def get_rule_recommendations(
    reaction_smiles: str,
    reaction_type: str,
) -> Tuple[str, List[List[Any]]]:
    """
    Get rule-based (SchemeConditionDB) recommendations.
    
    Args:
        reaction_smiles: Reaction SMILES string
        reaction_type: Reaction type name OR "Auto-detect"
                      Valid types: "C-N Coupling (Pd)", "Suzuki Coupling", etc.
    
    Returns:
        (summary_markdown, table_data)
    """
```

### 2. Core ML API (chemtools/recommend.py)

```python
def recommend_conditions_structured(
    reaction: str,
    reaction_type: Optional[str] = None,
    *,
    k: int = 50,
    limit: int = 5,
    relax: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Structured condition recommendations for API/UI consumers.
    
    Args:
        reaction: Reaction SMILES string
        reaction_type: ML family name (e.g., "C_N_Coupling_Pd", "Suzuki_CC")
                      If None, will auto-detect using rxn-insight or rules
        k: Number of nearest neighbors to consider (default: 50)
        limit: Maximum number of recommendations (default: 5)
        relax: Settings for relaxation/optimization
        constraints: Constraint rules for filtering
    
    Returns:
        Dict with keys: meta, input, detection, recommendations, alternatives, 
                       precedents, filters, role_featurization, rule_features
    """
```

### 3. Rule-Based API (scdb_matcher/matcher.py)

```python
def match(
    db: RuleDB,
    rxn_smiles: str | None = None,
    *,
    features: Mapping[str, Any] | None = None,
) -> MatchResult:
    """
    Match a reaction against a rule database.
    
    Args:
        db: Loaded rule database (scheme or selector)
        rxn_smiles: Reaction SMILES required for SMARTS-driven scheme databases
        features: Nested mapping of features required for selector databases
    
    Note: The reaction type is implicit in which database you load.
          Each database file represents a specific reaction type.
    """
```

## Parameter Flow Diagram

```
User Input
    ├─ reaction_smiles: "Brc1ccccc1.Nc1ccccc1>>product"
    └─ reaction_type: "Auto-detect" OR "C-N Coupling (Pd)"
         ↓
detect_and_map_reaction_type()
    ├─ If "Auto-detect":
    │   ├─ Try rxn-insight ML detection
    │   └─ Fallback to router rule-based detection
    └─ If manual: Use as-is
         ↓
    Returns:
        ├─ detected_family: "Buchwald_CN"
        ├─ ml_family: "C_N_Coupling_Pd"
        └─ rule_db_name: "C-N Coupling (Pd)"
         ↓
         ├─────────────────┬─────────────────┐
         ↓                 ↓                 ↓
    ML Path          Rule Path        Display
         ↓                 ↓                 ↓
recommend_conditions_    scdb_load_db()   Show detection
structured(              + scdb_match()    info to user
  reaction=smiles,
  reaction_type="C_N_Coupling_Pd"  # ML family name
)
```

## Usage Examples

### Example 1: Auto-Detection (Recommended)

```python
from app.ui_simple import get_ml_recommendations

# Let the system detect reaction type
summary, table = get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    reaction_type="Auto-detect",  # ← System will detect Suzuki
    top_k=3
)
```

**What happens:**
1. `detect_and_map_reaction_type()` analyzes SMILES
2. rxn-insight detects: "C-C Coupling" → "Suzuki coupling"
3. Maps to ML family: `Suzuki_CC`
4. Passes `reaction_type="Suzuki_CC"` to ML system
5. Returns recommendations

### Example 2: Manual Specification

```python
from app.ui_simple import get_ml_recommendations

# Explicitly specify reaction type
summary, table = get_ml_recommendations(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_type="C-N Coupling (Pd)",  # ← Manual specification
    top_k=3
)
```

**What happens:**
1. `detect_and_map_reaction_type()` sees manual type
2. Skips auto-detection
3. Maps "C-N Coupling (Pd)" → ML family `C_N_Coupling_Pd`
4. Passes to ML system directly
5. Returns recommendations

### Example 3: Rule-Based with Auto-Detection

```python
from app.ui_simple import get_rule_recommendations

# Auto-detect for rule-based system
summary, table = get_rule_recommendations(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    reaction_type="Auto-detect"  # ← System will detect and load DB
)
```

**What happens:**
1. `detect_and_map_reaction_type()` analyzes SMILES
2. Maps to rule DB name: "Suzuki Coupling"
3. Loads: `data/conditionDB/suzuki_db.json`
4. Matches reaction against rules
5. Returns recommendations

### Example 4: Direct API Call (Advanced)

```python
from chemtools import recommend

# Direct ML API call with explicit reaction type
data = recommend.recommend_conditions_structured(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    reaction_type="C_N_Coupling_Pd",  # ← ML family name
    k=10,
    limit=3,
    relax={
        "use_drfp": False,
        "use_rxn_insight": True,
    }
)

# Returns full structured data
recommendations = data["recommendations"]
detection = data["detection"]
precedents = data["precedents"]
```

## Reaction Type Naming Conventions

The system uses **three different naming conventions** depending on the layer:

### UI Layer (User-Facing)

Used in Gradio dropdown and user messages:

```python
"C-N Coupling (Pd)"
"C-N Coupling (Cu)"
"C-N Coupling (Ni)"
"Amide Formation"
"Suzuki Coupling"
"Auto-detect"
```

### ML Layer (Internal API)

Used in ML recommendation system:

```python
"C_N_Coupling_Pd"
"C_N_Coupling_Cu"
"C_N_Coupling_Ni"
"Amide_Coupling"
"Suzuki_CC"
None  # For auto-detect
```

### Rule-Based Layer (Database Files)

Used for rule database file selection:

```python
"C-N Coupling (Pd)"  → data/conditionDB/cn_coupling_pd_db.json
"C-N Coupling (Cu)"  → data/conditionDB/cn_coupling_cu_db.json
"Suzuki Coupling"    → data/conditionDB/suzuki_db.json
```

## Mapping Configuration

### ML Family Mapping (app/ui_simple.py)

```python
ML_FAMILY_MAP = {
    "Auto-detect": None,
    "C-N Coupling (Cu)": "C_N_Coupling_Cu",
    "C-N Coupling (Pd)": "C_N_Coupling_Pd",
    "C-N Coupling (Ni)": "C_N_Coupling_Ni",
    "Amide Formation": "Amide_Coupling",
    "Suzuki Coupling": "Suzuki_CC",
}
```

### Rule Database Mapping (app/ui_simple.py)

```python
RULE_DATABASES = {
    "C-N Coupling (Cu)": str(SCDB_DIR / "cn_coupling_cu_db.json"),
    "C-N Coupling (Pd)": str(SCDB_DIR / "cn_coupling_pd_db.json"),
    "C-N Coupling (Ni)": str(SCDB_DIR / "cn_coupling_ni.json"),
    "Amide Formation": str(SCDB_DIR / "amide_formation_db.json"),
    "Suzuki Coupling": str(SCDB_DIR / "suzuki_db.json"),
}
```

## Auto-Detection Flow

### Detection Function

```python
def detect_and_map_reaction_type(
    reaction_smiles: str, 
    requested_type: str
) -> Dict[str, Any]:
    """
    Detect reaction type and map to both ML and rule-based conventions.
    
    Returns:
        {
            'detected_family': str,      # Internal name (e.g., "Buchwald_CN")
            'ml_family': str | None,     # ML API name (e.g., "C_N_Coupling_Pd")
            'rule_db_name': str | None,  # Rule DB name (e.g., "C-N Coupling (Pd)")
            'confidence': float,         # Detection confidence 0-1
            'method': str,               # 'rxn_insight' or 'router'
            'success': bool,             # Whether detection succeeded
            'message': str,              # Human-readable message
        }
    """
```

### Detection Priority

1. **Manual Type** (if not "Auto-detect")
   - Return user's selection directly
   - Skip auto-detection

2. **rxn-insight ML Detector** (if available)
   - Use ML-based classification
   - Higher accuracy (~80% confidence)
   - Extracts catalyst information

3. **Router Rule-Based** (fallback)
   - Use functional group pattern matching
   - Always available
   - Lower confidence

## Integration Points

### From Gradio UI

```python
# In app/ui_simple.py Gradio interface setup
with gr.Row():
    reaction_type = gr.Dropdown(
        choices=[
            "Auto-detect",
            "C-N Coupling (Pd)",
            "C-N Coupling (Cu)",
            "Suzuki Coupling",
            # ...
        ],
        value="Auto-detect",
        label="Reaction Type"
    )

ml_button.click(
    fn=get_ml_recommendations,
    inputs=[reaction_input, reaction_type, top_k_slider],
    outputs=[ml_output, ml_table]
)
```

### From Python Script

```python
# Direct call example
from app.ui_simple import get_ml_recommendations

# With auto-detection
ml_summary, ml_table = get_ml_recommendations(
    "Brc1ccccc1.Nc1ccccc1>>product",
    "Auto-detect",
    top_k=3
)

# With manual type
ml_summary, ml_table = get_ml_recommendations(
    "Brc1ccccc1.Nc1ccccc1>>product",
    "C-N Coupling (Pd)",
    top_k=3
)
```

### From API (FastAPI endpoint)

```python
# In app/main.py (example)
from fastapi import FastAPI
from chemtools import recommend

app = FastAPI()

@app.post("/recommend/ml")
def recommend_ml(
    reaction_smiles: str,
    reaction_type: str | None = None,  # ← Optional parameter
    top_k: int = 3
):
    data = recommend.recommend_conditions_structured(
        reaction=reaction_smiles,
        reaction_type=reaction_type,  # Pass through
        k=10,
        limit=top_k
    )
    return data
```

## Error Handling

The system handles reaction type errors gracefully:

### Unsupported Type

```python
# If user selects/system detects unsupported type
# Returns helpful error:
"""
**Cannot Proceed with ML Recommendations**

❌ **No ML model available** for: `Esterification`

**Available ML reaction types:**
  - C-N Coupling (Pd) (`C_N_Coupling_Pd`)
  - Suzuki Coupling (`Suzuki_CC`)
  ...

**What to do:**
1. ✅ Try rule-based recommendations instead
2. 🔄 Manually select a supported reaction type
"""
```

### Detection Failure

```python
# If auto-detection fails
# Returns error with guidance:
"""
**Auto-detection failed:** [Error details]

**What to do:**
1. Manually select reaction type from dropdown
2. Verify SMILES format is correct
"""
```

## Best Practices

### When to Use Auto-Detect

✅ **Good for:**
- Standard, well-defined reactions
- When catalyst is in SMILES (`>Pd>`)
- Exploring unknown reactions
- Quick prototyping

❌ **Not ideal for:**
- Unusual or hybrid reactions
- When you know the exact type
- Low-confidence scenarios (<50%)

### When to Specify Manually

✅ **Good for:**
- You know the exact reaction type
- Auto-detection confidence is low
- Testing specific models
- Reproducible results

### Hybrid Approach (Recommended)

```python
# Check detection confidence first
from app.ui_simple import detect_and_map_reaction_type

result = detect_and_map_reaction_type(smiles, "Auto-detect")

if result['confidence'] < 0.6:
    print(f"⚠ Low confidence ({result['confidence']:.1%})")
    print("Consider manual selection")
    # Prompt user or use fallback
else:
    # Proceed with auto-detected type
    ml_summary, _ = get_ml_recommendations(smiles, "Auto-detect")
```

## Summary

✅ **Yes, both systems accept reaction type parameters!**

| System | Parameter | Format | Auto-Detect |
|--------|-----------|--------|-------------|
| **ML** | `reaction_type` | UI names or ML families | ✅ Yes |
| **Rule-Based** | `reaction_type` | UI names | ✅ Yes |
| **Core ML API** | `reaction_type` | ML family names | ✅ Yes (None) |
| **Core Rule API** | Implicit in DB | Database file path | ❌ No |

**Key Features:**
- ✅ Flexible input (auto-detect or manual)
- ✅ Intelligent mapping across naming conventions
- ✅ Catalyst-aware detection
- ✅ Graceful error handling
- ✅ Confidence scoring
- ✅ Always suggest alternatives

The system provides maximum flexibility while maintaining a clean, user-friendly API!
