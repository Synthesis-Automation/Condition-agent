# LLM Workflow Implementation Guide

## Overview

This document describes the **pure LLM workflow** for reagent taxonomy generation, implemented in response to user feedback that the mixed deterministic+LLM approach produced "very messy output."

## Architecture

### New Modules

1. **`llmtools/reagent_classifier.py`** - Core LLM classification functions
2. **`llmtools/prompts.py`** - Three new prompt templates (added)
3. **`reagent_taxonomy_qt.py`** - New `generate_taxonomy_entry_llm()` function

### Workflow Steps

```
Input: CAS number
  ↓
┌─────────────────────────────────────────────────────┐
│ Step 1: Resolve Identity (PubChem API)             │
│ → name, SMILES, formula, synonyms                  │
└─────────────────────────────────────────────────────┘
  ↓
┌─────────────────────────────────────────────────────┐
│ Step 2: LLM Role Classification                    │
│ → Classify into 13 role types (base, ligand, etc.) │
│ → Prompt: REAGENT_ROLE_CLASSIFICATION              │
│ → Output: {role, confidence, reasoning}            │
└─────────────────────────────────────────────────────┘
  ↓
┌─────────────────────────────────────────────────────┐
│ Step 3: LLM Field Assignment                       │
│ → Pick family within role                          │
│ → Assign all required fields                       │
│ → Prompt: REAGENT_FIELD_ASSIGNMENT                 │
│ → Output: {family, fields, abbreviations, ...}     │
└─────────────────────────────────────────────────────┘
  ↓
┌─────────────────────────────────────────────────────┐
│ Step 4: LLM Verification                           │
│ → Check for chemical accuracy errors               │
│ → Validate field consistency                       │
│ → Prompt: REAGENT_ENTRY_VERIFICATION               │
│ → Output: {approved, issues, suggestions}          │
└─────────────────────────────────────────────────────┘
  ↓
Output: {status, workflow, entry}
```

## API Reference

### Function: `generate_taxonomy_entry_llm()`

**Location**: `data-processor/reagent_taxonomy_qt.py`

**Purpose**: Pure LLM workflow replacing mixed deterministic+LLM approach

**Parameters**:
```python
def generate_taxonomy_entry_llm(
    *,
    cas: str,                           # CAS registry number
    registry_dir: Path,                 # Path to reagent registry
    llm_client: Any,                    # LLMClient from llmtools.clients
    name_override: Optional[str] = None,
    smiles_override: Optional[str] = None,
    resolver_timeout: float = DEFAULT_RESOLVER_TIMEOUT,
) -> Dict[str, Any]:
```

**Returns**:
```json
{
  "status": "ready_to_save" | "needs_review" | "error",
  "workflow": {
    "step1_identity": {
      "status": "success",
      "identity": {
        "name": "Triethylamine",
        "cas": "121-44-8",
        "smiles": "CCN(CC)CC",
        "molecular_formula": "C6H15N",
        "synonyms": ["TEA", "N,N-Diethylethanamine", ...]
      }
    },
    "step2_role": {
      "status": "success",
      "role": "base",
      "confidence": 0.95,
      "reasoning": "Tertiary amine with basicity",
      "model": "deepseek-v3",
      "tokens": 250,
      "latency_ms": 1200
    },
    "step3_fields": {
      "status": "success",
      "family": "tertiary_amines_aliphatic",
      "fields": {
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      },
      "abbreviations": ["TEA"],
      "additional_synonyms": [],
      "confidence": 0.92,
      "reasoning": "Aliphatic tertiary amine with strong basicity",
      "model": "deepseek-v3",
      "tokens": 450,
      "latency_ms": 1500
    },
    "step4_verification": {
      "status": "success",
      "approved": true,
      "issues": [],
      "suggestions": ["Consider adding pKa value to metadata"],
      "model": "deepseek-v3",
      "tokens": 350,
      "latency_ms": 1100
    }
  },
  "entry": {
    "name": "Triethylamine",
    "cas": "121-44-8",
    "molecular_formula": "C6H15N",
    "smiles": "CCN(CC)CC",
    "roles": {
      "base": {
        "families": ["tertiary_amines_aliphatic"],
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    },
    "abbreviations": ["TEA"],
    "synonyms": ["TEA", "N,N-Diethylethanamine", ...]
  },
  "message": "Entry successfully generated and verified by LLM. Ready to save!"
}
```

**Status Values**:
- **`ready_to_save`**: All steps succeeded, LLM approved entry, no errors
- **`needs_review`**: Entry generated but has warnings/errors from verification
- **`error`**: Workflow failed at some step (see `error` field)

### Module: `llmtools/reagent_classifier.py`

Three core functions for LLM-based classification:

#### 1. `classify_role()`

Classify reagent into one of 13 roles.

```python
def classify_role(
    identity: Dict[str, Any],
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 300,
) -> Dict[str, Any]:
```

**Prompt**: `REAGENT_ROLE_CLASSIFICATION`

**Output**:
```json
{
  "status": "success",
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary amine with strong basicity",
  "model": "deepseek-v3",
  "tokens": 250,
  "latency_ms": 1200
}
```

**Valid Roles**:
- `ligand`
- `metal_precursor`
- `preformed_metal_catalyst`
- `base`
- `acid`
- `condensation_agent`
- `oxidant`
- `reductant`
- `additive`
- `solvent`
- `organo_catalyst`
- `enzyme`
- `other_reagent` (last resort only)

#### 2. `assign_fields()`

Assign family and populate role-specific fields.

```python
def assign_fields(
    identity: Dict[str, Any],
    role: str,
    registry_dir: Path,
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 600,
) -> Dict[str, Any]:
```

**Prompt**: `REAGENT_FIELD_ASSIGNMENT`

**Output**:
```json
{
  "status": "success",
  "family": "tertiary_amines_aliphatic",
  "fields": {
    "basicity": "strong",
    "nucleophilicity": "moderate",
    "sterics": "unhindered"
  },
  "abbreviations": ["TEA"],
  "additional_synonyms": [],
  "confidence": 0.92,
  "reasoning": "Aliphatic tertiary amine...",
  "model": "deepseek-v3",
  "tokens": 450,
  "latency_ms": 1500
}
```

**Role-Specific Fields**:

| Role | Required Fields |
|------|----------------|
| `base` | `basicity`, `nucleophilicity`, `sterics` |
| `metal_precursor` | `metal`, `oxidation_states` |
| `preformed_metal_catalyst` | `metal`, `oxidation_states`, `ligand_type` |
| `solvent` | `proticity`, `polarity`, `coordination` |
| `ligand` | `donors`, `denticity` |
| `oxidant` | `strength_band` |
| `reductant` | `strength_band` |
| `condensation_agent` | `strength_band` |
| `acid` | `acidity` |
| `organo_catalyst` | `activation_mode`, `chirality` |
| `enzyme` | `source`, `cofactor_requirement` |

#### 3. `verify_entry()`

Quality control check for final entry.

```python
def verify_entry(
    entry: Dict[str, Any],
    llm_client: LLMClient,
    temperature: float = 0.0,
    max_tokens: int = 500,
) -> Dict[str, Any]:
```

**Prompt**: `REAGENT_ENTRY_VERIFICATION`

**Output**:
```json
{
  "status": "success",
  "approved": false,
  "issues": [
    {
      "severity": "error",
      "field": "roles.base.basicity",
      "message": "Tertiary amines typically have moderate-to-strong basicity, not weak"
    },
    {
      "severity": "warning",
      "field": "abbreviations",
      "message": "Consider adding common abbreviation 'Et3N'"
    }
  ],
  "suggestions": [
    "Add pKa value to metadata",
    "Include CAS cross-references"
  ],
  "model": "deepseek-v3",
  "tokens": 350,
  "latency_ms": 1100
}
```

**Issue Severities**:
- **`error`**: Chemical accuracy issue, must fix before saving
- **`warning`**: Suggestion for improvement, non-blocking

## Prompt Templates

### 1. REAGENT_ROLE_CLASSIFICATION

**Purpose**: Classify reagent into 13 role types

**Input Variables**:
- `{name}` - Chemical name
- `{cas}` - CAS number
- `{smiles}` - SMILES string
- `{molecular_formula}` - Molecular formula
- `{synonyms}` - Comma-separated synonyms

**Output Schema**:
```json
{
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary amine with strong basicity"
}
```

**Key Instructions**:
- Use chemical knowledge (functional groups, structure, reactivity)
- "other_reagent" is LAST RESORT only
- Provide confidence score (0.0-1.0)
- Explain reasoning in 1-2 sentences

### 2. REAGENT_FIELD_ASSIGNMENT

**Purpose**: Assign family and populate role-specific fields

**Input Variables**:
- `{name}`, `{cas}`, `{smiles}`, `{molecular_formula}` - Identity
- `{role}` - Reagent role (from step 2)
- `{families_description}` - Available families for this role
- `{fields_schema}` - Required fields and allowed values
- `{examples}` - Example entries from registry

**Output Schema**:
```json
{
  "family": "tertiary_amines_aliphatic",
  "fields": {
    "basicity": "strong",
    "nucleophilicity": "moderate",
    "sterics": "unhindered"
  },
  "abbreviations": ["TEA"],
  "additional_synonyms": [],
  "confidence": 0.92,
  "reasoning": "Aliphatic tertiary amine..."
}
```

**Key Instructions**:
- Pick ONE family from available list
- ALL required fields must be assigned
- Use EXACT allowed values (don't invent new values)
- Suggest common abbreviations
- High confidence if structure clearly matches family

### 3. REAGENT_ENTRY_VERIFICATION

**Purpose**: Quality control check for errors

**Input Variables**:
- `{entry_json}` - Complete entry JSON

**Output Schema**:
```json
{
  "approved": false,
  "issues": [
    {
      "severity": "error",
      "field": "roles.base.basicity",
      "message": "..."
    }
  ],
  "suggestions": ["...", "..."]
}
```

**Key Checks**:
- Chemical accuracy (basicity of amines, oxidation states of metals)
- Field value consistency (no contradictions)
- Required field completeness
- SMILES validity (if possible)

**Approval Logic**:
- `approved = false` if ANY issue has `severity = "error"`
- `approved = true` if only warnings or no issues

## Usage Examples

### Example 1: Triethylamine (Base)

**Input**:
```python
from pathlib import Path
from llmtools.clients import LLMClient
from reagent_taxonomy_qt import generate_taxonomy_entry_llm

client = LLMClient(provider="aliyun", model="deepseek-v3")
registry_dir = Path("data/reagents")

result = generate_taxonomy_entry_llm(
    cas="121-44-8",
    registry_dir=registry_dir,
    llm_client=client,
)
```

**Output**:
```json
{
  "status": "ready_to_save",
  "workflow": {
    "step1_identity": {
      "status": "success",
      "identity": {
        "name": "Triethylamine",
        "cas": "121-44-8",
        "smiles": "CCN(CC)CC",
        "molecular_formula": "C6H15N"
      }
    },
    "step2_role": {
      "status": "success",
      "role": "base",
      "confidence": 0.95
    },
    "step3_fields": {
      "status": "success",
      "family": "tertiary_amines_aliphatic",
      "fields": {
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    },
    "step4_verification": {
      "status": "success",
      "approved": true,
      "issues": []
    }
  },
  "entry": {
    "name": "Triethylamine",
    "cas": "121-44-8",
    "molecular_formula": "C6H15N",
    "smiles": "CCN(CC)CC",
    "roles": {
      "base": {
        "families": ["tertiary_amines_aliphatic"],
        "basicity": "strong",
        "nucleophilicity": "moderate",
        "sterics": "unhindered"
      }
    },
    "abbreviations": ["TEA"]
  },
  "message": "Entry successfully generated and verified by LLM. Ready to save!"
}
```

### Example 2: Pd(OAc)₂ (Metal Precursor)

**Input**:
```python
result = generate_taxonomy_entry_llm(
    cas="14024-61-4",
    registry_dir=registry_dir,
    llm_client=client,
)
```

**Output**:
```json
{
  "status": "ready_to_save",
  "workflow": {
    "step2_role": {
      "role": "metal_precursor",
      "confidence": 0.98
    },
    "step3_fields": {
      "family": "palladium_ii_carboxylates",
      "fields": {
        "metal": "Pd",
        "oxidation_states": [2]
      },
      "abbreviations": ["Pd(OAc)2"]
    }
  },
  "entry": {
    "name": "Palladium(II) acetate",
    "cas": "14024-61-4",
    "roles": {
      "metal_precursor": {
        "families": ["palladium_ii_carboxylates"],
        "metal": "Pd",
        "oxidation_states": [2]
      }
    }
  }
}
```

### Example 3: Error Handling

**Scenario**: Invalid CAS number

**Input**:
```python
result = generate_taxonomy_entry_llm(
    cas="999-99-9",  # Non-existent
    registry_dir=registry_dir,
    llm_client=client,
)
```

**Output**:
```json
{
  "status": "error",
  "error": "Failed to resolve CAS 999-99-9. No data from PubChem."
}
```

**Scenario**: Verification finds errors

**Input**:
```python
# Assume LLM assigns wrong basicity
result = generate_taxonomy_entry_llm(...)
```

**Output**:
```json
{
  "status": "needs_review",
  "workflow": {
    "step4_verification": {
      "approved": false,
      "issues": [
        {
          "severity": "error",
          "field": "roles.base.basicity",
          "message": "Tertiary amines should be 'moderate' or 'strong', not 'weak'"
        }
      ]
    }
  },
  "entry": {...},
  "message": "Entry has 1 error(s) and 0 warning(s). Please review before saving."
}
```

## Comparison: Old vs New

### Old Workflow (Mixed)

**Problems**:
- ❌ 7+ JSON keys (messy output)
- ❌ Unclear which entry is final (`entry_original` vs `llm_adjusted_entry`)
- ❌ No auto-detect role option
- ❌ Mixed deterministic token matching + LLM review
- ❌ Field suggestions not automatically applied

**Output Structure**:
```json
{
  "entry_original": {...},
  "llm_adjusted_entry": {...},
  "review_summary": {...},
  "changes_applied": [...],
  "field_suggestions": {...},
  "needs_attention": [...],
  "success": true
}
```

### New Workflow (Pure LLM)

**Benefits**:
- ✅ Clean 3-key output (`status`, `workflow`, `entry`)
- ✅ Clear final entry in `entry` key
- ✅ Auto-detect role via LLM
- ✅ Pure LLM pipeline (no deterministic matching)
- ✅ All fields assigned by LLM
- ✅ Built-in quality control (verification step)

**Output Structure**:
```json
{
  "status": "ready_to_save",
  "workflow": {
    "step1_identity": {...},
    "step2_role": {...},
    "step3_fields": {...},
    "step4_verification": {...}
  },
  "entry": {...}
}
```

## Integration with UI

### Planned UI Updates

1. **Mode Toggle**: "LLM-powered" vs "Deterministic (legacy)"
2. **Auto-detect Role**: New option "Auto-detect (LLM)" in role dropdown
3. **Clean Output Display**:
   - Show workflow steps with checkmarks
   - Display confidence scores
   - Highlight issues/warnings from verification
4. **Save Button Logic**:
   - Enable if `status == "ready_to_save"`
   - Show warning if `status == "needs_review"`
   - Disable if `status == "error"`

### Example UI Display

```
┌─────────────────────────────────────────────────┐
│ Workflow Progress                               │
├─────────────────────────────────────────────────┤
│ ✓ Step 1: Identity Resolved                    │
│   Triethylamine (CAS 121-44-8)                  │
│                                                 │
│ ✓ Step 2: Role Classification                  │
│   Role: base (confidence: 95%)                  │
│   Reasoning: "Tertiary amine with basicity"     │
│                                                 │
│ ✓ Step 3: Field Assignment                     │
│   Family: tertiary_amines_aliphatic             │
│   Fields: basicity=strong, nucleophilicity=...  │
│   (confidence: 92%)                             │
│                                                 │
│ ✓ Step 4: Verification                         │
│   Status: APPROVED ✓                            │
│   0 errors, 0 warnings                          │
│                                                 │
│ [Save Entry] [Review JSON] [Cancel]             │
└─────────────────────────────────────────────────┘
```

## Testing Plan

### Test Cases

1. **Common Bases**: Triethylamine, DIPEA, DBU, Cs₂CO₃
2. **Metal Precursors**: Pd(OAc)₂, Ni(acac)₂, CuI, FeCl₃
3. **Solvents**: DMF, THF, Toluene, Water
4. **Ligands**: PPh₃, dppf, BINAP, phenanthroline
5. **Edge Cases**: Multi-role reagents, unknown structures

### Validation

- ✅ All steps return `status = "success"`
- ✅ Entry matches registry schema
- ✅ Fields have valid values (no invented strings)
- ✅ Verification catches obvious errors
- ✅ Confidence scores are reasonable (>0.8 for common reagents)

## Next Steps

1. **UI Integration** (next):
   - Add mode toggle
   - Add "Auto-detect (LLM)" role option
   - Update output display for clean format
   
2. **Testing**:
   - Test with 10+ examples
   - Compare LLM workflow vs deterministic workflow
   - Validate field accuracy

3. **Documentation**:
   - Update main README with LLM workflow guide
   - Add troubleshooting section

## Files Modified

1. ✅ `llmtools/prompts.py` - Added 3 new templates
2. ✅ `llmtools/reagent_classifier.py` - Created (new file)
3. ✅ `reagent_taxonomy_qt.py` - Added `generate_taxonomy_entry_llm()`
4. ⏳ `reagent_taxonomy_qt.py` - UI updates (pending)

## Status

**Phase 1: Infrastructure** ✅ COMPLETE
- Prompt templates created
- Classifier module implemented
- Main workflow function added

**Phase 2: UI Integration** ⏳ NEXT
- Add mode toggle
- Update output display
- Connect new function to UI

**Phase 3: Testing** ⏳ PENDING
- Test with examples
- Validate accuracy
- Performance benchmarking
