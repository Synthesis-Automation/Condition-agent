# Proposed Streamlined LLM Workflow

## Problem Analysis

Current workflow (with LLM) has issues:
- ❌ **Messy output**: Original deterministic entry + LLM review + adjusted entry = confusing
- ❌ **Redundant processing**: Deterministic role inference, then LLM tries to fix it
- ❌ **Complex logic**: Two parallel paths (deterministic + LLM) that need merging
- ❌ **Unclear precedence**: Which fields come from deterministic vs. LLM?

## Proposed New Workflow (LLM Mode Only)

When **LLM is enabled**, skip deterministic inference entirely and use a clean 4-step LLM-driven pipeline:

---

### **Step 1: Resolve Chemical Identity** 
**Source**: PubChem / CAS resolver (existing `resolve_identity_from_cas`)

**Gather**:
```json
{
  "cas": "121-44-8",
  "name": "Triethylamine",
  "synonyms": ["TEA", "N,N-Diethylethanamine", ...],
  "smiles": "CCN(CC)CC",
  "inchi_key": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
  "molecular_formula": "C6H15N",
  "source": "pubchem"
}
```

**Output**: Chemical identity bundle

---

### **Step 2: LLM Role Classification** 
**When**: User selects "other_reagent" OR requests LLM auto-classification

**LLM Prompt**:
```
You are a chemical reagent classifier. Given this reagent information:

Name: {name}
CAS: {cas}
SMILES: {smiles}
Synonyms: {synonyms}

Available roles:
- ligand: Phosphines, NHCs, diimines, donor ligands
- metal_precursor: Metal salts/complexes (Pd, Ni, Cu, etc.)
- preformed_metal_catalyst: Precatalysts with ligands
- base: Bronsted/Lewis bases, amides, alkoxides, carbonates
- acid: Mineral, sulfonic, Lewis acids
- condensation_agent: Carbodiimides, uronium/phosphonium activators
- oxidant: Peroxides, hypervalent iodine, O2, oxidizers
- reductant: Hydrides, silanes, metal powders
- additive: Phase-transfer agents, halide scavengers
- solvent: Reaction media
- organo_catalyst: Small-molecule organocatalysts
- enzyme: Biocatalysts
- other_reagent: Generic (use as fallback only)

Respond ONLY with JSON:
{
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary aliphatic amine with strong basicity"
}
```

**Output**: 
```json
{
  "role": "base",
  "confidence": 0.95,
  "reasoning": "Tertiary aliphatic amine with strong basicity"
}
```

**Validation**:
- ✅ Check `role` is in valid role list
- ✅ Check `confidence` >= 0.7 (threshold)
- ❌ Reject if confidence too low → prompt user to manually select

---

### **Step 3: LLM Family & Field Assignment**
**When**: After role is determined

**LLM Prompt** (role-specific):
```
You are a chemical database curator. Classify this reagent into the correct family and assign all required fields.

Reagent:
- Name: {name}
- CAS: {cas}
- SMILES: {smiles}
- Role: {role}

Available families for role '{role}': {available_families_with_definitions}

Required fields for role '{role}': {field_schema_with_allowed_values}

Examples from database: {2-3_examples_from_existing_registry}

Respond ONLY with JSON matching this schema:
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
  "reasoning": "Unhindered tertiary amine, pKa ~10.7"
}
```

**Field schemas by role**:

**Base**:
```json
{
  "basicity": ["weak", "moderate", "strong", "superbase"],
  "nucleophilicity": ["weak", "moderate", "strong"],
  "sterics": ["unhindered", "moderate", "hindered"]
}
```

**Metal Precursor**:
```json
{
  "metal": "Pd|Ni|Cu|Fe|...",
  "oxidation_states": [0, 1, 2, ...],
  "ligand_type": "string (optional)"
}
```

**Solvent**:
```json
{
  "proticity": ["aprotic", "protic"],
  "polarity": ["polar", "nonpolar", "intermediate"],
  "coordination": ["coordinating", "weakly_coordinating", "non_coordinating"]
}
```

**Ligand**:
```json
{
  "donors": ["P", "N", "O", ...],
  "denticity": 1|2|3|4|...
}
```

**Output**:
```json
{
  "family": "tertiary_amines_aliphatic",
  "fields": {
    "basicity": "strong",
    "nucleophilicity": "moderate",
    "sterics": "unhindered"
  },
  "abbreviations": ["TEA"],
  "additional_synonyms": ["Triethyl amine"],
  "confidence": 0.92,
  "reasoning": "Clear tertiary amine structure"
}
```

**Validation**:
- ✅ Check `family` exists in role's family registry
- ✅ Check all `fields` match allowed values (enums)
- ✅ Check required fields are present
- ⚠️ Warn if confidence < 0.8

---

### **Step 4: Final Verification & Assembly**
**When**: Before saving

**LLM Verification Prompt**:
```
You are a quality control reviewer. Check this proposed reagent entry for errors.

Proposed entry:
{full_entry_json}

Check for:
1. Obvious mistakes (wrong metal, impossible oxidation state, etc.)
2. Field consistency (e.g., "superbase" + "weak nucleophilicity" is suspicious)
3. Missing critical information
4. SMILES vs. field mismatch (e.g., SMILES shows Pd but metal="Cu")

Respond ONLY with JSON:
{
  "approved": true|false,
  "issues": [
    {"severity": "error|warning", "field": "field_name", "message": "description"}
  ],
  "suggestions": ["Consider adding...", ...]
}
```

**Output**:
```json
{
  "approved": true,
  "issues": [],
  "suggestions": ["Consider adding pKa value to notes"]
}
```

**Actions**:
- If `approved = false` → Show issues, let user edit
- If `approved = true` → Proceed to save

**Final entry format**:
```json
{
  "id": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
  "name": "Triethylamine",
  "abbreviation": ["TEA"],
  "aliases": ["N,N-Diethylethanamine", "Triethyl amine"],
  "cas": "121-44-8",
  "inchi_key": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
  "smiles": "CCN(CC)CC",
  "roles": {
    "base": {
      "families": ["tertiary_amines_aliphatic"],
      "basicity": "strong",
      "nucleophilicity": "moderate",
      "sterics": "unhindered"
    }
  },
  "embedding_text": "type: BASE | family: Tertiary amines (aliphatic) ..."
}
```

---

## Clean Output Structure

### Success Case:
```json
{
  "status": "ready_to_save",
  "cas": "121-44-8",
  "workflow": {
    "step1_identity": {
      "source": "pubchem",
      "name": "Triethylamine",
      "smiles": "CCN(CC)CC"
    },
    "step2_role": {
      "role": "base",
      "confidence": 0.95,
      "method": "llm_classification"
    },
    "step3_fields": {
      "family": "tertiary_amines_aliphatic",
      "fields": {"basicity": "strong", ...},
      "confidence": 0.92
    },
    "step4_verification": {
      "approved": true,
      "issues": [],
      "suggestions": []
    }
  },
  "entry": {
    "id": "ZMANZCXQSJIPKH-UHFFFAOYSA-N",
    "name": "Triethylamine",
    ...
  }
}
```

### Needs Review Case:
```json
{
  "status": "needs_review",
  "cas": "121-44-8",
  "workflow": {
    "step2_role": {
      "role": "base",
      "confidence": 0.65,
      "warning": "Confidence below threshold (0.7)"
    },
    "step4_verification": {
      "approved": false,
      "issues": [
        {
          "severity": "error",
          "field": "oxidation_states",
          "message": "Oxidation state [5] not valid for Pd"
        }
      ]
    }
  },
  "entry": {...}
}
```

---

## Comparison: Old vs. New

| Aspect | Old (Deterministic + LLM) | New (Pure LLM) |
|--------|---------------------------|----------------|
| **Clarity** | ❌ 3 versions of entry (original, llm_review, adjusted) | ✅ 1 clean entry with workflow metadata |
| **Speed** | ❌ Deterministic inference + LLM (2 passes) | ✅ Direct LLM (1 pass per step) |
| **Accuracy** | ⚠️ LLM tries to fix deterministic errors | ✅ LLM builds from scratch |
| **Output** | ❌ Messy JSON with changes_applied, llm_auto_upgrade, etc. | ✅ Clean entry + workflow steps |
| **User experience** | ❌ Confusing which fields are final | ✅ Clear: entry is final, workflow is audit trail |
| **Validation** | ⚠️ Partial (LLM suggests, system validates) | ✅ Full LLM verification step |

---

## Implementation Plan

### Phase 1: New LLM-only function
```python
def generate_taxonomy_entry_llm(
    *,
    cas: str,
    registry_dir: Path,
    llm_client: LLMClient,
    user_role_hint: Optional[str] = None,  # "other_reagent" or specific role
) -> Dict[str, Any]:
    """Pure LLM workflow - no deterministic inference."""
    
    # Step 1: Resolve identity
    identity = resolve_identity_from_cas(cas)
    
    # Step 2: LLM role classification (if needed)
    if user_role_hint == "other_reagent" or user_role_hint is None:
        role_result = llm_classify_role(identity, llm_client)
    else:
        role_result = {"role": user_role_hint, "confidence": 1.0, "method": "user_specified"}
    
    # Step 3: LLM family & field assignment
    fields_result = llm_assign_fields(identity, role_result["role"], llm_client, registry_dir)
    
    # Step 4: LLM verification
    entry = build_entry(identity, role_result, fields_result)
    verification = llm_verify_entry(entry, llm_client)
    
    # Return clean result
    return {
        "status": "ready_to_save" if verification["approved"] else "needs_review",
        "cas": cas,
        "workflow": {
            "step1_identity": identity,
            "step2_role": role_result,
            "step3_fields": fields_result,
            "step4_verification": verification,
        },
        "entry": entry,
    }
```

### Phase 2: UI Changes
- **Mode selector**: "Deterministic" vs. "LLM-powered"
- **When LLM mode**: Role selector includes "Auto-detect (LLM)"
- **Output panel**: Show workflow steps + final entry
- **Review mode**: Let user edit fields before saving

### Phase 3: Prompt Templates
Create 3 new prompts in `llmtools/prompts.py`:
1. `ROLE_CLASSIFICATION` - Step 2
2. `FIELD_ASSIGNMENT` (role-specific variants) - Step 3
3. `ENTRY_VERIFICATION` - Step 4

---

## Benefits

1. ✅ **Cleaner output**: Single entry + workflow audit trail
2. ✅ **Better accuracy**: LLM builds from scratch, not fixing deterministic errors
3. ✅ **Validation built-in**: Step 4 catches mistakes
4. ✅ **Transparent**: Workflow shows each decision with confidence
5. ✅ **Flexible**: User can override role, or let LLM decide
6. ✅ **Matches database format**: Direct mapping to existing JSON structure

---

## Next Steps

If you approve this design:

1. Create new prompts in `llmtools/prompts.py`
2. Create `llmtools/reagent_classifier.py` with 3 LLM functions
3. Create `generate_taxonomy_entry_llm()` in `reagent_taxonomy_qt.py`
4. Update UI to toggle between modes
5. Test with your existing CAS numbers

Should I proceed with implementation?
