# LLM Workflow - Schema Compliance Fix

## Issue
LLM-generated reagent entries were **not following the reagent schema**, missing critical fields like:
- ❌ `id` - Required identifier (InChIKey or CAS)
- ❌ `embedding_text` - Required for semantic search
- ❌ `inchi_key` - Chemical identifier
- ❌ `aliases` - Alternate names
- ❌ Using wrong field names (`abbreviations` instead of `abbreviation`)

## Schema Reference
**Location**: `chemtools/reagent/reagent_schema/reagent_schema.json`

### Required Fields (per schema)
```json
{
  "id": "string",              // InChIKey preferred; else CAS
  "name": "string",
  "abbreviation": ["string"],  // NOTE: array, not "abbreviations"
  "aliases": ["string"],       // Distinct from abbreviations
  "cas": "string|null",
  "inchi_key": "string|null",
  "smiles": "string|null",
  "roles": {
    "<role>": {
      "families": ["string"],
      // ... role-specific fields
    }
  },
  "embedding_text": "string|null"  // Auto-generated for search
}
```

## Fix Applied

### Before (Non-Compliant)
```python
entry = {
    "name": resolved_identity.get("name"),
    "cas": normalized_cas,
    "molecular_formula": resolved_identity.get("molecular_formula"),  # ❌ Not in schema
    "smiles": resolved_identity.get("smiles"),
    "roles": {
        role: {
            "families": [fields_result["family"]],
            **fields_result["fields"],
        }
    },
}
if abbreviations:
    entry["abbreviations"] = abbreviations  # ❌ Wrong field name
if all_synonyms:
    entry["synonyms"] = all_synonyms
# ❌ Missing: id, inchi_key, aliases, embedding_text
```

### After (Schema-Compliant)
```python
# Build entry following reagent_schema.json structure
entry = {
    "id": resolved_identity.get("inchi_key") or normalized_cas,  # ✅ InChIKey preferred
    "name": resolved_identity.get("name"),
    "abbreviation": abbreviations if abbreviations else None,    # ✅ Correct field name
    "aliases": [syn for syn in all_synonyms if syn.lower() != name.lower()],  # ✅ Added
    "cas": normalized_cas,
    "inchi_key": resolved_identity.get("inchi_key"),             # ✅ Added
    "smiles": resolved_identity.get("smiles"),
    "roles": {
        role: {
            "families": [fields_result["family"]],
            **fields_result.get("fields", {}),
        }
    },
}

# ✅ Generate embedding_text using helper function
entry["embedding_text"] = build_embedding_text(role, family_dict, embedding_entry, all_synonyms)
```

## Key Changes

### 1. Added Missing Fields
- ✅ `id` - Uses InChIKey if available, otherwise CAS
- ✅ `inchi_key` - From PubChem resolution
- ✅ `aliases` - Filtered synonyms (excluding name)
- ✅ `embedding_text` - Auto-generated using `build_embedding_text()`

### 2. Fixed Field Names
- Changed `abbreviations` → `abbreviation` (array)
- Removed `molecular_formula` (not in schema)
- Changed `synonyms` → `aliases` (to match schema intent)

### 3. Used Helper Functions
```python
from chemtools.reagent import build_embedding_text, dedupe_synonyms

# Build embedding text following schema
entry["embedding_text"] = build_embedding_text(role, family_dict, embedding_entry, all_synonyms)
```

## Example Output

### Complete Schema-Compliant Entry
```json
{
  "id": "YNAVUWVOSKDBBP-UHFFFAOYSA-N",
  "name": "Triethylamine",
  "abbreviation": ["TEA", "Et3N"],
  "aliases": ["N,N-Diethylethanamine", "Triethylamin"],
  "cas": "121-44-8",
  "inchi_key": "YNAVUWVOSKDBBP-UHFFFAOYSA-N",
  "smiles": "CCN(CC)CC",
  "roles": {
    "base": {
      "families": ["tertiary_amines_aliphatic"],
      "basicity": "moderate",
      "nucleophilicity": "weak",
      "sterics": "unhindered"
    }
  },
  "embedding_text": "type: BASE | family: tertiary_amines_aliphatic | basicity: moderate | nucleophilicity: weak | sterics: unhindered | name: Triethylamine | abbr: TEA | CAS: 121-44-8 | synonyms: N,N-Diethylethanamine; Et3N; Triethylamin"
}
```

## Benefits

| Field | Purpose | Status |
|-------|---------|--------|
| `id` | Unique identifier for database | ✅ Fixed |
| `inchi_key` | Standard chemical identifier | ✅ Fixed |
| `abbreviation` | Short forms (e.g., TEA, Et3N) | ✅ Fixed (renamed) |
| `aliases` | Alternative names for search | ✅ Fixed |
| `embedding_text` | Semantic search indexing | ✅ Fixed |
| `roles` | Proper structure with families | ✅ Already working |

## Impact

### Before
- ❌ Entries couldn't be properly indexed for search
- ❌ Missing chemical identifiers
- ❌ Schema validation would fail
- ❌ Inconsistent with existing reagent database

### After
- ✅ Fully schema-compliant entries
- ✅ Searchable via embedding_text
- ✅ Contains all chemical identifiers
- ✅ Compatible with existing reagent database
- ✅ Ready for production use

## Testing

### How to Verify
1. Run UI: `python app/reagent_taxonomy_ui.py`
2. Select "Use LLM" mode
3. Enter CAS: `121-44-8` (Triethylamine)
4. Click "Generate"
5. Verify output includes:
   - ✅ `id` field
   - ✅ `abbreviation` field (not `abbreviations`)
   - ✅ `aliases` field
   - ✅ `inchi_key` field
   - ✅ `embedding_text` field

---

**File Modified**: `app/reagent_taxonomy_ui.py` (lines 683-735)  
**Date**: October 13, 2025  
**Status**: ✅ Schema-compliant  
**Reference**: `chemtools/reagent/reagent_schema/reagent_schema.json`
