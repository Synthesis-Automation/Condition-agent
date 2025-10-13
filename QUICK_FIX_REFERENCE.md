# Quick Fix Reference - Schema Path Resolution

## What Was Fixed

**Error**: `Families registry not found: C:\...\data\reagents\reagent_schema\families_registry.json`

**Solution**: Added automatic fallback to package-bundled schemas in `chemtools/reagent/reagent_schema/`

## Changed Files

- `llmtools/reagent_classifier.py` (2 functions updated)

## Three-Tier Path Resolution

```
1. registry_dir/reagent_schema/     ← User custom (if exists)
2. chemtools/reagent/reagent_schema/ ← Package bundled (NEW fallback)
3. Hardcoded defaults                ← Minimal fallback
```

## Before/After

### Before
```python
families_path = registry_dir / "reagent_schema" / "families_registry.json"
if not families_path.exists():
    return {"status": "error", "error": "..."}  # ❌ Failed immediately
```

### After
```python
families_path = registry_dir / "reagent_schema" / "families_registry.json"

if not families_path.exists():
    # Fallback to package
    import chemtools.reagent
    package_path = Path(chemtools.reagent.__file__).parent / "reagent_schema" / "families_registry.json"
    if package_path.exists():
        families_path = package_path  # ✅ Use package schemas
    else:
        return {"status": "error", "error": "..."}
```

## Testing

```bash
# Verify package schemas exist
python -c "from pathlib import Path; import chemtools.reagent; \
    print((Path(chemtools.reagent.__file__).parent / 'reagent_schema').exists())"
# Output: True ✅

# Run UI
python app/reagent_taxonomy_ui.py
```

## Expected Behavior

**LLM Workflow** (Use LLM mode):
- Step 1: ✅ Identity Resolution (PubChem)
- Step 2: ✅ Role Classification (LLM)
- Step 3: ✅ Field Assignment (LLM) ← **FIXED**
- Step 4: ✅ Entry Verification (LLM)

**Output** (Simplified):
```json
{
  "status": "ready_to_save",
  "entry": {
    "id": "...",
    "name": "...",
    "abbreviation": [...],  // ✅ Correct field name
    "aliases": [...],       // ✅ Added
    "cas": "...",
    "inchi_key": "...",     // ✅ Added
    "smiles": "...",
    "roles": {...},
    "embedding_text": "..." // ✅ Added
  }
}
```

## Complete Fix Chain

1. **Previous Session**: Entry building schema compliance
   - Fixed field names and structure
   - Added missing fields (id, embedding_text, etc.)

2. **Current Session**: Schema path resolution
   - Fixed schema file discovery
   - Zero-configuration setup

## Documentation

- `LLM_WORKFLOW_FIX_COMPLETE.md` - Complete summary
- `llmtools/SCHEMA_PATH_FIX.md` - Detailed technical docs
- `app/SCHEMA_COMPLIANCE_FIX.md` - Entry building fix

---

**Status**: ✅ Ready to test  
**Test Case**: DavePhos (CAS: 213697-53-1)
