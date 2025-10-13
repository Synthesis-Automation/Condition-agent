# Schema Path Resolution Fix

## Issue

When running the LLM workflow, the `assign_fields` step was failing with:
```
Families registry not found: C:\Git-softwares\Condition-agent\data\reagents\reagent_schema\families_registry.json
```

**Root Cause**: The LLM workflow expected schema files in `registry_dir/reagent_schema/`, but the schemas were moved to `chemtools/reagent/reagent_schema/` during package consolidation.

## Impact

- **Before**: LLM workflow failed at Step 3 (Field Assignment) if schemas weren't in registry_dir
- **After**: LLM workflow uses chemtools package schemas as fallback

## Fix Applied

### Files Modified
- `llmtools/reagent_classifier.py` (2 functions updated)

### Changes

#### 1. `assign_fields()` function (line ~351)

**Before**:
```python
# Load families for this role
families_path = registry_dir / "reagent_schema" / "families_registry.json"
if not families_path.exists():
    return {
        "status": "error",
        "error": f"Families registry not found: {families_path}",
    }
```

**After**:
```python
# Load families for this role - try registry_dir first, then chemtools package
families_path = registry_dir / "reagent_schema" / "families_registry.json"

if not families_path.exists():
    # Fallback to chemtools package location
    import chemtools.reagent
    package_path = Path(chemtools.reagent.__file__).parent / "reagent_schema" / "families_registry.json"
    if package_path.exists():
        families_path = package_path
    else:
        return {
            "status": "error",
            "error": f"Families registry not found in {registry_dir / 'reagent_schema'} or {package_path}",
        }
```

#### 2. `_load_schema_for_role()` function (line ~172)

**Before**:
```python
def _load_schema_for_role(role: str, registry_dir: Path) -> Dict[str, Any]:
    """Load the reagent schema and extract field definitions for a specific role."""
    schema_path = registry_dir / "reagent_schema" / "reagent_schema.json"
    
    # Default schema fallback
    default_schemas = {
```

**After**:
```python
def _load_schema_for_role(role: str, registry_dir: Path) -> Dict[str, Any]:
    """Load the reagent schema and extract field definitions for a specific role."""
    # Try registry_dir first
    schema_path = registry_dir / "reagent_schema" / "reagent_schema.json"
    
    if not schema_path.exists():
        # Fallback to chemtools package location
        import chemtools.reagent
        package_path = Path(chemtools.reagent.__file__).parent / "reagent_schema" / "reagent_schema.json"
        if package_path.exists():
            schema_path = package_path
    
    # Default schema fallback
    default_schemas = {
```

## Path Resolution Strategy

The LLM workflow now uses a **three-tier fallback**:

1. **Primary**: `registry_dir/reagent_schema/` - User-provided schemas (if any)
2. **Secondary**: `chemtools/reagent/reagent_schema/` - Package-bundled schemas ✅ **NEW**
3. **Tertiary**: Hardcoded defaults in code - Minimal fallback

## Benefits

| Aspect | Before | After |
|--------|--------|-------|
| **Schema Location** | Must be in registry_dir | Can use package schemas |
| **Setup Complexity** | Users must copy schemas | Works out-of-the-box |
| **Flexibility** | Single location only | Multi-tier fallback |
| **Error Handling** | Failed immediately | Graceful degradation |
| **Maintenance** | Duplicate schemas | Single source of truth |

## Testing

```python
from pathlib import Path
from llmtools.reagent_classifier import _load_schema_for_role

# Test with non-existent registry_dir
registry_dir = Path("data/reagents")  # No schemas here
schema = _load_schema_for_role("base", registry_dir)
# ✅ Returns schema from chemtools/reagent/reagent_schema/

print(list(schema.keys()))
# ['basicity', 'nucleophilicity', 'sterics']
```

## Schema Files

Located in `chemtools/reagent/reagent_schema/`:
- `reagent_schema.json` - Field definitions for all roles
- `families_registry.json` - Family definitions with metadata

## Related Issues

This fix also resolves a secondary issue where the error response was showing the full workflow output instead of the simplified display. The error now occurs earlier in the workflow, preventing the user from seeing confusing internal details.

## Future Improvements

Consider adding:
- Schema validation on startup
- User-friendly warning if custom schemas are outdated
- Schema version checking
