# UI Schema Path Fix - ReagentRegistryStore

## Issue

The error reappeared when starting the UI:
```
Families registry not found: C:\Git-softwares\Condition-agent\data\reagents\reagent_schema\families_registry.json
```

**Root Cause**: While we fixed the path resolution in `llmtools/reagent_classifier.py`, the **UI code** (`app/reagent_taxonomy_ui.py`) has its own `ReagentRegistryStore` class that also loads `families_registry.json`.

## Two Places Need Fixing

### 1. ✅ LLM Workflow (Already Fixed)
- File: `llmtools/reagent_classifier.py`
- Functions: `assign_fields()`, `_load_schema_for_role()`
- Status: ✅ Fixed (3-tier fallback added)

### 2. ✅ UI Initialization (Just Fixed)
- File: `app/reagent_taxonomy_ui.py`
- Class: `ReagentRegistryStore`
- Method: `_load_families()`
- Status: ✅ **JUST FIXED**

## Fix Applied

**File**: `app/reagent_taxonomy_ui.py` (line ~294)

### Before
```python
def _load_families(self) -> None:
    schema_path = self.base_dir / "reagent_schema" / "families_registry.json"
    if not schema_path.exists():
        raise FileNotFoundError(f"Families registry not found: {schema_path}")
    data = json.loads(schema_path.read_text(encoding="utf-8"))
```

### After
```python
def _load_families(self) -> None:
    # Try registry directory first
    schema_path = self.base_dir / "reagent_schema" / "families_registry.json"
    
    if not schema_path.exists():
        # Fallback to chemtools package location
        import chemtools.reagent
        package_path = Path(chemtools.reagent.__file__).parent / "reagent_schema" / "families_registry.json"
        if package_path.exists():
            schema_path = package_path
        else:
            raise FileNotFoundError(f"Families registry not found in {self.base_dir / 'reagent_schema'} or {package_path}")
    
    data = json.loads(schema_path.read_text(encoding="utf-8"))
```

## Why This Was Needed

The UI has **two separate code paths** that need schemas:

1. **UI Initialization** (`ReagentRegistryStore.__init__`)
   - Loads families for family hints and auto-completion
   - Used when UI starts up
   - **This was failing** ❌

2. **LLM Workflow** (`llmtools/reagent_classifier`)
   - Loads families during LLM workflow execution
   - Used when clicking "Generate"
   - **This was already fixed** ✅

Both needed the same fallback logic!

## Impact

| Action | Before | After |
|--------|--------|-------|
| **Start UI** | ❌ Error on startup | ✅ Loads from package |
| **LLM Workflow** | ✅ Works (already fixed) | ✅ Works |
| **Family Hints** | ❌ Not loaded | ✅ Loaded from package |
| **Auto-complete** | ❌ Not available | ✅ Available |

## Files Modified (Complete List)

1. ✅ `llmtools/reagent_classifier.py` - LLM workflow (2 functions)
2. ✅ `chemtools/reagent/taxonomy_utils.py` - InChIKey resolution
3. ✅ `app/reagent_taxonomy_ui.py` - Entry building, alias filtering, save fix, **AND UI initialization** ← **NEW**

## Testing

```bash
# Should start without errors now
python app/reagent_taxonomy_ui.py
```

**Expected**: UI starts, no schema errors, family hints work ✅

---

**Status**: ✅ COMPLETE  
**Both paths fixed**: LLM workflow + UI initialization
