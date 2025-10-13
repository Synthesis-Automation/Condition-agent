# Save Error Fix - Missing Family Argument

## Issue

When clicking "Save" button in UI, error occurred:
```
TypeError: add_entry() missing 1 required positional argument: 'entry'
```

## Root Cause

The `ReagentRegistryStore.add_entry()` method signature:
```python
def add_entry(self, role: str, family_id: str, entry: Dict[str, Any]) -> None:
```

But the UI was calling it with only 2 arguments:
```python
store.add_entry(role_for_save, entry_to_save)  # ❌ Missing family_id
```

## Fix Applied

**File**: `app/reagent_taxonomy_ui.py`  
**Line**: ~1681

### Before
```python
role_for_save = next(iter(roles_payload.keys()))
role_details = roles_payload.get(role_for_save) or {}
families = role_details.get("families")
family_for_save: Optional[str] = families[0] if isinstance(families, list) and families else None
result_context = self._last_result or {}
try:
    store.add_entry(role_for_save, entry_to_save)  # ❌ Missing argument
```

### After
```python
role_for_save = next(iter(roles_payload.keys()))
role_details = roles_payload.get(role_for_save) or {}
families = role_details.get("families")
family_for_save: Optional[str] = families[0] if isinstance(families, list) and families else None

if not family_for_save:
    self.show_error(f"Entry must include a 'families' list in the '{role_for_save}' role.")
    return

result_context = self._last_result or {}
try:
    store.add_entry(role_for_save, family_for_save, entry_to_save)  # ✅ All 3 arguments
```

## Changes

1. ✅ Added `family_for_save` as second argument to `add_entry()`
2. ✅ Added validation check for missing family
3. ✅ Clear error message if family is missing

## Testing

The save operation now works correctly:
```python
# Entry structure
entry = {
    "roles": {
        "ligand": {
            "families": ["phosphine_biphenyl"],  # ← This is extracted
            ...
        }
    }
}

# Correct call
store.add_entry("ligand", "phosphine_biphenyl", entry)  # ✅
```

## Impact

- **Before**: Save button clicked → Error dialog → No save
- **After**: Save button clicked → Entry saved to `ligand.json` → Success ✅

---

**Status**: ✅ FIXED  
**Part of**: LLM Workflow Complete Fix (Session 2)
