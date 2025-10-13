# Save Function Fix - Corrected Parameter Count

## Error
```
ReagentRegistryStore.add_entry() takes 3 positional arguments but 4 were given
```

## Root Cause

There are **TWO different `ReagentRegistryStore` classes** with **different signatures**:

### 1. chemtools Package Version
**File**: `chemtools/reagent/taxonomy_store.py`
```python
def add_entry(self, role: str, family_id: str, entry: Dict[str, Any]) -> None:
    # 3 parameters: role, family_id, entry
```

### 2. UI Local Version
**File**: `app/reagent_taxonomy_ui.py`
```python
def add_entry(self, role: str, entry: Dict[str, Any]) -> None:
    # 2 parameters: role, entry (NO family_id)
```

## The Mistake

I incorrectly assumed both classes had the same signature and "fixed" the call to:
```python
store.add_entry(role_for_save, family_for_save, entry_to_save)  # ❌ 3 args
```

But the UI version only needs:
```python
store.add_entry(role_for_save, entry_to_save)  # ✅ 2 args
```

## Fix Applied

**File**: `app/reagent_taxonomy_ui.py` (line ~1695)

### Reverted To
```python
role_for_save = next(iter(roles_payload.keys()))
role_details = roles_payload.get(role_for_save) or {}
families = role_details.get("families")
family_for_save: Optional[str] = families[0] if isinstance(families, list) and families else None
result_context = self._last_result or {}
try:
    store.add_entry(role_for_save, entry_to_save)  # ✅ 2 args only
```

## Why Two Different Classes?

The UI has a **simplified local version** of `ReagentRegistryStore` that:
- Doesn't need family-level organization
- Just appends entries to role files
- Simpler interface for UI needs

The chemtools package has a **full-featured version** that:
- Organizes entries by family
- Maintains family-level indexes
- Used by taxonomy generation tools

## Impact

| Action | Before | After |
|--------|--------|-------|
| **Click Save** | ❌ TypeError (4 args) | ✅ Works (2 args) |
| **Entry saved** | ❌ Failed | ✅ Appended to role file |

---

**Status**: ✅ FIXED  
**Correct call**: `store.add_entry(role, entry)` - 2 parameters only
