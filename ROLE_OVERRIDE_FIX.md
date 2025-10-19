# Role Override Fix - Summary

**Issue**: User-selected role was being ignored in Pure LLM workflow  
**Date Fixed**: October 19, 2025  
**Status**: ✅ FIXED

## Problem Description

When using the "Pure LLM Workflow" in the Reagent Taxonomy UI:
- User manually selects a role (e.g., "Metal precursor")
- System ignores the selection and asks LLM to classify the role
- LLM sometimes misclassifies (e.g., Copper → "reductant" instead of "metal_precursor")
- User's selection is completely ignored

### Example Issue

```json
{
  "name": "Copper",
  "roles": {
    "reductant": {  // ❌ Wrong! User selected "metal_precursor"
      "families": ["metal_powders"]
    }
  }
}
```

User explicitly selected "Metal precursor" but got "reductant" instead.

## Root Cause

The `generate_taxonomy_entry_llm()` function had no way to receive the user's role selection:

```python
def generate_taxonomy_entry_llm(
    cas: str,
    registry_dir: Path,
    llm_client: Any,
    name_override: Optional[str] = None,
    # ❌ Missing: role_override parameter
)
```

The function always called `classify_role()`, even when the user had already selected a role.

## Solution

### 1. Added `role_override` Parameter

```python
def generate_taxonomy_entry_llm(
    cas: str,
    registry_dir: Path,
    llm_client: Any,
    name_override: Optional[str] = None,
    role_override: Optional[str] = None,  # ✅ NEW
)
```

### 2. Updated Role Classification Logic

```python
# Step 2: Classify role (or use override)
if role_override:
    # ✅ Use the role provided by the user
    role = role_override
    workflow["step2_role"] = {
        "status": "success",
        "role": role,
        "source": "user_override",
    }
else:
    # LLM classifies the role (auto-detect)
    role_result = classify_role(resolved_identity, llm_client)
    role = role_result["role"]
```

### 3. Updated UI to Pass Role

```python
# Pass role if user selected one (skip if auto-detect)
role_to_pass = role if role and role != "__auto_detect__" else None

params = {
    "workflow_mode": "llm_workflow",
    "role_override": role_to_pass,  # ✅ NEW
    # ... other params
}
```

### 4. Updated Worker

```python
result = generate_taxonomy_entry_llm(
    cas=cas,
    registry_dir=registry_dir,
    llm_client=llm_client,
    name_override=name_override,
    role_override=role_override,  # ✅ NEW
)
```

## Workflow Comparison

### Before Fix ❌

```
User selects "Metal precursor"
    ↓
UI ignores selection
    ↓
Function calls LLM classify_role()
    ↓
LLM returns "reductant" (wrong!)
    ↓
Entry saved with wrong role
```

### After Fix ✅

```
User selects "Metal precursor"
    ↓
UI passes role_override="metal_precursor"
    ↓
Function skips LLM classification
    ↓
Uses role="metal_precursor" directly
    ↓
Entry saved with correct role
```

### Auto-Detect Still Works ✅

```
User selects "🤖 Auto-detect (LLM)"
    ↓
UI passes role_override=None
    ↓
Function calls LLM classify_role()
    ↓
LLM determines appropriate role
    ↓
Entry saved with LLM-selected role
```

## Files Modified

1. **`app/reagent_taxonomy_ui.py`**
   - Added `role_override` parameter to `generate_taxonomy_entry_llm()`
   - Updated role classification logic to check for override first
   - Updated UI to pass role from dropdown
   - Updated worker to pass role_override

## Testing

Created test script `test_role_override_fix.py` that validates:
- ✅ User-selected roles are passed correctly
- ✅ Auto-detect still triggers LLM classification
- ✅ None/empty selections are handled properly

All tests pass.

## Impact

- ✅ User selections are now respected
- ✅ No breaking changes to existing functionality
- ✅ Auto-detect still works when selected
- ✅ Backward compatible (role_override is optional)

## How to Use

### Option 1: Manual Role Selection (Recommended for known reagents)

1. Open Reagent Taxonomy UI
2. Enter CAS number
3. Select specific role from dropdown (e.g., "Metal precursor")
4. Click "Generate with LLM"
5. ✅ Entry will use your selected role

### Option 2: Auto-Detect (For unknown reagents)

1. Open Reagent Taxonomy UI
2. Enter CAS number
3. Select "🤖 Auto-detect (LLM)" from dropdown
4. Click "Generate with LLM"
5. ✅ LLM will classify the role

## Example: Copper as Metal Precursor

**Before Fix**:
```json
{
  "name": "Copper",
  "roles": {
    "reductant": {  // ❌ Wrong!
      "families": ["metal_powders"]
    }
  }
}
```

**After Fix** (with "Metal precursor" selected):
```json
{
  "name": "Copper",
  "roles": {
    "metal_precursor": {  // ✅ Correct!
      "metal": "Cu",
      "oxidation_states": [0, 1, 2]
    }
  }
}
```

## Related Issues

This fix resolves:
- User selections being ignored in Pure LLM workflow
- Misclassification of reagents when user knows the correct role
- Confusion about whether role selection has any effect

## Future Enhancements

Consider:
- Add visual feedback when role_override is used (show "Using your selection" message)
- Log source of role in workflow output ("user_override" vs "llm_classified")
- Add option to "Review LLM suggestion" before accepting auto-detect result
