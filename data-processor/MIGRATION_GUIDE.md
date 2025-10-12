# Migration Guide: Old → New LLM Workflow

## Overview

This guide helps you transition from the old mixed workflow to the new pure LLM workflow.

## Side-by-Side Comparison

### Old Workflow (Still Available)

```python
from reagent_taxonomy_qt import generate_taxonomy_entry

result = generate_taxonomy_entry(
    cas="121-44-8",
    role="base",                    # ⚠️ Manual role selection required
    registry_dir=Path("data/reagents"),
    allow_default_family=True,
    dry_run=False,
    llm_options={
        "provider": "aliyun",
        "model": "deepseek-v3",
        "enable_llm": True,
    },
)

# Output: 7+ keys (messy)
# {
#   "entry_original": {...},
#   "llm_adjusted_entry": {...},
#   "review_summary": {...},
#   "changes_applied": [...],
#   "field_suggestions": {...},
#   "needs_attention": [...],
#   "success": true
# }
```

### New Workflow (Recommended)

```python
from llmtools.clients import LLMClient
from reagent_taxonomy_qt import generate_taxonomy_entry_llm

client = LLMClient(provider="aliyun", model="deepseek-v3")

result = generate_taxonomy_entry_llm(
    cas="121-44-8",
    # ✅ No role parameter - LLM auto-detects!
    registry_dir=Path("data/reagents"),
    llm_client=client,
)

# Output: 3 keys (clean)
# {
#   "status": "ready_to_save",
#   "workflow": {...},
#   "entry": {...}
# }
```

## Key Differences

| Aspect | Old Workflow | New Workflow |
|--------|-------------|--------------|
| **Role Selection** | Manual (`role="base"`) | Auto-detect (LLM) |
| **Output Keys** | 7+ keys | 3 keys |
| **Final Entry** | Unclear (`entry_original` vs `llm_adjusted_entry`) | Clear (`entry` key) |
| **LLM Usage** | Optional, review-only | Required, full pipeline |
| **Quality Check** | Manual | Built-in (Step 4) |
| **Confidence** | ❌ Not provided | ✅ Per-step scores |
| **Workflow Trace** | ❌ Hidden | ✅ Full trace |
| **Status** | `success: true/false` | `ready_to_save` / `needs_review` / `error` |

## Migration Steps

### Step 1: Update Imports

**Before**:
```python
from reagent_taxonomy_qt import generate_taxonomy_entry
```

**After**:
```python
from llmtools.clients import LLMClient
from reagent_taxonomy_qt import generate_taxonomy_entry_llm
```

### Step 2: Initialize LLM Client

**Before**:
```python
llm_options = {
    "provider": "aliyun",
    "model": "deepseek-v3",
    "enable_llm": True,
}
```

**After**:
```python
client = LLMClient(provider="aliyun", model="deepseek-v3")
```

### Step 3: Call New Function

**Before**:
```python
result = generate_taxonomy_entry(
    cas=cas_number,
    role=role,                      # Manual
    registry_dir=registry_dir,
    allow_default_family=True,
    dry_run=False,
    llm_options=llm_options,
)
```

**After**:
```python
result = generate_taxonomy_entry_llm(
    cas=cas_number,
    # No role parameter!
    registry_dir=registry_dir,
    llm_client=client,
)
```

### Step 4: Update Result Handling

**Before**:
```python
if result.get("success"):
    entry = result.get("llm_adjusted_entry") or result.get("entry_original")
    # ⚠️ Which entry is correct?
    print(f"Generated: {entry['name']}")
else:
    print("Failed")
```

**After**:
```python
status = result.get("status")

if status == "ready_to_save":
    entry = result["entry"]        # ✅ Always this one
    print(f"✅ Generated: {entry['name']}")
    # Save to registry
    
elif status == "needs_review":
    entry = result["entry"]
    issues = result["workflow"]["step4_verification"]["issues"]
    print(f"⚠️ Generated but has {len(issues)} issue(s)")
    # Show issues, allow manual review
    
else:  # "error"
    error = result.get("error")
    print(f"❌ Failed: {error}")
```

## Complete Example

### Old Code

```python
from pathlib import Path
from reagent_taxonomy_qt import generate_taxonomy_entry

registry_dir = Path("data/reagents")

# Process triethylamine
result = generate_taxonomy_entry(
    cas="121-44-8",
    role="base",                    # Had to know this!
    registry_dir=registry_dir,
    allow_default_family=True,
    dry_run=False,
    llm_options={
        "provider": "aliyun",
        "model": "deepseek-v3",
        "enable_llm": True,
    },
)

# Extract final entry (confusing)
if result.get("success"):
    # Which one???
    entry = result.get("llm_adjusted_entry") or result.get("entry_original")
    
    # Check if LLM made changes
    changes = result.get("changes_applied", [])
    if changes:
        print(f"LLM made {len(changes)} changes")
    
    # Check for issues
    needs_attention = result.get("needs_attention", [])
    if needs_attention:
        print(f"⚠️ Issues: {needs_attention}")
    
    # Save entry (if looks good)
    save_to_registry(entry)
```

### New Code

```python
from pathlib import Path
from llmtools.clients import LLMClient
from reagent_taxonomy_qt import generate_taxonomy_entry_llm

registry_dir = Path("data/reagents")
client = LLMClient(provider="aliyun", model="deepseek-v3")

# Process triethylamine
result = generate_taxonomy_entry_llm(
    cas="121-44-8",
    # No role parameter - LLM detects it!
    registry_dir=registry_dir,
    llm_client=client,
)

# Handle result (clear)
if result["status"] == "ready_to_save":
    entry = result["entry"]        # Always this one!
    
    # Get workflow details
    role = result["workflow"]["step2_role"]["role"]
    confidence = result["workflow"]["step2_role"]["confidence"]
    print(f"✅ Detected role: {role} (confidence: {confidence:.0%})")
    
    # Save entry directly
    save_to_registry(entry)
    
elif result["status"] == "needs_review":
    entry = result["entry"]
    issues = result["workflow"]["step4_verification"]["issues"]
    
    print(f"⚠️ Generated but needs review:")
    for issue in issues:
        print(f"  [{issue['severity']}] {issue['field']}: {issue['message']}")
    
    # Allow user to review and decide
    if user_approves(entry, issues):
        save_to_registry(entry)
        
else:  # "error"
    print(f"❌ Failed: {result['error']}")
```

## Backward Compatibility

**Good news**: Old workflow still works!

You can migrate gradually:
1. Keep old code for existing scripts
2. Use new workflow for new features
3. Migrate old code when convenient

Both functions coexist:
- `generate_taxonomy_entry()` - Old workflow (mixed)
- `generate_taxonomy_entry_llm()` - New workflow (pure LLM)

## Batch Processing Migration

### Old Code

```python
for cas in cas_list:
    result = generate_taxonomy_entry(
        cas=cas,
        role=role_map.get(cas, "other_reagent"),  # Manual mapping
        registry_dir=registry_dir,
        allow_default_family=True,
        dry_run=False,
        llm_options=llm_options,
    )
    
    if result.get("success"):
        entry = result.get("llm_adjusted_entry") or result.get("entry_original")
        entries.append(entry)
```

### New Code

```python
client = LLMClient(provider="aliyun", model="deepseek-v3")

for cas in cas_list:
    result = generate_taxonomy_entry_llm(
        cas=cas,
        # No role needed - LLM auto-detects!
        registry_dir=registry_dir,
        llm_client=client,
    )
    
    if result["status"] in ["ready_to_save", "needs_review"]:
        entry = result["entry"]
        entries.append(entry)
        
        # Log if needs review
        if result["status"] == "needs_review":
            log_for_review(cas, result["workflow"]["step4_verification"]["issues"])
```

## UI Integration Example

### Old UI Code

```python
def on_generate_clicked(self):
    cas = self.cas_input.text()
    role = self.role_combo.currentText()  # User selects
    
    result = generate_taxonomy_entry(
        cas=cas,
        role=role,
        registry_dir=self.registry_dir,
        allow_default_family=True,
        dry_run=False,
        llm_options=self.get_llm_options(),
    )
    
    # Display messy output
    self.output_text.setPlainText(json.dumps(result, indent=2))
```

### New UI Code

```python
def on_generate_clicked(self):
    cas = self.cas_input.text()
    # No role combo needed if using LLM mode!
    
    client = LLMClient(
        provider=self.provider_combo.currentText(),
        model=self.model_combo.currentText(),
    )
    
    result = generate_taxonomy_entry_llm(
        cas=cas,
        registry_dir=self.registry_dir,
        llm_client=client,
    )
    
    # Display clean output with workflow steps
    self.display_workflow(result)
    
    # Enable/disable save button based on status
    self.save_button.setEnabled(result["status"] in ["ready_to_save", "needs_review"])
```

## Common Pitfalls

### ❌ Pitfall 1: Checking Wrong Status

**Old habit**:
```python
if result.get("success"):
    # ...
```

**New way**:
```python
if result["status"] == "ready_to_save":
    # ...
```

### ❌ Pitfall 2: Using Wrong Entry Key

**Old habit**:
```python
entry = result.get("llm_adjusted_entry") or result.get("entry_original")
```

**New way**:
```python
entry = result["entry"]  # Always this one!
```

### ❌ Pitfall 3: Requiring Role Parameter

**Old habit**:
```python
result = generate_taxonomy_entry_llm(
    cas=cas,
    role="base",  # ❌ This parameter doesn't exist!
    ...
)
```

**New way**:
```python
result = generate_taxonomy_entry_llm(
    cas=cas,
    # No role parameter - LLM detects it!
    ...
)

# Get detected role from workflow
role = result["workflow"]["step2_role"]["role"]
```

## Performance Considerations

**Old Workflow**:
- 1 LLM call (optional review step)
- Fast if LLM disabled

**New Workflow**:
- 3 LLM calls (role, fields, verification)
- Total time: ~3-5 seconds
- More expensive but more accurate

**Tip**: For batch processing, consider:
- Parallel processing (multiple workers)
- Caching for duplicate CAS numbers
- Using faster models (e.g., deepseek-r1-distill-qwen-7b)

## Troubleshooting

### "LLM support not available"

**Old workflow**: Works without LLM
**New workflow**: Requires LLM

**Solution**:
```bash
pip install openai  # or your LLM provider
export ALIYUN_API_KEY=your-key
```

### "Missing role parameter"

**Old workflow**: `role` is required
**New workflow**: `role` is auto-detected

**Solution**: Remove `role` parameter from call

### "Can't find final entry"

**Old workflow**: Check `llm_adjusted_entry` or `entry_original`
**New workflow**: Always use `result["entry"]`

## Migration Checklist

- [ ] Install LLM dependencies (`pip install openai`)
- [ ] Set API key environment variable
- [ ] Update imports (add `LLMClient`, `generate_taxonomy_entry_llm`)
- [ ] Replace function calls
- [ ] Remove manual role selection
- [ ] Update result handling (status-based)
- [ ] Update entry extraction (`result["entry"]`)
- [ ] Test with example CAS numbers
- [ ] Update UI (if applicable)
- [ ] Update documentation/comments

## Need Help?

Just ask! I can:
- Help migrate specific code
- Explain any differences
- Add features to new workflow
- Create custom migration scripts

## Summary

**Why migrate?**
- ✅ Cleaner output (3 keys vs 7+)
- ✅ Auto-detect role (no manual selection)
- ✅ Built-in quality check
- ✅ Full workflow transparency
- ✅ Future-proof design

**When NOT to migrate?**
- LLM API not available
- Cost is a concern (3 LLM calls)
- Need offline/deterministic mode
- Existing code works fine

**Bottom line**: New workflow is cleaner and more powerful, but old workflow still works for backward compatibility.
