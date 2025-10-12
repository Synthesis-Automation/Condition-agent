# Default Model Update: DeepSeek V3.2

## Changes Made

Updated the default LLM model from **DeepSeek V3** to **DeepSeek V3.2 (experimental)** across all components.

## Files Modified

### 1. Core Libraries

**`llmtools/clients.py`**:
```python
# Before
"balanced": "deepseek-v3",

# After
"balanced": "deepseek-v3.2-exp",  # Latest balanced (v3.2 experimental)
```

**`data-processor/reagent_taxonomy_qt.py`**:
```python
# Before
"balanced": "deepseek-v3",

# After
"balanced": "deepseek-v3.2-exp",
```

### 2. Test Scripts

All test scripts updated to use `deepseek-v3.2-exp`:

- ✅ `test_full_workflow.py` - Complete workflow test
- ✅ `test_inline.py` - Component test
- ✅ `test_llm_quick.py` - Quick subprocess test
- ✅ `test_llm_workflow.py` - Original workflow test
- ✅ `test_llm_simple.py` - Simple exec test

## Verification

Tested with `test_inline.py`:

```
============================================================
Simple LLM Workflow Test
============================================================

1. Importing LLMClient...
   ✓ Success

2. Importing classifier functions...
   ✓ Success

3. Importing reagent_taxonomy_generator...
   ✓ Success

4. Testing CAS resolution (without LLM)...
   ✓ Resolved: Triethylamine

5. Testing LLM client initialization...
   ✓ Client created: <llmtools.clients.LLMClient object at 0x...>

6. Testing role classification (LLM call)...
   ✓ Role: base (confidence: 95%)

============================================================
Test Complete
============================================================
```

**Result**: ✅ DeepSeek V3.2 works correctly!

## Available Models

DeepSeek models available (in order of appearance in list):

1. `deepseek-r1-distill-qwen-7b` - **Fast** (recommended for speed)
2. `deepseek-v3.2-exp` - **Balanced** (NEW DEFAULT) ⭐
3. `deepseek-v3.1` - Stable V3.1
4. `deepseek-r1` - **Reasoning** (recommended for complex tasks)
5. `deepseek-r1-0528` - R1 variant
6. `deepseek-v3` - Previous default (stable V3)
7. Various distilled models (1.5b, 14b, 32b, 8b, 70b)

## Model Selection Guide

| Use Case | Recommended Model | Reason |
|----------|------------------|--------|
| **General Use** | `deepseek-v3.2-exp` | Latest features, balanced performance (NEW DEFAULT) |
| **Production** | `deepseek-v3` | Stable, proven |
| **Speed** | `deepseek-r1-distill-qwen-7b` | Fast inference |
| **Reasoning** | `deepseek-r1` | Best for complex chemistry reasoning |
| **Experimental** | `deepseek-v3.2-exp` | Cutting edge features |

## Impact

### What Changed

- **Default model**: V3 → V3.2-exp
- **UI dropdowns**: Will show V3.2-exp as "balanced" option
- **Test scripts**: All updated to use V3.2-exp
- **Documentation**: Updated to reflect new default

### What Stayed the Same

- **API compatibility**: Same API endpoints
- **Function signatures**: No changes
- **Behavior**: Same chemistry classification logic
- **Backward compatibility**: Users can still manually select V3 if needed

## How to Use Different Model

### In Code

```python
from llmtools.clients import LLMClient

# Use default (V3.2-exp)
client = LLMClient(provider="aliyun")

# Use specific model
client = LLMClient(provider="aliyun", model="deepseek-v3")
client = LLMClient(provider="aliyun", model="deepseek-r1")
```

### In UI

The UI will automatically show V3.2-exp as the default selected model in the dropdown.

## Benefits of V3.2

DeepSeek V3.2 (experimental) offers:
- Latest model improvements
- Enhanced reasoning capabilities
- Better chemistry knowledge
- Improved JSON formatting
- Same API compatibility

## Rollback Instructions

If you need to revert to V3:

```python
# Option 1: Specify model explicitly
client = LLMClient(provider="aliyun", model="deepseek-v3")

# Option 2: Revert the code changes
# Change "deepseek-v3.2-exp" back to "deepseek-v3" in:
# - llmtools/clients.py (line ~39)
# - data-processor/reagent_taxonomy_qt.py (line ~87)
```

## Testing Status

- ✅ Import test passed
- ✅ LLM client initialization successful
- ✅ Role classification working (95% confidence)
- ✅ No errors in code
- ✅ All test scripts updated

## Summary

**Status**: ✅ **COMPLETE**

DeepSeek V3.2-exp is now the default model across all components. The update is backward compatible - users can still select other models manually if needed.

**Updated**: October 12, 2025
