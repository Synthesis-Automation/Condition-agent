# Fix: DeepSeek Markdown Code Fence Issue

## Problem

DeepSeek (Aliyun provider) and some other LLMs wrap their JSON responses in markdown code fences:

```
```json
{
  "status": "reject",
  "proposed_role": "ligand",
  ...
}
```
```

This causes `json.loads()` to fail with:
```
JSONDecodeError: Expecting value: line 1 column 1 (char 0)
```

## Solution

Added `_strip_markdown_fences()` helper function in `llmtools/reagent_review.py` that:

1. Detects markdown code fences (` ```json ... ``` `)
2. Strips both opening and closing fences
3. Handles variations: `json`, `JSON`, no language specifier
4. Preserves the actual JSON content

## Code Changes

### File: `llmtools/reagent_review.py`

**Added helper function**:
```python
def _strip_markdown_fences(content: str) -> str:
    """
    Remove markdown code fences from LLM response.
    
    Some LLM providers (especially DeepSeek) wrap JSON in ```json ... ``` fences.
    This function extracts the content between fences.
    """
    content = content.strip()
    
    # Check for markdown code fence with language specifier
    if content.startswith("```json") or content.startswith("```JSON"):
        # Remove opening fence
        content = content[7:].lstrip()
    elif content.startswith("```"):
        # Remove opening fence without language
        content = content[3:].lstrip()
    
    # Remove closing fence
    if content.endswith("```"):
        content = content[:-3].rstrip()
    
    return content.strip()
```

**Updated parsing logic** (line ~155):
```python
if raw_content:
    # Strip markdown code fences that some LLMs add
    cleaned_content = _strip_markdown_fences(raw_content)
    try:
        parsed = json.loads(cleaned_content)
    except json.JSONDecodeError as exc:
        parse_error = f"{exc.__class__.__name__}: {exc}"
```

## Testing

Verified with your actual DeepSeek response - see `test_markdown_fix.py`:

```
✅ SUCCESS! JSON parsed correctly

Parsed data:
  Status: reject
  Proposed role: ligand
  Proposed family: phosphine
  Confidence: 0.9
  Field suggestions: {'SMILES': '...', 'molecular_formula': 'C15H15P'}
```

Also tested edge cases:
- ✅ Generic fence (` ``` `)
- ✅ Uppercase JSON (` ```JSON `)
- ✅ No fence (passthrough)
- ✅ Extra whitespace

## Impact

**Before**: DeepSeek responses failed with parse_error
**After**: All LLM providers work correctly, including DeepSeek

## Compatibility

- ✅ Backward compatible (handles responses without fences)
- ✅ Works with OpenAI (no fences)
- ✅ Works with Aliyun/DeepSeek (with fences)
- ✅ No performance impact (simple string operations)

## Example Usage

Now your DeepSeek LLM review will work:

```json
{
  "llm_review": {
    "enabled": true,
    "status": "ok",  // ← Changed from "parse_error"
    "provider": "aliyun",
    "model": "deepseek-v3",
    "analysis": {
      "status": "reject",
      "proposed_role": "ligand",
      "proposed_family": "phosphine",
      "confidence": 0.9,
      "field_suggestions": {
        "SMILES": "C=CCP(c1ccccc1)c2ccccc2",
        "molecular_formula": "C15H15P"
      }
    }
  }
}
```

## Testing Your Fix

Run the test:
```powershell
python data-processor\test_markdown_fix.py
```

Or test in the UI:
1. Enable LLM with Aliyun/DeepSeek provider
2. Generate a reagent entry
3. Verify `llm_review.status == "ok"` (not "parse_error")
4. Check that `llm_review.analysis` is populated

## Related Files

- **Fix**: `llmtools/reagent_review.py`
- **Test**: `data-processor/test_markdown_fix.py`
- **Docs**: This file

Fixed! 🎉
