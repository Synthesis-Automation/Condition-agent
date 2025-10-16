# Web CLI Sync Complete

## Summary
Updated `scripts/web_recommendation_cli.py` to match the UI and behavior of `scripts/local_recommendation_cli.py`.

## Changes Applied

### 1. **Import Updates**
- ✅ Added `choose_catalyst` import
- ✅ Added `summarize_llm_synthesis` import
- ❌ Removed `summarize_fusion` import (deprecated)
- ❌ Removed `FUSION_VARIANTS_DEFAULT` import (deprecated)

### 2. **Command-Line Arguments**
- ✅ Added `--catalyst` with choices: `["None", "Pd", "Cu", "Ni", "other"]`
- ✅ Added `--llm-timeout` (default: 180 seconds)
- ✅ Updated `--family` choices to: `["Suzuki", "C_N_Coupling", "Amide_formation"]`
  - Removed deprecated catalyst-specific variants (C_N_Coupling_Cu, Buchwald_CN, etc.)
- ✅ Updated `--strategy` choices to: `["all", "rule", "ml", "protocol", "llm_synthesis"]`
  - Removed deprecated `fusion` strategy
- ❌ Removed `--fusion-variants` argument (deprecated)

### 3. **Main Workflow Logic**
- ✅ **Auto-Detection**: Calls `detect_family_from_reaction()` when user selects "Auto-detect"
- ✅ **Mandatory Catalyst Selection**: Always prompts for catalyst (5 options)
  - Options: "Catalyst - optional; any or none" (default), Pd, Cu, Ni, other
- ✅ **Catalyst Parameter Passing**: Passes `catalyst_preference` to all API calls via `relax` parameter

### 4. **API Functions Updated**

#### `call_rule_based()`
```python
# Added parameter:
catalyst_preference: Optional[str] = None

# Passes catalyst via relax:
if catalyst_preference:
    payload["relax"] = {"catalyst_class": catalyst_preference}
```

#### `call_ml_recommendation()`
```python
# Added parameter:
catalyst_preference: Optional[str] = None

# Passes catalyst via relax:
if catalyst_preference:
    payload["relax"]["catalyst_class"] = catalyst_preference
```

#### `call_llm_synthesis()` - NEW
```python
def call_llm_synthesis(
    base_url: str,
    reaction: str,
    family: Optional[str] = None,
    catalyst_preference: Optional[str] = None,
    timeout: int = 180,
) -> Dict[str, Any]:
    """Execute the LLM synthesis recommendation endpoint."""
    # Calls: POST /api/v1/recommend/llm_synthesis
    # Uses configurable timeout (default 180s)
```

#### `call_fusion_recommendation()` - REMOVED
- Deprecated function removed entirely
- Users should use `call_ml_recommendation()` with `rerank_strategy='rule'`

### 5. **Execution Flow**
```
1. Prompt for reaction SMILES (or use --rxn)
2. Prompt for reaction type (or use --family)
3. Auto-detect if type is None
   └─> detect_family_from_reaction(reaction)
4. Prompt for catalyst (MANDATORY, or use --catalyst)
5. Run selected strategies:
   - rule: call_rule_based() with catalyst_preference
   - ml: call_ml_recommendation() with catalyst_preference
   - protocol: call_protocol_recommendation()
   - llm_synthesis: call_llm_synthesis() with catalyst_preference + timeout
6. Summarize and save results
```

### 6. **Usage Examples Updated**
```bash
# Interactive mode (prompts for input):
python scripts/web_recommendation_cli.py

# Provide reaction, auto-detect type, specify catalyst:
python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --catalyst Cu

# Specify reaction type and catalyst:
python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --family C_N_Coupling --catalyst Cu

# Run only LLM synthesis with longer timeout:
python scripts/web_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --strategy llm_synthesis --llm-timeout 300

# Custom server URL:
python scripts/web_recommendation_cli.py --url http://myserver:8080
```

## Key Improvements

### ✅ Consistent UX
- Web CLI now matches local CLI behavior exactly
- Same prompts, same flow, same options

### ✅ Auto-Detection
- Properly detects reaction type before catalyst prompt
- Uses `chemtools.router.detect_family_from_reaction()`

### ✅ Mandatory Catalyst Selection
- Always prompts for catalyst preference
- 5 options with clear default ("Catalyst - optional; any or none")

### ✅ LLM Timeout Fix
- Configurable timeout (default 180s vs old 60s)
- `--llm-timeout` argument for flexibility

### ✅ Simplified Choices
- Removed deprecated catalyst-specific reaction types
- Cleaner, more maintainable interface

## Compatibility Notes

### Breaking Changes
- ❌ `--fusion-variants` removed (use ML with rerank instead)
- ❌ Deprecated family names removed (C_N_Coupling_Cu → C_N_Coupling + --catalyst Cu)
- ❌ `fusion` strategy removed (use `ml` with `--rerank rule`)

### Migration Guide
```bash
# OLD (deprecated):
--family C_N_Coupling_Cu

# NEW (current):
--family C_N_Coupling --catalyst Cu

# OLD (deprecated):
--strategy fusion --fusion-variants 5

# NEW (current):
--strategy ml --rerank rule
```

## Testing Checklist
- [ ] Test interactive mode (prompts for all inputs)
- [ ] Test auto-detection with --rxn only
- [ ] Test catalyst selection (all 5 options)
- [ ] Test --strategy llm_synthesis with --llm-timeout
- [ ] Test --catalyst argument (None/Pd/Cu/Ni/other)
- [ ] Test --rerank rule vs analytics
- [ ] Verify API calls pass catalyst_preference correctly
- [ ] Confirm no errors when running all strategies

## Files Modified
1. `scripts/web_recommendation_cli.py` (541 lines, previously 518 lines)
   - Updated imports
   - Added/removed arguments
   - Updated main workflow
   - Replaced fusion with llm_synthesis
   - Added catalyst support to API calls

## Related Files (No Changes Needed)
- `scripts/local_recommendation_cli.py` - Already updated (reference implementation)
- `scripts/recommendation_cli_utils.py` - Already updated (shared utilities)
- `chemtools/router.py` - Already supports auto-detection

## Status
✅ **COMPLETE** - Web CLI now fully synced with local CLI
