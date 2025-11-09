# CLI Update: Default top_k Changed to 1

## Change Summary

**Date**: November 9, 2025  
**Changed**: Default `top_k` value from **5** to **1**

## Modified Files

- `app/unified_rule_protocol_interactive_cli.py`
  - Line ~518: `initial_k: int = 1` (was 5)
  - Line ~806: `default=1` in argparse (was 5)
  - Line ~808: Help text updated to "default: 1" (was "default: 5")

## Behavior

### Before (top_k = 5)
```
Found 5 recommendation(s)
📖 [1] Suzuki-Miyaura Coupling (0.157)
📖 [2] Buchwald-Hartwig C-N Coupling (0.076)
📖 [3] Ullmann-Goldberg Cu C-N Coupling (0.056)
📋 [4] Protocol X (0.045)
📋 [5] Protocol Y (0.032)
```

### After (top_k = 1)
```
Found 1 recommendation(s)
📖 [1] Suzuki-Miyaura Coupling (0.157)
```

## Usage Examples

### Default (top 1 result)
```bash
python app/unified_rule_protocol_interactive_cli.py
reaction> Brc1ccccc1.B(O)Oc1ccccc1>>c1ccc(-c2ccccc2)cc1
# Shows only top 1 match
```

### Show more results
```bash
reaction> /k 3
Set top_k = 3
reaction> <reaction_smiles>
# Shows top 3 matches
```

### Filter by type (each shows top 1)
```bash
reaction> /type rule
Set source_type = rule
reaction> <reaction_smiles>
# Shows top 1 rule

reaction> /type protocol
Set source_type = protocol
reaction> <reaction_smiles>
# Shows top 1 protocol
```

### Split mode (top 1 rule + top 1 protocol)
```bash
reaction> /split on
Split view enabled
reaction> <reaction_smiles>
━━━ TOP RULE ━━━
📖 [1] Rule name...

━━━ TOP PROTOCOL ━━━
📋 [1] Protocol name...
```

## Rationale

### Why Change to top_k = 1?

1. **Focus on best match**: Most users want the single best recommendation
2. **Reduced noise**: Fewer irrelevant results (e.g., Buchwald-Hartwig false positives)
3. **Cleaner output**: Easier to read and process
4. **Automation-friendly**: When using automation format, typically only need the top match
5. **Easy to increase**: Users can still use `/k 5` if they want more results

### Benefits

✅ **Better UX**: Top result is usually sufficient  
✅ **Less confusion**: Avoids showing marginally relevant rules  
✅ **Faster**: Less output to process  
✅ **Flexible**: Can increase with `/k` command anytime  

### When to Use More Results

Use `/k 3` or `/k 5` when:
- Exploring different reaction options
- Comparing multiple approaches
- Low similarity scores (may need alternatives)
- Research/analysis (not automation)

## Testing

Run tests:
```bash
# Test default behavior
python test_top_k_default.py

# Demo with different filters
python demo_top_k_1.py
```

## Backward Compatibility

✅ **CLI commands unchanged**: `/k` command still works  
✅ **API unchanged**: `recommender.recommend(top_k=5)` still works  
✅ **Only default changed**: Users can override anytime  

## Documentation Updates

Updated files:
- `app/unified_rule_protocol_interactive_cli.py` (code + help text)
- `test_top_k_default.py` (test script)
- `demo_top_k_1.py` (demo script)
- `CLI_TOP_K_CHANGE.md` (this document)
