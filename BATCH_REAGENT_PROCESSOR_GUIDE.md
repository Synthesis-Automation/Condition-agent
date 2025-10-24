# Batch Reagent Registry Processor - Quick Start Guide

## Overview

The `batch_add_reagents_to_registry.py` script automatically processes reagents from `reagents_mapped_to_registry_roles.md` and adds them to the reagent registry system. It handles:

- ✅ CAS number lookup via PubChem API
- ✅ Duplicate detection (skips existing reagents)
- ✅ Automatic role assignment
- ✅ Optional LLM enhancement for better metadata
- ✅ Dry-run mode for testing

## Prerequisites

1. **Reagent mapping file** must exist: `reagents_mapped_to_registry_roles.md`
2. **Registry directory** must exist: `data/reagents/`
3. **Internet connection** for PubChem CAS lookup
4. **(Optional)** LLM API keys for enhanced processing

## Basic Usage

### 1. Dry Run (Preview Only - Recommended First)

```powershell
python batch_add_reagents_to_registry.py --dry-run --max-reagents 10
```

This will:
- Process first 10 reagents
- Search for CAS numbers
- Show what would be added
- **NOT save** to registry

### 2. Add High-Frequency Reagents (Conservative)

```powershell
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5
```

This will:
- Only process reagents with ≥100 occurrences
- Add 1.5 second delay between API calls (rate limiting)
- Save to registry

### 3. Add All Reagents (Full Batch)

```powershell
python batch_add_reagents_to_registry.py --min-occurrences 1 --delay 1.0 --output batch_results.json
```

This will:
- Process all 552 reagents
- Save results to JSON file
- Take ~10 minutes (with rate limiting)

### 4. Skip Certain Roles

```powershell
python batch_add_reagents_to_registry.py --skip-roles solvent additive --min-occurrences 50
```

This will:
- Skip solvents and additives
- Only process bases, ligands, catalysts, condensation agents
- Minimum 50 occurrences

### 5. Use LLM for Enhanced Processing

```powershell
python batch_add_reagents_to_registry.py --use-llm --llm-provider aliyun --llm-model deepseek-r1-distill-qwen-7b --max-reagents 20
```

This will:
- Use LLM to classify families and assign properties
- Process first 20 reagents
- Requires LLM API keys configured

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--md-file` | `reagents_mapped_to_registry_roles.md` | Input markdown file |
| `--registry-dir` | `data/reagents/` | Registry directory |
| `--dry-run` | False | Preview only, don't save |
| `--max-reagents` | None (all) | Limit number to process |
| `--min-occurrences` | 10 | Minimum usage count |
| `--skip-roles` | None | Roles to skip (space-separated) |
| `--delay` | 1.0 | Seconds between API calls |
| `--use-llm` | False | Enable LLM workflow |
| `--llm-provider` | aliyun | LLM provider |
| `--llm-model` | deepseek-r1-distill-qwen-7b | LLM model |
| `--output` | None | Save results to JSON |

## Workflow

For each reagent in the markdown file:

1. **Check if already exists** (by name)
   - If yes → Skip

2. **Search for CAS number** via PubChem API
   - If not found → Skip
   - If found → Continue

3. **Check if exists** (by CAS)
   - If yes → Skip

4. **Add to registry**
   - Deterministic mode: Use default family, create basic entry
   - LLM mode: Use LLM to classify family and assign properties

5. **Save to registry file** (e.g., `data/reagents/condensation_agent.json`)

## Recommended Strategy

### Phase 1: High-Confidence Reagents (Dry Run)

```powershell
# Preview top 50 most-used reagents
python batch_add_reagents_to_registry.py --dry-run --min-occurrences 100 --max-reagents 50
```

### Phase 2: Add High-Confidence Reagents

```powershell
# Add reagents with ≥100 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5 --output phase2_results.json
```

**Expected:** ~50-100 reagents added

### Phase 3: Medium-Frequency Reagents

```powershell
# Add reagents with 20-99 occurrences
python batch_add_reagents_to_registry.py --min-occurrences 20 --delay 1.5 --output phase3_results.json
```

**Expected:** ~150-200 reagents added

### Phase 4: Low-Frequency Reagents (Optional)

```powershell
# Add all remaining reagents
python batch_add_reagents_to_registry.py --min-occurrences 1 --delay 2.0 --output phase4_results.json
```

**Expected:** ~200-300 reagents added

### Phase 5: Review and LLM Enhancement (Optional)

```powershell
# Use LLM to enhance specific roles
python batch_add_reagents_to_registry.py --use-llm --skip-roles solvent additive --min-occurrences 50 --max-reagents 100
```

## Output Format

The script prints progress for each reagent:

```
[1/552] Processing: DIC (role: condensation_agent, occurrences: 1718)
  Searching for CAS number...
  ✓ Found CAS: 693-13-0
  Adding with deterministic workflow...
  ✓ Added successfully (deterministic)

[2/552] Processing: K3PO4 (role: base, occurrences: 10296)
  ✓ Already exists in registry (by name), skipping
```

## Summary Statistics

After completion:

```
================================================================================
BATCH PROCESSING SUMMARY
================================================================================

Total processed:      100
Already exists:       25
CAS not found:        10
Added successfully:   60
Failed:               3
Skipped:              2

================================================================================
```

## Results JSON Structure

If `--output` is specified, saves detailed results:

```json
{
  "total_processed": 100,
  "already_exists": 25,
  "cas_not_found": 10,
  "added_successfully": 60,
  "failed": 3,
  "skipped": 2,
  "details": [
    {
      "name": "DIC",
      "role": "condensation_agent",
      "cas": "693-13-0",
      "family": "carbodiimides",
      "status": "added",
      "method": "deterministic"
    },
    {
      "name": "K3PO4",
      "role": "base",
      "status": "already_exists",
      "reason": "Found by name in registry"
    },
    ...
  ]
}
```

## Troubleshooting

### Issue: "CAS number not found"

**Cause:** Reagent name doesn't match PubChem database

**Solutions:**
- Check reagent name spelling in markdown file
- Manually look up CAS number and add via UI
- Try abbreviation instead of full name

### Issue: Rate limiting / 429 errors

**Cause:** Too many API requests to PubChem

**Solution:** Increase `--delay` parameter:
```powershell
python batch_add_reagents_to_registry.py --delay 2.0
```

### Issue: LLM errors

**Cause:** API keys not configured or quota exceeded

**Solutions:**
- Check environment variables (ALIYUN_API_KEY, OPENAI_API_KEY)
- Use deterministic mode (remove `--use-llm`)
- Switch provider/model

### Issue: Family not found

**Cause:** Role doesn't have default family in ROLE_CONFIG

**Solution:** 
- Check `app/reagent_taxonomy_ui.py` ROLE_CONFIG
- Ensure `default_family` is set for the role
- Or use LLM mode to infer family

## Files Modified

The script modifies these registry files:

- `data/reagents/additive.json`
- `data/reagents/base.json`
- `data/reagents/condensation_agent.json`
- `data/reagents/ligand.json`
- `data/reagents/metal_precursor.json`
- `data/reagents/preformed_metal_catalyst.json`
- `data/reagents/solvent.json`

**⚠️ Recommendation:** Commit changes to git before running to enable rollback if needed.

## Next Steps After Batch Processing

1. **Review added entries** in registry JSON files
2. **Verify families** are appropriate
3. **Add missing metadata:**
   - Refine abbreviations
   - Clean up synonyms
   - Add SMILES/InChI if missing
4. **Use reagent_taxonomy_ui.py** for manual corrections
5. **Run validation** to ensure schema compliance

## Example: Conservative Production Run

```powershell
# Step 1: Dry run first
python batch_add_reagents_to_registry.py --dry-run --min-occurrences 50 --max-reagents 10

# Step 2: Process in batches
python batch_add_reagents_to_registry.py --min-occurrences 100 --delay 1.5 --output results_100plus.json
python batch_add_reagents_to_registry.py --min-occurrences 50 --delay 1.5 --output results_50plus.json
python batch_add_reagents_to_registry.py --min-occurrences 20 --delay 2.0 --output results_20plus.json

# Step 3: Review results
cat results_100plus.json | jq '.added_successfully'
```

---

## Support

For issues or questions:
1. Check the summary statistics for error types
2. Review `--output` JSON for detailed error messages
3. Use `--dry-run` to preview before committing changes
4. Manually add problematic reagents via `reagent_taxonomy_ui.py`
