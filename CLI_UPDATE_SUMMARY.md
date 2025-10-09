# CLI Update Summary - local_recommendation_cli.py

## Changes Made

Updated `scripts/local_recommendation_cli.py` to support command-line arguments while maintaining backward compatibility with interactive mode.

## Key Features

### 1. **Smart Input Mode**
- **Non-interactive**: When `--rxn` and `--family` are provided, runs without prompting
- **Interactive**: When arguments are missing, prompts user for input
- **Hybrid**: Can provide some arguments (e.g., reaction) and be prompted for others (e.g., type)

### 2. **New Command-Line Arguments**

```bash
--rxn, --reaction       # Reaction SMILES (required if non-interactive)
--family, --type        # Reaction type (optional - will auto-detect or prompt)
--k                     # Number of precedents for ML (default: 50)
--limit                 # Number of ML recommendations (default: 5)
--fusion-variants       # Number of fusion variants (default: 5)
--save-dir             # Output directory (default: "results")
--strategy             # Which to run: all/rule/ml/fusion (default: all)
--pretty               # Pretty-print JSON (flag)
```

### 3. **Flexible Strategy Selection**

Can now run only specific recommendation strategies:
- `--strategy rule` - Only rule-based matching
- `--strategy ml` - Only ML precedent search
- `--strategy fusion` - Only fusion recommendation (fastest)
- `--strategy all` - All three (default)

### 4. **Custom Output Directory**

Results can be saved to any directory via `--save-dir`:
```bash
python scripts/local_recommendation_cli.py \
  --rxn "..." \
  --save-dir my_experiment_1
```

## Usage Examples

### Example 1: Fully automated (CI/CD, scripts)
```bash
python scripts/local_recommendation_cli.py \
  --rxn "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1" \
  --family Ullmann_CN \
  --strategy fusion \
  --save-dir automated_tests
```
✅ No prompts, runs directly

### Example 2: Quick test with auto-detection
```bash
python scripts/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccccc1-c1ccccc1" \
  --strategy fusion
```
✅ Auto-detects as Suzuki, no type prompt

### Example 3: Interactive exploration
```bash
python scripts/local_recommendation_cli.py
```
✅ Prompts for reaction SMILES and type selection

### Example 4: Custom parameters
```bash
python scripts/local_recommendation_cli.py \
  --rxn "..." \
  --family Buchwald_CN \
  --k 100 \
  --fusion-variants 10 \
  --save-dir high_diversity_test
```
✅ Uses 100 precedents, generates 10 variants

## Backward Compatibility

✅ **Fully backward compatible**
- Old usage: `python scripts/local_recommendation_cli.py` still works
- Still prompts interactively when no arguments provided
- All previous functionality preserved

## Benefits

1. **Automation-friendly**: Can now be used in batch scripts, CI/CD pipelines
2. **Faster workflow**: Skip interactive prompts when you know what you want
3. **Flexible**: Mix command-line args with interactive prompts as needed
4. **Better organization**: Custom output directories for different experiments
5. **Efficient**: Run only the strategies you need

## Testing Results

### Test 1: Ullmann C-N Coupling (Non-interactive)
```bash
python scripts/local_recommendation_cli.py \
  --rxn "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1" \
  --family Ullmann_CN \
  --strategy fusion \
  --save-dir test_output
```

**Result**: ✅ 
- No prompts shown
- Correctly identified as Ullmann reaction
- Cu-based catalyst recommended (Copper(I) iodide)
- Saved to `test_output/` directory

### Test 2: Suzuki Coupling with Custom Parameters
```bash
python scripts/local_recommendation_cli.py \
  --rxn "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccccc1-c1ccccc1" \
  --family Suzuki \
  --strategy fusion \
  --k 30 \
  --save-dir suzuki_test
```

**Result**: ✅
- No prompts shown
- Used 30 precedents (custom k value)
- Pd-based catalyst recommended (Tetrakis(triphenylphosphine)palladium(0))
- Saved to `suzuki_test/` directory

### Test 3: Interactive Mode
```bash
python scripts/local_recommendation_cli.py --strategy fusion
```

**Result**: ✅
- Prompted for reaction SMILES
- Prompted for reaction type selection
- Ran fusion strategy only
- Saved to default `results/` directory

## Code Changes

1. Added `import argparse` 
2. Created comprehensive argument parser with help text and examples
3. Modified `main()` to check for arguments before prompting
4. Added conditional strategy execution
5. Implemented custom output directory handling
6. Maintained all existing functionality

## Files Modified

- `scripts/local_recommendation_cli.py` - Main CLI script (added argument parsing)

## Files Created

- `TEST_CLI_MODES.md` - Documentation of CLI modes
- `test_ullmann_fix.py` - Test script for verifying Ullmann reaction fix
