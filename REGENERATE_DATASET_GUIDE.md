# How to Regenerate Dataset with DRFP Precomputation

## Why Regenerate?

Your current Buchwald dataset was created before DRFP precomputation was implemented. Regenerating it will:
- Add precomputed DRFP fingerprints to each reaction
- Make precedent searches **18x faster** (from ~9 seconds to <0.5 seconds)
- Eliminate the slow startup time you're experiencing

## Steps to Regenerate

### 1. Open the RDF Processor
```powershell
python data-processor/Scifinder_rdf_processer.py
```

### 2. In the GUI:
- Click "Select Folder" 
- Navigate to your RDF files folder containing Buchwald reactions
- Click "Process Folder"

### 3. Watch for Preprocessing Statistics
At the end of processing, you should see:
```
============================================================
PREPROCESSING STATISTICS
============================================================
Total reactions:           1343
Successfully preprocessed: 1343 (100.0%)
Failed:                    0 (0.0%)
Skipped:                   0 (0.0%)
DRFP fingerprints:         1343 (100.0%)  <-- IMPORTANT!

Family Detection Sources:
  SciFinder exact match:   1200 (89.4%)
  SciFinder partial match: 100 (7.4%)
  SciFinder unmapped:      43 (3.2%)
  Unknown/Missing:         0 (0.0%)
============================================================
```

**Key thing to check**: `DRFP fingerprints: 1343 (100.0%)`
- If you see 100%, all reactions have precomputed DRFP ✓
- If you see <100%, some reactions failed DRFP computation

### 4. Output Location
The new JSONL file will be in:
```
data/reaction_dataset/C_N_Coupling_Pd.jsonl
```

## Verification

After regenerating, test the performance:

```powershell
python test_new_drfp_code.py
```

Expected output:
```
[1] Testing with DRFP enabled...
   ✓ Completed in 0.4-0.5s     <-- Should be fast!
   ✓ Found 5 precedents

[2] Testing second search (cached)...
   ✓ Completed in 0.2-0.3s
   ✓ Found 5 precedents
   ✓ Speedup: 1.5x faster

RESULT:
✓ PASS: Precedent search is fast (<2s)
  First search: 0.45s
  Second search: 0.28s
```

## What Changed in the Dataset?

### Before (Old Dataset):
```json
{
  "reaction_id": "...",
  "smiles": {"reactants": "...", "products": "..."},
  "precomputed": {
    "reaction_smiles": "...",
    "features": {...}
  }
}
```

### After (New Dataset):
```json
{
  "reaction_id": "...",
  "smiles": {"reactants": "...", "products": "..."},
  "precomputed": {
    "reaction_smiles": "...",
    "features": {...},
    "drfp_fp": [0, 1, 0, 1, ...],     <-- NEW! 4096 elements
    "drfp_n_bits": 4096,              <-- NEW!
    "drfp_radius": 3                  <-- NEW!
  }
}
```

## File Size Increase

- **Old dataset**: ~2-3 MB
- **New dataset**: ~15-18 MB
- **Increase**: ~5-6x larger, but **18x faster** searches!

## Troubleshooting

### If DRFP fingerprints: 0 (0.0%)
This means DRFP package isn't installed. Install it:
```powershell
pip install drfp
```

### If only some reactions have DRFP (e.g., 90%)
Some reactions have invalid SMILES that DRFP can't process. This is OK - the system will compute DRFP on-demand for those reactions.

### If processing seems slow
DRFP computation adds ~8-10 seconds to dataset creation for 1343 reactions. This is a one-time cost that saves 9 seconds on every future precedent search!

## Performance Expectations

| Operation | Old Dataset | New Dataset | Improvement |
|-----------|-------------|-------------|-------------|
| Dataset creation | 5-10s | 15-20s | Slower (one-time) |
| First precedent search | ~9s | ~0.5s | **18x faster** |
| Subsequent searches | ~9s | ~0.3s | **30x faster** |
| Total for 10 searches | ~90s | ~3s | **30x faster** |

**Amortization Point**: After just **2 precedent searches**, the new dataset has already paid for itself!

## Next Steps

After regenerating the dataset:
1. ✓ Verify fast precedent searches with `test_new_drfp_code.py`
2. ✓ Test the full UI with `python app/ui_simple.py`
3. ✓ Regenerate other datasets (Suzuki, Ullmann, etc.) if you have them

## Documentation

See `DRFP_PRECOMPUTATION_OPTIMIZATION.md` for technical details on:
- Performance profiling results
- Implementation details
- Benchmark comparisons
- Storage trade-offs
