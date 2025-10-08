# Demo Files - Quick Reference

**Location**: `tests/` directory  
**Date**: January 2025

---

## Available Demos

### 1. Basic Tools Demo
**File**: `tests/demo_basic_tools.py`

**Run from project root**:
```bash
python tests/demo_basic_tools.py
```

**Run from tests directory**:
```bash
cd tests
python demo_basic_tools.py
```

**Features Demonstrated**:
- ✅ SMILES normalization (canonicalization)
- ✅ Reaction family detection (with/without catalysts)
- ✅ Molecular featurization (leaving groups, nucleophiles)
- ✅ Property lookup (full reagent database)
- ✅ DRFP similarity calculations

**Duration**: ~5 seconds

---

### 2. Recommendations Demo
**File**: `tests/demo_recommendations.py`

**Run from project root**:
```bash
python tests/demo_recommendations.py
```

**Run from tests directory**:
```bash
cd tests
python demo_recommendations.py
```

**Features Demonstrated**:
- ✅ Precedent search (k-NN with DRFP)
- ✅ Condition recommendations (NEW modular API)
- ✅ ML-enhanced recommendations (hybrid approach)
- ✅ Plate design for high-throughput screening
- ✅ Structured condition outputs

**Duration**: ~10-15 seconds

---

## Quick Test Commands

### Run Both Demos
```powershell
# PowerShell
python tests/demo_basic_tools.py
python tests/demo_recommendations.py
```

```bash
# Bash/Linux
python tests/demo_basic_tools.py && python tests/demo_recommendations.py
```

---

### Run from Tests Directory
```powershell
# PowerShell
cd tests
python demo_basic_tools.py
python demo_recommendations.py
```

---

## Key Features

### Path Handling ✅
Both demos now include automatic path setup:
```python
import sys
from pathlib import Path

# Add parent directory to path for chemtools import
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))
```

This allows the demos to work from:
- ✅ Project root: `python tests/demo_basic_tools.py`
- ✅ Tests directory: `cd tests && python demo_basic_tools.py`
- ✅ Any subdirectory: Works automatically

---

## What Changed

### Before (Root Directory)
```
Condition-agent/
  ├── demo_basic_tools.py        ❌ Old location
  ├── demo_recommendations.py    ❌ Old location
  └── chemtools/
```

### After (Tests Directory)
```
Condition-agent/
  ├── tests/
  │   ├── demo_basic_tools.py    ✅ New location
  │   └── demo_recommendations.py ✅ New location
  └── chemtools/
```

**Benefits**:
- ✅ Better organization (demos with tests)
- ✅ Cleaner root directory
- ✅ Follows standard Python project structure

---

## Expected Output

### Basic Tools Demo
```
======================================================================
  ChemTools Basic Tools Demo
  October 2025 - Core Utilities
======================================================================

✅ 1. SMILES Normalization
   → Normalized: Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1

✅ 2. Reaction Family Detection
   → Family: Ullmann_CN
   → Confidence: 0.90

✅ 3. Molecular Featurization
   → LG: Br
   → Nucleophile class: amine_primary

✅ 4. Property Lookup (Full Reagent Database)
   → DMF found: N,N-Dimethylformamide
   → BINAP found: BINAP (CAS: 2250-01-3)

✅ 5. DRFP Similarity
   → Similarity: 0.375 (Br vs Cl)

======================================================================
  ✅ Basic Tools Demo Complete!
======================================================================

📖 Next: python tests/demo_recommendations.py
🎨 Or try: python app/ui_gradio.py
```

### Recommendations Demo
```
======================================================================
  ChemTools Recommendations Demo
  October 2025 - Condition Recommendations & ML
======================================================================

✅ 1. Precedent Search (k-NN)
   → Found: 5 precedents
   → Top condition: Pd(OAc)2, XPhos, K3PO4, Toluene

✅ 2. Condition Recommendations (NEW API)
   → Family: C_N_Coupling_Pd
   → Conditions: 3 recommendations

✅ 3. Structured Condition Outputs
   → Structured format with full details

✅ 4. ML-Enhanced Recommendations
   → Hybrid strategy: rule-based + ML
   → Top predictions with confidence scores

✅ 5. Plate Design (High-Throughput Screening)
   → Plate layout: 96 wells
   → 5 cores × 3 variants each

======================================================================
  ✅ Recommendations Demo Complete!
======================================================================

📖 See also: python tests/demo_basic_tools.py
🎨 Try interactive UI: python app/ui_gradio.py
```

---

## Troubleshooting

### Import Errors
If you see:
```
ModuleNotFoundError: No module named 'chemtools'
```

**Solution**: Make sure you're running from the project root:
```bash
cd /path/to/Condition-agent
python tests/demo_basic_tools.py
```

Or ensure the parent directory path setup is present in the demo file.

---

### DRFP Not Installed
If you see:
```
⚠️  DRFP not installed
```

**Solution** (optional):
```bash
pip install drfp
```

DRFP is optional. The demos will skip similarity calculations if not installed.

---

### ML Models Not Found
If you see:
```
⚠️  ML models not available
```

This is normal. ML features are optional and require:
- Model files in `models/` directory
- ML dependencies: `pip install -r requirements-ml.txt`

The demos will skip ML sections if models aren't available.

---

## Related Files

### Other Demo Files
- `demo_get_all_reagents.py` - Reagent database exploration
- `test_get_all_reagents.py` - Comprehensive reagent tests
- `demo_chemtools_complete.py` - Legacy comprehensive demo
- `demo_chemtools_quick.py` - Legacy quick demo

### Documentation
- `DEMO_SPLIT_GUIDE.md` - Demo organization guide
- `DEMO_QUICK_REF.md` - Quick reference card
- `GET_ALL_REAGENTS_FEATURE.md` - Reagent lookup feature guide
- `PROPERTIES_REMOVAL_COMPLETE.md` - Properties.py removal summary

---

## Integration with Other Tools

### Run with Gradio UI
```bash
python app/ui_gradio.py
```

### Run with Simple UI
```bash
python app/ui_simple.py
```

### Run Tests
```bash
pytest tests/ -v
```

### Run Registry CLI
```bash
python -m chemtools.scdb_matcher.cli --query "Toluene"
```

---

## Tips

### Quick Demo Run
```bash
# Just see basic tools
python tests/demo_basic_tools.py

# Just see recommendations
python tests/demo_recommendations.py

# Both in sequence
python tests/demo_basic_tools.py && python tests/demo_recommendations.py
```

### Save Output to File
```powershell
# PowerShell
python tests/demo_basic_tools.py | Tee-Object -FilePath demo_output.txt
```

```bash
# Bash
python tests/demo_basic_tools.py | tee demo_output.txt
```

### Run Specific Section
Both demos are structured with individual test functions. You can modify `main()` to run only specific tests:

```python
# In demo_basic_tools.py
tests = [
    test_1_smiles_normalization,  # Only run this one
    # test_2_family_detection,    # Comment out others
    # test_3_molecular_featurization,
    # ...
]
```

---

## Summary

✅ **Demos moved to `tests/` directory**  
✅ **Automatic path setup for flexible execution**  
✅ **Works from project root or tests directory**  
✅ **Clean, organized project structure**  

**Run now**:
```bash
python tests/demo_basic_tools.py
python tests/demo_recommendations.py
```

🎉 **All working!**
