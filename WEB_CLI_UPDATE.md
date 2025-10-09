# Web CLI Updated to Match Local CLI ✅

## Question: Does web version work the same as local version?

**Before**: ❌ NO - The web version was missing several key features  
**After**: ✅ YES - Both versions now have feature parity

## Changes Made to `scripts/web_recommendation_cli.py`

### 1. Added New Parameters to `call_ml_recommendation()` ✅

**Before**:
```python
def call_ml_recommendation(
    base_url: str,
    reaction: str,
    reaction_type: Optional[str],
    k: int,
    limit: int,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "reaction_type": reaction_type or None,
        "k": k,
        "limit": limit,
        "relax": {},
        "constraints": {},
    }
```

**After**:
```python
def call_ml_recommendation(
    base_url: str,
    reaction: str,
    reaction_type: Optional[str],
    k: int,
    limit: int,
    rerank_strategy: str = 'rule',           # NEW
    filter_unknown_reagents: bool = False,   # NEW
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "reaction": reaction,
        "reaction_type": reaction_type or None,
        "k": k,
        "limit": limit,
        "relax": {},
        "constraints": {},
        "rerank_strategy": rerank_strategy,        # NEW
        "filter_unknown_reagents": filter_unknown_reagents,  # NEW
    }
```

### 2. Added Deprecation Warning to `call_fusion_recommendation()` ✅

```python
def call_fusion_recommendation(...):
    """
    Execute the fusion recommendation endpoint.
    
    NOTE: The /api/v1/recommend/fusion endpoint is DEPRECATED.
    Use /api/v1/recommend/conditions with rerank_strategy='rule' instead.
    """
    import warnings
    warnings.warn(
        "Fusion endpoint is deprecated. Use call_ml_recommendation() "
        "with rerank_strategy='rule' instead.",
        DeprecationWarning,
        stacklevel=2
    )
    ...
```

### 3. Added Command-Line Arguments ✅

**Before**: Always interactive (no CLI args)  
**After**: Full argument support matching local version

```bash
# New command-line arguments:
--url, --base-url         # Server URL (default: http://localhost:8000)
--rxn, --reaction         # Reaction SMILES
--family, --type          # Reaction family
--k                       # Number of precedents
--limit                   # Number of recommendations
--fusion-variants         # Fusion variants count
--save-dir               # Output directory
--strategy               # Which strategies to run (all/rule/ml/fusion)
--rerank                 # Reranking strategy (none/rule/analytics)
--filter-unknown         # Filter unknown reagents
```

### 4. Added Selective Strategy Execution ✅

**Before**: Always runs all 3 strategies (rule, ML, fusion)  
**After**: Can run specific strategies with `--strategy`

```python
# Run selected strategies
run_rule = args.strategy in ["all", "rule"]
run_ml = args.strategy in ["all", "ml"]
run_fusion = args.strategy in ["all", "fusion"]

if run_rule:
    rule_result = call_rule_based(...)
if run_ml:
    ml_result = call_ml_recommendation(...)
if run_fusion:
    fusion_result = call_fusion_recommendation(...)
```

### 5. Added Custom Output Directory ✅

**Before**: Hardcoded `results/` directory  
**After**: Configurable with `--save-dir`

```python
# Update output directory from args
output_dir = Path(args.save_dir)
output_dir.mkdir(parents=True, exist_ok=True)

def save_to_dir(data: Dict[str, Any], filename: str) -> Path:
    output_path = output_dir / filename
    ...
```

## Feature Parity Comparison

| Feature | Local CLI | Web CLI (Before) | Web CLI (After) |
|---------|-----------|------------------|-----------------|
| **CLI Arguments** | ✅ | ❌ | ✅ |
| **`--rxn` argument** | ✅ | ❌ | ✅ |
| **`--family` argument** | ✅ | ❌ | ✅ |
| **`--k` argument** | ✅ | ❌ | ✅ |
| **`--limit` argument** | ✅ | ❌ | ✅ |
| **`--strategy` selection** | ✅ | ❌ | ✅ |
| **`--rerank` parameter** | ✅ | ❌ | ✅ |
| **`--filter-unknown` parameter** | ✅ | ❌ | ✅ |
| **`--save-dir` parameter** | ✅ | ❌ | ✅ |
| **Deprecation warnings** | ✅ | ❌ | ✅ |
| **Selective execution** | ✅ | ❌ | ✅ |
| **Interactive mode** | ✅ | ✅ | ✅ |

## Usage Examples

### Interactive Mode (Same as Before)
```bash
python scripts/web_recommendation_cli.py
# Prompts for reaction SMILES and family
```

### Command-Line Mode (NEW!)
```bash
# Basic usage with reaction and family
python scripts/web_recommendation_cli.py \
  --rxn "C/C=C/Br.C=CB(O)O>>C/C=C/C=C" \
  --family Suzuki

# With reranking and filtering
python scripts/web_recommendation_cli.py \
  --rxn "C/C=C/Br.C=CB(O)O>>C/C=C/C=C" \
  --family Suzuki \
  --rerank rule \
  --filter-unknown

# Run only ML strategy
python scripts/web_recommendation_cli.py \
  --rxn "C/C=C/Br.C=CB(O)O>>C/C=C/C=C" \
  --strategy ml \
  --k 100 \
  --limit 5

# Custom server and output directory
python scripts/web_recommendation_cli.py \
  --url http://production-server:8080 \
  --rxn "..." \
  --save-dir my_results
```

## Benefits of Updates

### 1. Consistency ✅
- Both CLIs now have identical functionality
- Same command-line arguments
- Same behavior patterns

### 2. Automation ✅
- Can script web API testing
- No need for interactive input
- Easier CI/CD integration

### 3. New Features ✅
- **Reranking strategies**: Can test rule vs analytics vs none
- **Filter unknown reagents**: Test reagent database coverage
- **Selective execution**: Test individual strategies faster
- **Custom output**: Organize results better

### 4. Better Testing ✅
- Test server API with same options as local
- Compare web vs local outputs
- Verify catalyst specificity improvement via API

## Backward Compatibility

✅ **Fully backward compatible**

- Interactive mode still works (no arguments = prompts)
- All old functionality preserved
- New features are opt-in via arguments
- No breaking changes

## Testing

**Test both modes work**:

```bash
# Test interactive mode (prompts)
python scripts/web_recommendation_cli.py

# Test command-line mode with new features
python scripts/web_recommendation_cli.py \
  --rxn "C/C=C/Br.C=CB(O)O>>C/C=C/C=C" \
  --family Suzuki \
  --strategy ml \
  --rerank rule \
  --limit 2
```

**Expected output**: Should show specific catalyst complexes (Pd(dppf)Cl2) instead of generic "Palladium", matching the local CLI improvement!

## Files Modified

1. ✅ `scripts/web_recommendation_cli.py`
   - Updated `call_ml_recommendation()` - added rerank/filter parameters
   - Updated `call_fusion_recommendation()` - added deprecation warning
   - Updated `main()` - added full argparse support
   - Added selective strategy execution
   - Added custom output directory support

## Conclusion

**YES, the web version now works the same as the local version!** ✅

Both CLIs now:
- ✅ Support the same command-line arguments
- ✅ Have the same reranking options
- ✅ Can filter unknown reagents
- ✅ Support selective strategy execution
- ✅ Show deprecation warnings for fusion
- ✅ Provide specific catalyst information (via the updated backend)

The web CLI is now a complete HTTP-based equivalent of the local CLI, making it easy to test the API server with the same features available for local testing.
