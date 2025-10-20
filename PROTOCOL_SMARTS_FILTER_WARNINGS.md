# Protocol SMARTS Filter Warnings - Implementation Summary

## Problem

The protocol recommender was silently filtering out protocols based on SMARTS pattern matching without informing the user, making it unclear why certain protocols were or weren't appearing in results.

## Solution

Added **explicit warnings** when SMARTS filtering removes candidates, helping users understand:
1. How many protocols were filtered out
2. Why they might not be seeing expected results
3. How to bypass SMARTS filtering if needed

## Implementation

### Warning System

The system now tracks and reports:
- Number of protocols before SMARTS filtering
- Number of protocols after SMARTS filtering  
- Difference (how many were removed)
- User-friendly warning message with suggestions

### Warning Messages

#### When Some Protocols Are Filtered

```
⚠️  SMARTS FILTER WARNING:
    SMARTS filtering removed 15 protocol(s) that did not match 
    the reaction structure. 1 of 16 protocols remain. 
    Use --no-smarts-filter to see all protocols ranked by similarity.
```

#### When All Protocols Are Filtered

```
⚠️  SMARTS FILTER WARNING:
    No protocols found matching the reaction SMARTS pattern. 
    All 16 candidate(s) were filtered out. 
    Try --no-smarts-filter to see protocols ranked by chemical 
    similarity only.
```

## Example: Alkyl-Aryl Suzuki Coupling

**Query:** `CCBr.c1ccccc1B(O)O>>CCc1ccccc1`
- Alkyl bromide (`CCBr`) + phenylboronic acid
- Most Suzuki protocols are for aryl-aryl couplings

### With SMARTS Filtering (Default)

```bash
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 3
```

**Result:**
- ⚠️ Warning shown: 15 of 16 protocols filtered out
- Only 1 protocol shown: Cu-catalyzed alkyl-aryl Suzuki (0.288 similarity)
- This is the ONLY protocol with SMARTS matching alkyl halides

### Without SMARTS Filtering

```bash
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 3 --no-smarts-filter
```

**Result:**
- No warning (filter disabled)
- 3 protocols shown ranked by DRFP similarity:
  1. Pd-catalyzed Suzuki (0.418 similarity) - **Better match by DRFP!**
  2. Ni-catalyzed Suzuki (0.295 similarity)
  3. Cu-catalyzed Suzuki (0.288 similarity)

## Benefits

### 1. **Transparency**
Users now understand why results may be limited or unexpected.

### 2. **Actionable Guidance**
Clear suggestion to use `--no-smarts-filter` for broader results.

### 3. **Better Decision Making**
Users can choose between:
- **Structural matching** (SMARTS) - precise but restrictive
- **Chemical similarity** (DRFP) - flexible but may include different reactions

### 4. **Debugging Aid**
Helps protocol database maintainers identify:
- Overly restrictive SMARTS patterns
- Missing protocol coverage
- SMARTS pattern errors

## Technical Details

### Code Changes

**`chemtools/protocol/recommend.py`:**

1. Track filtering statistics:
```python
num_before_smarts = len(candidates)
if use_smarts_filter:
    candidates = self._filter_by_smarts(reaction_smiles, candidates)
    num_after_smarts = len(candidates)
```

2. Generate warning messages:
```python
if num_after_smarts < num_before_smarts:
    logger.warning(f"SMARTS filtering removed {removed} protocol(s)")
    smarts_filter_warning = "SMARTS filtering removed..."
```

3. Include in response:
```python
'extras': {
    'smarts_filter_warning': smarts_filter_warning,
    'num_before_smarts_filter': num_before_smarts,
    'num_after_smarts_filter': num_after_smarts
}
```

**`chemtools/protocol/recommend_cli.py`:**

Display warnings prominently:
```python
if extras.get('smarts_filter_warning'):
    print("⚠️  SMARTS FILTER WARNING:")
    print(f"    {extras['smarts_filter_warning']}")
```

### Logging Levels

- **WARNING**: Logged when protocols are filtered out
- **INFO**: Normal operation messages
- **DEBUG**: Detailed filtering information

## Usage Examples

### Check if SMARTS is Limiting Results

```bash
# Try with SMARTS filtering (default)
python -m chemtools.protocol.recommend_cli "your-reaction" --k 5

# If warning appears, try without SMARTS
python -m chemtools.protocol.recommend_cli "your-reaction" --k 5 --no-smarts-filter
```

### Compare Structural vs Chemical Similarity

```bash
# Structural match (SMARTS)
python -m chemtools.protocol.recommend_cli "CCBr.ArB>>CCAr"

# Chemical similarity (DRFP only)
python -m chemtools.protocol.recommend_cli "CCBr.ArB>>CCAr" --no-smarts-filter
```

## When to Use Each Mode

### Use SMARTS Filtering (Default) When:
- ✅ You want protocols for the **exact reaction type**
- ✅ Structural features are critical (sp2 vs sp3, halide type, etc.)
- ✅ You prefer fewer, more precise matches

### Use `--no-smarts-filter` When:
- ✅ SMARTS filtering gives **too few results**
- ✅ You want to explore **similar reaction types**
- ✅ Protocol SMARTS patterns may be too restrictive
- ✅ You need inspiration from **chemically similar** reactions

## Future Enhancements

Potential improvements:
- [ ] Show which specific SMARTS patterns caused filtering
- [ ] Add "partial match" mode (fuzzy SMARTS matching)
- [ ] Display sample SMARTS patterns in warning messages
- [ ] Track and report common filtering patterns
- [ ] Auto-suggest alternative search terms
- [ ] Explain SMARTS mismatch (e.g., "sp3 vs sp2 halide")

## Related Documentation

- `chemtools/protocol/readme.md` - Protocol CLI usage guide
- `PROTOCOL_DRFP_OPTIMIZATION_SUMMARY.md` - DRFP storage details
- Protocol database SMARTS best practices in protocol README
