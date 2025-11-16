# Z-Score Based Ranking Implementation

## Summary

The HTE condition recommendation system now uses **z-score as the primary ranking metric** for conditions. The z-score measures how successful a condition is relative to all experiments in the database.

## Changes Made

### 1. Updated Data Model (`chemtools/HTE/recommender.py`)

Added `avg_z_score` field to `ConditionRecommendation`:

```python
@dataclass
class ConditionRecommendation:
    """Single condition recommendation with metadata
    
    Recommendations are ranked primarily by avg_z_score, which measures
    the success of a condition relative to all experiments in the database.
    """
    # Statistics
    success_rate: float = 0.0
    avg_yield: float = 0.0
    median_yield: float = 0.0
    num_experiments: int = 0
    avg_z_score: float = 0.0  # Average z-score (PRIMARY ranking metric)
    confidence_score: float = 0.0  # Secondary score
```

### 2. Updated Confidence Score Calculation

The confidence score now prioritizes z-score:

```python
def _calculate_confidence_score(
    self,
    avg_z_score: float,
    num_experiments: int,
    avg_yield: float
) -> float:
    """
    Calculate confidence score combining multiple factors.
    
    Formula: weighted combination of z-score (primary), sample size, and avg yield
    Z-score is the primary metric as it measures the success of a condition.
    """
    # Normalize z-score: typical range is -3 to +3, scale to 0-1
    z_score_normalized = max(0.0, min(1.0, (avg_z_score + 3.0) / 6.0))
    
    # Combined score with weights
    confidence = (
        0.6 * z_score_normalized +  # 60% from z-score (primary)
        0.25 * sample_score +        # 25% from sample size
        0.15 * yield_weight          # 15% from avg yield
    ) * 100.0
```

**Previous formula:**
- 50% success rate
- 30% sample size
- 20% avg yield

**New formula:**
- 60% z-score (primary)
- 25% sample size
- 15% avg yield

### 3. Updated Ranking Logic

Conditions are now sorted primarily by average z-score:

```python
# Sort by average z-score (primary metric), then confidence score
recommendations.sort(key=lambda x: (x.avg_z_score, x.confidence_score), reverse=True)
```

### 4. Updated Output Formatting

#### Terminal Output
```
================================================================================
Recommendation #1
================================================================================
⭐ Avg Z-Score: 2.61 (Primary Ranking Metric)
Confidence Score: 63.8/100
Success Rate: 0.0% (1 experiments)
Avg Yield: 49.8% | Median: 49.8%

🧪 CONDITIONS:
  Catalyst: CuI
  Ligand: PPBO
  Base: NaOtBu
  Solvent: tAmOH

📊 STATISTICS:
  Reaction Type: C_N_Coupling
  Z-Score: Avg=2.61, Range=[2.61, 2.61]
  Reactant Types: ArBr + ArNH2
```

#### Compact Format
```
TOP RECOMMENDATION (Z-Score: 2.61, Confidence: 63.8/100)
  Catalyst: CuI
  Ligand: PPBO
  Base: NaOtBu
  Solvent: tAmOH
```

#### JSON Format
```json
{
  "recommendations": [
    {
      "rank": 1,
      "avg_z_score": 2.61,
      "confidence_score": 63.8,
      "success_rate": 0.0,
      "avg_yield": 49.8,
      ...
    }
  ]
}
```

#### Demo Script Output
```
Rank Z-Score Catalyst Ligand    Base    Solvent Experiments Avg_Yield_% Confidence
1    2.61    CuI      PPBO      NaOtBu  tAmOH   1           49.8        63.8
2    2.59    CuI      DMPAO     NaOtBu  tAmOH   1           49.4        63.6
3    2.21    CuI      2-Isobu.. NaOtBu  tAmOH   1           43.9        58.9
```

## Usage Examples

### CLI with Z-Score Ranking

```bash
# Get top conditions ranked by z-score
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --catalyst copper -k 10 --min-exp 1

# Compact output showing z-score
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --catalyst copper --compact

# JSON output with z-score
python -m chemtools.HTE.cli -a "c1ccc(Br)cc1" -b "c1ccc(N)cc1" \
    --catalyst copper --json -o output.json
```

### Python API

```python
from chemtools.HTE import HTERecommender

recommender = HTERecommender()
result = recommender.recommend(
    reactant_a_smiles="c1ccc(Br)cc1",
    reactant_b_smiles="c1ccc(N)cc1",
    catalyst_filter="Cu",
    top_k=10,
    min_experiments=1
)

# Results are sorted by avg_z_score
for rec in result.recommendations:
    print(f"Z-Score: {rec.avg_z_score:.2f}")
    print(f"  {rec.catalyst}/{rec.ligand}/{rec.base}/{rec.solvent}")
    print(f"  Yield: {rec.avg_yield:.1f}%, Confidence: {rec.confidence_score:.1f}")
```

## Benefits of Z-Score Ranking

1. **Relative Performance**: Z-score measures how much better/worse a condition is compared to the database average
2. **Statistical Significance**: Higher z-scores indicate conditions that are statistically significantly better
3. **Standardized Metric**: Z-scores are normalized across different reaction types
4. **Identifies Outliers**: High z-scores (>2.0) indicate exceptional conditions

## Z-Score Interpretation

- **Z > 2.0**: Excellent condition (top ~2.5%)
- **Z > 1.0**: Good condition (better than ~84%)
- **Z = 0.0**: Average condition
- **Z < 0.0**: Below average condition
- **Z < -2.0**: Poor condition (bottom ~2.5%)

## Impact on Copper C-N Coupling Example

**Before (success rate ranking):**
- Conditions ranked by success rate and yield
- No differentiation for conditions with similar yields

**After (z-score ranking):**
- Top condition: CuI/PPBO/NaOtBu/tAmOH (Z=2.61, 49.8% yield)
- Second: CuI/DMPAO/NaOtBu/tAmOH (Z=2.59, 49.4% yield)
- Clear ranking based on relative performance vs. database

## Files Modified

1. `chemtools/HTE/recommender.py` - Core ranking logic
2. `chemtools/HTE/cli.py` - CLI output formatting
3. `demo_hte_condition_summary.py` - Demo script table

## Testing

All functionality tested with copper-catalyzed C-N coupling example:
- ✅ CLI standard output shows z-score prominently
- ✅ Compact format includes z-score
- ✅ JSON output includes avg_z_score field
- ✅ Demo script shows z-score in table
- ✅ Conditions correctly ranked by z-score

## Notes

- Z-score calculation uses the existing `z-Score` column from HTE database
- Confidence score still considers sample size and yield, but z-score dominates
- All existing functionality preserved, just re-weighted for ranking
