# Recommendation API Overview

The recommendation system is unified around one index that spans datasets, protocols, and HTE screens.

## Main entry points

- `chemtools.recommend.UnifiedRecommender`
- `chemtools.recommend.recommend_from_reaction(...)`
- `chemtools.recommend.recommend_conditions_structured(...)`

## Example

```python
from chemtools.recommend import recommend_from_reaction

result = recommend_from_reaction(
    reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    k=5,
)

for rec in result["recommended_conditions"]:
    print(rec["rank"], rec["source_type"], rec["similarity"])
```

## Index build

```bash
python data-processor/build_unified_recommendation_index.py
```
