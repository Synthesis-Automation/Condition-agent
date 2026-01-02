# Unified Recommender (Dataset + Protocol + HTE)

The unified recommender searches across three sources in one system:

- Reaction datasets (`data/reaction_dataset/*.jsonl`)
- Literature protocols (`data/protocol_db_v2/*.json`)
- HTE screens (`data/HTE_db/*.csv`)

It combines:
- DRFP similarity for entries with reaction SMILES
- Feature-tag similarity for entries without reaction SMILES (HTE)

## Build the index

```bash
python data-processor/build_unified_recommendation_index.py
```

Default output:

```
build/unified_recommendation_index/
  index.jsonl
  fingerprints.npz
  stats.json
```

## Basic usage

```python
from chemtools.recommend import UnifiedRecommender

recommender = UnifiedRecommender()
results = recommender.recommend(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    top_k=5,
)

for r in results:
    print(r.rank, r.source_type, r.family, r.similarity)
```

## Output fields

Each `RecommendationResult` includes:

- `source_type`: `dataset`, `protocol`, or `hte`
- `family`: canonical reaction family (best effort)
- `similarity`: combined DRFP + feature score
- `drfp_similarity`: DRFP component (if available)
- `feature_similarity`: tag similarity (if available)

Use `include_details=True` in `recommender.recommend()` to load the full
record for the top hits.
