# Fusion API - Quick Reference

## Start Server

```powershell
uvicorn app.main:app --reload --port 8000
```

## Test Endpoints

```powershell
# Run simple test suite
python test_fusion_api_simple.py

# Run interactive test suite
python test_fusion_api.py
```

## Endpoints

### Fusion Recommendations (NEW)

```bash
POST http://localhost:8000/api/v1/recommend/fusion

{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 50,
  "max_variants": 5
}
```

**Returns:** Recommendations + `fusion_meta` with adaptive weights

### Baseline Recommendations

```bash
POST http://localhost:8000/api/v1/recommend

{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 50
}
```

**Returns:** Recommendations (no fusion_meta)

## Fusion Metadata Structure

```json
{
  "fusion_meta": {
    "adaptive_weights": {
      "α": 0.329,  // Precedent
      "β": 0.503,  // Analytics
      "γ": 0.168,  // Rules
      "δ": 0.000   // ML (not yet integrated)
    },
    "evidence_summary": {
      "precedents": 10,
      "diversity": 0.087,
      "dataset_size": 5552
    },
    "reasoning": [
      "Low diversity (0.09) → precedents may be biased",
      "Low similarity (0.50) → precedents less relevant"
    ]
  }
}
```

## Quick Test

```python
import requests

r = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1', 'k': 10, 'max_variants': 3}
)

meta = r.json()['fusion_meta']
print(f"Weights: α={meta['adaptive_weights']['α']:.2f}, β={meta['adaptive_weights']['β']:.2f}")
```

## Status

✅ **COMPLETE** - All tests passing (2/2)

## Documentation

- Full docs: `FUSION_API_COMPLETE.md`
- API docs: http://localhost:8000/docs
