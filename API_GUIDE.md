# ChemTools API Guide

Complete guide to all available API endpoints in the ChemTools FastAPI application.

## Quick Start

### Start the Server

```powershell
# Development mode (auto-reload)
uvicorn app.main:app --reload --port 8000

# Or use Makefile
make run
```

### Access Interactive Documentation

Once the server is running:
- **Swagger UI (Interactive):** http://localhost:8000/docs
- **ReDoc (Clean Docs):** http://localhost:8000/redoc

### Test All Endpoints

```powershell
# Run all API tests
pytest tests/

# Test specific endpoints
python test_fusion_api_simple.py
```

---

## 📋 Endpoint Categories

1. [Health & Monitoring](#health--monitoring)
2. [SMILES Operations](#smiles-operations)
3. [Reaction Detection](#reaction-detection)
4. [Featurization](#featurization)
5. [Precedent Search](#precedent-search)
6. [Recommendations](#recommendations) ⭐ Most Used
7. [Condition Core](#condition-core)
8. [Constraints & Filtering](#constraints--filtering)
9. [Explanations](#explanations)
10. [Plate Design](#plate-design)

---

## Health & Monitoring

### GET `/health`

Check if the API server is running.

**Request:**
```bash
curl http://localhost:8000/health
```

**Response:**
```json
{"ok": true}
```

**Python:**
```python
import requests
r = requests.get('http://localhost:8000/health')
print(r.json())  # {'ok': True}
```

---

### GET `/metrics`

Get server metrics and performance statistics.

**Request:**
```bash
curl http://localhost:8000/metrics
```

---

## SMILES Operations

### POST `/api/v1/smiles/normalize`

Normalize and canonicalize SMILES strings.

**Request:**
```json
{
  "smiles": "c1ccccc1"
}
```

**Response:**
```json
{
  "normalized": "c1ccccc1",
  "canonical": "c1ccccc1"
}
```

**Python:**
```python
import requests

r = requests.post(
    'http://localhost:8000/api/v1/smiles/normalize',
    json={'smiles': 'c1ccccc1'}
)
print(r.json())
```

---

## Reaction Detection

### POST `/api/v1/router/detect-family`

Detect the reaction family (e.g., C-N coupling, Suzuki, etc.).

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
}
```

**Response:**
```json
{
  "family": "C_N_Coupling_Cu",
  "confidence": 0.95,
  "detected_by": "SMARTS_match"
}
```

---

### POST `/api/v1/reaction/detect-type`

Detect detailed reaction type using RxnMapper/RxnInsight.

**Request:**
```json
{
  "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
}
```

**Response:**
```json
{
  "reaction_type": "Suzuki coupling",
  "confidence": 0.98
}
```

---

## Featurization

### POST `/api/v1/featurize/ullmann`

Generate Ullmann-specific features for C-N coupling reactions.

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
}
```

**Response:**
```json
{
  "features": {
    "aryl_halide": {...},
    "amine": {...},
    "coupling_partner_descriptors": [...]
  }
}
```

---

### POST `/api/v1/featurize/role-aware/molecule`

Featurize a molecule with role-aware descriptors.

**Request:**
```json
{
  "smiles": "c1ccccc1",
  "roles": ["aromatic", "solvent"]
}
```

---

### POST `/api/v1/featurize/role-aware/reaction`

Featurize a complete reaction with role-aware features.

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "roles": {
    "reactant_1": ["aryl_halide"],
    "reactant_2": ["amine"]
  }
}
```

---

### GET `/api/v1/featurize/role-aware/fields`

Get available role-aware feature fields.

**Request:**
```bash
curl http://localhost:8000/api/v1/featurize/role-aware/fields
```

**Response:**
```json
{
  "fields": ["aromatic", "aliphatic", "electron_rich", "electron_poor", ...]
}
```

---

## Precedent Search

### GET `/api/v1/cores`

List available catalytic cores.

**Request:**
```bash
curl "http://localhost:8000/api/v1/cores?family=C_N_Coupling_Pd&limit=10"
```

**Query Parameters:**
- `family` (optional): Filter by reaction family
- `limit` (default: 200): Maximum number of cores
- `counts` (default: true): Include usage counts

**Response:**
```json
{
  "cores": [
    {"core": "Pd", "count": 5234},
    {"core": "Cu", "count": 3421},
    ...
  ]
}
```

---

### POST `/api/v1/precedent/knn`

Find k-nearest neighbor precedents.

**Request:**
```json
{
  "family": "C_N_Coupling_Pd",
  "features": {"aryl_halide": {...}, "amine": {...}},
  "k": 50,
  "relax": {}
}
```

**Response:**
```json
{
  "precedents": [
    {
      "core": "Pd",
      "base": "K3PO4",
      "solvent": "toluene",
      "similarity": 0.95,
      "reaction_id": "..."
    },
    ...
  ]
}
```

---

### POST `/api/v1/core/search`

Search for cores matching specific criteria.

**Request:**
```json
{
  "family": "C_N_Coupling_Pd",
  "query": "Pd"
}
```

---

## Recommendations

### POST `/api/v1/recommend` ⭐

**Get reaction condition recommendations (baseline k-NN method).**

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 50,
  "relax": {},
  "constraints": {}
}
```

**Response:**
```json
{
  "input": {
    "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "family": "C_N_Coupling_Cu"
  },
  "recommendation": {
    "core": "Cu",
    "base": "Tripotassium phosphate",
    "solvent": "Sulfolane"
  },
  "formatted": {
    "recommended_conditions": [
      {
        "rank": 1,
        "summary": {
          "core": "Cu",
          "base": "Tripotassium phosphate",
          "solvent": "Sulfolane",
          "confidence": "high"
        },
        "details": {...}
      },
      ...
    ]
  }
}
```

**Python Example:**
```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/recommend',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50
    }
)

result = response.json()
top_rec = result['formatted']['recommended_conditions'][0]
print(f"Top recommendation: {top_rec['summary']}")
```

---

### POST `/api/v1/recommend/fusion` ⭐ NEW

**Get fusion recommendations with adaptive evidence weighting.**

Combines multiple evidence sources (precedent, analytics, rules, ML) with intelligent adaptive weighting based on data quality.

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "k": 50,
  "max_variants": 5,
  "relax": {},
  "constraints": {}
}
```

**Response:**
```json
{
  "input": {...},
  "recommendation": {...},
  "formatted": {
    "recommended_conditions": [...]
  },
  "fusion_meta": {
    "adaptive_weights": {
      "α": 0.329,  // Precedent weight
      "β": 0.503,  // Analytics weight
      "γ": 0.168,  // Rule weight
      "δ": 0.000   // ML weight (not yet integrated)
    },
    "evidence_summary": {
      "precedents": 10,
      "diversity": 0.087,
      "dataset_size": 5552
    },
    "reasoning": [
      "Low diversity (0.09) → precedents may be biased",
      "Low similarity (0.50) → precedents less relevant",
      "No strong rule match → low rule weight"
    ]
  }
}
```

**Key Differences from Baseline:**
- ✅ Includes `fusion_meta` with adaptive weights
- ✅ Combines 4 evidence sources (vs 1 in baseline)
- ✅ Provides reasoning for recommendations
- ✅ Detects and adjusts for data quality issues
- ✅ More robust for novel/unusual reactions

**Python Example:**
```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50,
        'max_variants': 5
    }
)

result = response.json()
meta = result['fusion_meta']

print("Adaptive Weights:")
print(f"  Precedent (α): {meta['adaptive_weights']['α']:.2f}")
print(f"  Analytics (β): {meta['adaptive_weights']['β']:.2f}")
print(f"  Rules (γ):     {meta['adaptive_weights']['γ']:.2f}")
print(f"  ML (δ):        {meta['adaptive_weights']['δ']:.2f}")

print("\nReasoning:")
for reason in meta['reasoning']:
    print(f"  - {reason}")
```

---

### POST `/api/v1/recommend/conditions`

Recommend conditions with structured output.

**Request:**
```json
{
  "family": "C_N_Coupling_Pd",
  "features": {...},
  "k": 50
}
```

---

## Plate Design

### POST `/api/v1/design_plate`

Design a plate of reaction conditions for high-throughput screening.

**Request:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
  "plate_size": 96,
  "relax": {},
  "constraints": {}
}
```

**Response:**
```json
{
  "plate": {
    "size": 96,
    "conditions": [
      {"well": "A1", "core": "Pd", "base": "K3PO4", "solvent": "toluene"},
      {"well": "A2", "core": "Pd", "base": "Cs2CO3", "solvent": "dioxane"},
      ...
    ]
  }
}
```

---

## Condition Core

### POST `/api/v1/condition-core/parse`

Parse reaction conditions from text or reagent list.

**Request:**
```json
{
  "reagents": ["Pd(OAc)2", "XPhos", "K3PO4"],
  "text": "Heated in toluene at 100°C"
}
```

**Response:**
```json
{
  "core": "Pd",
  "base": "K3PO4",
  "ligand": "XPhos",
  "solvent": "toluene",
  "temperature": "100°C"
}
```

---

### POST `/api/v1/condition-core/validate-dataset`

Validate condition cores against dataset schema.

**Request:**
```json
{
  "family": "C_N_Coupling_Pd",
  "cores": ["Pd", "Cu", "Ni"]
}
```

---

## Constraints & Filtering

### POST `/api/v1/constraints/filter`

Filter candidates based on constraint rules.

**Request:**
```json
{
  "candidates": [
    {"core": "Pd", "base": "K3PO4", "solvent": "toluene"},
    {"core": "Cu", "base": "Cs2CO3", "solvent": "DMF"}
  ],
  "rules": {
    "exclude_solvents": ["DMF"],
    "require_core": ["Pd"]
  }
}
```

**Response:**
```json
{
  "filtered": [
    {"core": "Pd", "base": "K3PO4", "solvent": "toluene"}
  ],
  "excluded_count": 1
}
```

---

## Explanations

### POST `/api/v1/explain/precedents`

Get explanations for why specific precedents were selected.

**Request:**
```json
{
  "pack": {
    "precedents": [...]
  },
  "features": {...}
}
```

**Response:**
```json
{
  "explanations": [
    {
      "precedent_id": "...",
      "similarity": 0.95,
      "reasons": [
        "Similar aryl halide structure",
        "Matching amine nucleophilicity"
      ]
    }
  ]
}
```

---

## Legacy Endpoints

### POST `/match`

Legacy SCDB matching endpoint (deprecated, use `/api/v1/router/detect-family`).

---

## Common Request Patterns

### Basic Recommendation Workflow

```python
import requests

BASE_URL = 'http://localhost:8000'

# 1. Detect reaction family
family_response = requests.post(
    f'{BASE_URL}/api/v1/router/detect-family',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'}
)
family = family_response.json()['family']
print(f"Detected family: {family}")

# 2. Get recommendations
rec_response = requests.post(
    f'{BASE_URL}/api/v1/recommend',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50
    }
)
recommendations = rec_response.json()['formatted']['recommended_conditions']

# 3. Display top 3
for rec in recommendations[:3]:
    summary = rec['summary']
    print(f"Rank {summary['rank']}: {summary['core']} + {summary['base']} in {summary['solvent']}")
```

### Fusion Recommendation with Constraints

```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50,
        'max_variants': 5,
        'constraints': {
            'exclude_solvents': ['DMF', 'DMSO'],
            'prefer_core': ['Pd']
        }
    }
)

result = response.json()

# Check why weights were chosen
print("Weight Reasoning:")
for reason in result['fusion_meta']['reasoning']:
    print(f"  - {reason}")

# Get recommendations
for rec in result['formatted']['recommended_conditions'][:3]:
    print(f"\n{rec['summary']}")
```

---

## Testing Endpoints

### Quick Health Check

```bash
curl http://localhost:8000/health
```

### Test Recommendation Endpoint

```bash
curl -X POST "http://localhost:8000/api/v1/recommend" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 10
  }'
```

### Test Fusion Endpoint

```bash
curl -X POST "http://localhost:8000/api/v1/recommend/fusion" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "k": 50,
    "max_variants": 5
  }'
```

### Run Test Suites

```powershell
# Fusion API tests
python test_fusion_api_simple.py

# All integration tests
pytest tests/ -v

# Specific test file
pytest tests/test_fusion_integration.py -v
```

---

## Error Handling

All endpoints return standard HTTP status codes:

- `200 OK` - Success
- `400 Bad Request` - Invalid input
- `404 Not Found` - Resource not found
- `500 Internal Server Error` - Server error

**Error Response Format:**
```json
{
  "detail": "Error message describing what went wrong"
}
```

**Example:**
```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/recommend',
    json={'reaction': 'invalid'}
)

if response.status_code != 200:
    error = response.json()
    print(f"Error: {error['detail']}")
```

---

## Performance Tips

1. **Cache Results**: Identical reactions return identical results
2. **Batch Requests**: Use plate design for multiple conditions
3. **Limit k**: Start with k=50, increase only if needed
4. **Use Constraints**: Filter early to reduce processing time

---

## Interactive Testing

The **easiest way** to test endpoints is using the built-in Swagger UI:

1. Start the server: `uvicorn app.main:app --reload --port 8000`
2. Open http://localhost:8000/docs in your browser
3. Click on any endpoint to expand it
4. Click "Try it out"
5. Fill in the request body
6. Click "Execute"
7. See the response immediately

**No code required!** 🎉

---

## Next Steps

- **Fusion System:** See `FUSION_API_COMPLETE.md` for detailed fusion documentation
- **ChemTools Library:** See `CHEMTOOLS_QUICKSTART.md` for Python library usage
- **Development:** See `AGENTS.md` for contribution guidelines

---

## Support

- **API Docs:** http://localhost:8000/docs
- **Issues:** File on GitHub repository
- **Documentation:** See `docs/` folder for more guides
