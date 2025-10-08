# ChemTools API Documentation Summary

Complete guide to starting, testing, and using all ChemTools API endpoints.

---

## 📚 Available Documentation

| Document | Purpose | Best For |
|----------|---------|----------|
| **API_QUICK_START.md** | Quick reference card | Getting started fast |
| **API_GUIDE.md** | Complete endpoint reference | Detailed API usage |
| **FUSION_API_COMPLETE.md** | Fusion system details | Understanding fusion recommendations |
| **FUSION_API_QUICK_REF.md** | Fusion quick reference | Using fusion endpoint |

---

## 🚀 Quick Start (60 seconds)

### 1. Start the Server

```powershell
uvicorn app.main:app --reload --port 8000
```

Expected output:
```
INFO:     Uvicorn running on http://127.0.0.1:8000 (Press CTRL+C to quit)
INFO:     Started server process
INFO:     Application startup complete.
```

### 2. Test It Works

**Option A: Browser** (Easiest)
- Open http://localhost:8000/docs
- Click on any endpoint
- Click "Try it out"
- Fill in parameters
- Click "Execute"

**Option B: Command Line**
```powershell
curl http://localhost:8000/health
# Should return: {"ok":true}
```

**Option C: Python Script**
```python
import requests
r = requests.get('http://localhost:8000/health')
print(r.json())  # {'ok': True}
```

### 3. Try a Real Recommendation

```powershell
python test_fusion_api_simple.py
```

Expected output:
```
✅ PASS: Basic Fusion Endpoint
✅ PASS: Fusion vs Baseline Comparison
Results: 2/2 tests passed
```

---

## 🎯 Most Common Use Cases

### Use Case 1: Get Reaction Recommendations

**Scenario:** You have a reaction and want condition recommendations.

**Endpoint:** `POST /api/v1/recommend`

**Code:**
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
top = result['formatted']['recommended_conditions'][0]
print(f"Use {top['summary']['core']} + {top['summary']['base']} in {top['summary']['solvent']}")
```

**Output:**
```
Use Cu + Tripotassium phosphate in Sulfolane
```

---

### Use Case 2: Get Fusion Recommendations (NEW)

**Scenario:** You want smarter recommendations with adaptive weighting.

**Endpoint:** `POST /api/v1/recommend/fusion`

**Code:**
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
print(f"  Precedent: {meta['adaptive_weights']['α']:.1%}")
print(f"  Analytics: {meta['adaptive_weights']['β']:.1%}")
print(f"  Rules: {meta['adaptive_weights']['γ']:.1%}")

print("\nWhy these weights?")
for reason in meta['reasoning']:
    print(f"  - {reason}")
```

**Output:**
```
Adaptive Weights:
  Precedent: 32.9%
  Analytics: 50.3%
  Rules: 16.8%

Why these weights?
  - Low diversity (0.09) → precedents may be biased
  - Low similarity (0.50) → precedents less relevant
  - No strong rule match → low rule weight
```

---

### Use Case 3: Detect Reaction Family

**Scenario:** You need to know what type of reaction you have.

**Endpoint:** `POST /api/v1/router/detect-family`

**Code:**
```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/router/detect-family',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1'}
)

family = response.json()['family']
print(f"Detected: {family}")
```

**Output:**
```
Detected: C_N_Coupling_Cu
```

---

### Use Case 4: Design a Screening Plate

**Scenario:** You want 96 diverse conditions for HTS.

**Endpoint:** `POST /api/v1/design_plate`

**Code:**
```python
import requests

response = requests.post(
    'http://localhost:8000/api/v1/design_plate',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'plate_size': 96
    }
)

plate = response.json()['plate']
print(f"Generated {len(plate['conditions'])} conditions")
for cond in plate['conditions'][:3]:
    print(f"  {cond['well']}: {cond['core']} + {cond['base']} in {cond['solvent']}")
```

---

## 🧪 Testing

### Run All Tests

```powershell
# Fusion API tests
python test_fusion_api_simple.py

# All integration tests
pytest tests/ -v

# Specific test file
pytest tests/test_fusion_integration.py -v
```

### Test Individual Endpoints

**Using cURL:**
```bash
# Health check
curl http://localhost:8000/health

# Get recommendations
curl -X POST "http://localhost:8000/api/v1/recommend" \
  -H "Content-Type: application/json" \
  -d '{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "k": 10}'

# Fusion recommendations
curl -X POST "http://localhost:8000/api/v1/recommend/fusion" \
  -H "Content-Type: application/json" \
  -d '{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "k": 50, "max_variants": 5}'
```

**Using Python:**
```python
import requests

# Test baseline
r1 = requests.post(
    'http://localhost:8000/api/v1/recommend',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1', 'k': 10}
)
print(f"Baseline: {r1.status_code}")

# Test fusion
r2 = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1', 'k': 50, 'max_variants': 5}
)
print(f"Fusion: {r2.status_code}")
print(f"Has fusion_meta: {'fusion_meta' in r2.json()}")
```

---

## 📋 Complete Endpoint List

### High-Level Endpoints (Most Used)

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/v1/recommend` | POST | Get condition recommendations (baseline) |
| `/api/v1/recommend/fusion` | POST | Get fusion recommendations (adaptive) ⭐ NEW |
| `/api/v1/router/detect-family` | POST | Detect reaction family |
| `/api/v1/design_plate` | POST | Design HTS plate |
| `/health` | GET | Check server health |

### Mid-Level Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/v1/cores` | GET | List catalytic cores |
| `/api/v1/precedent/knn` | POST | Find similar precedents |
| `/api/v1/constraints/filter` | POST | Filter by constraints |
| `/api/v1/explain/precedents` | POST | Explain precedent selection |

### Low-Level Endpoints (Advanced)

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/v1/smiles/normalize` | POST | Normalize SMILES |
| `/api/v1/featurize/ullmann` | POST | Generate Ullmann features |
| `/api/v1/featurize/role-aware/molecule` | POST | Featurize with role info |
| `/api/v1/condition-core/parse` | POST | Parse condition text |

**See `API_GUIDE.md` for complete details on all endpoints.**

---

## 🔧 Configuration

### Environment Variables

```bash
# Disable RDKit if not installed
export CHEMTOOLS_DISABLE_RDKIT=1

# Precomputed DRFP fingerprints (optional)
export CHEMTOOLS_DRFPPATH=/path/to/drfp_bundle.npz
```

### Server Options

```powershell
# Development (auto-reload)
uvicorn app.main:app --reload --port 8000

# Production (with workers)
uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4

# Custom port
uvicorn app.main:app --reload --port 9000
```

---

## 📊 Baseline vs Fusion Comparison

| Feature | Baseline (`/recommend`) | Fusion (`/recommend/fusion`) |
|---------|------------------------|------------------------------|
| **Method** | k-NN only | Multi-evidence fusion |
| **Evidence Sources** | Precedents only | Precedents + Analytics + Rules + ML |
| **Weighting** | Fixed | Adaptive (data-quality aware) |
| **Metadata** | Basic | Rich (weights, reasoning, quality) |
| **Novel Reactions** | May struggle | More robust |
| **Batch Effects** | Not detected | Detected and adjusted |
| **Explainability** | Limited | High (reasoning provided) |
| **Speed** | Fast | Slightly slower |
| **Best For** | Quick lookups | Complex/novel reactions |

---

## 💡 Tips & Best Practices

### Performance

1. **Start with k=50**: Good balance of speed and coverage
2. **Use constraints**: Filter early to reduce processing
3. **Cache results**: Same reaction = same output
4. **Batch with plate design**: More efficient than individual requests

### Accuracy

1. **Use fusion for novel reactions**: Better handling of unusual cases
2. **Check fusion reasoning**: Understand why weights were chosen
3. **Watch diversity scores**: Low diversity = potential batch effects
4. **Compare baseline vs fusion**: See if fusion adds value

### Debugging

1. **Check `/health` first**: Verify server is running
2. **Use Swagger UI**: Interactive testing in browser
3. **Check status codes**: 200=success, 400=bad input, 500=server error
4. **Read error messages**: `detail` field explains what went wrong

---

## 🐛 Troubleshooting

### Server won't start

```powershell
# Check if port 8000 is in use
netstat -ano | findstr :8000

# Try a different port
uvicorn app.main:app --reload --port 9000
```

### Endpoint returns 500 error

```python
# Check the response for details
import requests
r = requests.post('http://localhost:8000/api/v1/recommend', json={'reaction': 'invalid'})
print(r.status_code)
print(r.json()['detail'])  # Error message
```

### Tests fail

```powershell
# Make sure server is running first
curl http://localhost:8000/health

# Check Python dependencies
pip install -r requirements.txt

# Run with verbose output
pytest tests/ -v -s
```

---

## 📖 Learning Path

### Beginner (First 10 minutes)

1. Start server: `uvicorn app.main:app --reload --port 8000`
2. Open browser: http://localhost:8000/docs
3. Try `/health` endpoint in Swagger UI
4. Try `/api/v1/recommend` with example reaction

### Intermediate (Next 30 minutes)

1. Read `API_QUICK_START.md`
2. Run `python test_fusion_api_simple.py`
3. Try fusion endpoint in Swagger UI
4. Compare baseline vs fusion results

### Advanced (Next hour)

1. Read `API_GUIDE.md` for all endpoints
2. Read `FUSION_API_COMPLETE.md` for fusion details
3. Write custom scripts using the API
4. Explore advanced features (constraints, plate design)

---

## 🔗 Quick Links

- **Interactive Docs:** http://localhost:8000/docs
- **Alternative Docs:** http://localhost:8000/redoc
- **Health Check:** http://localhost:8000/health

---

## 📞 Next Steps

- **Try It Now:** Start server and open http://localhost:8000/docs
- **Read More:** See `API_GUIDE.md` for complete endpoint reference
- **Understand Fusion:** See `FUSION_API_COMPLETE.md` for fusion system details
- **Use Python Library:** See `CHEMTOOLS_QUICKSTART.md` for direct library usage

---

**Ready to start?**

```powershell
uvicorn app.main:app --reload --port 8000
```

Then open http://localhost:8000/docs in your browser! 🚀
