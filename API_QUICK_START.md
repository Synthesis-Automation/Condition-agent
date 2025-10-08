# API Quick Reference

## 🚀 Start Server

```powershell
uvicorn app.main:app --reload --port 8000
```

## 📖 Interactive Docs

- **Swagger UI:** http://localhost:8000/docs (Try endpoints in browser!)
- **ReDoc:** http://localhost:8000/redoc

## 🧪 Most Common Endpoints

### Get Recommendations (Baseline)
```bash
POST /api/v1/recommend
{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "k": 50}
```

### Get Recommendations (Fusion - NEW)
```bash
POST /api/v1/recommend/fusion
{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", "k": 50, "max_variants": 5}
```

### Detect Reaction Family
```bash
POST /api/v1/router/detect-family
{"reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"}
```

### List Available Cores
```bash
GET /api/v1/cores?family=C_N_Coupling_Pd&limit=10
```

### Health Check
```bash
GET /health
```

## 🐍 Python Quick Test

```python
import requests

# Get fusion recommendations
r = requests.post(
    'http://localhost:8000/api/v1/recommend/fusion',
    json={
        'reaction': 'Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1',
        'k': 50,
        'max_variants': 5
    }
)

result = r.json()
print(result['fusion_meta']['adaptive_weights'])
```

## ✅ Run Tests

```powershell
# Fusion API tests
python test_fusion_api_simple.py

# All integration tests
pytest tests/ -v
```

## 📚 Full Documentation

- **All Endpoints:** `API_GUIDE.md`
- **Fusion Details:** `FUSION_API_COMPLETE.md`
- **ChemTools Library:** `CHEMTOOLS_QUICKSTART.md`
