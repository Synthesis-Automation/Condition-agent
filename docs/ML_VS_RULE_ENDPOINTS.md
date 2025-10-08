# Recommendation Endpoints - ML vs Rule-Based

## Summary

The ChemTools API provides **BOTH** ML-based and Rule-based recommendation endpoints:

| Endpoint | Type | Description |
|----------|------|-------------|
| `POST /api/v1/recommend/conditions` | **ML-based** | Precedent KNN recommendations |
| `POST /api/v1/recommend` | **ML-based** | Simplified ML recommendations |
| `POST /match` | **Rule-based** | SchemeConditionDB pattern matching |

---

## 1. ML-Based Recommendations (Primary)

### Endpoint: `POST /api/v1/recommend/conditions`

**Type:** Machine Learning (Precedent KNN)

**What it does:**
- Searches a database of precedent reactions
- Uses k-nearest neighbors (KNN) to find similar reactions
- Returns top conditions based on structural similarity
- Provides confidence scores based on precedent support

**Input:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
  "reaction_type": "C_N_Coupling_Pd",  // Optional
  "k": 50,                              // Number of neighbors
  "limit": 3                            // Number of recommendations
}
```

**Output Format:** Compact, UI-formatted (3-5 KB)
```json
{
  "meta": {
    "generated_at": "2025-10-06T14:30:00Z",
    "model": "ML-precedent-knn",
    "status": "success",
    "processing_time_ms": 95.5
  },
  "input": {
    "reaction_smiles": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    "requested_reaction_type": "C_N_Coupling_Pd"
  },
  "detection": {
    "detected_reaction_type": "C-N Coupling (Pd/Buchwald)",
    "confidence": 0.95,
    "method": "rxn-insight-ml"
  },
  "recommendations": [
    {
      "rank": 1,
      "confidence": 0.475,
      "reagents": [
        {"role": "base", "name": "Potassium tert-butoxide", "cas": "865-48-5"},
        {"role": "solvent", "name": "Toluene", "cas": "108-88-3"}
      ],
      "conditions": {
        "temperature": {"value": 110, "unit": "°C"},
        "time": {"value": 24, "unit": "hours"}
      },
      "precedent_count": 19
    }
  ]
}
```

**Key Features:**
- ✅ Auto-detects reaction type if not specified
- ✅ Returns multiple ranked recommendations
- ✅ Provides confidence scores (0-1)
- ✅ Includes precedent counts (how many similar reactions)
- ✅ Compact output (~3-5 KB)
- ✅ CAS numbers for all reagents

**Use Cases:**
- Novel reactions without exact rule matches
- Optimization based on literature precedents
- Getting multiple alternative conditions
- When you want confidence scores

---

## 2. ML-Based Recommendations (Alternative)

### Endpoint: `POST /api/v1/recommend`

**Type:** Machine Learning (Raw output)

**What it does:**
- Same underlying ML engine as `/api/v1/recommend/conditions`
- Returns raw, unformatted output
- Includes full feature vectors and internal data

**Input:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
  "k": 50,
  "relax": {},
  "constraints": {}
}
```

**Output Format:** Large, raw output (~200-300 KB)
- Includes all internal ML features
- Contains feature vectors, fields, masks
- More detailed but less robot-friendly

**Use Cases:**
- Debugging ML model behavior
- Research/analysis of recommendation logic
- When you need full internal details

---

## 3. Rule-Based Recommendations

### Endpoint: `POST /match`

**Type:** Rule-based (SchemeConditionDB)

**What it does:**
- Pattern matching against curated rule databases
- Exact or fuzzy structural matching
- Returns deterministic, expert-curated conditions
- No confidence scores (rules are binary: match or no match)

**Input:**
```json
{
  "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
  "db": "data/conditionDB/cn_coupling_pd_db.json",  // Optional
  "include_trace": false
}
```

**Available Databases:**
- `cn_coupling_pd_db.json` - C-N Coupling (Pd/Buchwald)
- `cn_coupling_cu_db.json` - C-N Coupling (Cu/Ullmann)
- `cn_coupling_ni.json` - C-N Coupling (Ni)
- `amide_formation_db.json` - Amide Formation
- `suzuki_db.json` - Suzuki Coupling

**Output Format:**
```json
{
  "match_type": "scheme",
  "entry_id": "SCDB-BH-ARBR-ANILINE-core",
  "scheme_name": "ArBr + Aniline",
  "conditions": {
    "core": "Pd/XPhos",
    "base": "K3PO4",
    "solvent": "Toluene",
    "temperature": "110°C",
    "time": "24h",
    "ligand": "XPhos",
    "catalyst": "Pd2(dba)3",
    "pd_source": ["Pd2(dba)3", "Pd(OAc)2"],
    "loadings": {
      "Pd_mol%": [0.5, 2.0],
      "ligand_mol%": [1.0, 4.0]
    }
  },
  "note": "Standard Buchwald-Hartwig conditions"
}
```

**Key Features:**
- ✅ Expert-curated conditions
- ✅ Multiple options per category (e.g., multiple bases)
- ✅ Loading ranges (mol%)
- ✅ Exact or fuzzy matching
- ✅ Fast (<50ms typically)
- ❌ No confidence scores (deterministic)
- ❌ Limited to curated reaction types

**Use Cases:**
- Well-known reaction types with established protocols
- When you want expert-recommended conditions
- Fast lookups without ML overhead
- Validating ML recommendations against rules

---

## Comparison Table

| Feature | ML (`/api/v1/recommend/conditions`) | Rule-Based (`/match`) |
|---------|-------------------------------------|----------------------|
| **Method** | Precedent KNN | Pattern matching |
| **Speed** | Slow (90-100s) | Fast (<1s) |
| **Coverage** | Any reaction in precedent DB | Only curated types |
| **Confidence** | Yes (0-1 score) | No (binary match) |
| **Flexibility** | Handles novel reactions | Exact patterns only |
| **Output Size** | Compact (3-5 KB) | Very compact (1-2 KB) |
| **Multiple Options** | Yes (ranked by confidence) | Yes (alternatives per category) |
| **Precedent Info** | Yes (literature refs) | No |
| **Auto-detection** | Yes | No |

---

## Recommended Workflow for Robots

### Option 1: ML-Only (Flexible)
```
1. Send reaction to /api/v1/recommend/conditions
2. Get 3-5 ranked recommendations
3. Use top recommendation (highest confidence)
4. Fall back to rank 2 if rank 1 fails
```

**Pros:** Works for any reaction, provides alternatives  
**Cons:** Slower (90s), may have lower confidence for novel reactions

---

### Option 2: Rule-First, ML-Fallback (Recommended)
```
1. Try /match first (fast, <1s)
2. If match found:
   → Use rule-based conditions (expert-curated)
3. If no match:
   → Fall back to /api/v1/recommend/conditions
   → Use ML recommendations
```

**Pros:** Fast for known reactions, ML backup for novel ones  
**Cons:** Requires two endpoints

---

### Option 3: Hybrid Validation
```
1. Get ML recommendations from /api/v1/recommend/conditions
2. Get rule recommendations from /match
3. Compare:
   - If they agree → high confidence, use either
   - If they differ → use rule-based (expert-curated)
   - If no rule match → use ML only
```

**Pros:** Maximum confidence, validation  
**Cons:** Slower, more complex

---

## Code Examples

### Python: ML Recommendations
```python
import requests

response = requests.post(
    "http://localhost:8000/api/v1/recommend/conditions",
    json={
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "reaction_type": "C_N_Coupling_Pd",
        "limit": 3
    }
)

result = response.json()
top_rec = result["recommendations"][0]

print(f"Confidence: {top_rec['confidence']:.2%}")
print(f"Reagents:")
for reagent in top_rec["reagents"]:
    if reagent["role"] not in ["starting_material"]:
        print(f"  - {reagent['name']} [{reagent['role']}]: CAS {reagent['cas']}")
```

### Python: Rule-Based
```python
import requests

response = requests.post(
    "http://localhost:8000/match",
    json={
        "reaction": "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
        "db": "data/conditionDB/cn_coupling_pd_db.json"
    }
)

result = response.json()
conditions = result["conditions"]

print(f"Match: {result['entry_id']}")
print(f"Core: {conditions['core']}")
print(f"Base options: {conditions.get('base', [])}")
print(f"Solvent options: {conditions.get('solvent', [])}")
```

### Python: Hybrid Approach
```python
import requests

BASE_URL = "http://localhost:8000"
reaction = "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1"

# Try rule-based first
try:
    rule_response = requests.post(
        f"{BASE_URL}/match",
        json={"reaction": reaction, "db": "data/conditionDB/cn_coupling_pd_db.json"}
    )
    if rule_response.status_code == 200:
        print("✅ Rule match found - using expert conditions")
        conditions = rule_response.json()["conditions"]
        # Use rule-based conditions
    else:
        raise ValueError("No rule match")
except:
    # Fall back to ML
    print("⚠️ No rule match - using ML recommendations")
    ml_response = requests.post(
        f"{BASE_URL}/api/v1/recommend/conditions",
        json={"reaction": reaction, "limit": 3}
    )
    recommendations = ml_response.json()["recommendations"]
    # Use ML recommendations
```

---

## Summary

**For Robot Integration:**

✅ **Use `/api/v1/recommend/conditions` (ML) as your PRIMARY endpoint**
- Compact output (3-5 KB)
- Works for any reaction
- Provides confidence scores
- CAS numbers for reagent dispensing
- Multiple ranked alternatives

✅ **Use `/match` (Rule-based) for VALIDATION or FAST LOOKUP**
- Expert-curated conditions
- Very fast (<1s)
- Good for common reaction types

❌ **Avoid `/api/v1/recommend` (raw ML)**
- Large output (200+ KB)
- Same ML engine, just unformatted
- Not robot-friendly

**Recommended Strategy:** Start with ML endpoint, optionally validate against rules for known reaction types.
