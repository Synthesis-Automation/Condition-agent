# Suzuki API Endpoint Testing Guide

## Overview

`test_suzuki_api.py` is a comprehensive test script that validates all three recommendation approaches for Suzuki coupling reactions via FastAPI endpoints.

## Quick Start

### 1. Start the FastAPI Server

```powershell
# Activate your virtual environment first
.\.venv\Scripts\Activate.ps1

# Start the server
uvicorn app.main:app --reload --port 8000
```

The server should be running at `http://localhost:8000`. You can verify by visiting `http://localhost:8000/docs` in your browser.

### 2. Run the Test Script

In a new terminal:

```powershell
# Activate virtual environment
.\.venv\Scripts\Activate.ps1

# Run tests
python test_suzuki_api.py
```

## What Gets Tested

The script tests **5 diverse Suzuki reactions** using **3 different approaches**:

### Tested Reactions

1. **Simple Ph-Ph Coupling** (benchmark)
   - `Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1`

2. **Electron-poor ArCl** (challenging)
   - `Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1`

3. **Electron-rich ArBr** (donor groups)
   - `Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1>>COc1ccc(-c2ccccc2)cc1`

4. **Heteroaryl Pyridine** (coordination issues)
   - `Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1`

5. **Ortho-substituted** (sterically hindered)
   - `Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C`

### Three Approaches

#### 1. Rule-Based (`/match`)
- **Method**: SMARTS pattern matching against rule database
- **Endpoint**: `POST /match`
- **Database**: `data/conditionDB/Suzuki_db.json`
- **Pros**: Fast, deterministic, chemistry knowledge-based
- **Cons**: Limited to known patterns in the database

**Example Request**:
```json
{
  "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
  "db": "data/conditionDB/Suzuki_db.json",
  "include_trace": true
}
```

#### 2. ML-Based (`/api/v1/recommend/conditions`)
- **Method**: DRFP-based k-NN precedent search
- **Endpoint**: `POST /api/v1/recommend/conditions`
- **Pros**: Data-driven, handles novel substrates, generalizes well
- **Cons**: Requires precedents, can be biased by batch effects

**Example Request**:
```json
{
  "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
  "reaction_type": null,
  "k": 50,
  "limit": 5,
  "relax": {},
  "constraints": {}
}
```

#### 3. Fusion (`/api/v1/recommend/fusion`)
- **Method**: Multi-source evidence with adaptive weighting
- **Endpoint**: `POST /api/v1/recommend/fusion`
- **Evidence Sources**:
  - **α (Precedent)**: DRFP-based k-NN with diversity scoring
  - **β (Analytics)**: Dataset statistics (catalytic systems, bases, solvents)
  - **γ (Rules)**: SMARTS-based scheme matching
  - **δ (ML)**: DRFP yield predictor (currently 0)
- **Pros**: Most robust, adaptive, transparent reasoning
- **Cons**: More complex, slower

**Example Request**:
```json
{
  "reaction": "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
  "k": 50,
  "max_variants": 5,
  "relax": {},
  "constraints": {}
}
```

## Output Examples

### Rule-Based Output
```
✓ Match Type: exact
  Entry: Ph-Br_standard_conditions

  Recommended Conditions:
    Catalyst: Pd(OAc)2 + XPhos
    Base: K3PO4
    Solvent: Dioxane
    Temperature: 100°C
    Time: 12 hours
```

### ML-Based Output
```
✓ Detected Type: Suzuki (confidence: 98.50%)
  Found 5 recommendations

  1. Confidence: 95.20% | Support: 47 precedents
     Catalyst: Pd(OAc)2/XPhos
     Base: K3PO4
     Solvent: 1,4-Dioxane
     Temperature: 100°C
     Time: 12 hours
```

### Fusion Output
```
✓ Fusion System Active
  Adaptive Weights:
    α (precedents) = 0.650
    β (analytics)  = 0.250
    γ (rules)      = 0.080
    δ (ML)         = 0.020

  Evidence Quality:
    Precedent count: 47
    Diversity score: 0.425 (OK)

  Reasoning:
    - High precedent count (47) supports strong evidence
    - Good diversity (0.425) suggests unbiased sampling
    - Rule match found, adding confidence

  1. Core: Pd/XPhos | Confidence: HIGH
     Base: K3PO4
     Solvent: Dioxane
     Scores: PS=0.892, AS=0.734, RS=0.900, MS=0.000
```

## Interpreting Fusion Scores

### Component Scores (PS, AS, RS, MS)
- **PS (Precedent Score)**: Quality from k-NN precedent search
- **AS (Analytics Score)**: Agreement with dataset statistics
- **RS (Rule Score)**: Match quality from SMARTS rules
- **MS (ML Score)**: Predicted yield from ML model (currently 0)

### Adaptive Weights (α, β, γ, δ)
- Automatically adjusted based on evidence quality
- Low diversity → reduces α (precedent weight)
- Strong rule match → increases γ (rule weight)
- Sum always equals 1.0

### Diversity Score
- **> 0.5**: Excellent diversity, reliable precedents
- **0.3-0.5**: Good diversity
- **< 0.3**: Low diversity, potential batch effect (precedent weight reduced)

## Troubleshooting

### Server Not Running
```
✗ Cannot connect to server at http://localhost:8000
  Error: [Connection Error]

  Please start the server with:
  uvicorn app.main:app --reload --port 8000
```
**Solution**: Start the FastAPI server in a separate terminal.

### Missing Database File
```
✗ Error: HTTP 404
  Database file not found: data/conditionDB/Suzuki_db.json
```
**Solution**: Ensure the database file exists in the correct location.

### Import Errors
```
ModuleNotFoundError: No module named 'chemtools'
```
**Solution**: Ensure you're in the correct directory and virtual environment is activated.

## Customization

### Test Different Reactions
Edit the `test_reactions` list in `main()`:

```python
test_reactions = [
    {
        "smiles": "YOUR_REACTION_SMILES",
        "description": "Description",
        "note": "Special notes"
    },
    # ... more reactions
]
```

### Test Single Approach
Modify the `test_reaction_comprehensive()` call:

```python
results = test_reaction_comprehensive(
    reaction=rxn_data["smiles"],
    description=rxn_data["description"],
    test_rule=True,      # Set to False to skip
    test_ml=True,        # Set to False to skip
    test_fusion=False    # Set to False to skip
)
```

### Change k Value
Modify the k parameter in test functions:

```python
test_ml_based(reaction, description, k=100)  # More precedents
test_fusion(reaction, description, k=25)     # Fewer precedents
```

## Related Files

- **API Server**: `app/main.py`
- **Recommendation Core**: `chemtools/recommend/core.py`
- **Fusion Engine**: `chemtools/ml/fusion_recommender.py`
- **Rule Matcher**: `chemtools/rule_scdb_matcher/`
- **Database**: `data/conditionDB/Suzuki_db.json`

## API Documentation

Full API documentation is available at:
- Interactive docs: http://localhost:8000/docs
- ReDoc: http://localhost:8000/redoc
