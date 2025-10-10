# Protocol Database - Quick Start Guide

## 🚀 Get Started in 2 Minutes

### Step 1: Build the Index

```bash
cd c:\Git-softwares\Condition-agent
python -m chemtools.protocol.cli build
```

Output:
```
✅ Index built successfully!
Indexed 16 protocols
```

### Step 2: Test It

```bash
python test_protocol_recommendation.py
```

Expected:
```
✅ All tests completed successfully!
```

### Step 3: Use It

```python
from chemtools.protocol import ProtocolRecommender

# Initialize (loads index)
recommender = ProtocolRecommender()

# Find similar protocols
results = recommender.recommend(
    reaction_smiles='BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1',
    k=5
)

# Print matches
for match in results['matches']:
    print(f"{match['similarity']:.3f}: {match['source_title']}")
```

Output:
```
0.802: Copper-Catalyzed Suzuki-Miyaura Coupling of Unactivated Alkyl Halides
0.293: Palladium-catalyzed Suzuki–Miyaura Cross-coupling Reaction
...
```

## 📖 Common Tasks

### View Statistics

```bash
python -m chemtools.protocol.cli stats
```

### List All Families

```bash
python -m chemtools.protocol.cli list-families
```

### Show Protocols for a Family

```bash
python -m chemtools.protocol.cli show-family "Suzuki_Cu_alkyl_halide+aryl_boron"
```

### Filter by Tags

```python
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3,
    tags=['suzuki', 'palladium']
)
```

### Get Detailed Conditions

```python
results = recommender.recommend_with_details(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

for match in results['matches']:
    cond = match['conditions']
    print(f"Catalyst: {cond['catalyst']}")
    print(f"Solvent: {cond['solvent']}")
    print(f"Temp: {cond['temperature_C']} °C")
```

## 🔧 Maintenance

### Add New Protocols

1. Add JSON file to `data/protocol_db/`
2. Rebuild index: `python -m chemtools.protocol.cli build`
3. Done! (incremental update, very fast)

### Force Rebuild

```bash
python -m chemtools.protocol.cli build --force
```

## 📚 Learn More

- **Full Documentation:** `docs/PROTOCOL_MODULE.md`
- **Implementation Details:** `docs/PROTOCOL_MODULE_SUMMARY.md`
- **Schema:** `data/protocol_db/examples/Structured_Output_schema.json`

## ✨ Key Features

- ✅ **DRFP similarity search** (same as ML recommendation)
- ✅ **Fast indexing** with incremental updates
- ✅ **Filter by family and tags**
- ✅ **Extract experimental conditions**
- ✅ **CLI tools** for management
- ✅ **Production-ready** API

## 🎯 Use Cases

1. **Find similar protocols** for your reaction
2. **Get experimental procedures** from Organic Syntheses
3. **Filter by catalyst type** (Pd, Cu, Ni, etc.)
4. **Compare multiple protocols** side-by-side
5. **Extract standard conditions** for common reactions

That's it! The protocol recommendation system is ready to use. 🎉
