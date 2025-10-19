# Quick Start: Using SMARTS-Based Protocol Matching

## Overview

The protocol recommendation system now uses `reaction_SMARTS` patterns to intelligently match your input reactions to the most structurally similar protocols.

## Basic Usage

```python
from chemtools.protocol.recommend import ProtocolRecommender

# Initialize the recommender
recommender = ProtocolRecommender()

# Get recommendations with SMARTS filtering (recommended)
results = recommender.recommend(
    reaction_smiles='O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1',
    k=5,
    use_smarts_filter=True  # Default: enabled
)

# Access recommendations
for rec in results['recommended_conditions']:
    print(f"Rank {rec['rank']}: {rec['protocol']['title']}")
    print(f"  Confidence: {rec['confidence']:.4f}")
    print(f"  Family: {rec['protocol']['reaction_family']}")
    print(f"  SMARTS patterns: {rec['protocol']['reaction_smarts']}")
```

## When to Use SMARTS Filtering

### ✅ Use SMARTS Filtering (Default) When:
- You want structurally similar protocols only
- You have a specific reaction type in mind
- You want high-precision recommendations
- Your protocols have well-defined SMARTS patterns

### ⚠️ Disable SMARTS Filtering When:
- You want broad exploratory results
- Protocols lack SMARTS patterns
- You want to see all similar reactions regardless of structure
- You're using purely fingerprint-based similarity

```python
# Disable SMARTS filtering for broader results
results = recommender.recommend(
    reaction_smiles='...',
    k=5,
    use_smarts_filter=False  # Only use DRFP similarity
)
```

## Example: Comparing Results

```python
reaction = "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>O=C(C)c1ccc(C2C(CCC2)=O)cc1"

# With SMARTS filtering
precise = recommender.recommend(reaction, k=5, use_smarts_filter=True)
print(f"Precise matches: {len(precise['recommended_conditions'])} protocols")

# Without SMARTS filtering
broad = recommender.recommend(reaction, k=5, use_smarts_filter=False)
print(f"Broad matches: {len(broad['recommended_conditions'])} protocols")
```

**Expected output:**
```
Precise matches: 1 protocols  # Only alpha-arylation
Broad matches: 3 protocols    # Includes other reactions with some similarity
```

## Adding SMARTS to Your Protocols

Edit your protocol JSON to include `reaction_SMARTS`:

```json
{
  "reaction": {
    "reaction_smiles": "O=C1CCCC1.Brc2ccc(C(C)=O)cc2>>O=C(C)c3ccc(C4C(CCC4)=O)cc3",
    "reaction_SMARTS": [
      "O=C1CCCC1.Br[a]>>[a]C2C(CCC2)=O"
    ],
    "family": "Pd/Enamine_α-Arylation_C(sp3)-C(sp2)",
    ...
  }
}
```

### SMARTS Pattern Tips:
- Use `[a]` for generic aromatic carbons
- Use `[A]` for any heavy atom
- Use specific patterns for key functional groups
- Include both reactant and product patterns
- Can include multiple patterns per protocol

After adding patterns, rebuild the index:
```bash
python -m chemtools.protocol.cli build
```

## API Endpoint (If exposed via FastAPI)

```bash
curl -X POST http://localhost:8000/api/v1/protocol/recommend \
  -H "Content-Type: application/json" \
  -d '{
    "reaction_smiles": "O=C1CCCC1.Brc1ccc(C(C)=O)cc1>>...",
    "k": 5,
    "use_smarts_filter": true
  }'
```

## Troubleshooting

### No matches returned with SMARTS filtering
**Problem**: SMARTS filtering returns 0 protocols

**Solutions**:
1. Check if protocols have `reaction_SMARTS` defined
2. Verify SMARTS patterns are correct
3. Try disabling SMARTS filtering temporarily
4. Use more generic SMARTS patterns

```python
# Debug: Check what SMARTS patterns exist
for filename, record in recommender.indexer.records.items():
    print(f"{filename}: {record.reaction_smarts}")
```

### Index needs rebuilding
**Problem**: New SMARTS patterns not recognized

**Solution**: Rebuild the protocol index
```bash
python -m chemtools.protocol.cli build --force
```

### RDKit not available
**Problem**: System falls back to permissive mode

**Note**: The system gracefully handles missing RDKit. Install RDKit for full functionality:
```bash
pip install rdkit
```

## Performance Considerations

- **SMARTS filtering is fast**: Typically < 10ms per protocol
- **Reduces DRFP computation**: Only computes similarity for matched protocols
- **Recommended for production**: Improves both speed and quality

## For More Information

See full documentation: `docs/PROTOCOL_SMARTS_ENHANCEMENT.md`
