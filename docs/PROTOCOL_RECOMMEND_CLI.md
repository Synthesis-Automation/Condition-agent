# Protocol Recommendation CLI - Quick Start

## Overview

The protocol recommendation CLI allows you to find matching experimental protocols for your reactions using DRFP (Differential Reaction Fingerprint) similarity and optional SMARTS structural filtering.

## Installation

The CLI is part of the `chemtools.protocol` package. No additional installation needed if you have the main project set up.

## Quick Start

### Basic Usage

```bash
# Find protocols for a Suzuki coupling
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
```

### Common Examples

```bash
# Get top 5 matches
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --k 5

# Filter by reaction family
python -m chemtools.protocol.recommend_cli "Brc1ccccc1.Nc1ccccc1>>" --family Buchwald

# Filter by tags
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --tags "Suzuki,Pd"

# Save results to JSON
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --output results.json --pretty

# Set minimum similarity threshold
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --min-similarity 0.5

# Disable SMARTS filtering (use DRFP only)
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --no-smarts-filter
```

## Output Format

### Human-Readable Format (Default)

```
======================================================================
Protocol Recommendations
======================================================================

Model: Protocol-DRFP
Status: success
Processing time: 184.6 ms

Detected type: Suzuki_Miyaura
Confidence: 0.295

Reaction: CCBr.c1ccccc1B(O)O>>CCc1ccccc1

Found 2 matching protocol(s):

======================================================================
Rank 1 - Similarity: 0.295
======================================================================
Title: Nickel-Catalyzed Suzuki-Miyaura Coupling...
Journal: Organic Syntheses (2016)
Family: Ni_Suzuki_ArylHalide+BoronicAcid_tAmOH

Conditions:
  Catalyst: NiCl2(PCy3)2
  Base: K3PO4
  Solvent: t-AmOH
  Temperature: 120 °C
  Time: 1.0 h

Source file: Ni_Suzuki_ArylHalide+BoronicAcid_tAmOH.json
```

### JSON Format (with --output)

```json
{
  "meta": {
    "model_type": "Protocol-DRFP",
    "status": "success",
    "processing_time_ms": 184.6
  },
  "input": {
    "reaction_smiles": "CCBr.c1ccccc1B(O)O>>CCc1ccccc1"
  },
  "detection": {
    "detected_type": "Suzuki_Miyaura",
    "confidence": 0.295
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.295,
      "conditions": {...},
      "protocol_metadata": {...}
    }
  ]
}
```

## Command-Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `reaction_smiles` | - | Reaction SMILES (required) | - |
| `--k` | `-n` | Number of results | 5 |
| `--family` | `-f` | Filter by reaction family | None |
| `--tags` | `-t` | Filter by tags (comma-separated) | None |
| `--min-similarity` | - | Minimum similarity threshold (0-1) | 0.0 |
| `--no-smarts-filter` | - | Disable SMARTS structural filtering | False |
| `--output` | `-o` | Save to JSON file | None |
| `--pretty` | `-p` | Pretty print JSON | False |
| `--verbose` | `-v` | Verbose logging | False |
| `--legacy-format` | - | Use legacy output format | False |

## Filtering Options

### By Reaction Family

```bash
# Find only Suzuki protocols
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>" --family Suzuki

# Find only Buchwald-Hartwig protocols
python -m chemtools.protocol.recommend_cli "Brc1ccccc1.Nc1ccccc1>>" --family Buchwald
```

Available families can be found by running:
```bash
python -m chemtools.protocol.cli list-families
```

### By Tags

```bash
# Find protocols tagged with "Pd" and "coupling"
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>" --tags "Pd,coupling"

# Find protocols with specific conditions
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>" --tags "aryl_halide,boronic_acid"
```

Available tags can be found by running:
```bash
python -m chemtools.protocol.cli list-tags
```

## How It Works

### Step-by-Step Process

1. **Load Protocol Index**: Loads precomputed DRFP fingerprints for all protocols
2. **Compute Query DRFP**: Generates DRFP fingerprint for your reaction
3. **Filter Candidates**: 
   - Optional family filter (e.g., only "Suzuki" protocols)
   - Optional tag filter (e.g., only protocols with "Pd" tag)
   - **SMARTS structural matching** (default: enabled)
4. **Rank by Similarity**: Computes cosine similarity between query and candidates
5. **Return Top-K**: Returns the most similar protocols

### SMARTS Filtering

By default, the CLI uses **both DRFP similarity AND SMARTS structural matching**:

- **DRFP**: Fuzzy similarity (87% similar)
- **SMARTS**: Exact structural requirements (e.g., aryl bromide + boronic acid)

You can disable SMARTS filtering with `--no-smarts-filter` to rely purely on DRFP similarity.

## Integration with Other Tools

### Python API

```python
from chemtools.protocol import ProtocolRecommender

recommender = ProtocolRecommender()
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=5,
    reaction_family='Suzuki',
    use_smarts_filter=True
)

for rec in results['recommended_conditions']:
    print(f"Rank {rec['rank']}: {rec['protocol_metadata']['title']}")
```

### FastAPI Endpoint (if integrated)

```bash
curl -X POST "http://localhost:8000/api/v1/protocol/recommend" \
  -H "Content-Type: application/json" \
  -d '{
    "reaction_smiles": "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    "k": 5
  }'
```

## Troubleshooting

### No Results Found

```bash
# Check if protocols are indexed
python -m chemtools.protocol.cli stats

# Try without SMARTS filtering
python -m chemtools.protocol.recommend_cli "YOUR_SMILES" --no-smarts-filter

# Lower similarity threshold
python -m chemtools.protocol.recommend_cli "YOUR_SMILES" --min-similarity 0.0
```

### Index Not Found

```bash
# Rebuild the index
python -m chemtools.protocol.cli build --force
```

### Invalid SMILES

Make sure your reaction SMILES:
- Uses `>>` to separate reactants and products
- Has valid SMILES on both sides
- Example: `CCBr.c1ccccc1B(O)O>>CCc1ccccc1`

## Examples by Reaction Type

### Suzuki Coupling

```bash
python -m chemtools.protocol.recommend_cli "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" --family Suzuki
```

### Buchwald-Hartwig Amination

```bash
python -m chemtools.protocol.recommend_cli "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --family Buchwald
```

### Sonogashira Coupling

```bash
python -m chemtools.protocol.recommend_cli "Ic1ccccc1.C#C>>c1ccccc1C#C" --family Sonogashira
```

### Cross-Electrophile Coupling

```bash
python -m chemtools.protocol.recommend_cli "CI.C(Cl)=O>>CC(=O)C" --tags "Ni,cross-electrophile"
```

## See Also

- [Protocol Database Tools README](readme.md) - Full documentation
- [SMARTS Fixes Summary](../../docs/SMARTS_FIXES_SUMMARY.md) - Common issues
- [Protocol Validation Tool](../../docs/PROTOCOL_VALIDATION_TOOL.md) - Validation guide

---

**Created**: October 20, 2025  
**Module**: `chemtools.protocol.recommend_cli`
