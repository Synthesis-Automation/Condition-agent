# Protocol Module - Standard JSON Output Format

## Overview

The Protocol Recommendation module now outputs **standard JSON format** consistent with other ChemTools recommendation systems (ML-based, Rule-based, Fusion).

This ensures:
- ✅ Uniform API responses across all recommendation modes
- ✅ Easy integration with existing ChemTools workflows
- ✅ Consistent client code for parsing results
- ✅ Compatible with `chemtools.output_formatter` utilities

## Standard Output Structure

```json
{
  "meta": { ... },
  "input": { ... },
  "detection": { ... },
  "recommended_conditions": [ ... ],
  "extras": { ... }
}
```

### 1. Meta Section

Model metadata and processing information:

```json
{
  "meta": {
    "generated_at": "2025-10-12T03:09:15.413034Z",
    "model": "Protocol-DRFP",
    "model_version": "1.0.0",
    "status": "success",
    "version": "2.0",
    "processing_time_ms": 1702.6
  }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `generated_at` | string | ISO 8601 timestamp |
| `model` | string | Model identifier ("Protocol-DRFP") |
| `model_version` | string | Version of recommendation engine |
| `status` | string | "success" or "error" |
| `version` | string | Schema version |
| `processing_time_ms` | float | Processing time in milliseconds |

### 2. Input Section

Query information and options:

```json
{
  "input": {
    "reaction_smiles": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    "requested_reaction_type": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "options": {
      "k": 3,
      "tags": ["suzuki", "palladium"]
    }
  }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `reaction_smiles` | string | Input reaction SMILES |
| `requested_reaction_type` | string | Optional requested family filter |
| `options` | object | Query options (k, tags, etc.) |

### 3. Detection Section

Reaction type detection results:

```json
{
  "detection": {
    "family": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "detected_reaction_type": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "method": "protocol-similarity",
    "confidence": 0.8018
  }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `family` | string | Detected reaction family |
| `detected_reaction_type` | string | Detected type (same as family) |
| `method` | string | Detection method ("protocol-similarity") |
| `confidence` | float | Confidence score (0.0-1.0), based on top match similarity |

### 4. Recommended Conditions (Main Results)

Array of protocol recommendations, ordered by similarity:

```json
{
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.8018,
      "protocol": {
        "filename": "Suzuki_Cu_C(sp3)-C(sp2).json",
        "title": "Copper-Catalyzed Suzuki-Miyaura Coupling...",
        "journal": "Organic Syntheses",
        "year": 2025,
        "doi": "10.15227/orgsyn.102.0086",
        "url": "https://www.orgsyn.org/demo.aspx?prep=v102p0086",
        "reaction_smiles": "BrC1CCCCC1.c1ccc(Cl)cc1B(...)>>...",
        "reaction_family": "Suzuki_Cu_alkyl_halide+aryl_boron",
        "tags": ["Suzuki", "Cu", "alkyl_halide", "Coupling"],
        "notes": "CuBr·SMe2 (5 mol%), bathophenanthroline (7.5 mol%)..."
      },
      "similarity": 0.8018,
      "source": "protocol_database"
    }
  ]
}
```

#### Recommendation Entry Fields

| Field | Type | Description |
|-------|------|-------------|
| `rank` | int | Rank position (1, 2, 3, ...) |
| `confidence` | float | Similarity score (same as `similarity`) |
| `protocol` | object | Full protocol metadata |
| `similarity` | float | DRFP similarity score (0.0-1.0) |
| `source` | string | Always "protocol_database" |

#### Protocol Object Fields

| Field | Type | Description |
|-------|------|-------------|
| `filename` | string | Protocol JSON filename |
| `title` | string | Full protocol title |
| `journal` | string | Source journal |
| `year` | int | Publication year |
| `doi` | string | DOI reference |
| `url` | string | Protocol URL |
| `reaction_smiles` | string | Protocol reaction SMILES |
| `reaction_family` | string | Reaction family classification |
| `tags` | array | Protocol tags |
| `notes` | string | Quick reference notes |

### 5. Extras Section

Additional search metadata:

```json
{
  "extras": {
    "num_candidates": 16,
    "num_total_protocols": 16,
    "num_matches": 3
  }
}
```

| Field | Type | Description |
|-------|------|-------------|
| `num_candidates` | int | Number of protocols after filtering |
| `num_total_protocols` | int | Total protocols in index |
| `num_matches` | int | Number of results returned |

## With Condition Extraction

When using `recommend_with_details()`, each recommendation gets a `conditions` field:

```json
{
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.8018,
      "protocol": { ... },
      "similarity": 0.8018,
      "source": "protocol_database",
      "conditions": {
        "catalyst": "CuBr·SMe2",
        "ligand": "Bathophenanthroline",
        "base": "NaOtBu",
        "solvent": "toluene",
        "additives": [],
        "temperature_C": 80.0,
        "time_h": 24.0,
        "atmosphere": "air"
      }
    }
  ]
}
```

## Comparison with Other Modes

| Feature | ML-based | Rule-based | Protocol-DRFP |
|---------|----------|------------|---------------|
| Model name | `ML-precedent-knn` | `Rule-based-SCDB` | `Protocol-DRFP` |
| Detection method | `rxn-insight-ml` | `pattern-match` | `protocol-similarity` |
| Recommendations | Condition sets from precedents | Condition sets from rules | Full protocols with extracted conditions |
| Confidence | ML model score | 1.0 (if matched) | DRFP similarity score |
| Source | Reaction database | Rule database | Protocol database |

## Python Usage Examples

### Basic Recommendation

```python
from chemtools.protocol import ProtocolRecommender

recommender = ProtocolRecommender()

results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

# Access standard fields
print(f"Model: {results['meta']['model']}")
print(f"Status: {results['meta']['status']}")
print(f"Processing time: {results['meta']['processing_time_ms']:.2f} ms")

print(f"\nDetected family: {results['detection']['family']}")
print(f"Confidence: {results['detection']['confidence']:.3f}")

print(f"\nRecommendations:")
for rec in results['recommended_conditions']:
    print(f"  Rank {rec['rank']}: {rec['protocol']['title']}")
    print(f"    Similarity: {rec['similarity']:.3f}")
    print(f"    DOI: {rec['protocol']['doi']}")
```

### With Condition Extraction

```python
results = recommender.recommend_with_details(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

for rec in results['recommended_conditions']:
    protocol = rec['protocol']
    conditions = rec.get('conditions', {})
    
    print(f"\nRank {rec['rank']}: {protocol['title']}")
    print(f"  Similarity: {rec['similarity']:.3f}")
    print(f"  Catalyst: {conditions.get('catalyst')}")
    print(f"  Ligand: {conditions.get('ligand')}")
    print(f"  Base: {conditions.get('base')}")
    print(f"  Solvent: {conditions.get('solvent')}")
    print(f"  Temperature: {conditions.get('temperature_C')} °C")
    print(f"  Time: {conditions.get('time_h')} h")
```

### Error Handling

```python
results = recommender.recommend(
    reaction_smiles='INVALID>>SMILES',
    k=3
)

if results['meta']['status'] == 'error':
    print(f"Error: {results['extras']['error']}")
elif len(results['recommended_conditions']) == 0:
    print("No matching protocols found")
    print(f"Searched {results['extras']['num_candidates']} candidates")
else:
    print(f"Found {results['extras']['num_matches']} matches")
```

## Legacy Format Support

For backward compatibility, the legacy format can still be accessed:

```python
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3,
    use_standard_format=False  # Returns legacy format
)

# Legacy format has 'matches' instead of 'recommended_conditions'
for match in results['matches']:
    print(f"{match['similarity']:.3f}: {match['source_title']}")
```

## Integration with ChemTools

The standard format integrates seamlessly with other ChemTools components:

```python
from chemtools.output_formatter import ensure_standard_output

# Protocol results are already in standard format
protocol_results = recommender.recommend(...)

# Can be used directly with formatters
print(protocol_results['meta']['model'])  # "Protocol-DRFP"
print(protocol_results['detection']['family'])
```

## Migration Guide

If you have code using the old format:

### Old Code (Legacy Format)

```python
results = recommender.recommend(...)

for match in results['matches']:
    print(match['similarity'])
    print(match['source_title'])
    print(match['source_doi'])
```

### New Code (Standard Format)

```python
results = recommender.recommend(...)

for rec in results['recommended_conditions']:
    protocol = rec['protocol']
    print(rec['similarity'])  # or rec['confidence']
    print(protocol['title'])
    print(protocol['doi'])
```

## Summary

✅ **Unified Format**: Protocol recommendations now match ML/Rule/Fusion output structure  
✅ **Easy Integration**: Compatible with existing ChemTools workflows and utilities  
✅ **Backward Compatible**: Legacy format still available with `use_standard_format=False`  
✅ **Rich Metadata**: Includes timing, detection confidence, and search statistics  
✅ **Clear Structure**: Well-defined sections (meta, input, detection, recommended_conditions, extras)
