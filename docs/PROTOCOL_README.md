# Protocol Recommendation Module - README

## Quick Summary

The **Protocol Recommendation Module** finds the most relevant experimental protocol for your reaction using DRFP similarity search.

### Key Features

- ✅ **DRFP-based similarity**: Uses molecular fingerprints to find similar reactions
- ✅ **Standard JSON output**: Compatible with ML-based and Rule-based recommendation outputs
- ✅ **Fast indexing**: Precomputes DRFP fingerprints for instant search
- ✅ **Flexible filtering**: Filter by reaction family and tags
- ✅ **Condition extraction**: Automatically extracts catalyst, solvent, temperature, time, etc.
- ✅ **CLI tools**: Command-line interface for index management and testing

### What's New: Standard JSON Output! 🎉

The protocol module now outputs **standard format** matching other ChemTools recommendation systems:

```json
{
  "meta": {
    "model": "Protocol-DRFP",
    "status": "success",
    "processing_time_ms": 1702.6,
    ...
  },
  "input": {
    "reaction_smiles": "...",
    "options": {"k": 3}
  },
  "detection": {
    "family": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "confidence": 0.8018,
    "method": "protocol-similarity"
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.8018,
      "protocol": {
        "filename": "...",
        "title": "...",
        "doi": "...",
        ...
      },
      "similarity": 0.8018,
      "source": "protocol_database"
    }
  ],
  "extras": {
    "num_candidates": 16,
    "num_total_protocols": 16
  }
}
```

**Benefits:**
- Same structure as ML-based and Rule-based recommendations
- Easy to integrate into existing workflows
- Compatible with `chemtools.output_formatter` utilities
- Clear separation of concerns (meta, input, detection, results)

## Quick Start

### 1. Build the Index

```bash
python -m chemtools.protocol.cli build
```

This scans `data/protocol_db/*.json` and builds `.protocol_index.json` with DRFP fingerprints.

### 2. Use in Python

```python
from chemtools.protocol import ProtocolRecommender

# Initialize
recommender = ProtocolRecommender()

# Get recommendations
results = recommender.recommend(
    reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
    k=3
)

# Standard format
for rec in results['recommended_conditions']:
    protocol = rec['protocol']
    print(f"Rank {rec['rank']}: {protocol['title']}")
    print(f"  Similarity: {rec['similarity']:.3f}")
    print(f"  DOI: {protocol['doi']}")
```

### 3. Test with Interactive CLI

```bash
# Windows PowerShell
.\run_protocol_cli.ps1

# Or with Python
python test_protocol_cli.py
```

Then enter reaction SMILES to get protocol recommendations.

## Documentation

- **[PROTOCOL_MODULE.md](PROTOCOL_MODULE.md)** - Complete technical documentation
- **[PROTOCOL_OUTPUT_FORMAT.md](PROTOCOL_OUTPUT_FORMAT.md)** - Detailed output format specification
- **[PROTOCOL_QUICKSTART.md](PROTOCOL_QUICKSTART.md)** - 2-minute quick start guide
- **[PROTOCOL_CLI_GUIDE.md](PROTOCOL_CLI_GUIDE.md)** - Interactive CLI user guide

## Architecture

### Components

1. **ProtocolIndexer** (`chemtools/protocol/indexer.py`)
   - Scans protocol JSON files
   - Computes DRFP fingerprints
   - Builds searchable index with MD5-based incremental updates

2. **ProtocolRecommender** (`chemtools/protocol/recommend.py`)
   - Loads precomputed index
   - Computes query DRFP
   - Finds top-k similar protocols using cosine similarity
   - Returns standard JSON format

3. **CLI** (`chemtools/protocol/cli.py`)
   - Build, stats, list-families, show-family commands
   - Index management

4. **Interactive Tester** (`test_protocol_cli.py`)
   - User-friendly testing interface
   - Single query mode and interactive mode

### Data Flow

```
Protocol JSON Files
        ↓
   [Indexer]
        ↓
   DRFP Index (.protocol_index.json)
        ↓
   [Recommender]
        ↓
   Standard JSON Output
```

## Example Output

### Basic Recommendation

```python
{
  "meta": {
    "generated_at": "2025-10-12T03:09:15.413034Z",
    "model": "Protocol-DRFP",
    "model_version": "1.0.0",
    "status": "success",
    "version": "2.0",
    "processing_time_ms": 1702.6
  },
  "input": {
    "reaction_smiles": "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1",
    "options": {"k": 3}
  },
  "detection": {
    "family": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "detected_reaction_type": "Suzuki_Cu_alkyl_halide+aryl_boron",
    "method": "protocol-similarity",
    "confidence": 0.8018
  },
  "recommended_conditions": [
    {
      "rank": 1,
      "confidence": 0.8018,
      "protocol": {
        "filename": "Suzuki_Cu_C(sp3)-C(sp2).json",
        "title": "Copper-Catalyzed Suzuki-Miyaura Coupling of Unactivated Alkyl Halides with Arylborons",
        "journal": "Organic Syntheses",
        "year": 2025,
        "doi": "10.15227/orgsyn.102.0086",
        "url": "https://www.orgsyn.org/demo.aspx?prep=v102p0086",
        "reaction_smiles": "BrC1CCCCC1.c1ccc(Cl)cc1B(OC(C)(C)C)OC(C)(C)C>>Clc1ccc(C2CCCCC2)cc1",
        "reaction_family": "Suzuki_Cu_alkyl_halide+aryl_boron",
        "tags": ["Suzuki", "Cu", "alkyl_halide", "Coupling"],
        "notes": "CuBr·SMe2 (5 mol%), bathophenanthroline (7.5 mol%), NaOtBu (1.5 equiv), toluene (20 mL), 80 °C, 24 h, air"
      },
      "similarity": 0.8018,
      "source": "protocol_database"
    }
  ],
  "extras": {
    "num_candidates": 16,
    "num_total_protocols": 16,
    "num_matches": 1
  }
}
```

### With Extracted Conditions

```python
results = recommender.recommend_with_details(...)

# Each recommendation now has 'conditions' field:
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
    "temperature_C": 80.0,
    "time_h": 24.0,
    "atmosphere": "air"
  }
}
```

## Comparison: Protocol vs ML vs Rule

| Feature | ML-based | Rule-based | Protocol-DRFP |
|---------|----------|------------|---------------|
| **Output** | Condition sets | Condition sets | Full protocols |
| **Source** | Reaction database | Rule database | Protocol database |
| **Model** | `ML-precedent-knn` | `Rule-based-SCDB` | `Protocol-DRFP` |
| **Method** | `rxn-insight-ml` | `pattern-match` | `protocol-similarity` |
| **Confidence** | ML score | 1.0 (if matched) | DRFP similarity |
| **Format** | ✅ Standard JSON | ✅ Standard JSON | ✅ Standard JSON |

**All three modes now use the same standard output format!**

## Status

- ✅ Indexing: Complete and tested
- ✅ Recommendation: Complete and tested
- ✅ Standard output: Implemented and documented
- ✅ CLI tools: Complete
- ✅ Interactive tester: Complete
- ✅ Documentation: Complete
- ✅ Tests: All passing

**Current database**: 16 protocols indexed

## Future Enhancements

- [ ] Add API endpoint: `POST /api/v1/protocol/recommend`
- [ ] Integrate with existing condition recommendation workflow
- [ ] Add more protocols to database
- [ ] Full-text search of procedure text
- [ ] Protocol validation and schema checking

## Related Documentation

- Main README: `README.md`
- API Documentation: `docs/API_DOCUMENTATION.md`
- Output Formatter: `chemtools/output_formatter.py`

---

**Questions?** See the detailed documentation files or run `python test_protocol_cli.py` to try it interactively!
