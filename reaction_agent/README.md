# Reaction SMILES Analysis Agent

**Self-contained POC module** for analyzing reaction SMILES using deterministic cheminformatics + LLM interpretation.

## Overview

This module combines:
- **Deterministic analysis** (RDKit + rxnmapper) for ground truth facts
- **LLM interpretation** (GPT-4) for mechanistic reasoning and classification

**Design principle**: Tools compute facts; GPT explains facts.

## Quick Start

```python
from llmtools.clients import LLMClient
from reaction_agent import ReactionSMILESAnalyzer

# Initialize
client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)

# Analyze
result = analyzer.analyze("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")

# Access results
print(result["interpretation"]["overall_class"])  # "cross_coupling"
print(result["interpretation"]["confidence"])      # 0.85
```

## Module Structure

```
reaction_agent/
├── __init__.py           # Main exports
├── core.py               # Deterministic analysis (RDKit, rxnmapper)
├── prompts.py            # LLM prompt templates
├── agent.py              # LLM integration and orchestration
├── tests/
│   └── test_core.py      # Unit tests (12 passing)
├── scripts/
│   └── demo.py           # Demo with example reactions
└── docs/
    ├── README.md         # Detailed documentation
    └── IMPLEMENTATION_SUMMARY.md
```

## Installation

```bash
# Core dependencies
pip install rdkit openai

# Optional (for full functionality)
pip install rxnmapper
```

## Usage

### Basic Analysis

```python
from reaction_agent import ReactionSMILESAnalyzer
from llmtools.clients import LLMClient

client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)

result = analyzer.analyze("CCBr>>CCN")
```

### Functional Interface

```python
from reaction_agent import analyze_reaction_smiles
from llmtools.clients import LLMClient

client = LLMClient(provider="openai", model="gpt-4o")

result = analyze_reaction_smiles(
    rxn_smiles="CCBr>>CCN",
    client=client,
    temperature=0.0
)
```

### Deterministic Only (No LLM)

```python
from reaction_agent import analyze_deterministic

result = analyze_deterministic("CCBr>>CCN")
# Returns: input, tool_facts (no LLM interpretation)
```

## Running the Demo

```bash
export OPENAI_API_KEY='your-key-here'
python reaction_agent/scripts/demo.py
```

## Running Tests

```bash
# All tests
pytest reaction_agent/tests/

# Specific test
pytest reaction_agent/tests/test_core.py -v

# Skip tests requiring rxnmapper
pytest reaction_agent/tests/ -m "not skipif"
```

## Output Schema

```json
{
  "schema_version": "reaction_analysis.v1",
  "input": {
    "rxn_smiles_raw": "...",
    "rxn_smiles_clean": "...",
    "spectators": [...],
    "parse_warnings": [...]
  },
  "tool_facts": {
    "mapped_rxn_smiles": "...",
    "mapping_qc": {"ok": true, "confidence": 0.95},
    "bond_changes": [
      {"id": "BC1", "change": "broken", "a1": 12, "a2": 57, "bond": "single"},
      {"id": "BC2", "change": "formed", "a1": 12, "a2": 89, "bond": "single"}
    ],
    "reaction_center_atoms": [12, 57, 89]
  },
  "interpretation": {
    "overall_class": "nucleophilic_substitution",
    "tags": ["SNAr"],
    "roles": {...},
    "events": [...],
    "mechanism_summary": [...],
    "confidence": 0.78
  },
  "metadata": {
    "model": "gpt-4o-mini",
    "total_tokens": 1245,
    "latency_ms": 2340
  }
}
```

## Key Features

### ✅ Deterministic Foundation
- RDKit canonicalization
- Spectator detection (Cl⁻, Na⁺, etc.)
- rxnmapper atom mapping
- Bond-level diff analysis

### ✅ LLM Integration
- Structured JSON output
- Confidence-based QC gating
- No hallucination (uses only tool facts)
- Detailed mechanistic interpretation

### ✅ Reliability
- Conservative cleanup (preserves raw input)
- Warnings over guesses
- Graceful degradation on failures
- Comprehensive metadata

### ✅ Self-Contained
- No modifications to existing codebase
- All code in `/reaction_agent` folder
- Independent tests and documentation

## Components

### `core.py` - Deterministic Analysis

Functions:
- `clean_reaction_smiles()` - Canonicalize and detect spectators
- `map_reaction()` - Atom mapping via rxnmapper
- `extract_bond_changes()` - Bond-level diff
- `analyze_deterministic()` - Full deterministic pipeline

Data classes:
- `CleanReport`, `MappingReport`, `BondChangeReport`
- `BondChange`

### `prompts.py` - LLM Templates

- `REACTION_SMILES_ANALYSIS` - Main prompt template
- `get_template()` - Template accessor

### `agent.py` - LLM Integration

- `ReactionSMILESAnalyzer` - Agent class
- `analyze_reaction_smiles()` - Main function
- QC gating and confidence adjustment

## Documentation

See `reaction_agent/docs/`:
- **README.md** - Complete usage guide and API reference
- **IMPLEMENTATION_SUMMARY.md** - Implementation details and metrics

Design specification: `docs/reaction_smiles_analysis_agent_simple_v1.md`

## Troubleshooting

### rxnmapper not installed
```bash
pip install rxnmapper
```
Analyzer will work without it but skip atom mapping.

### API key issues
```bash
export OPENAI_API_KEY='sk-...'
```

### Import errors
Make sure project root is in your Python path:
```python
import sys
sys.path.insert(0, '/path/to/Condition-agent')
```

## Limitations (POC)

- No condition prediction (only SMILES analysis)
- Simple spectator list (expandable)
- Single-step focus (multi-step gets basic analysis)
- No stereochemistry interpretation
- rxnmapper optional (graceful fallback)

## Future Enhancements

- Extended spectator dictionary (OTf, TsO, etc.)
- Multi-stage analysis for tandem reactions
- OpenAI structured output API (GPT-5+)
- Batch processing optimizations
- API endpoints (FastAPI)

## License

Part of the Condition-agent project.

## Version

v0.1.0 (POC)
