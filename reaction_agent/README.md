# Reaction SMILES Analysis Agent

**Self-contained POC module** for analyzing reaction SMILES using deterministic cheminformatics + LLM interpretation.

## Overview

This module combines:
- **Deterministic analysis** (RDKit + rxnmapper) for ground truth facts
- **LLM interpretation** (GPT-4/GPT-5) for mechanistic reasoning and classification

**Design principle**: Tools compute facts; GPT explains facts.

## Quick Start

### Interactive CLI (Recommended)

```bash
# Set API key
export OPENAI_API_KEY='sk-...'

# Start interactive mode
python reaction_agent/cli.py

# Or analyze directly
python reaction_agent/cli.py --reaction "CCBr>>CCN"

# Batch processing
python reaction_agent/cli.py --batch reactions.txt

# No LLM (deterministic only)
python reaction_agent/cli.py --reaction "CCBr>>CCN" --no-llm
```

**See [CLI_GUIDE.md](docs/CLI_GUIDE.md) for complete CLI documentation.**

### Quantitative Validation (NEW!)

Validate reliability beyond LLM self-reported confidence:

```python
from reaction_agent.scripts.quantitative_validation import validate_reaction

# Test with multiple models for cross-validation
validation = validate_reaction(
    rxn_smiles="Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1",
    models=['gpt-4o-mini', 'gpt-4o', 'o3-mini']
)

print(f"Overall Validated Score: {validation['overall_score']:.3f}")
print(f"Reliability: {validation['reliability']}")  # HIGH/MEDIUM/LOW/VERY_LOW
print(f"Recommendation: {validation['recommendation']}")
```

**Key Features**:
- ✅ **Deterministic quality score** (objective, tool-based)
- ✅ **Cross-model consistency** (multiple models agree?)
- ✅ **Specificity analysis** (detailed = more reliable)
- ✅ **Ensemble confidence** (weighted voting)
- ✅ **Comprehensive reliability rating**

**Quick Check** (single metric - atom mapping quality):
```bash
python reaction_agent/cli.py --reaction "..." --no-llm
```
Mapping confidence >0.8 = reliable, <0.6 = questionable.

**See [VALIDATION_QUICK_REFERENCE.md](docs/VALIDATION_QUICK_REFERENCE.md) for practical guide.**

### LLM-Assisted Mapping (NEW! Experimental)

When rxnmapper fails or gives low confidence, use LLM reasoning to improve reliability:

```python
from reaction_agent.scripts.llm_assisted_mapping import hybrid_mapping_workflow

# Hybrid workflow: rxnmapper → LLM validation/analysis if needed
result = hybrid_mapping_workflow(
    rxn_smiles="complex_tandem_reaction",
    confidence_threshold=0.6
)

print(f"Final confidence: {result['final_confidence']:.3f}")
if 'llm_analysis' in result:
    # LLM identified mechanism and mapping errors
    analysis = result['llm_analysis']['llm_analysis']
    print(f"Reaction type: {analysis['reaction_analysis']['type']}")
    print(f"Stages: {len(analysis['reaction_analysis']['stages'])}")
```

**Key Benefits**:
- ✅ **Validates rxnmapper** when confidence is borderline (0.4-0.7)
- ✅ **Analyzes complex reactions** when rxnmapper fails (<0.4)
- ✅ **Identifies mapping errors** with mechanistic reasoning
- ✅ **Suggests corrections** for manual review

**Test Results**:
- Simple reactions: rxnmapper only (fast, $0)
- Complex tandem: LLM correctly identified 2-stage mechanism + 4 mapping errors
- Cost: ~$0.003-0.006 per borderline/complex reaction

**See [LLM_MAPPING_SUMMARY.md](docs/LLM_MAPPING_SUMMARY.md) for details and examples.**

### Python API

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
├── README.md             # This file
├── __init__.py           # Main exports
├── cli.py                # ⭐ Interactive CLI application
├── core.py               # Deterministic analysis (RDKit, rxnmapper)
├── prompts.py            # LLM prompt templates
├── agent.py              # LLM integration and orchestration
├── examples/
│   ├── sample_reactions.txt  # Sample reactions for testing
│   └── llm_test_rxn.csv      # Complex reactions for validation
├── tests/
│   └── test_core.py      # Unit tests (✅ 12/12 passing)
├── scripts/
│   ├── demo.py                       # Demo with example reactions
│   ├── compare_models.py             # Model comparison for reactions
│   ├── test_workflows.py             # Workflow testing framework
│   ├── quantitative_validation.py    # ⭐ Validation framework
│   ├── test_validation_comparison.py # Compare simple vs complex reactions
│   └── llm_assisted_mapping.py       # ⭐ NEW: LLM-assisted atom mapping
├── results/
│   ├── WORKFLOW_TESTING_RESULTS.md  # Comprehensive testing results
│   ├── QUICK_SUMMARY.md             # Visual summary
│   └── validation_example.json      # Example validation output
└── docs/
    ├── CLI_GUIDE.md                    # ⭐ Complete CLI documentation
    ├── README.md                       # Detailed API documentation
    ├── IMPLEMENTATION_SUMMARY.md       # Implementation details
    ├── VALIDATION_GUIDE.md             # ⭐ Validation metrics explained
    ├── VALIDATION_QUICK_REFERENCE.md   # ⭐ Quick validation guide
    ├── LLM_ASSISTED_MAPPING.md         # ⭐ NEW: LLM mapping guide
    └── LLM_MAPPING_SUMMARY.md          # ⭐ NEW: Quick mapping summary
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

### ✅ Quantitative Validation (NEW!)
- **Multiple independent metrics** beyond LLM confidence
- **Deterministic quality score** (35% weight) - most objective
- **Cross-model consistency** (25% weight) - multiple models agree?
- **Specificity analysis** (20% weight) - detailed = reliable
- **Ensemble confidence** - weighted voting from models
- **Reliability ratings**: HIGH/MEDIUM/LOW/VERY_LOW
- **Practical recommendations** for production use

### ✅ Model Testing & Optimization
- **Workflow comparison** across 5 model categories
- **Performance benchmarks** for complex reactions
- **Cost-quality-speed analysis**
- **Model recommendations** by reaction type

### ✅ LLM-Assisted Mapping (NEW! Experimental)
- **Validates rxnmapper** when confidence is borderline
- **Analyzes complex reactions** using LLM reasoning (o3, o3-mini)
- **Identifies mapping errors** via mechanistic understanding
- **Suggests corrections** for manual review
- **Hybrid workflow**: rxnmapper first, LLM helps when needed
- Successfully identified 2-stage mechanism + 4 mapping errors in test

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
- **CLI_GUIDE.md** - Complete interactive CLI guide
- **README.md** - Complete usage guide and API reference
- **IMPLEMENTATION_SUMMARY.md** - Implementation details and metrics
- **VALIDATION_GUIDE.md** - ⭐ Detailed validation metrics explanation
- **VALIDATION_QUICK_REFERENCE.md** - ⭐ Quick reference for practical validation
- **LLM_ASSISTED_MAPPING.md** - ⭐ **NEW: LLM-assisted atom mapping guide**
- **LLM_MAPPING_SUMMARY.md** - ⭐ **NEW: Quick mapping summary with examples**

See `reaction_agent/results/`:
- **WORKFLOW_TESTING_RESULTS.md** - Comprehensive model comparison on complex reactions
- **QUICK_SUMMARY.md** - Visual summary of model performance

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
