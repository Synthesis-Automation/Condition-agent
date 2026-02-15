# Reaction SMILES Analysis Agent - POC

**Status**: Proof of Concept (v1)
**Design**: Based on `docs/reaction_smiles_analysis_agent_simple_v1.md`

## Overview

This POC demonstrates a hybrid approach to reaction analysis that combines:

1. **Deterministic cheminformatics** (RDKit + rxnmapper) for "ground truth" facts
2. **LLM interpretation** (GPT-4) for mechanistic reasoning and classification

**Design principle**: Tools compute facts; GPT explains facts. If tools can't compute reliable facts, we return warnings and low confidence rather than inventing.

## Architecture

```
┌─────────────────────────┐
│  Reaction SMILES        │
│  (reactants>>products)  │
└───────────┬─────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────┐
│  chemtools.reaction_smiles_analyzer (Deterministic)         │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐  │
│  │   Cleaning   │→ │   Mapping    │→ │  Bond Changes    │  │
│  │  (RDKit)     │  │ (rxnmapper)  │  │  (diff analysis) │  │
│  └──────────────┘  └──────────────┘  └──────────────────┘  │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           │ tool_facts
                           ▼
┌─────────────────────────────────────────────────────────────┐
│  llmtools.reaction_smiles_agent (LLM Interpretation)        │
│  ┌──────────────────┐  ┌──────────────────────────────┐    │
│  │  Format Prompt   │→ │  GPT Call (Structured JSON)  │    │
│  │  (with facts)    │  │  + QC Gating                 │    │
│  └──────────────────┘  └──────────────────────────────┘    │
└──────────────────────────┬──────────────────────────────────┘
                           │
                           ▼
                    ┌────────────┐
                    │   Result   │
                    │   (JSON)   │
                    └────────────┘
```

## Components

### 1. Core Deterministic Module

**File**: `chemtools/reaction_smiles_analyzer.py`

Functions:
- `clean_reaction_smiles()` - Canonicalize and detect spectators
- `map_reaction()` - Atom mapping via rxnmapper
- `extract_bond_changes()` - Bond-level diff from mapped reaction
- `analyze_reaction_smiles()` - Full deterministic pipeline

Key features:
- Conservative cleanup (preserves raw input)
- Spectator detection (common ions/salts)
- QC gating (warns if mapping fails)
- No guessing - returns warnings instead

### 2. LLM Integration

**File**: `llmtools/reaction_smiles_agent.py`

Classes/Functions:
- `ReactionSMILESAnalyzer` - Agent class
- `analyze_reaction_smiles()` - Function interface

Key features:
- Structured JSON output from LLM
- Markdown fence stripping
- QC-based confidence gating
- Detailed metadata tracking

### 3. Prompt Template

**File**: `llmtools/prompts.py`

Template: `REACTION_SMILES_ANALYSIS`

Key constraints:
- "Use ONLY provided tool facts"
- "Do NOT invent reagents or conditions"
- "Reference bond changes by ID"
- Strict JSON schema with enums

## Installation

### Requirements

```bash
# Core dependencies (should already be installed)
pip install rdkit openai

# Optional but recommended for full functionality
pip install rxnmapper
```

### Environment Setup

```bash
# Set OpenAI API key
export OPENAI_API_KEY='your-key-here'
```

## Usage

### Simple Example

```python
from llmtools import LLMClient, ReactionSMILESAnalyzer

# Initialize
client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)

# Analyze
result = analyzer.analyze("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")

# Access results
print(result["interpretation"]["overall_class"])  # e.g., "cross_coupling"
print(result["interpretation"]["tags"])            # e.g., ["Buchwald-Hartwig"]
print(result["interpretation"]["confidence"])      # e.g., 0.85
```

### Functional Interface

```python
from llmtools import LLMClient, analyze_reaction_smiles

client = LLMClient(provider="openai", model="gpt-4o-mini")

result = analyze_reaction_smiles(
    rxn_smiles="CCBr>>CCN",
    client=client,
    drop_spectators=True,
    temperature=0.0
)
```

### Batch Processing

```python
reactions = [
    "CCBr>>CCN",
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1",
    # ... more reactions
]

results = analyzer.analyze_batch(reactions)
```

### Test Script

Run the included test script:

```bash
python scripts/test_reaction_smiles_agent.py
```

## Output Schema

```json
{
  "schema_version": "reaction_analysis.v1",
  "input": {
    "rxn_smiles_raw": "...",
    "rxn_smiles_clean": "...",
    "reactants_clean": [...],
    "products_clean": [...],
    "spectators": [...],
    "parse_warnings": [...],
    "standardization_actions": [...]
  },
  "tool_facts": {
    "mapped_rxn_smiles": "...",
    "mapping_qc": {
      "ok": true,
      "confidence": 0.95,
      "notes": []
    },
    "bond_changes": [
      {
        "id": "BC1",
        "change": "broken",
        "a1": 12,
        "a2": 57,
        "bond": "single"
      },
      {
        "id": "BC2",
        "change": "formed",
        "a1": 12,
        "a2": 89,
        "bond": "single"
      }
    ],
    "reaction_center_atoms": [12, 57, 89]
  },
  "interpretation": {
    "overall_class": "nucleophilic_substitution",
    "tags": ["SNAr"],
    "roles": {
      "electrophile": "aromatic chloride",
      "nucleophile": "dimethylamine",
      "leaving_group": "chloride"
    },
    "events": [
      {
        "event_id": "E1",
        "event_type": "substitution",
        "bond_change_refs": ["BC1", "BC2"],
        "short_rationale": "C–Cl cleavage followed by C–N formation",
        "confidence": 0.85
      }
    ],
    "mechanism_summary": [
      "Nucleophilic attack by amine on electron-deficient aromatic carbon",
      "Displacement of chloride leaving group"
    ],
    "warnings": [],
    "confidence": 0.78
  },
  "metadata": {
    "model": "gpt-4o-mini",
    "provider": "openai",
    "total_tokens": 1245,
    "latency_ms": 2340
  }
}
```

## Reliability Features

### 1. Conservative Cleanup
- Keeps raw input unchanged
- Minimal standardization (canonicalization only)
- Logs all modifications

### 2. Spectator Handling
- Detects common ions/salts (Cl⁻, Na⁺, etc.)
- Removes from mapping for cleaner analysis
- Preserves in report for transparency

### 3. QC Gating
Never overclaims confidence:
- If mapping fails → confidence capped at 0.3
- If mapping weak → warnings added
- If no bond changes → degraded analysis

### 4. No Hallucination
- LLM prompted to use ONLY tool facts
- No invention of reagents/conditions
- Returns warnings over guesses

## Limitations (POC)

1. **No condition prediction** - Only analyzes SMILES, not reagents/conditions
2. **Simple spectator list** - Could be expanded
3. **Single-step focus** - Multi-step reactions get basic analysis
4. **No stereochemistry analysis** - Preserved but not interpreted
5. **rxnmapper required** - Mapping fails gracefully but reduces quality

## Future Enhancements

See design doc section 9 for planned improvements:
- Fallback mode for parse/mapping failures
- Extended spectator dictionary (OTf, TsO, etc.)
- Multi-stage analysis for tandem reactions
- Structured output via OpenAI API schema enforcement (for GPT-5+)

## Testing

Run the test script with different models:

```python
# Fast and cheap
client = LLMClient(provider="openai", model="gpt-4o-mini")

# More capable
client = LLMClient(provider="openai", model="gpt-4o")

# For complex reactions (if needed)
client = LLMClient(provider="openai", model="gpt-4-turbo")
```

## Troubleshooting

### rxnmapper not installed
```
pip install rxnmapper
```

If installation fails, the analyzer will still work but skip atom mapping.

### API key issues
```bash
# Check if key is set
echo $OPENAI_API_KEY

# Set it
export OPENAI_API_KEY='sk-...'
```

### JSON parsing errors
The agent handles common issues:
- Strips markdown fences (```json ... ```)
- Reports parsing errors instead of crashing
- Returns raw response for debugging

### Low confidence warnings
This is expected when:
- Mapping confidence < 0.5
- Complex/unusual reactions
- Insufficient structural info

Lower confidence = more honest uncertainty!

## Development

### Adding new spectators

Edit `chemtools/reaction_smiles_analyzer.py`:

```python
SPECTATORS = {
    "Cl", "Br", "I", "F",
    "[OTf-]",  # Add triflate
    "[BF4-]",  # Add tetrafluoroborate
    # ... etc
}
```

### Customizing prompt

Edit `llmtools/prompts.py` → `REACTION_SMILES_ANALYSIS`

### Extending output schema

1. Update prompt template with new fields
2. Document in this README
3. Update test assertions

## References

- Design doc: `docs/reaction_smiles_analysis_agent_simple_v1.md`
- Repository guidelines: `AGENTS.md`
- RDKit docs: https://www.rdkit.org/docs/
- rxnmapper: https://github.com/rxn4chemistry/rxnmapper
