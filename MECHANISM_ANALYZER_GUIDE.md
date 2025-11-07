# Reaction Mechanism Analyzer - Implementation Guide

## Overview

To add a **reaction mechanism analyzer** tool to the agent, you'll need to integrate or create capabilities for analyzing reaction mechanisms, including electron flow, intermediate steps, and transition states.

## Agent Integration Objectives

- **Primary user stories**
  - _Mechanism QA_: Chemist asks the agent why a product forms; the agent should call the analyzer, extract salient mechanism facts, and respond with natural-language reasoning that references concrete evidence.
  - _Feasibility checks_: During planning, the agent runs the analyzer to detect implausible steps (e.g., missing leaving groups) before suggesting routes.
  - _Teaching mode_: Agent can walk a learner through each elementary step, rendering optional visuals.
- **Contract-level requirements**
  - Input must accept `reaction_smiles`, optional `context` (temperature, solvent, catalyst), and `detail_level`.
  - Output should return `mechanism_type`, `confidence`, ordered `steps` (with intermediates + electron movement descriptors), and `warnings`.
  - Responses must include `evidence_refs` that point back to internal detectors (bond changes, functional groups, precedent lookup) so the agent can cite sources.
- **System integration**
  - Tool exposed in `chemtools_wrapper` with clear schema + JSONSchema metadata for LLM agents.
  - Wrapper should orchestrate: SMILES normalization → reaction family detection → mechanism classifier → electron flow + intermediates → natural-language summary template.
  - Add idempotent cache hooks so repeated calls with identical inputs reuse expensive computations (RXNMapper, DRFP similarities).

---

## 🎯 What a Mechanism Analyzer Should Do

### Core Capabilities

1. **Identify Reaction Center** - Detect which atoms/bonds change ✅ (Already available!)
2. **Electron Flow Analysis** - Track electron movement (arrow pushing)
3. **Intermediate Identification** - Predict reaction intermediates
4. **Transition State Estimation** - Estimate activation energy/TS structure
5. **Mechanism Classification** - Identify mechanism type (SN2, E2, radical, etc.)
6. **Step-by-Step Breakdown** - Break mechanism into elementary steps

---

## ✅ What We Already Have

### 1. Bond Analysis Tool (Already Integrated!)

Location: `chemtools._atom_mapping` and `chem_assistant/chemtools_wrapper.py`

**Tool**: `analyze_bond_changes_tool`

**Capabilities**:

- Identifies broken bonds
- Identifies formed bonds
- Detects leaving groups
- Provides confidence scores
- Uses RXNMapper + MCS hybrid approach

**Example Usage**:

```python
from chem_assistant.chemtools_wrapper import analyze_bond_changes_tool

result = analyze_bond_changes_tool.invoke({
    "reaction_smiles": "CCBr.c1ccccc1B(O)O>>CCc1ccccc1",
    "use_hybrid": True
})

# Returns:
# {
#   "broken_bonds": [{"atoms": [1, 7], "bond_type": "SINGLE"}],
#   "formed_bonds": [{"atoms": [1, 8], "bond_type": "SINGLE"}],
#   "leaving_groups": ["Br"],
#   "interpretation": "C-Br bond broken, C-C bond formed..."
# }
```

This provides the **foundation** for mechanism analysis!

---

## ❌ What's Missing for Full Mechanism Analysis

### 1. Electron Flow Tracking ⚠️ HIGH PRIORITY

**What**: Track how electrons move during the reaction (arrow pushing)

**Implementation Options**:

#### Option A: Rule-Based (Simpler)

Create SMARTS patterns for common mechanisms:

- Nucleophilic substitution (SN1/SN2)
- Elimination (E1/E2)
- Addition reactions
- Oxidation/reduction
- Radical mechanisms

**Example Pattern**:

```python
# SN2 mechanism
{
    "name": "SN2",
    "pattern": {
        "nucleophile": "[N,O,S,C-]",
        "electrophile": "[CX4][Cl,Br,I]",
        "leaving_group": "[Cl,Br,I]"
    },
    "electron_flow": [
        {"from": "nucleophile_lone_pair", "to": "C_sigma_bond"},
        {"from": "C_X_bond", "to": "leaving_group"}
    ]
}
```

#### Option B: ML-Based (Advanced)

Use machine learning to predict electron flow:

- Train on reaction databases with known mechanisms
- Use graph neural networks (GNNs) to predict electron paths
- Requires labeled training data

**Recommended**: Start with Option A (rule-based)

---

### 2. Intermediate Prediction ⚠️ MEDIUM PRIORITY

**What**: Identify reaction intermediates (carbocations, carbanions, radicals)

**Implementation**:

```python
def predict_intermediates(
    reaction_smiles: str,
    mechanism_type: str
) -> List[Dict[str, Any]]:
    """
    Predict intermediates based on mechanism type.

    Returns:
        List of intermediate structures with stability info
    """
    # For SN1: carbocation intermediate
    # For E1: carbocation intermediate
    # For radical: radical intermediate
    # etc.
```

**Dependencies**:

- RDKit for structure manipulation
- Stability scoring (resonance, hyperconjugation)

---

### 3. Mechanism Classification ⚠️ MEDIUM PRIORITY

**What**: Classify mechanism type (SN1, SN2, E1, E2, radical, etc.)

**Implementation**:

```python
def classify_mechanism(
    reaction_smiles: str,
    bond_changes: Dict[str, Any],
    functional_groups: Dict[str, Any]
) -> Dict[str, Any]:
    """
    Classify reaction mechanism.

    Returns:
        {
            "mechanism_type": "SN2",
            "confidence": 0.85,
            "evidence": [
                "Primary alkyl halide",
                "Strong nucleophile",
                "One-step bond formation/breaking"
            ]
        }
    """
```

**Logic**:

- Check substrate type (1°, 2°, 3°)
- Check nucleophile strength
- Check leaving group ability
- Check solvent polarity
- Check stereochemistry changes

---

### 4. Transition State Estimation 🔬 ADVANCED (Optional)

**What**: Estimate transition state structure and energy

**Requirements**:

- Quantum chemistry calculations (Gaussian, ORCA, Psi4)
- Significant computation time
- Expert-level chemistry knowledge

**Not Recommended** for initial implementation (too complex)

---

### 5. Agent-Facing Narrative & Evidence Layer ⚠️ MEDIUM PRIORITY

**What**: Convert raw mechanism predictions into structured, citeable explanations the agent can safely quote.

**Implementation Ideas**:

- Build a `mechanism_renderer.py` that maps classifier output + electron flow traces into:
  - `steps`: ordered list with `title`, `description`, `key_atoms`, `electron_flow`, `evidence_refs`.
  - `narrative`: templated paragraph tuned for LLM responses (short sentences, explicit caveats).
- Add `explain.py` helper functions to stitch together bond analysis interpretation, functional group matches, and precedent snippets.
- Provide severity-graded warnings (e.g., "Mechanism unknown", "Multiple plausible pathways") so the agent can hedge naturally.

---

### 6. Telemetry, Guardrails & Cost Controls ⚠️ MEDIUM PRIORITY

**What**: Ensure the agent can diagnose failures and avoid costly or unsafe calls.

**Implementation Ideas**:

- Emit structured metrics (`mechanism.analysis.latency_ms`, `rxnmapper.cache_hit`) via existing logging hooks.
- Capture failure snapshots (input SMILES + error) for offline triage, with automatic retries on transient RDKit/RXNMapper issues.
- Enforce budget-aware switches: `detail_level` controls whether to run electron flow, intermediate prediction, or TS estimation.
- Add validation to block unreasonably large reactions (e.g., >80 atoms) unless `allow_large_systems=True`.

---

---

## 🛠️ Recommended Implementation Plan

### Phase 1: Basic Mechanism Analyzer (1-2 weeks)

**Goal**: Add rule-based mechanism analysis using existing bond change data

**Steps**:

1. **Create Mechanism Classifier** (`chemtools/mechanism/classifier.py`)

```python
"""
Mechanism classification based on structural features.
"""

def classify_mechanism(
    reactants: List[str],
    products: List[str],
    bond_changes: Dict[str, Any],
    reaction_family: str
) -> Dict[str, Any]:
    """
    Classify reaction mechanism type.

    Uses:
    - Bond change analysis (already available)
    - Functional group detection (already available)
    - Reaction family (already available)

    Returns:
        {
            "mechanism_type": "SN2",
            "electron_flow": [...],
            "key_steps": [...],
            "confidence": 0.85
        }
    """
```

2. **Create Electron Flow Rules** (`chemtools/mechanism/electron_flow.py`)

```python
"""
Rule-based electron flow prediction.
"""

MECHANISM_RULES = {
    "SN2": {
        "description": "Nucleophilic substitution (bimolecular)",
        "requirements": {
            "substrate": ["primary", "secondary"],
            "nucleophile": "strong",
            "leaving_group": "good"
        },
        "steps": [
            {
                "step": 1,
                "description": "Nucleophile attacks electrophilic carbon",
                "electron_flow": "Nu: → C-X → X:",
                "energy": "single_transition_state"
            }
        ]
    },
    # Add more mechanisms...
}

def predict_electron_flow(
    mechanism_type: str,
    reactants: List[str],
    products: List[str]
) -> Dict[str, Any]:
    """Generate electron flow diagram."""
```

3. **Create Agent Tool Wrapper** (`chem_assistant/chemtools_wrapper.py`)

```python
@tool(args_schema=MechanismAnalysisInput)
def analyze_mechanism_tool(
    reaction_smiles: str,
    include_electron_flow: bool = True,
    include_intermediates: bool = True
) -> Dict[str, Any]:
    """
    Analyze reaction mechanism including electron flow and intermediates.

    Provides:
    - Mechanism type (SN1, SN2, E1, E2, etc.)
    - Electron flow description (arrow pushing)
    - Predicted intermediates
    - Step-by-step mechanism breakdown
    - Rate-determining step identification

    This complements bond change analysis by providing mechanistic insights.

    Args:
        reaction_smiles: Reaction SMILES
        include_electron_flow: Include electron flow arrows
        include_intermediates: Predict intermediate structures

    Returns:
        Dict with mechanism analysis

    Example:
        >>> analyze_mechanism_tool(
        ...     "CCBr.NCC>>CCNCC",
        ...     include_electron_flow=True
        ... )
    """
    try:
        from chemtools.mechanism import classify_mechanism, predict_electron_flow
        from chemtools import analyze_bond_changes

        # Get bond changes first
        bond_analysis = analyze_bond_changes(reaction_smiles)

        # Classify mechanism
        mechanism = classify_mechanism(
            reaction_smiles=reaction_smiles,
            bond_changes=bond_analysis
        )

        result = {
            "mechanism_type": mechanism.get("type"),
            "confidence": mechanism.get("confidence"),
            "description": mechanism.get("description"),
            "key_steps": mechanism.get("steps", [])
        }

        if include_electron_flow:
            flow = predict_electron_flow(
                mechanism["type"],
                reaction_smiles
            )
            result["electron_flow"] = flow

        if include_intermediates:
            result["intermediates"] = mechanism.get("intermediates", [])

        return _success_response(result)

    except Exception as e:
        return _error_response(str(e))
```

4. **Add to CHEMTOOLS_TOOLS list**

```python
CHEMTOOLS_TOOLS = [
    # ... existing tools ...
    analyze_bond_changes_tool,
    analyze_mechanism_tool,  # NEW
    # ...
]
```

---

### Phase 2: Advanced Features (2-4 weeks)

1. **Add More Mechanism Types**

   - Pericyclic reactions (Diels-Alder, Cope rearrangement)
   - Oxidation/reduction mechanisms
   - Organometallic mechanisms (Pd-catalyzed cross-coupling)
   - Radical mechanisms

2. **Intermediate Stability Scoring**

   - Resonance stabilization
   - Hyperconjugation
   - Inductive effects
   - Steric factors

3. **Visual Mechanism Diagrams**
   - Generate arrow-pushing diagrams
   - Energy profile diagrams
   - 3D structure visualization

---

## 📦 Required Dependencies

### Already Installed ✅

- RDKit (for structure manipulation)
- ChemTools bond analysis

### May Need to Add

```bash
# For advanced visualization
pip install py3Dmol         # 3D molecule viewer
pip install cairosvg        # SVG rendering
pip install pillow          # Image manipulation

# For energy calculations (optional, advanced)
pip install psi4            # Quantum chemistry (if needed)
```

---

## 🎯 Quick Start (Minimal Implementation)

If you want to **add basic mechanism analysis quickly**, here's a minimal approach:

### 1. Create `chemtools/mechanism/__init__.py`

```python
"""
Reaction mechanism analysis module.
"""

from typing import Dict, Any, List

def classify_mechanism_simple(
    reaction_smiles: str,
    reaction_family: str
) -> Dict[str, Any]:
    """
    Simple mechanism classification based on reaction family.
    """
    # Map reaction families to mechanism types
    FAMILY_TO_MECHANISM = {
        "Buchwald_CN": "oxidative_addition_reductive_elimination",
        "Suzuki": "transmetalation_coupling",
        "Amide_formation": "nucleophilic_acyl_substitution",
        "C_N_Coupling": "oxidative_addition_reductive_elimination",
        "SN2_Alkylation": "SN2",
        "Williamson_ether": "SN2",
    }

    mechanism = FAMILY_TO_MECHANISM.get(
        reaction_family,
        "unknown"
    )

    # Get descriptions
    DESCRIPTIONS = {
        "SN2": "Nucleophilic substitution (bimolecular) - backside attack, inversion",
        "oxidative_addition_reductive_elimination": "Pd-catalyzed: Pd(0) → Pd(II) → Pd(0)",
        "transmetalation_coupling": "Boronic acid transfers aryl to Pd catalyst",
        "nucleophilic_acyl_substitution": "Nucleophile attacks carbonyl, tetrahedral intermediate",
    }

    return {
        "mechanism_type": mechanism,
        "description": DESCRIPTIONS.get(mechanism, "Mechanism details not available"),
        "confidence": 0.7 if mechanism != "unknown" else 0.3
    }

__all__ = ["classify_mechanism_simple"]
```

### 2. Add Tool to Wrapper

```python
@tool
def analyze_mechanism_tool(reaction_smiles: str) -> Dict[str, Any]:
    """Analyze reaction mechanism."""
    from chemtools import detect_reaction, analyze_bond_changes
    from chemtools.mechanism import classify_mechanism_simple

    # Detect family
    detection = detect_reaction(reaction_smiles)
    family = detection.get("family", "unknown")

    # Get bond changes
    bonds = analyze_bond_changes(reaction_smiles)

    # Classify mechanism
    mechanism = classify_mechanism_simple(reaction_smiles, family)

    return {
        "success": True,
        "mechanism_type": mechanism["mechanism_type"],
        "description": mechanism["description"],
        "bond_changes": {
            "broken": bonds.get("broken_bonds", []),
            "formed": bonds.get("formed_bonds", [])
        },
        "reaction_family": family
    }
```

This gives you a **functional mechanism analyzer in ~1 hour** that provides basic insights!

---

## 📚 Recommended Resources

### For Learning Mechanisms

1. **Reaction Mechanism Databases**

   - name-reactions.com
   - ReactionFlash (mechanism database)
   - Organic Chemistry Portal

2. **SMARTS Patterns**

   - RDKit SMARTS tutorial
   - Daylight Chemical Information Systems

3. **Mechanism Classification**
   - "March's Advanced Organic Chemistry" (textbook)
   - "Strategic Applications of Named Reactions" (book)

---

## 🎯 Summary & Recommendation

### ✅ You Already Have (Use These!)

- **Bond change analysis** - Core foundation ✅
- **Functional group detection** - 80+ groups ✅
- **Reaction family detection** - Classifies reactions ✅
- **SMILES normalization** - Data preprocessing ✅

### 🚀 Recommended Next Steps

#### Phase Roadmap Snapshot

| Phase               | Duration   | Capabilities                                                            | Agent Outcome                                                               | Exit Criteria                                                                        |
| ------------------- | ---------- | ----------------------------------------------------------------------- | --------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ |
| Option 1 – MVP      | 1-2 hours  | Family→mechanism mapping, bond-change summary, tool wiring              | Agent can name plausible mechanism + cite bond evidence                     | ✅ Structured tool output (schema documented) <br> ✅ Unit tests for classifier stub |
| Option 2 – Enhanced | 1-2 weeks  | Rule-based electron flow, intermediate prediction, renderer, safeguards | Agent narrates multi-step mechanism, highlights risks, toggles detail level | ✅ ≥10 golden reactions validated <br> ✅ Telemetry + caching hooks in wrapper       |
| Option 3 – Research | 1-2 months | ML/quantum augmentations, TS estimation, visuals                        | Agent provides research-grade justification w/ energy profile               | ✅ Benchmarked vs literature set <br> ✅ Safety review for long-running jobs         |

**Option 1: Quick & Simple** (1-2 hours)

1. Create minimal `chemtools/mechanism` module
2. Map reaction families → mechanism types
3. Add simple `analyze_mechanism_tool`
4. Return: mechanism type, description, bond changes

**Option 2: Comprehensive** (1-2 weeks)

1. Build rule-based mechanism classifier
2. Add electron flow prediction
3. Predict intermediates
4. Generate step-by-step mechanisms
5. Add visual diagrams

**Option 3: Research-Grade** (1-2 months)

1. Everything from Option 2
2. ML-based mechanism prediction
3. Quantum chemistry integration
4. Transition state calculations
5. Energy profiles

### 💡 My Recommendation

Start with **Option 1** to get something working quickly, then incrementally add features from Option 2 as needed. The existing bond analysis tool gives you 80% of what users need for mechanism understanding!

---

## 🔧 Files to Create/Modify

### Create New Files:

```
chemtools/mechanism/
├── __init__.py (mechanism classification)
├── classifier.py (rule-based classifier)
├── electron_flow.py (electron flow rules)
└── intermediates.py (intermediate prediction)
```

### Modify Existing Files:

```
chem_assistant/chemtools_wrapper.py
  - Add MechanismAnalysisInput schema
  - Add analyze_mechanism_tool function
  - Add to CHEMTOOLS_TOOLS list
  - Update module docstring

chem_assistant/chemtools_agent.py
  - Update system prompt with mechanism analysis instructions
```

### Test Files:

```
tests/test_mechanism_analysis.py
test_mechanism_tool.py
```

## Agent Interaction Flow & Tool Contract

1. **Normalize & validate input** – sanitize SMILES, ensure atom count limits, pull context (temperature, solvent, catalyst) if provided.
2. **Synchronous analysis pipeline** – run bond analysis → functional group detection → reaction family classifier → mechanism classifier/electron flow/intermediate predictors.
3. **Evidence assembly** – collect references to precedent reactions, SMARTS matches, and atom indices so the agent can cite concrete facts.
4. **Narrative rendering** – convert structured steps into a short explanation plus optional detail-level blocks (summary, step-by-step, caveats).
5. **Response shaping for LLM** – include confidence, warnings, and cost flags; log telemetry + cache key for observability.

### Sample Tool Payload

```json
{
  "reaction_smiles": "CCBr.NH2CC>>CCNH2CC",
  "context": {
    "solvent": "DMF",
    "temperature_c": 80,
    "catalyst": "Pd(PPh3)4"
  },
  "detail_level": "high",
  "include_visuals": false
}
```

### Output Contract Highlights

- `mechanism_type` (string) and `confidence` (0-1)
- `steps`: ordered list with `description`, `electron_flow`, `key_intermediates`, `evidence_refs`
- `warnings`: e.g., "Ambiguous electrophile", "Multiple pathways plausible"
- `precedents`: optional list of similar reactions (source, similarity score)
- `metrics`: timings + cache hits for tracing

## Validation & Evaluation Plan

- **Unit tests** for classifier rules, electron-flow templates, renderer edge cases (use fixtures in `tests/conftest.py`).
- **Golden dataset** of ~25 curated reaction SMILES (SN1/SN2/E1/E2, Pd-couplings, radical) with expected mechanism outputs.
- **Integration tests** exercising the agent tool chain via `chemtools_wrapper` to catch schema regressions.
- **LLM-in-the-loop regression**: scripted prompt set (`scripts/run_mechanism_eval.py`) that ensures the agent references evidence and respects warnings.
- **Performance budget checks**: benchmark typical calls (<2 s) and add alerts if latency exceeds thresholds.

## Risks & Open Questions

1. **Reaction coverage** – rule-based approach may miss pericyclic/photoredox pathways; mitigate by allowing "unknown" with clear warnings.
2. **RDKit/RXNMapper brittleness** – long/charged SMILES can fail; add retries + graceful degradation to family-based heuristics.
3. **Intermediate explosion** – complex catalytic cycles could generate many intermediates; cap at `max_steps` and surface truncation notice.
4. **LLM hallucination** – ensure tool always returns machine-readable evidence so prompts can force the agent to quote data, not invent it.
5. **Data provenance** – precedent sourcing must avoid proprietary datasets; stick to bundled samples or user-provided libraries.

---

**Total Estimated Effort**:

- Minimal (Option 1): 1-2 hours
- Basic (Option 2): 1-2 weeks
- Advanced (Option 3): 1-2 months

**Dependencies**: Mostly already satisfied! Just need RDKit (already installed)

**Recommendation**: Start with Option 1 (minimal) to validate user interest, then expand based on feedback.
