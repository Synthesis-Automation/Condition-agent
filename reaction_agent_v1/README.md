# Reaction Agent Foundation (v1)

This folder is the first implementation of a general-purpose reaction analysis agent loop:

- `gateway` -> session lane orchestration
- `planner` -> dynamic next-step selection
- `tools` via `tool_registry` -> deterministic analysis, validation, coverage advice
- `memory` -> persistent in-process session state
- `evidence` -> canonical `ReactionEvidence` object shared by all tools

## Current toolchain

1. `reaction_diff`: canonical evidence extraction (principal pair, MCS, formula deltas, motifs, candidates)
2. `fallback_candidate_retrieval`: heuristic fallback candidate recovery from reaction diff evidence
3. `validate_decision`: agent-level confidence/candidate gate
4. `confidence_calibrator`: placeholder contract (not implemented yet)
5. `llm_rerank_constrained`: placeholder contract (not implemented yet)
6. `precedent_lookup`: placeholder contract (not implemented yet)
7. `coverage_advice`: taxonomy/tool coverage expansion suggestions for unknown cases

## Workflow (evidence-first foundation mode)

Deterministic sequence in current foundation implementation:

1. `reaction_diff`
2. `fallback_candidates` (only when deterministic candidates are empty)
3. `validate`
4. `confidence_calibration` (placeholder)
5. `llm_rerank` (placeholder)
6. `precedent_lookup` (placeholder)
7. `coverage` (only when final decision is `unknown`)
8. `finalize`

## Run

```bash
python -m reaction_agent_v1.cli "O=C(O)c1ccccc1.NCc1ccccc1>>O=C(NCc1ccccc1)c1ccccc1"
```

Fallback recovery example (SNAr-like):

```bash
python -m reaction_agent_v1.cli "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1" --json
```

Unknown case (with coverage suggestions):

```bash
python -m reaction_agent_v1.cli "CCO.CN>>CCN" --json
```
