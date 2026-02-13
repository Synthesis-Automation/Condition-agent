# Reaction Agent Foundation (v1)

This folder is the first implementation of a general-purpose reaction analysis agent loop:

- `gateway` -> session lane orchestration
- `planner` -> dynamic next-step selection
- `tools` via `tool_registry` -> deterministic analysis, validation, coverage advice
- `memory` -> persistent in-process session state

## Current toolchain

1. `analyze_reaction`: wraps `poc_gpt52_reaction_v2` deterministic taxonomy analyzer
2. `fallback_candidate_retrieval`: placeholder contract (not implemented yet)
3. `validate_decision`: agent-level confidence/candidate gate
4. `confidence_calibrator`: placeholder contract (not implemented yet)
5. `llm_rerank_constrained`: placeholder contract (not implemented yet)
6. `precedent_lookup`: placeholder contract (not implemented yet)
7. `coverage_advice`: taxonomy/tool coverage expansion suggestions for unknown cases

## Workflow (foundation mode)

Deterministic sequence in current foundation implementation:

1. `analyze`
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

Unknown case (with coverage suggestions):

```bash
python -m reaction_agent_v1.cli "Clc1ncc(-c2ccccc2)cn1.NN>>NN=c1ncc(-c2ccccc2)c[nH]1" --json
```
