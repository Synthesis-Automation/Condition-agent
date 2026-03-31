---
id: forward_prediction
name: Forward Prediction
version: 1
summary: Product prediction and forward synthesis evaluation for reactant-first queries.
category: chemistry
tool_policy: forward_specialist
workflow_targets:
  - forward_synthesis
triggers:
  keywords:
    - predict the product
    - what is the product
    - expected product
  requires_any:
    - reaction_smiles
tool_allowlist:
  - forward_synthesis_step
  - plan_forward_route
  - validate_synthesis_proposal
eligibility:
  rdkit: true
  python_modules:
    - rdkit
  data_files: []
  env_vars: []
  taxonomy_ids: []
priority: 75
---

## When to use

Use this skill for reactant-first forward prediction or when ranking plausible product outcomes.

## Required behavior

1. Distinguish product prediction from explanation or retrosynthesis.
2. Prefer validated proposals over unsupported product guesses.
3. Say when multiple outcomes remain plausible.
