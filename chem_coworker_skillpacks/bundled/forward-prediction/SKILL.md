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
prompting:
  inject_mode: on_demand
  max_tokens_hint: 650
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
priority: 75
---

## When to use

Use this skill for reactant-first forward prediction or when ranking plausible product outcomes.

## Required behavior

1. Distinguish prediction from explanation.
2. Prefer validated proposals over unsupported product guesses.
