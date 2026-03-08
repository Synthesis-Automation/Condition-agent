---
id: retrosynthesis_route_planning
name: Retrosynthesis Route Planning
version: 1
summary: Route planning and proposal validation for target-oriented synthesis questions.
category: chemistry
tool_policy: retro_specialist
workflow_targets:
  - retrosynthesis
triggers:
  keywords:
    - retrosynthesis
    - plan a route
    - how to make
    - route to
  requires_any:
    - molecule_smiles
tool_allowlist:
  - retrosynthesis_step
  - plan_route
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
  max_tokens_hint: 800
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
priority: 85
---

## When to use

Use this skill for target-first route planning, disconnection analysis, and precursor proposal evaluation.

## Required behavior

1. Keep route proposals grounded in explicit step or validation tool evidence.
2. Surface weak points, protecting group risk, and uncertain disconnections.
