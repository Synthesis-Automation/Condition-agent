---
id: retrosynthesis_route_planning
name: Retrosynthesis Route Planning
version: 2
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
    - synthesize
    - synthesis of
    - how would you prepare
    - starting materials
    - disconnection
    - precursors for
    - synthetic route
  requires_any:
    - molecule_smiles
tool_allowlist:
  - retrosynthesis_step
  - plan_route
  - plan_route_candidates
  - validate_synthesis_proposal
  - identify_retrons
  - generate_disconnections
  - search_hte_precedent
  - search_by_product_similarity
  - apply_hte_templates
  - inspect_target
tool_default_args:
  retrosynthesis_step:
    condition_strategy: full
    include_precedent: true
    include_conditions: true
  plan_route:
    max_depth: 4
    max_branches: 2
  plan_route_candidates:
    top_k: 5
    beam_width: 8
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: false
  must_surface_warnings: true
  forbid_knowledge_only_conditions: true
eligibility:
  rdkit: true
  python_modules:
    - rdkit
  data_files:
    - data/HTE_db
  env_vars: []
  taxonomy_ids: []
priority: 85
---

## When to use

Use this skill for target-first route planning, disconnection analysis, and precursor proposal evaluation.

## Required behavior

1. For full-route tasks, prefer `plan_route_candidates` over a single greedy route, then compare the top candidates against the user's stated strategy.
2. Keep route proposals grounded in explicit step or validation tool evidence.
3. Surface weak points, protecting group risk, and uncertain disconnections.
4. Avoid presenting unverified routes as settled.
5. Never fabricate routes from general knowledge alone — all disconnections and conditions must be backed by tool output.
