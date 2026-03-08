---
id: condition_recommendation
name: Condition Recommendation
version: 1
summary: HTE-backed condition selection with explicit evidence and caveats.
category: chemistry
tool_policy: conditions_specialist
workflow_targets:
  - forward_chemistry
  - forward_synthesis
triggers:
  keywords:
    - recommend condition
    - suggest condition
    - what catalyst
    - what solvent
  requires_any:
    - reaction_smiles
tool_allowlist:
  - analyze_reaction
  - recommend_reaction_conditions
  - reagent_assistant
tool_preferred_order:
  - analyze_reaction
  - recommend_reaction_conditions
provides_context:
  - reaction_type
  - recommendations
requires_context:
  - reaction_type
eligibility:
  rdkit: true
  python_modules:
    - rdkit
  data_files:
    - data/HTE_db
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 700
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
  forbid_knowledge_only_conditions: true
priority: 90
---

## When to use

Use this skill for catalyst, ligand, solvent, base, additive, and temperature recommendations for a defined transformation.

## Required behavior

1. Establish reaction identity from taxonomy-backed evidence before recommending conditions.
2. Prefer HTE-backed output over unaudited model knowledge.
3. If evidence is thin or missing, say so and lower confidence.
