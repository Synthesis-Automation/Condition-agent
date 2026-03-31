---
id: condition_recommendation
name: Condition Recommendation
version: 2
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
    - condition recommendation
    - what catalyst
    - what solvent
    - what base
    - what ligand
    - reaction condition
    - optimal condition
  requires_any:
    - reaction_smiles
tool_allowlist:
  - analyze_reaction
  - recommend_reaction_conditions
  - reagent_assistant
tool_default_args:
  recommend_reaction_conditions:
    condition_strategy: full
    top_k: 5
  analyze_reaction:
    include_conditions: true
    condition_strategy: full
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
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
priority: 90
---

## When to use

Use this skill for catalyst, ligand, solvent, base, additive, and temperature recommendations for a defined transformation.

## Required behavior

1. Confirm taxonomy-backed reaction identity before recommending conditions.
2. Prefer tool-backed condition evidence over model-only suggestions.
3. Surface uncertainty when evidence is weak, sparse, or absent.
4. Never fabricate conditions from general knowledge alone — all recommendations must be backed by tool output.
