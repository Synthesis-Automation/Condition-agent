---
id: reaction_analysis
name: Reaction Analysis
version: 1
summary: Taxonomy-backed reaction interpretation and evidence gathering.
category: chemistry
tool_policy: general_chemistry
workflow_targets:
  - forward_chemistry
triggers:
  keywords:
    - analyze reaction
    - what reaction
    - identify reaction
  requires_any:
    - reaction_smiles
tool_allowlist:
  - analyze_reaction
  - featurize_molecule
eligibility:
  rdkit: true
  python_modules:
    - rdkit
  data_files: []
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 500
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
priority: 70
---

## When to use

Use this skill when the user asks what a reaction is doing, which reaction family it belongs to,
or what the substrate and product changes imply mechanistically.

## Required behavior

1. Prefer taxonomy-backed analysis before free-form explanation.
2. Surface ambiguity when reaction identity is weak or mixed.
3. Keep the answer tied to observed bond changes and detected reaction family.
