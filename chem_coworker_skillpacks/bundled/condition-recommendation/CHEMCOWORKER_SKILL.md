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
