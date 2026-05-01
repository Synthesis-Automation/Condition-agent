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
  - build_reaction_context
  - detect_reaction_type
  - inspect_functional_groups
  - featurize_molecule
  - get_literature_condition_evidence
  - get_motif_condition_evidence
  - get_similarity_condition_evidence
  - get_rule_condition_evidence
  - compose_condition_candidates
  - score_condition_candidates
  - recommend_reaction_conditions
  - reagent_assistant
tool_default_args:
  compose_condition_candidates:
    top_k: 5
  score_condition_candidates:
    top_k: 5
  get_literature_condition_evidence:
    top_k: 5
  get_motif_condition_evidence:
    top_k: 5
  get_similarity_condition_evidence:
    top_k: 5
  get_rule_condition_evidence:
    top_k: 5
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

1. Start with a short heuristic chemistry analysis before tool use.
2. Build reaction context first, then gather source-specific evidence, compose candidates, then score candidates against evidence and the user's stated constraints.
3. Prefer atomic condition tools over the legacy one-shot condition facade.
4. Use `score_condition_candidates` for mild, cheap, metal-free, green, scale-up, HTE-only, literature-only, or substrate-compatibility preferences.
5. Surface uncertainty when evidence is weak, sparse, or absent.
6. Never fabricate conditions from general knowledge alone — all recommendations must be backed by tool output.
