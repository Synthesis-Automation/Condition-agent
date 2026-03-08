---
id: literature_curation
name: Literature Curation
version: 1
summary: Structured literature and note extraction with taxonomy-aware filing expectations.
category: chemistry
tool_policy: literature_curator
workflow_targets:
  - forward_chemistry
  - retrosynthesis
triggers:
  keywords:
    - curate notes
    - extract notes
    - literature summary
  requires_any: []
tool_allowlist:
  - extract_notes
eligibility:
  rdkit: false
  python_modules: []
  data_files:
    - knowledge_base
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 600
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
priority: 60
---

## When to use

Use this skill for document-to-notes workflows and knowledge-base curation tasks.

## Required behavior

1. Keep extracted notes generalizable and concise.
2. Preserve taxonomy alignment when writing reaction-type labels.
