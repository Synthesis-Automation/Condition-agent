# Condition investigation adviser v1

This is an optional menu, not a required workflow. Skip, reorder, repeat, or replace
suggestions according to the question and available evidence. Scientific contracts,
recorded provenance, and honest uncertainty still apply.

Questions worth resolving:

- What exact graph transformation, chemical forms, stereochemistry, scale, and user
  constraints must the proposed experiment address?
- Which precedents match the reactive sites and surrounding substrate, and which
  differences could change selectivity or compatibility?
- Does the source describe a complete experiment, or only a product, generic scope,
  or incomplete indexed record?

Useful recorded calls (with an open workspace `w` and the supplied reaction):

```python
check = w.run("analyze_reaction", {"reaction_smiles": reaction_smiles})
hits = w.run("recommend_conditions", {"reaction_smiles": reaction_smiles, "top_k": 3})
```

Inspect the returned precedent IDs with `inspect_condition_precedents` and
`get_procedures` when useful. Compare whole observed recipes: component identities,
roles and amounts, activation, order of addition, solvent, temperature, time,
atmosphere, workup and isolation. Preserve missing values. Do not combine individually
attractive ingredients from different experiments and describe the mixture as observed.
Use `resolve_recipe` and `assess_recipe` for identity and compatibility;
`propose_condition_adaptation` records a complete proposed recipe, changed fields,
reasons, assumptions and risks separately from its source observation.

If direct source retrieval fails, inspect available primary sources with the browser.
Preserve actual visible text with `w.capture_source(text, url=url, locator=locator)`
and exact supporting passages with `w.record_source_excerpt`. This is agent-supplied
text, not independently verified retrieval. Retain access failures for debugging;
do not turn snippets, missing procedures, or an inaccessible document into evidence.

Stop when the user's decision has sufficient support or the remaining evidence is
unavailable. Explain the best-supported complete recipe, transfer risks, unresolved
fields, and useful alternatives. Retrieval ranking and passed compatibility checks
do not establish experimental yield or success.
