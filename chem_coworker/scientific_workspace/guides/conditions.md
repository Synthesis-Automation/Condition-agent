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

For a diverse screening panel, use the built-in recorded operation:

```python
panel = w.run("generate_weak_label_screening_array", {
    "reaction_smiles": reaction_smiles, "array_size": 96,
})
```

It returns up to the requested number of distinct intact recipes (default 24,
currently maximum 250) after graph-derived query and compatibility checks. It
requires `weak_label_records`; baseline capture also pins the sibling recipe
catalog. It does not require either structural condition index. An optional
`source_reaction_type_hint` can narrow supported labels, but contradictory labels
are rejected. Unsupported transformations and missing datasets remain explicit.

These are weak-label screening suggestions: source reaction structures are not
verified, and historical yields do not predict this query's yield. Cite recipe
IDs, source row numbers and the call artifact; preserve warnings, ambiguous source
types, missing values and complete recipes. Do not pass weak-label row numbers to
`get_precedents` or `get_procedures`, which query the structural corpus. Report
shortfalls when fewer compatible recipes exist. The full array remains under
`result.recommendations` in the call artifact; summaries show only a preview.
For requested JSON/CSV exports, read the full artifact, save inside the
investigation and attach the files with the generating call as evidence.

If direct source retrieval fails, inspect available primary sources with the browser.
Preserve actual visible text with `w.capture_source(text, url=url, locator=locator)`
and exact supporting passages with `w.record_source_excerpt`. This is agent-supplied
text, not independently verified retrieval. Retain access failures for debugging;
do not turn snippets, missing procedures, or an inaccessible document into evidence.

Stop when the user's decision has sufficient support or the remaining evidence is
unavailable. Explain the best-supported complete recipe, transfer risks, unresolved
fields, and useful alternatives. Retrieval ranking and passed compatibility checks
do not establish experimental yield or success.
