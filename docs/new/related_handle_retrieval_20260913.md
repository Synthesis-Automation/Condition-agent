# Related-handle retrieval and simpler ranking — 2026-09-13

## Problem and resulting behavior

The Workbench exposed nine overlapping ranking factors. Changing these weights
could not retrieve Ar-Br precedents for an Ar-I query because structural lookup
keys and departing-fragment gates operated before ranking. Progressive search
could also report its target reached based on distinct recipes without enough
independent evidence.

Search scope is now separate from ranking preferences. Automatic broadening
checks a bounded related-handle tier when the close same-handle pool misses the
independent-support target (default 2). Explicit broader scope checks that tier
regardless of close support. Same-handle mode omits it and stops signed search
after the bond-edit tier, before broad core/edit-graph neighbours. A supported
related-handle tier can return fewer than `top_k` recipes; it does not force a
potentially expensive broader search to fill the display limit.

The main ranking choices are Balanced, Closest chemistry and Strongest
supporting evidence. Existing detailed weights remain under Advanced. Default
weights have not changed; the new closest-chemistry preset is an explicit,
uncalibrated prior. Historical yield and functional-group evidence labels no
longer imply predicted yield or demonstrated general tolerance.

## Chemistry contract

`reactive_taxonomy.handle_variants` owns the graph operation. It requires one
observed substitution event with a broken aromatic C-I/C-Br bond and one formed
bond at that same carbon. Multiple events, retained leaving atoms, conflicting
indices/maps, unsupported handles and order changes do not produce variants.
Only the departing atom is changed in a hypothetical lookup reaction; products
and agents are retained. The original signature, query and atom correspondence
are not changed or reported as if an analogue were observed.

`related_handle_retrieval.v1@1.0` permits I-to-Br and Br-to-I for single-bond
formation to C, N, O or S. The recommender featurizes the lookup hypothesis and
uses its existing exact facet key. Candidate departing-fragment multisets must
match in full, including any partner leaving fragment. Compatibility and user
condition constraints are assessed against the original query before support
counting and ranking. A different cyanide source is not automatically treated
as equivalent. Missing signatures/core evidence cannot enable this tier.

This is a candidate-transfer policy requiring chemical review, not a claim of
experimental interchangeability or a universal ordering of halogen reactivity.
No named reaction family is required or forced.

## Evidence and API changes

`precedent_match_levels.v1@1.0` defines ordinal distances:

1. Same normalized reactants and products.
2. Close precedent with matching reaction features; molecules may differ.
3. Explicitly permitted related handle, with the query/precedent elements shown.
4. Broader structural analogue requiring inspection.

These levels precede numerical ranking within each returned tier. Canonical
identity also takes priority when choosing a representative condition variant.
Recipe cores remain deduplicated; closer evidence is not displaced by a broader
recipe with a higher historical yield or score. `top_k` remains a maximum, so a
full set of closer recipes can hide retrieved lower-tier suggestions. JSON
traces retain the attempt, independent support and broadening reason.

- Recommendation result schema: 4.0 -> 4.1 (additive match fields, search scope,
  and retrieval-trace broadening provenance).
- Chemist ranking profile definition: 1.0 -> 1.1; preferences schema stays 1.0.
- New request field: `search_scope`, default `automatic`; invalid values receive
  HTTP 422. Python entry points accept the same independent keyword.
- Exact facet explanations now say matching features, rather than suggesting
  identical molecules. Empty progressive searches retain their traces; their
  no-candidate level is `no_chemically_compatible_precedent`.
- An externally mapped core remains available as fallback when its signature
  does not retrieve a usable pool; mapping provenance and cautions are retained.

The three primary packages retain their ownership boundaries. No application
chemistry rule, legacy path, alternative conversion pipeline or new runtime
dependency was introduced.

## Migration and validation

Baseline: clean repository at commit `0f6e7b2b`. Existing taxonomy signatures,
facet definitions, converted records, admission tiers and SQLite index schema
6.5 are unchanged. The new tier queries existing keys, so neither the Full nor
Compact dataset needs rebuilding. Browser assets need rebuilding and the server
needs restarting to load the new controls and Python code.

No dataset snapshots or untouched evaluation sets were changed. Coverage can
increase at query time through explicitly marked Level 3 evidence, while the
same-handle scope deliberately excludes the broad core/edit-graph tiers.
Statistical calibration and independent chemist release gates remain pending.

Full-index smoke check used all 571,157 trusted indexed records and query
`Ic1n[nH]c2ccccc12.N#C[Cu]>>N#Cc1n[nH]c2ccccc12`. Automatic/broader scope returned
three Br-based Level 3 precedents: `US07166621B2:767236_0`,
`US05801183:417039_0`, and `US08450363B2:1326394_0`. These are structural
analogues, not the identical reaction. Results and traces are local artifacts
under `results/related_handle_validation/full_*.json`.

The browser regression uses a one-bromide SQLite fixture. It checks the three
primary profiles, advanced functional-group label, request scope/profile,
visible Level 3 explanation, and absence of the analogue in same-handle mode.
Run it against a Workbench started with that fixture using
`RELATED_HANDLE_FIXTURE=1`, `WEBUI_TEST_PROFILE=research_workbench` and the
appropriate `WEBUI_TEST_URL`; the test otherwise skips the fixture-only check.

Final automated validation results are recorded below after completion.
