# Automatic broadening: honor the requested recipe count

The Workbench passed `top_k` correctly, but shared-core retrieval stopped at the
first direct tier with sufficient independent evidence. For the reported query,
automatic mode returned five detailed-local recipes even with `top_k=10`, while
Broader analogues returned ten:

```text
Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1
```

Automatic mode now continues through permitted direct tiers until it has both
the minimum independent support and the requested recipe count. If direct tiers
are exhausted first, it consults the existing auxiliary channels. The count uses
the same match-level, source-context and canonical recipe-core grouping as final
ranking, after graph and condition compatibility checks. Repeated observations
of one recipe cannot fill the target by themselves.

Closer matches remain first. Same reactive handle and Broader analogues retain
their scope behavior. Candidate budgets, chemistry gates, evidence independence,
ranking weights and source requirements are unchanged. The requested count is
still an upper bound: a library or bounded search may provide fewer recipes.

The retrieval definition is now `shared_core_retrieval.v3@3.1`, with an explicit,
validated automatic stopping policy. API routes and request/response schemas,
reaction definitions, index identities and stored projections are unchanged.
No corpus conversion or index rebuild is required. Restart the web server to
load the backend change; rebuild the frontend for the updated scope explanation.

Regression coverage includes the supplied Suzuki query in memory and SQLite,
preservation of closer results, duplicate-recipe aggregation, independent support,
auxiliary retrieval after a recipe shortfall, condition constraints, and all three
search scopes. This is a retrieval-control fix, not a new chemistry validation or
release-gate claim.

## Local verification

The Workbench HTTP endpoint was exercised against the existing Full and Compact
libraries. Full automatic retrieval with `top_k=10` changed from five L0 recipes
to ten recipes: the original five L0 recipe IDs and precedent IDs in their
original order, then five L1 recipes. Full `top_k=5` and Same reactive handle
still return the original five. Full automatic retrieval returns ten with
RXNMapper enabled or disabled; Compact returns ten (two L0 and eight L1).
Responses are saved locally under `results/automatic_broadening_fix/`.

The 35 focused shared-core tests, complete `pytest -q` suite (1,541 passed in
394.67 seconds), and frontend production build pass. Definition validation rejects
unsupported schema, version, and stopping-policy values; Python name checks pass.
Interactive browser verification was unavailable because no browser was connected;
the rebuilt Workbench entry point and HTTP recommendation endpoint were verified.
