# Condition Coworker

`chem_coworker` is the small application layer for structure-first reaction-
condition recommendation. It contains no chemistry rules and no tool registry.

The request flow is deliberately short:

```text
ConditionRequest
  -> optional reaction-source completion confirmation
  -> condition_recommender.GenericConditionRecommender
  -> GenericRecommendationResult
  -> deterministic text rendering
```

`reactive_taxonomy` remains the molecular source of truth,
`condition_registry` owns condition identities and roles, and
`condition_recommender` owns admission, compatibility, retrieval, scoring, and
recipe aggregation. A named reaction family is optional.

Start the interactive CLI with:

```powershell
python -m chem_coworker
```

Enter reaction SMILES repeatedly at the `reaction>` prompt. The app guides
missing-source confirmation and supports `/help`, `/settings`, `/top-k`,
`/minimum`, `/profile`, `/json`, `/last`, `/save`, and `/quit`.

Run one non-interactive recommendation with:

```powershell
python -m chem_coworker "reactants>>product"
```

Use `--json` for the full, untruncated typed result. RXNMapper is optional and
is enabled explicitly with `--use-rxnmapper`.
