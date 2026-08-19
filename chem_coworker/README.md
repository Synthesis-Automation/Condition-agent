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
python -m pip install -r requirements-cli.txt
python -m chem_coworker
```

Startup checks the persisted index contract. It prefers the full index when it
is current and automatically falls back to the compact index when the full
artifact is stale or unavailable. An explicit `--index` never falls back.

Enter reaction SMILES repeatedly at the `reaction>` prompt. The app guides
missing-source confirmation and supports `/help`, `/settings`, `/top-k`,
`/minimum`, `/profile`, `/json`, `/last`, `/save`, `/clear`, and `/quit`.

On an interactive terminal, the enhanced interface adds:

- persistent history and up-arrow recall;
- tab completion for slash commands and ranking profiles;
- live analysis/retrieval status indicators;
- structured result tables and evidence/caution panels;
- multiline editing with `Alt+Enter` and submission with `Enter`.

Reaction history is stored in `~/.chemcoworker/history`. Use `--no-history`
for confidential sessions, or `--plain` for the dependency-free interface.

Run one non-interactive recommendation with:

```powershell
python -m chem_coworker "reactants>>product"
```

Use `--json` for the full, untruncated typed result. RXNMapper is optional and
is enabled explicitly with `--use-rxnmapper`.
