# Condition Coworker

`chem_coworker` is the small application layer for structure-first reaction-
condition recommendation. It contains no chemistry rules and no tool registry.

The request flow is deliberately short, with the LLM used only where contextual
judgment is useful:

```text
ConditionRequest
  -> optional reaction-source completion confirmation
  -> condition_recommender.GenericConditionRecommender
  -> GenericRecommendationResult
  -> optional bounded LLM condition review
  -> deterministic rendering of the result and advisory review
```

`reactive_taxonomy` remains the molecular source of truth,
`condition_registry` owns condition identities and roles, and
`condition_recommender` owns admission, compatibility, retrieval, scoring, and
recipe aggregation. A named reaction family is optional.

The deterministic result remains the source of truth. The LLM reviews only the
top existing recipes using a typed evidence packet. It can mark a recipe as
`keep`, `downrank`, `flag`, or `needs_information`, and it can change the display
order when that setting is enabled. It cannot modify a recipe or score, invent a
precedent, restore a chemistry-filtered candidate, or replace the original
ranking stored under `result` in JSON output. Model output with unknown recipe or
evidence IDs is rejected as a failed review.

Start the interactive CLI with:

```powershell
python -m pip install -r requirements-cli.txt
$env:OPENAI_API_KEY = "your-key"
python -m chem_coworker
```

Startup checks the persisted index contract. It prefers the full index when it
is current and automatically falls back to the compact index when the full
artifact is stale or unavailable. An explicit `--index` never falls back.

Enter reaction SMILES repeatedly at the `reaction>` prompt. The app guides
missing-source confirmation and supports `/help`, `/settings`, `/top-k`,
`/minimum`, `/profile`, `/model`, `/provider`, `/review`, `/reasoning`,
`/review-limit`, `/max-tokens`, `/review-order`, `/json`, `/last`, `/save`,
`/clear`, and `/quit`.

The main LLM controls are:

```text
/model                         list available models
/model gpt-5.6-terra           select model and infer its provider
/review off|auto|always        select when review runs
/reasoning none|low|medium|high|xhigh|max
/review-limit 5                review at most the first five recipes
/max-tokens 2000               cap model output
/review-order on|off           apply or disable advisory display ordering
```

`auto` reviews results with explicit uncertainty signals such as warnings,
cautions, weak support, broad fallback retrieval, or closely scored leading
recipes. `always` reviews every valid result that contains recommendations.
`off` performs deterministic recommendation only. Settings changed in the
interactive app are saved in `~/.chemcoworker/config.json`.

OpenAI review uses the Responses API with schema-constrained output. Aliyun is
also available through its OpenAI-compatible chat endpoint; set
`ALIYUN_API_KEY` before selecting an Aliyun model. If a key is missing or a
provider call fails, the deterministic recommendation is still returned and the
review is visibly marked unavailable.

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

For example, run a deterministic-only request or force review with a selected
model:

```powershell
python -m chem_coworker --review off "reactants>>product"
python -m chem_coworker --review always --model gpt-5.6-sol --reasoning-effort high "reactants>>product"
```

Use `--json` for the full, untruncated typed result. RXNMapper is optional and
is enabled explicitly with `--use-rxnmapper`.

## Windows and IDE launch command

From any PowerShell working directory, use:

```powershell
cmd /d /c "cd /d C:\Git-softwares\Condition-agent && C:\Users\xubar\AppData\Local\Programs\Python\Python312\python.exe -m chem_coworker"
```

For an IDE run configuration, use:

```text
Program: C:\Windows\System32\cmd.exe
Arguments: /d /c "cd /d C:\Git-softwares\Condition-agent && C:\Users\xubar\AppData\Local\Programs\Python\Python312\python.exe -m chem_coworker"
```

The value after `-m` must be the importable module name `chem_coworker`, not a
filesystem path.
