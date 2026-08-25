# Chem Coworker

`chem_coworker` is the small application layer for structure-first condition
recommendation and chemistry-first one- and multistep retrosynthesis. It contains no
chemistry rules and no tool registry.

The canonical optional LLM path is the bounded assistance controller in
`chem_coworker.assistance`. It uses one immutable evidence ledger and one shared
schema-constrained transport for conditions, one-step retrosynthesis, and
multistep retrosynthesis. Chemistry, recommendation, and route-search packages
remain authoritative; model claims are advisory and cannot mutate their typed
results.

Assistance defaults to `off`. The supported rollout sequence is `off`, `shadow`,
`advisory`, and—only after the documented chemistry and human-review gates—
`canonical_advisory`. The older single-pass review contracts remain a gated
compatibility surface until those external activation gates are complete, but
their provider and retry implementation now uses the same shared transport.

The deterministic request flow remains deliberately short:

```text
ConditionRequest
  -> optional reaction-source completion confirmation
  -> condition_recommender.GenericConditionRecommender
  -> GenericRecommendationResult
  -> optional bounded LLM condition review
  -> deterministic rendering of the result and advisory review

RetrosynthesisRequest
  -> core_retrosynthesis generic operator retrieval
  -> reverse graph application and hard forward-signature validation
  -> deterministic STRAT1 strategy grouping and ranking
  -> optional condition_recommender evidence for each proposed reaction
  -> optional bounded LLM feasibility and chemoselectivity review
  -> deterministic rendering of the strategies and advisory review

MultistepRetrosynthesisRequest
  -> bounded core_retrosynthesis search over validated one-step operators
  -> explicit supplier-stock/literature or molecular-weight terminal decisions
  -> deterministic beam search, route scoring, and solved/partial separation
  -> optional deterministic condition evidence for every retained step
  -> optional assistance inspection of typed route issues
  -> at most one policy-bounded deterministic route refinement search
  -> optional bounded LLM review of fixed route candidates as complete sequences
  -> deterministic rendering with original and advisory display ranks
```

Run an explicit one-shot advisory session with:

```powershell
python -m chem_coworker --assistance advisory --review off "reactants>>product"
python -m chem_coworker --mode retro --assistance advisory --review off "target_smiles"
python -m chem_coworker --mode multistep --assistance advisory --review off "target_smiles"
```

Use `--assistance shadow` to exercise the controller without displaying its
answer. Traces are not persisted by default; `--save-assistance-trace PATH` is
the explicit opt-in. If the controller asks a condition clarification in
advisory mode, the CLI accepts exactly one answer, resolves it through
`condition_registry`, and recomputes through `condition_recommender`.

Hard condition restrictions can also be supplied directly:

```powershell
python -m chem_coworker --exclude-condition "Pd(PPh3)4" --maximum-temperature-c 60 "reactants>>product"
```

The local API exposes `/api/v1/experimental/assistance` and
`/api/v1/experimental/assistance/confirm-condition` only when an
`AssistanceApplicationService` is explicitly injected. An unconfigured server
returns `ASSISTANCE_NOT_CONFIGURED`; ordinary deterministic endpoints are
unchanged.

`reactive_taxonomy` remains the molecular source of truth,
`condition_registry` owns condition identities and roles, and
`condition_recommender` owns admission, compatibility, retrieval, scoring, and
recipe aggregation. A named reaction family is optional.

The deterministic result remains the source of truth. When review is enabled,
`top-k` is the requested number of final condition strategies, not the size of
the initial recipe pool. The coworker retrieves at least twice `top-k` (capped at
50, and at least the configured review limit), then sends the leading recipes to
the LLM in a typed evidence packet. It can mark a recipe as
`keep`, `downrank`, `flag`, or `needs_information`, and it can change the display
order when that setting is enabled.

The LLM must also place every reviewed recipe into exactly one condition-strategy
group. Minor base, solvent, concentration, temperature, or workup differences
may be grouped when the chemistry-defining catalyst/ligand or activation system
is the same. Materially different catalyst systems, redox regimes, coupling
reagents, or compatibility risks remain separate. For example, several Suzuki
recipes using tetrakis(triphenylphosphine)palladium, carbonate base, and aqueous
ethereal solvent can be shown as one strategy, while a Pd(dppf)Cl2 recipe remains
distinct. The code chooses one group representative from the LLM verdict and the
original deterministic rank.

The LLM cannot modify a recipe or score, invent a precedent, restore a
chemistry-filtered candidate, or replace the original ranking stored under
`result` in JSON output. Model output with unknown recipe/evidence IDs, duplicate
group membership, or incomplete grouping is rejected as a failed review.

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

Switch to one-step retrosynthesis with `/mode retro`, then enter a single target
SMILES at the `target>` prompt. Return with `/mode conditions`. Retrosynthesis
uses the validated generic operator library in
`results/operator_retrosynthesis_poc/full_scale_v3/compact/operator_library_v3.json.gz`
by default. It returns distinct `STRAT1` strategies, not a list inflated by
leaving-group or protecting-group variants.

Only forward-validated graph proposals reach the LLM. The review can annotate
and reorder those strategies for functional-group compatibility,
chemoselectivity, competing reactive sites, precursor plausibility, protecting-
group needs, precedent transfer, and condition feasibility. It cannot generate
or edit precursors, restore a rejected proposal, or change the deterministic
ranking stored under `strategies` in JSON.

The model-facing packet uses short request-local aliases such as `strategy-1`
and `evidence.strategy-1.conditions`. Code maps those aliases back to the full
deterministic `STRAT1` identities after schema validation, avoiding copy errors
in long hash IDs without accepting fuzzy or invented evidence references.

The retrosynthesis controls are:

```text
/mode retro                    switch the prompt to target SMILES
/top-k 5                       request five final strategies
/realizations 3                retain up to three concrete variants per strategy
/validate-limit 100            bound expensive forward-validation attempts
/retro-conditions on|off       attach deterministic condition recommendations
/review off|auto|always        control advisory LLM route review
```

Switch to multistep retrosynthesis with `/mode multistep`. The deterministic
planner searches only to depth 2 or 3 and stops at an explicit expansion budget.
It never treats a partial route as solved. Supplier-stock provenance is preferred
for terminal starting materials; the literature molecule index is a fallback,
and low molecular weight remains a clearly labeled heuristic.

The older single-pass route-review contract remains outside the search loop and
can only annotate fixed candidates. The assistance controller adds a separate,
policy-bounded agent loop. It may inspect application-derived route issues and
request one `refine_route` action, but it still cannot add a structure, edit
SMILES, select atom edits, author conditions, mark a partial route solved, or
override deterministic chemistry.

`refine_route` accepts only an existing route alias, step index, cited typed
issue evidence, a closed objective, and one of two initial methods:

- `alternate_disconnection` excludes the selected step's exact `STRAT1`
  strategy for that product expansion;
- `alternate_realization` excludes only its exact concrete realization while
  retaining other realizations of the same strategy.

The deterministic planner derives the molecular exclusion from the stored
route, reruns bounded search, recomputes step condition evidence, and compares
the new routes with the source using objective-specific issue counts. The source
route is never mutated or discarded. Refinement evidence records parent result
and route IDs, the excluded candidate scope, accepted candidate IDs, remaining
issue counts, and whether a verified deterministic improvement was found.

Initial typed issue kinds are precursor compatibility, selectivity, insufficient
condition evidence, unresolved leaves, and tactical functional-state changes.
Resolving one of these signals does not establish experimental feasibility.
Short aliases such as `route-1`, `route-1.step-2`, and `route-1.issue-1` prevent
long internal IDs from being copied incorrectly.

Multistep interactive controls are:

```text
/mode multistep               switch the prompt to a route target
/depth 2|3                    set maximum search depth
/per-step 5                   retain five actions per molecule expansion
/beam 20                      retain twenty frontier route states
/expansions 12                cap expanded route states
/guidance TEXT                set advisory route-ranking preferences
/guidance clear               remove route guidance
/retro-conditions on|off      attach deterministic conditions to every step
```

Guidance can express preferences such as “avoid palladium”, “prefer a convergent
route”, or “minimize protecting groups”. It is untrusted preference data and only
affects the LLM presentation order; it does not change route generation, route
validity, terminal status, or the stored deterministic ranking.

`auto` reviews a multi-strategy result and also responds to selectivity
warnings, precursor-compatibility cautions, weak support, or generalized
operator fallback. Condition evidence is advisory: lack of a matching recipe
does not invalidate a graph-verified route.

The main LLM controls are:

```text
/model                         list available models
/model gpt-5.6-terra           select model and infer its provider
/review off|auto|always        select when review runs
/reasoning none|low|medium|high|xhigh|max
/review-limit 10               review at most the first ten recipes
/max-tokens 8000               cap model output
/review-order on|off           apply or disable advisory display ordering
```

`auto` reviews results with explicit uncertainty signals such as warnings,
cautions, weak support, broad fallback retrieval, or closely scored leading
recipes. Multiple retrieved condition variants also trigger grouping in `auto`
mode. `always` reviews every valid result that contains recommendations.
`off` performs deterministic recommendation only. Settings changed in the
interactive app are saved in `~/.chemcoworker/config.json`.

OpenAI review uses the Responses API with schema-constrained output. Aliyun is
also available through its OpenAI-compatible chat endpoint; set
`ALIYUN_API_KEY` before selecting an Aliyun model. If a key is missing or a
provider call fails, the deterministic recommendation is still returned and the
review is visibly marked unavailable.

Structured condition review defaults to 8,000 output tokens because reviewing
and grouping ten recipes can require substantially more output than a simple
answer. If OpenAI or Aliyun returns incomplete, empty, or truncated JSON, the
coworker retries once with a larger budget (up to 20,000) and reduced reasoning
effort. A successful retry is reported in the review metadata. If both attempts
fail, the warning reports the provider completion state instead of exposing a
raw Pydantic end-of-file validation error.

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

Run one non-interactive one-step retrosynthesis with:

```powershell
python -m chem_coworker --mode retro "target_smiles"
python -m chem_coworker --mode retro --top-k 5 --review always --model gpt-5.6-terra "target_smiles"
```

Run a bounded multistep search with:

```powershell
python -m chem_coworker --mode multistep --max-depth 3 --top-k 5 "target_smiles"
python -m chem_coworker --mode multistep --review always --guidance "prefer convergence and avoid protecting-group churn" "target_smiles"
```

Additional one-shot bounds are `--per-step-top-k`, `--beam-width`, and
`--max-expansions`. Use `--stock-index PATH` to select a compatible supplier
portfolio or literature molecule index. `--json` includes every route, partial
route, leaf decision, step-level condition result, diagnostic counter, and LLM
evidence citation.

Retro and multistep commands show a compact elapsed-time status while loading
the operator library and running the bounded search. Interactive terminals
refresh one line; redirected runs emit only start and completion lines. Progress
is written to stderr, so `--json` stdout remains valid JSON. Use `--no-progress`
to suppress it for quiet automation.

Use `--no-retro-conditions` for a faster structure-only search. Override the
operator artifact with `--retrosynthesis-library PATH`; use
`--max-realizations`, `--max-templates-to-apply`,
`--max-candidates-to-validate`, `--no-include-l0`, or `--no-retro-context` for
bounded search experiments.

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
