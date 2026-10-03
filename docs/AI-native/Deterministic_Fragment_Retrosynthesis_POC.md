# Deterministic fragment-guided retrosynthesis POC

Status: development investigation, 2026-10-03. This is a single-step experiment,
not an independent chemistry review, untouched evaluation, or complete route planner.

The runnable example is
[`fragment_guided_retrosynthesis_poc.py`](../../examples/ai_native/fragment_guided_retrosynthesis_poc.py).
It uses existing domain packages, without an LLM, web search, or a second
reaction converter. Its bounded selection and execution settings are validated
from [`fragment_guided_retro_policy.v1.json`](../../core_retrosynthesis/definitions/fragment_guided_retro_policy.v1.json).

## Regular research web UI

Start or restart with `python -m app.web_api --workbench --build --port 8000`.
Select **Fragment-guided retro**. The default **Assisted fragment research**
workflow keeps the full target separate from the editable search query. Paste or
draw a connected target, optionally suggest strategic regions, then choose or
edit a core and click **Search chosen fragment**. Inspect source reactions,
construction/retention evidence and linked procedures. Revise the query explicitly
or request smaller target-derived suggestions; searches never broaden themselves.
Record why you changed the query, choose up to six observations and click
**Assess selected precedents on target**. Discovery needs the prepared fragment
index; assessment also needs the selected compact/full operator library.

The assessment rechecks query membership and selected observation IDs against the
current index. IDs outside the fresh returned-hit sample are rejected. Complete
searches, validated mapping and resolved internal formed-bond evidence are
required for guidance. Partial results remain inspectable. All bounded target
embeddings are preserved as alternative hypotheses; symmetry is not silently
resolved. Truncated alignment enumeration or excessive projection work cannot
seed guidance. These limits are validated from
[`fragment_transfer.v1.json`](../../core_retrosynthesis/definitions/fragment_transfer.v1.json)
(`fragment_transfer_policy.v1`, definition `fragment_transfer.v1@1.0`): at most
64 target alignments and 2,048 projected witness combinations. Focus bonds still
rank by independent reference support, with atom IDs breaking ties.

Selected-source operators and witness-directed library proposals use the same
compiler admission and forward signature/formed-bond gates as the automatic POC.
Source-product/target comparisons retain unmatched regions, attachments, stereo
and alignment uncertainty in exported evidence. They do not validate conditions
or experimental feasibility. SMILES queries permit extra substitution unless
constrained: omitting a substituent stops requiring it, rather than requiring its
absence. SMARTS edits are explicit and are not automatically certified as logical
generalizations. The tool uses the local index and its procedures, without a CAS
search integration or automatic literature-method adaptation.

The bounded unrestricted baseline is optional, disabled by default for manual
assessment. Without it, no additional-coverage claim is made. Browser history
retains up to 20 completed search/assessment attempts or errors, parent links,
queries, notes and source choices. Restore earlier attempts or export research
JSON before refreshing/closing the page. This is browser memory, not durable
storage or a pinned investigation.

For the existing automatic experiment, choose **Automatic POC comparison**, then
enter/draw a target or use **Fragment retro example**, and click **Evaluate
fragment-guided retro**. It requires both index and operator library.
`--fragment-index` or `FRAGMENT_PRECEDENT_INDEX` chooses the discovery index;
`CORE_RETROSYNTHESIS_LIBRARY_ROOT` chooses the directory containing `compact/` and `full/` libraries.

The UI shows selected fragment regions, complete/partial search status, source
compiler admissions, baseline and guided proposals, extra precursor sets and
actual search work. Controls bound fragment queries, construction bonds and
proposals per arm. Source admission and forward signature gates cannot be disabled.
**Export JSON** preserves source records, matches, exclusions, policy versions,
resource file metadata and all candidate diagnostics. Interactive calls execute
once and do not claim repeated-run determinism or independent evaluation.
They do not create a scientific-workspace investigation; use the CLI below for
recorded, pinned, replayable investigations. Start a new investigation after this
shared domain-code change; prior saved artifacts remain readable.

The graph projection and transfer evaluation now live in
[`core_retrosynthesis.fragment_guidance`](../../core_retrosynthesis/fragment_guidance.py).
The CLI and HTTP runtime compose that same implementation; no duplicate chemistry
path or corpus rebuild was introduced.

### HTTP contracts

`POST /api/v1/fragments/search` accepts an optional `target_smiles` alongside
the existing `query`, `query_format`, `topology`, `limit` and `timeout_seconds`.
Supplying a target requires the compiled query to match it before index search.
Existing standalone searches remain supported without a target.

`POST /api/v1/retrosynthesis/fragment-transfer` requires that target and query,
plus `selected_observation_ids` (one to six distinct IDs). It accepts
`library_mode`, `max_focus_bonds`, `top_k` and `include_baseline` (default false).
The server loads authoritative source records; client-supplied source chemistry
is rejected. Responses retain `query_search`, projection policy/alignment
evidence, `source_comparisons`, source admissions and transfer arms. Both routes
are research-workbench capabilities; the focused condition-only profile hides
the transfer route.

### Authored manual development check

For the example target above, an explicit query revision from
`Fc1cnc2cn[nH]c2c1` to `c1cnc2cn[nH]c2c1` stops requiring fluorine while retaining
the fused core. In the current local index, complete searches found 5 versus 55
observations, including 5 versus 14 with construction evidence. Groups overlap;
more hits alone do not establish more useful or transferable precedents.

Selecting `US09073918B2:1639626_0` and `US09029536B2:1615846_0` from the relaxed
query produced four distinct verified precursor sets through witness-directed
library searches. Both complete source reactions failed operator admission with
`materialized_core_not_verified`, so direct-source transfer returned no proposals.
No baseline was requested. The recorded calculation repeated identically; the
live HTTP/browser assessment returned the same precursor count and rejections.
This is one authored development case, without experimental or benchmark
effectiveness claims. Local records are under
`results/ai_native/fragment_research_20261003/`; browser evidence and exported
session are under `results/ai_native/fragment_research_web_20261003/`.

## Run

From the repository root, using the environment with RDKit and RDChiral:

```powershell
python -m examples.ai_native.fragment_guided_retrosynthesis_poc --output results/ai_native/my_fragment_retro_poc
```

The default target is `Fc(cn1)cc2c1c(c3ccccc3OC)n[nH]2`. Use `--target` for another
connected target, `--index` for an existing fragment index, and `--library` for an
existing generic operator library. Each run requires a **new output directory**.
The example never rebuilds the corpus or its indexes. It locally compiles at most
six selected source observations into a small experimental operator library.

Initialization pins the scientific code, definitions, environment, selected data,
and POC policy. Fragment calls and subsequent custom computations are recorded in
the scientific workspace. The runner snapshots its complete Python source into
the investigation before execution; custom execution records preserve that source,
its hash, parameters, input evidence, errors and output.

## What the code explores

1. Generate target-derived queries and choose the first three in the existing
   deterministic suggestion order. These are overlapping search regions, not a
   partition or concrete precursors.
2. Search each query with target-membership validation. Incomplete searches,
   unresolved or truncated embeddings, retention, and boundary-only changes
   cannot seed construction-site hypotheses.
3. Project internal formed-bond witnesses through the query alignment onto
   canonical target atom IDs. Preserve each supporting source and embedding;
   this is an analogue hypothesis, not observed target atom mapping.
4. Compile eligible **complete source reactions** through `build_generic_library`
   in data-driven mode, retaining the existing `pass_only` chemistry policy and
   source round-trip admission. Keep rejected observations and their diagnostics.
5. Compare three kinds of output: an unrestricted prepared-library baseline;
   focused applications of admitted search-hit operators; and focused search of
   the prepared library guided by the projected bond witnesses.
6. Require returned proposals to have verified signatures and focused proposals
   to have final formed-bond verification. Preserve other edits, compatibility
   assessments, contextual relaxation levels, precursor structures and source IDs.

A witness-directed prepared-library proposal can use a different source operator
and change other bonds. It must not be attributed as a direct transfer of the
fragment hit. Whole-source compilation may fail even when a local witness is valid.
Those failures remain evidence gaps, not declarations of synthetic impossibility.

## Outputs and reproducibility

`report.json` contains the full comparison and saved artifact references.
`report.html` draws proposed transformations and selected source reactions, and
includes compiler rejection details. The custom evaluation execution directory
also retains `selected_source_operators.json`, including an empty library if all
selected sources fail admission.

Every completed fragment search is repeated; comparison excludes only execution
measurements. Precursor searches are run twice and compared with their stage
diagnostics. Search deadlines remain wall-clock limits: equality on a completed
run does not guarantee equality under timeout, on other machines, or after a
version change. Incomplete searches cannot seed this experiment.

The arms have different total work budgets. Candidate counts are descriptive,
not an equal-cost effectiveness benchmark. No conditions, stock availability,
selectivity, precursor preparation or complete-route feasibility are established.
The main system's chemistry-review and untouched-evaluation release gates remain.

Regression coverage includes admitted source transfer, whole-source rejection
despite a local witness, unresolved and incomplete evidence, symmetric alignment
hypotheses, conflicting projected bonds, policy gates and repeated results.

## Findings for the default target

The recorded local run is
[`fragment_guided_retro_20261003_final`](../../results/ai_native/fragment_guided_retro_20261003_final/report.html).
The three searches completed and found respectively 14, 5 and 4 observations
with construction evidence. These groups overlap; their counts must not be added
to estimate distinct support. The returned sample yielded two observations with
fully resolved embeddings, supporting three target-bond construction hypotheses.

Witness-directed prepared-library searches returned four distinct precursor sets,
including a ketone plus hydrazine and a protected hydrazone precursor. Three
sets were absent from the unrestricted baseline's top three. All returned
proposals passed signature validation, and focused proposals additionally passed
formed-bond verification. Search and precursor-generation repeat checks matched.
The different total budgets prevent interpreting this as an effectiveness gain.

Direct source transfer produced no admitted operators. Both selected reactions
have valid observations and verified product completeness, but whole-core quality
is `review`: `not_all_edits_graph_checked`. Their checked-edit fractions are
6/9 and 7/9. The compiler's existing `pass_only` gate therefore rejects them as
`materialized_core_not_verified`. This is a validator-coverage gap, not evidence
that the reported chemistry is impossible. Additional recorded inspections are
saved in `source_gate_diagnostics.json` inside that investigation.

Integration also required restoring long source SMILES from the discovery
index's chunked-text representation. The POC checks offsets, length and the
source SHA-256 before compilation; truncated or conflicting text cannot be
treated as a complete source structure. No chemistry definitions, core admission
policies, corpus records or public application APIs were changed.
