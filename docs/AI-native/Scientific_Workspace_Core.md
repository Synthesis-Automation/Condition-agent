# Scientific workspace core

Status: core separation implemented; scientific release gates remain unchanged.
Date: 2026-10-02

## Purpose and boundary

The core connects a capable agent to a prepared, recorded scientific environment.
The agent chooses its investigation strategy, uses domain tools, inspects sources,
and writes custom analysis inside the investigation. Chemistry remains in the
owning scientific packages. Task playbooks supply optional advice; presentation
profiles supply default writing and display instructions.

| Responsibility | Owner |
| --- | --- |
| Scientific call recording, execution outcomes, evidence links and replay | `workspace.py`, `core/store.py`, `core/execution.py` |
| Scientific code, definitions, dependency and dataset identity | `core/baseline.py` |
| Capability discovery and explicit invocation contract | `core/operation_contracts.py` |
| Domain-specific adapters, prerequisites and replay projections | `adapters/operations.py` and `adapters/*.py` |
| Turn lifecycle, cancellation, runtime continuation and answer publication | `runtime/conversation.py`, `runtime/agent_runtime.py` |
| Per-turn application identity and immutable instruction snapshots | `agent_context/context.py` |
| Prompt composition | `agent_context/prompts.py` |
| Shared investigator rules and workspace usage | `agent_instructions/core.md`, `agent_instructions/workspace_usage.md` |
| Optional conditions and retrosynthesis strategy | `task_playbooks/*.md` |
| Optional procedural memory | `agent_context/learning.py`, `agent_instructions/learning.md` |
| Answer authoring contract and current default presentation | `agent_instructions/answer_authoring.md`, `presentation/default.md` |
| Scientific attribution validation and browser rendering | `answers/*.py`, saved projections in `views/*.py`, and `app/web_api/scientific_presentation.py` |

The runtime chooses actions. The core does not classify a question into a reaction
family, enforce a task checklist, or implement scientific scoring. Observation,
interpretation, compatibility, admission, and uncertainty keep their domain meaning.
Execution completion is recorded separately from chemical validity.

## Package organization

The package root retains `ScientificWorkspace`, the public exports and the
`python -m chem_coworker.scientific_workspace` CLI. Internal modules have one
canonical location; the removed flat modules have no forwarding wrappers.

| Subpackage | Responsibility |
| --- | --- |
| `core/` | Baseline identity, recording, artifact storage, execution and capability contracts |
| `adapters/` | Domain calls, pinned catalogs, exact saved-step selection and source capture |
| `runtime/` | Agent processes, conversation lifecycle, activity and runtime settings |
| `agent_context/` | Per-turn context, prompt composition, frozen guides and procedural memory |
| `answers/` | Typed answer contracts, handoff, finalization and evidence attribution validation |
| `views/` | Shared bounded projections, detailed/brief summaries, artifact inspection and saved precedent views |

The names distinguish three different responsibilities:

| Previous directory | Current directory | Contents and purpose |
| --- | --- | --- |
| `guidance/` | `agent_context/` | Python code that snapshots turn context, composes prompts and recalls procedural lessons; see `context.py`, `prompts.py` and `learning.py` |
| `guides/` | `task_playbooks/` | Optional Markdown strategies for particular investigations, such as conditions, screening and retrosynthesis |
| `instructions/` | `agent_instructions/` | Shared Markdown rules for agent behavior, workspace usage, answer authoring and procedural memory |

For example, `agent_context/prompts.py` composes a prompt from the shared rules
in `agent_instructions/core.md`, optionally supplemented by the investigation
strategy in `task_playbooks/conditions.md`. The instructions apply across tasks;
the playbook supplies optional questions and stopping advice. Neither owns the
chemistry rules enforced by the scientific packages.

`presentation/` remains the home of default writing/display preferences.
The `task_guide()` API and versioned `guides`/`guidance` data fields retain their
names; serialized contracts are independent of directory names. Recorded old
instruction snapshots remain readable without live-file fallback, and v1
verification remains strict for original `guides/` snapshots as well as renamed
playbooks. Current source directories and imports use only the new names.
`paths.py` supplies installed resource and repository locations so moving adapters
does not change worker imports, default repository selection or resource fallback.

Shared child-process options and process-tree termination live in
`core/process_utils.py`. Scientific scripts and workers no longer import the
agent harness for these helpers. This module is covered by scientific identity;
changing execution behavior requires a new investigation. Runtime and guidance
ownership remains an explicit file inventory, including the new nested paths;
an unknown Python module is still scientific regardless of its folder.

Step precedent lookup lives in `adapters/step_precedents.py`, exact saved-call
selection in `adapters/step_selection.py`, and attribution/mandatory support
validation in `answers/step_precedents.py`. Condition attribution is in
`answers/condition_precedents.py`. Display recovery and source-observation views
are in `views/precedents.py` and `views/condition_precedents.py` and do not invoke
scientific tools. Detailed and brief summaries keep their existing schemas and
bounds; `views/artifact_inspection.py` keeps paged evidence inspection separate.

This reorganization changes code paths and scientific baseline identity. Restart
the server and create a new investigation to run or replay scientific calls with
this code. Historical evidence, answers and conversation rendering remain
readable; recorded baselines and artifacts are never rewritten or automatically
migrated. Domain schemas, chemistry definitions and dataset artifacts do not
change or require rebuilding.

## Stable capability contract

`OperationProvider` exposes three methods: `definition(name)`, `catalog()`, and
`invoke(name, arguments)`. The default provider is `ScientificOperations`.
`ScientificWorkspace(root, operations=provider)` accepts another explicitly supplied
Python provider using the same recording and replay interface.

Immutable `OperationDefinition` declarations specify a name, contract version,
required artifacts, evidence-reference arguments, optional usage policy, optional
execution-status field, and scientific replay projection. Domain schemas remain
owned by the domain packages; this phase adds contract identity and argument
discovery rather than a second schema system.

The catalogue preserves `name`, `signature`, and `description`, and adds declaration
metadata. Prerequisites describe required inputs; presence does not establish their
validity. Conditional requirements remain checked by the invoked adapter.
Executable behavior is registered Python code, never imported from arbitrary JSON.
New adapters belong in baseline-identified source files. A changed implementation
or contract must not be introduced silently through an external provider.

The workspace uses declarations to link verified input evidence, retain incomplete
execution outcomes, and compare replay results. These behaviors contain no tool-name
branches. Forward and fragment telemetry exclusions live with those adapters and
preserve the previous comparison scope. Unknown/private invocation is rejected.
Declared execution status must be `completed`, `error`, `timed_out` or `cancelled`.
A missing or invalid declared status is recorded as an error with the returned
payload retained; it cannot become completed computation evidence.

Calls record `operation_contract_version` and the available `scientific_identity`.
Replay rejects a changed declared contract before invocation. Historical successful
calls without contract metadata retain their existing replay behavior under their
original baseline checks; saved artifacts are never rewritten.

## Scientific identity and application context

New investigations use `scientific_baseline.v2`. Scientific identity covers code,
definitions, environment versions and selected data identities. Explicitly owned
runtime, guide, instruction and presentation files are recorded independently.
Unknown new Python adapters remain scientific by default. Execution, storage and
evidence-validation changes still invalidate the scientific baseline.

At each conversation turn, `scientific_application_context.v1` records the current
application layers, exact Markdown resources, runtime settings and optional advice.
The turn and saved answer reference this immutable event. Each runtime attempt saves
its exact composed input as `investigation-prompt.txt`; the runtime adapter may add
its transport instructions in its own `prompt.txt`. Attempt files are hashed in the
saved answer trace. Application-context records are configuration, never scientific
evidence or citations.

Guide changes take effect at the next explicitly recorded turn. Within a turn,
composition and `task_guide` use saved resources. Procedural lessons stay pinned to
the investigation; lessons published during a turn are available to later
investigations. Headless workspaces retain their initial frozen guidance until
application context is explicitly recorded.

Changed application layers or runtime settings start a new agent thread over saved
investigation history. Unchanged context resumes the existing thread. This preserves
evidence while preventing a previous instruction set from silently carrying forward.
Scientific code/data changes require a new investigation; linked baseline upgrades
and automatic scientific-environment restoration are outside this phase.

Explicit `scientific_baseline.v1` records retain their original strict code/guide
verification scope. They are not silently migrated to v2. Historical answers and
artifacts remain readable; continuing or replaying still requires their environment.

## Prompt ownership

The composer joins common instructions, available capabilities, workspace usage,
optional selected guides, optional memory instructions, answer authoring and a
presentation profile. It discovers guide names from recorded context; it does not
route by question keywords. Task guides can be read lazily with `w.task_guide(name)`.
An explicit client can include guides through `task_names=(...)`.

Every substantive instruction has one owner. Task-specific strategy belongs in its
guide, invocation constraints belong in declared capability policy, and formatting
belongs in presentation. The [presentation layer](Scientific_Workspace_Presentation.md)
owns the chemist-facing cards and default output style. Optional attributed rationale
and saved condition-inspection links extend the answer view without adding a new
scientific API or changing the agent execution loop.

## Task playbooks

The task layer consists of two main Markdown guides and four supporting guides.
Each uses the same four sections: purpose and context, scientific questions,
evidence distinctions, and stopping criteria. Main guides describe the decision;
supporting guides explain a specialized investigation only when relevant.

| Main task | Optional supporting advice |
| --- | --- |
| [Conditions](../../chem_coworker/scientific_workspace/task_playbooks/conditions.md) | [Screening](../../chem_coworker/scientific_workspace/task_playbooks/conditions_screening.md) |
| [Retrosynthesis](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis.md) | [Fragment precedents](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_fragments.md), [route revision](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_revision.md), [forward challenges](../../chem_coworker/scientific_workspace/task_playbooks/retrosynthesis_forward.md) |

The existing discovery mechanism advertises all guide names; text is read through
`w.task_guide(name)` or explicitly included with `task_names`. A main guide names
its supporting guides, without automatically loading their text. All text is part
of the recorded context. No task registry, keyword router or workflow engine is added.
Procedural lesson categories remain `general`, `conditions` and `retrosynthesis`.

Task advice owns questions, transfer judgments, optional heuristics and stopping
decisions. Tool signatures, prerequisites and limits remain in the capability
catalogue and [tool reference](readme.md#available-operations-and-ownership).
Mandatory chemistry and attribution checks remain in domain validators and the
shared answer contract. Presentation profiles own default answer formatting.
Avoid copying arguments, numerical limits, answer JSON or runtime instructions
into task guides. Evolve each owner directly so the core can remain stable.

The task cleanup changes Markdown guidance and documentation only. Scientific
implementation, contracts, evidence validators and UI layout retain their behavior.
V2 investigations can receive these guides at the next recorded application turn;
headless and legacy baseline behavior remains as described above.

## Extending the system

1. Add scientific behavior to its owning package with applicable chemistry tests.
2. Add a thin adapter and an explicit capability declaration. Version breaking API
   changes and substantive scientific behavior; identify changed data snapshots.
3. Add or revise a task guide when useful. A new guide needs no core task switch.
4. Update presentation independently; retain attributed underlying evidence.
5. Preserve historical artifacts. Remove obsolete runtime paths after parity is
   established, with clear removal criteria for any temporary compatibility.

Focused regressions cover new capability extension, contract drift, partial outcomes,
scientific versus application drift, frozen context, prompt ownership and continuation.
Run the complete suite with `pytest -q` before handing off changes. These architecture
checks do not satisfy independent chemistry-review or untouched-evaluation gates.

## Verification on 2026-10-02

The original core-separation `pytest -q` run finished with **2,271 passed, 2 skipped and 2 failed**.
At that point, the failures were `test_validate_reports_current_registry_state` and
`test_registry_audit_reconciles_all_rows` in `tests/condition_registry/`: both
reported two undeclared ambiguous identifiers. Both failures were reproduced using
the original packages and tests from revision
`2c49cc620a0ef6f152c07fb145c4abb497aff626` in an isolated snapshot. The core changes
do not change those registry definitions or relax validation.

Core, conversation, activity, scientific-domain and application regressions passed
apart from those existing registry checks. Ruff checks passed for the changed
Python modules and tests. Historical evidence remains readable; start a new
investigation after this scientific-core upgrade to use the v2 baseline.

The registry defects were subsequently resolved by the separate
[BINAP curation correction](../new/binap_registry_identity_correction_20261002.md).
Both audit tests remain in place; see that correction for current validation results.

### Package organization verification on 2026-10-02

After the module moves and responsibility splits, the complete `pytest -q` run
finished with **2,298 passed and 2 skipped**. The skipped PDF-parser tests require
the optional `pypdf` dependency. Workspace lint checks and whitespace checks passed.
Both root CLIs load; the research application boots and its scientific UI,
configuration and conversation-list endpoints return HTTP 200.

New regressions require process-control changes to invalidate scientific identity,
retain conservative identity for unclassified modules in nested application
folders, and prevent scientific execution from importing the agent runtime.
Existing summary, attribution, saved-answer, worker, baseline and chemistry
regressions remain in place. These engineering checks do not satisfy independent
chemistry-review or untouched-evaluation release gates.

The directory-naming follow-up (`agent_context/`, `task_playbooks/`, and
`agent_instructions/`) passed the complete suite with **2,301 passed and 2 skipped**
on 2026-10-02. The same optional PDF dependency accounts for the skips. Added
regressions preserve original instruction snapshots, reject mixing old resources
into incomplete current snapshots, and retain strict v1 verification of original
guide paths. Workspace lint checks passed and all 21 workspace documentation links
resolve after the renames.
