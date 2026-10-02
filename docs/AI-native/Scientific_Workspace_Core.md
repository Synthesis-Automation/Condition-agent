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
| Scientific call recording, execution outcomes, evidence links and replay | `workspace.py`, `store.py`, `execution.py` |
| Scientific code, definitions, dependency and dataset identity | `baseline.py` |
| Capability discovery and explicit invocation contract | `operation_contracts.py` |
| Domain-specific adapters, prerequisites and replay projections | `operations.py` and its existing adapter modules |
| Turn lifecycle, cancellation, runtime continuation and answer publication | `conversation.py`, `agent_runtime.py` |
| Per-turn application identity and immutable instruction snapshots | `context.py` |
| Prompt composition | `prompts.py` |
| Shared investigator rules and workspace usage | `instructions/core.md`, `instructions/workspace_usage.md` |
| Optional conditions and retrosynthesis strategy | `guides/*.md` |
| Optional procedural memory | `learning.py`, `instructions/learning.md` |
| Answer authoring contract and current default presentation | `instructions/answer_authoring.md`, `presentation/default.md` |
| Scientific attribution validation and browser rendering | Existing answer validators and `app/web_api/scientific_presentation.py` |

The runtime chooses actions. The core does not classify a question into a reaction
family, enforce a task checklist, or implement scientific scoring. Observation,
interpretation, compatibility, admission, and uncertainty keep their domain meaning.
Execution completion is recorded separately from chemical validity.

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
belongs in presentation. The existing answer schema, validators and browser layout
are preserved in this phase. User-controlled UI layout and simpler task views are
the next presentation phase, not a new scientific API.

## Task playbooks

The task layer consists of two main Markdown guides and four supporting guides.
Each uses the same four sections: purpose and context, scientific questions,
evidence distinctions, and stopping criteria. Main guides describe the decision;
supporting guides explain a specialized investigation only when relevant.

| Main task | Optional supporting advice |
| --- | --- |
| [Conditions](../../chem_coworker/scientific_workspace/guides/conditions.md) | [Screening](../../chem_coworker/scientific_workspace/guides/conditions_screening.md) |
| [Retrosynthesis](../../chem_coworker/scientific_workspace/guides/retrosynthesis.md) | [Fragment precedents](../../chem_coworker/scientific_workspace/guides/retrosynthesis_fragments.md), [route revision](../../chem_coworker/scientific_workspace/guides/retrosynthesis_revision.md), [forward challenges](../../chem_coworker/scientific_workspace/guides/retrosynthesis_forward.md) |

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
