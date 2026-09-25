# AI-Native Scientific Tools: Phased Implementation Plan

Status: local workspace, development pilots, and user-testable agent chat implemented; scientific release gates pending
Revision: 2

Date: 2026-09-23

Implementation progress, validation, and outstanding gates are tracked in
[Scientific_Workspace_Implementation_Status.md](Scientific_Workspace_Implementation_Status.md).
Use the [workspace quickstart](Scientific_Workspace_Quickstart.md) for the
implemented Python/CLI interface and optional local browser conversations. Phase
acceptance gates below are not all met. The user-testable development milestone
now supports natural-language questions through the installed Codex runtime,
saved evidence, and follow-ups; this does not waive Phase 6 or release gates.

The next development slice, structured scientific answers, is implemented as
`scientific_answer.v2`: explicit molecules, steps, routes, attributed conditions
and yields, and captured source references. Reaction and dependency drawings
are views of agent-authored objects, not newly verified scientific records.
The condition-investigation development phase now adds structured precedent
inspection, attributed adaptation records, corrected direct-assessment semantics,
recorded custom Python execution and one bounded answer correction. Transfer
validity still needs source-supported trials and independent review. Broader live
route revision and comparative/untouched evaluation remain subsequent priorities.

## 1. Objective and governing design

Implement the [AI-native scientific tools design](AI_Native_Scientific_Tools_Design.md)
as a reusable scientific environment for advanced agents. The initial scope is
reaction analysis, dataset access, condition recommendation, and retrosynthesis.

The first useful milestone is a recorded scientific investigation using the
current Codex environment, existing Python and CLI operations, and custom
analysis where useful. Prepare and exercise that environment before building a
broad tool interface. This revision moves the pilot ahead of contract expansion
and MCP implementation; later phases formalize what the pilot demonstrates.

The [current chemistry implementation roadmap](../new/type_agnostic_reaction_recommendation_implementation.md)
remains authoritative for package boundaries, chemistry validation, dataset
admission, and production release. This plan does not waive its release gates.
Hardware execution is out of scope.

## 2. Architectural commitments

- Keep chemistry rules in their owning packages. Application adapters expose
  operations and manage artifacts without creating a second scientific path.
- Reuse `reactive_taxonomy`, `condition_registry`, `condition_recommender`, and
  the existing `core_retrosynthesis` operations. Keep `chem_coworker` thin.
- Support both structured operations and a programmable Python/shell workspace.
  Start with existing APIs and CLI access, then harden demonstrated needs and
  add thin MCP access for clients that need it. Do not require an internal LLM
  call before an external agent can use a scientific operation.
- Treat Codex as the first integration and evaluation client, not a dependency
  of scientific packages. Keep provider configuration in optional adapters.
- Keep evidence, hypotheses, proposed adaptations, and observed records distinct.
  Agent proposals never bypass admission or overwrite contradictory evidence.
- Preserve canonical recommendation behavior, including compatibility before
  ranking, resolved recipes, provenance, fallback reporting, and abstention.
- Allow unsupported chemistry to remain investigable without labeling it
  verified or experimentally feasible.
- Save custom code and its inputs, outputs, and assumptions. Changes to
  scientific rules belong in separately versioned development work; an agent
  must not change its own scientific baseline during fixed-baseline evaluation.
- Use one shared service layer for future web and MCP interfaces. A browser
  application may call that layer directly; MCP is not required internally.

## 3. Phase overview and dependencies

| Phase | Outcome | Dependency |
| --- | --- | --- |
| 0 | Recorded baseline, prepared workspace, capability inventory, and evaluation protocol | Existing code and artifacts |
| 1 | Real investigations using the existing agent and environment | Phase 0 development readiness |
| 2 | Hardened scientific contracts and thin MCP access where needed | Phase 1 findings |
| 3 | Persistent, resumable investigations | Phase 2 |
| 4 | Evidence-backed condition investigations | Phase 3 |
| 5 | Iterative retrosynthesis with route revisions | Phase 3 |
| 6 | Measured comparison and untouched evaluation | Phases 4 and 5 |
| 7 | Consolidated interfaces and deployment readiness | Phase 6 and roadmap release gates |

Phases 4 and 5 can proceed independently and mature the early pilot tasks.
Initial artifact recording starts in Phase 1; Phase 3 formalizes storage and
resumption. A pilot or interface prototype is a development milestone, not a
production chemistry-validation claim. Advance by acceptance criteria rather
than calendar estimates.

## 4. Phase 0: Baseline and workspace preparation

### Work

1. Record code revision, schemas, definition hashes, dataset checksums, index
   identity, and the validation status of each artifact.
2. Run the full test suite and taxonomy, registry, conversion, and index checks
   applicable to the selected baseline.
3. Investigate the two registry failures observed during the preceding design
   review: `test_validate_reports_current_registry_state` and
   `test_registry_audit_reconciles_all_rows` reported two undeclared ambiguous
   identifiers. Recheck current state; this observation is historical and must
   not be assumed to describe a later checkout.
4. Classify existing operations as directly reusable, coupled to the assistance
   controller, or missing. Record owners and evidence requirements.
5. Select development cases spanning straightforward chemistry, ambiguity,
   conflicting evidence, unsupported transformations, condition-transfer
   problems, and multistep route revision.
6. Reserve untouched evaluation cases before development. Follow the roadmap's
   reference/reaction-connected partitioning and leakage controls.
7. Verify the Python/RDKit environment and runnable entry points for the selected
   tasks. Record unavailable datasets, model adapters, services, and credentials
   as limitations. Do not assume every module is a usable capability.
8. Write a short entry guide with working invocations, dataset locations and
   coverage, result schemas, and evidence rules. Create a per-investigation
   output directory and a minimal manifest/report convention. Reuse the current
   environment; containerization is optional at this stage.

### Deliverables

- A baseline manifest and validation report.
- A capability inventory mapping operations to existing code and contracts.
- Development tasks and a predefined evaluation rubric.
- A separately controlled untouched evaluation partition.
- A working scientific environment, entry guide, and pilot artifact convention.

### Acceptance gate

For development readiness, a fresh session can execute the selected operations,
locate the intended data, and save inspectable results. Known defects and affected
operations are recorded; affected outputs cannot be presented as validated.
Development cases and untouched cases are separated.

Baseline defects must be understood and resolved before freezing a release
candidate. A pilot on a documented development snapshot does not satisfy that
release gate. No expected snapshot is changed merely to match current output.

## 5. Phase 1: Investigate with the existing agent

### Work

1. Run the canonical workflow on a development case with a condition-transfer
   mismatch. Record the baseline output and unresolved scientific question.
2. Ask the current agent to investigate using existing APIs, full precedent
   records, source inspection, and custom Python analysis as needed. Record
   external research separately from local-corpus evidence.
3. Attempt a second task requiring reconsideration of an earlier retrosynthetic
   choice after a downstream problem. Preserve unsupported steps and unsuccessful
   attempts as explicit findings.
4. Save inputs, calls, scripts, versions, outputs, evidence links, hypotheses,
   decisions, and unresolved questions from the start. Use files before building
   a general investigation database.
5. Record friction: API discovery, missing procedure fields, misleading statuses,
   unavailable assessments, slow queries, and repeated custom calculations.
6. Have a chemist assess whether the investigation added supported insight.
   Use this feedback to prioritize the next phase; do not tune against untouched
   cases or require a positive result to preserve the pilot record.

The agent may choose its sequence of actions. This phase needs no new UI,
internal LLM controller, or MCP server. If a task reveals a scientific defect,
record it and handle the fix as separate versioned work; preserve the original
run. Do not edit baseline chemistry or admission rules to obtain a favorable
pilot outcome.

### Deliverables

- Two investigation records, including useful findings and documented limits.
- Saved custom analyses with reproducible inputs and outputs.
- A prioritized backlog linking proposed contracts to observed scientific needs.

### Acceptance gate

The records show what the agent could investigate, which conclusions are
supported, and where science or access was insufficient. Deterministic
calculations can be replayed, and observations remain distinct from proposals.
This is the first reviewable development milestone; it does not claim a
validated improvement over the baseline.

## 6. Phase 2: Harden scientific contracts and portable access

### Work

Define typed, independently callable operations around existing domain APIs,
prioritized by Phase 1 findings. Reuse useful code in
`chem_coworker/assistance/` without requiring its controller or provider transport.

| Area | Initial operation scope |
| --- | --- |
| Analysis | Analyze a molecule; observe a reaction; inspect transformation evidence. |
| Datasets | Search precedents; retrieve complete records and available procedures. |
| Conditions | Resolve a recipe; recommend conditions; assess a supplied recipe. |
| Retrosynthesis | Generate candidates; retrieve step precedents; verify a step or route. |

Each operation specifies schemas, units, assumptions, provenance, versions,
limitations, and an inspectable result reference. Support partial observations
when interpretation fails. Preserve convenient complete workflows and composable
operations over the same implementation, alongside programmable workspace access.

Audit status semantics. In the inspected code, `assess_reaction_recipe()`
classifies a missing verified signature as a hard conflict. Determine whether
this represents inability to assess or an actual chemistry violation. Correct
the owning contract where necessary, with migration and regression coverage;
do not silently reinterpret the result in MCP. Distinguish pass, fail, unknown,
unsupported, and execution error through explicit domain-status mappings.

For a selected client that needs portable tool access, implement one thin MCP
server over the useful contracts. Supply clear descriptions, structured failures,
pagination, limits, concise summaries, and full-evidence access. Verify parity
against direct calls on the pilot cases. Keep client-specific configuration in
optional adapters. If no client yet needs MCP, record its deferral and validate
the Python/CLI contracts directly; it is not a prerequisite for further science.

Begin locally. A later remote transport must preserve scientific contracts.
Custom server-side execution, if needed by a remote client, is a separate
capability with per-investigation workers; publishing ordinary scientific MCP
tools does not automatically grant Python or filesystem access.

### Deliverables

- Public operation contracts, examples, and an updated environment guide.
- Contract, architecture, and direct-call parity tests.
- A thin MCP adapter and integration results if required by the selected client,
  otherwise a documented deployment deferral.

### Acceptance gate

An external agent can use the selected operations without the built-in
assistance loop. Interfaces preserve authoritative behavior, versions,
provenance, admission status, and abstention. Mandatory chemistry filters remain
inside canonical workflows. Domain changes have positive, negative, ambiguous,
conflict, and invariance regressions where applicable.

## 7. Phase 3: Persistent investigation state

### Work

Add application-owned, revisioned records for investigations, hypotheses,
candidate recipes, route alternatives, assessments, decisions, and unresolved
questions. Reference existing domain objects instead of copying them into a
universal reaction dictionary.

Evolve these contracts from the pilot's saved files. Preserve and migrate useful
pilot evidence instead of requiring investigations to restart in a new format.

Persist tool inputs and result references, dataset and definition versions, and
model metadata when relevant. Store concise decision rationales and evidence
links; private model reasoning is not required.

Support saving, resuming, and retrieving an investigation independently of a
provider's conversation format. Preserve cancellation and explicit completion,
insufficient-evidence, failure, and budget-exhaustion states. Add resumable job
references for expensive operations when required.

### Deliverables

- Investigation contracts, storage, and artifact-access operations.
- Tests for serialization, revisions, resumption, and deterministic replay.
- A development case resumed by a fresh agent session.

### Acceptance gate

A fresh session can recover the objective, constraints, alternatives, evidence,
and open questions. Hypotheses remain distinct from source observations. Recorded
deterministic computations are reproducible; fresh model decisions need not be.

## 8. Phase 4: Condition investigations

### Work

Use a substrate mismatch between a query and its nearest precedents as the first
task. Enable the agent to inspect complete recipes and available procedures,
compare structural environments and outcomes, investigate contextual component
roles, identify counterevidence, and propose an explicitly labeled adaptation.

Submit proposed recipes to available checks without promoting them to observed
recipes. Add missing domain assessments only when required by the task and
supported by explicit definitions and regressions.

Preserve whole-recipe co-occurrence, operating variants, independent support,
and missingness. Do not invent quantities, temperature, concentration, time, or
addition order. Do not treat unreported outcomes as failures or association as
proof of mechanism.

### Deliverables

- Recipe comparison and targeted assessment capabilities.
- A condition investigation artifact containing observations, proposed changes,
  supporting and conflicting evidence, assumptions, and unresolved concerns.
- Blind-review development examples.

### Acceptance gate

Independent review finds useful condition-transfer reasoning. Proposed changes
are traceable and labeled; incompatibilities and unknowns are preserved. A
compatibility assessment is not represented as an outcome prediction.

## 9. Phase 5: Iterative retrosynthesis

### Work

Reuse existing generation, external-proposal assessment, route verification, and
repair operations. Add missing interfaces for inspecting an intermediate,
assessing a proposed step, and revising an individual branch.

Let the agent maintain alternative routes, investigate weak steps, compare step
orders, and reconsider an earlier disconnection. Retain unsupported proposals as
hypotheses outside verified admission.

Keep structural consistency, experimental support, selectivity, operating
conditions, and route practicality as distinct assessments. Route revisions must
preserve history and trigger the relevant checks on affected dependencies.

### Deliverables

- Versioned route graphs with alternatives and step assessments.
- Targeted inspection and revision operations.
- A development investigation that revises an earlier choice after discovering
  a downstream problem.

### Acceptance gate

The revision addresses an evidenced problem and remains structurally consistent
within the available checks. The final artifact exposes unresolved feasibility
and selectivity questions rather than hiding them behind a single score.

## 10. Phase 6: Comparative evaluation

### Work

Compare three systems:

1. The deterministic workflow.
2. The deterministic workflow with LLM review.
3. An advanced agent independently using the scientific toolkit.

Within the agent condition, compare structured tools alone against the same
tools plus programmable workspace access. Separate a controlled comparison
with fixed scientific code, data, and source evidence from an open-web condition.
Record external sources and check for leakage of evaluation answers. Additional
evidence or compute must not be misattributed to runtime quality.

Use the same underlying model where possible, record tool/data access, compare
matched resource budgets, and repeat agent runs to measure variation. Evaluate
another compatible client to test portability; do not assume equivalent quality
from protocol compatibility alone.

Blind chemistry review should assess plausible and supported improvements,
missed contradictions, unsupported claims, appropriate abstention, and whether
citations substantiate the associated claim. Report latency, tokens, cost,
execution failures, and recoverability alongside scientific outcomes.

Resolve development disagreements before running untouched evaluation. Follow
the governing roadmap's gate order and predeclared acceptance criteria.

### Deliverables

- Recorded comparative runs and numerical reports.
- Blind-review packets and adjudication records.
- An untouched evaluation report tied to exact versions and data partitions.

### Acceptance gate

The evidence demonstrates an improvement on predefined scientific criteria
without violating domain boundaries or the roadmap's incompatibility gates.
Longer answers and more tool calls do not count as improvement by themselves.

## 11. Phase 7: Consolidation and deployment readiness

### Work

- Make internal assistance, CLI, web, and MCP adapters use the same scientific
  capabilities. Remove temporary duplicate adapters once parity is established.
- Publish stable contracts, examples, supported limitations, and upgrade rules.
- Complete full-corpus and persisted-artifact validation required by the roadmap.
- Keep optional runtime integrations replaceable and provider credentials outside
  scientific packages and browser code.
- For remote use, add authenticated access, project isolation, background jobs,
  cancellation, and persistent investigation storage.
- Package the prepared environment reproducibly, with pinned dependencies and
  identified code/data snapshots. Use separate writable investigation workers
  and shared read-only validated data where appropriate.
- For hosted ChatGPT access, expose supported remote tools. Do not describe this
  as replacing ChatGPT's built-in sandbox. Programmable execution on our server
  requires an explicit service or a compatible runtime integration; see the
  design's deployment choices and official documentation links.
- If a dedicated website is selected, expose the shared services through a web
  API and provide structure, route, recipe, and evidence views. Run its selected
  agent runtime on the backend. The website and remote MCP endpoint share domain
  implementations.

### Deliverables

- A consolidated toolkit and migration/removal notes.
- Release validation and operational documentation.
- Optional remote MCP and web deployments after their requirements are selected.

### Acceptance gate

Required machine and human release gates pass. Public interfaces reach the
canonical workflows, saved investigations remain interpretable, and deployment
does not create a second chemistry implementation.

## 12. First implementation slice and verification policy

The first slice should prepare the existing environment and run one recorded
condition-transfer investigation with the current agent. Use existing analysis,
retrieval, full-record inspection, and recommendation operations directly. Save
any custom scripts and evidence, then attempt the route-revision pilot. Use
observed obstacles to choose contract improvements. A new MCP server, a custom
agent runtime, and a dedicated website are not prerequisites for this slice.

For implementation changes, run focused tests and the complete `pytest -q` suite
before handoff. Validate loaders and schemas when definitions change. Document
chemistry and dataset impacts explicitly; never update snapshots merely because
outputs changed. Keep new experiments on development data until the untouched
evaluation gate.

This plan defines proposed work. Passing an interface or deployment test does
not establish scientific accuracy or complete a chemistry release gate.
