# AI-Native Scientific Tools: Phased Implementation Plan

Status: proposed; phases below are not marked as implemented  
Date: 2026-09-22

## 1. Objective and governing design

Implement the [AI-native scientific tools design](AI_Native_Scientific_Tools_Design.md)
as a reusable scientific environment for advanced agents. The initial scope is
reaction analysis, dataset access, condition recommendation, and retrosynthesis.

The first useful milestone is an external agent independently investigating a
reaction through existing scientific capabilities. Persistent investigations,
condition adaptation, and iterative route revision follow.

The [current chemistry implementation roadmap](../new/type_agnostic_reaction_recommendation_implementation.md)
remains authoritative for package boundaries, chemistry validation, dataset
admission, and production release. This plan does not waive its release gates.
Hardware execution is out of scope.

## 2. Architectural commitments

- Keep chemistry rules in their owning packages. Application adapters expose
  operations and manage artifacts without creating a second scientific path.
- Reuse `reactive_taxonomy`, `condition_registry`, `condition_recommender`, and
  the existing `core_retrosynthesis` operations. Keep `chem_coworker` thin.
- Expose a runtime-independent Python interface and a thin MCP interface.
  Preserve CLI access where useful. Do not require an internal LLM call before
  an external agent can use a scientific operation.
- Treat Codex as the first integration and evaluation client, not a dependency
  of scientific packages. Keep provider configuration in optional adapters.
- Keep evidence, hypotheses, proposed adaptations, and observed records distinct.
  Agent proposals never bypass admission or overwrite contradictory evidence.
- Preserve canonical recommendation behavior, including compatibility before
  ranking, resolved recipes, provenance, fallback reporting, and abstention.
- Allow unsupported chemistry to remain investigable without labeling it
  verified or experimentally feasible.
- Use one shared service layer for future web and MCP interfaces. A browser
  application may call that layer directly; MCP is not required internally.

## 3. Phase overview and dependencies

| Phase | Outcome | Dependency |
| --- | --- | --- |
| 0 | Reproducible baseline, capability inventory, and evaluation protocol | Existing code and artifacts |
| 1 | Independent scientific operation contracts | Phase 0 |
| 2 | External agent accesses and composes the tools | Phase 1 |
| 3 | Persistent, resumable investigations | Phase 2 |
| 4 | Evidence-backed condition investigations | Phase 3 |
| 5 | Iterative retrosynthesis with route revisions | Phase 3 |
| 6 | Measured comparison and untouched evaluation | Phases 4 and 5 |
| 7 | Consolidated interfaces and deployment readiness | Phase 6 and roadmap release gates |

Phases 4 and 5 can proceed independently. An interface prototype is a development
milestone, not a production chemistry-validation claim. Advance by acceptance
criteria rather than calendar estimates.

## 4. Phase 0: Baseline, inventory, and evaluation tasks

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

### Deliverables

- A baseline manifest and validation report.
- A capability inventory mapping operations to existing code and contracts.
- Development tasks and a predefined evaluation rubric.
- A separately controlled untouched evaluation partition.

### Acceptance gate

Baseline defects are understood and resolved before freezing a release
candidate. Artifacts are version-consistent and reproducible. Development cases
and untouched cases are separated. No expected snapshot is changed merely to
match current output.

## 5. Phase 1: Independent scientific contracts

### Work

Define typed, independently callable operations around existing domain APIs.
Reuse useful code in `chem_coworker/assistance/` without requiring its controller
or provider transport.

| Area | Initial operation scope |
| --- | --- |
| Analysis | Analyze a molecule; observe a reaction; inspect transformation evidence. |
| Datasets | Search precedents; retrieve complete records and available procedures. |
| Conditions | Resolve a recipe; recommend conditions; assess a supplied recipe. |
| Retrosynthesis | Generate candidates; retrieve step precedents; verify a step or route. |

Each operation must specify input schema, units, assumptions, result schema,
provenance, versions, limitations, and an inspectable result reference. Support
partial observations when interpretation fails.

Audit status semantics before exposing them. In the inspected code,
`assess_reaction_recipe()` classifies a missing verified signature as a hard
conflict. Determine whether this represents inability to assess or an actual
chemistry violation. Correct the owning contract where necessary, with migration
and regression coverage; do not silently reinterpret the result in MCP.

Distinguish check outcomes such as pass, fail, unknown, unsupported, and execution
error. Preserve and explicitly map existing domain statuses. None of these
outcomes alone establishes experimental success.

### Deliverables

- Public scientific operation contracts and thin application adapters.
- Focused contract, architecture, and parity tests.
- Examples covering success, ambiguity, contradiction, unsupported cases, and
  missing data.

### Acceptance gate

Operations run without an internal LLM controller and match the authoritative
package behavior. Domain changes have positive, negative, ambiguous, conflict,
and invariance regressions where applicable. Mandatory chemistry filters remain
inside the canonical workflow.

## 6. Phase 2: Direct external-agent access

### Work

1. Implement one thin MCP server exposing the Phase 1 operations.
2. Supply clear tool descriptions, structured failures, pagination, explicit
   limits, concise summaries, and access to complete records and evidence.
3. Provide a short chemistry usage guide documenting assumptions and suitable
   investigation patterns without prescribing every reasoning step.
4. Connect the first external agent and run a complete development investigation.
5. Preserve high-level recommendation and planning operations as conveniences
   built on the same implementations as the composable tools.
6. Document protocol/transport requirements and optional client configuration.
   Keep Codex-specific settings separate from the server and domain packages.

Begin with local deployment. Keep the adapter structured so a remote HTTP
transport can be added without changing scientific contracts. Do not build a new
agent runtime or split every capability into a separate deployed service.

### Deliverables

- A working MCP toolkit and setup instructions.
- Tool contract tests and recorded integration smoke results.
- One complete external-agent investigation with inspectable evidence.

### Acceptance gate

An external agent can analyze a reaction, search and inspect precedents, obtain
conditions, and explain limitations without invoking the built-in assistance
loop. Tool outputs preserve versions, provenance, admission status, and
abstention. This is the first development release.

## 7. Phase 3: Persistent investigation state

### Work

Add application-owned, revisioned records for investigations, hypotheses,
candidate recipes, route alternatives, assessments, decisions, and unresolved
questions. Reference existing domain objects instead of copying them into a
universal reaction dictionary.

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

The first slice should contain baseline cleanup, reaction analysis, precedent
retrieval and full-record inspection, and condition recommendation through one
independent interface. Demonstrate one evidence-backed investigation before
expanding the tool catalog or building a larger application.

For implementation changes, run focused tests and the complete `pytest -q` suite
before handoff. Validate loaders and schemas when definitions change. Document
chemistry and dataset impacts explicitly; never update snapshots merely because
outputs changed. Keep new experiments on development data until the untouched
evaluation gate.

This plan defines proposed work. Passing an interface or deployment test does
not establish scientific accuracy or complete a chemistry release gate.
