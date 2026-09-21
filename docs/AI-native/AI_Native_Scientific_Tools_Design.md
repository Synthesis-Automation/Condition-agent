# AI-Native Scientific Tools and Agent-Assisted Chemistry

Status: design proposal, not an implemented capability or validation claim  
Date: 2026-09-22

## 1. Purpose

Make reaction analysis, datasets, condition recommendation, and retrosynthesis
available as modular scientific capabilities that advanced agents can compose
into an investigation.

The intended improvement is the ability to choose useful questions, gather
evidence, propose hypotheses, compare alternatives, and revise a plan as new
information arrives. Adding an LLM explanation to a fixed pipeline does not
provide all of these capabilities.

This document develops the intent in [my_idea.md](my_idea.md) and refines the
information-system portion of
[AI_Chemistry_MCP_System_Design.md](AI_Chemistry_MCP_System_Design.md).
It does not replace the package boundaries or release gates in the
[current implementation roadmap](../new/type_agnostic_reaction_recommendation_implementation.md).
Hardware control and autonomous laboratory execution are outside this proposal.

## 2. Model, agent runtime, and scientific environment

These are separate architectural responsibilities:

| Part | Responsibility |
| --- | --- |
| Model | Reason about the available context and propose the next action or response. |
| Agent runtime | Manage the investigation loop, context, tool execution, progress, failures, stopping, and continuation. |
| Scientific environment | Supply domain operations, datasets, evidence, persistent artifacts, and explicit checks. |

The same model can behave differently when provided with a different runtime,
tools, context, and execution budget. An advanced agent can decide what evidence
to acquire and revisit earlier decisions; a fixed sequence of LLM calls exposes
only the reasoning opportunities anticipated by its author.

Codex is a candidate runtime, not a required dependency of the chemistry
packages. OpenAI describes its harness as managing context, tools, conversation
state, and ongoing execution, and documents integration through its SDK and
app-server. See [Codex as a platform](https://developers.openai.com/blog/codex-as-a-platform).
This supports the integration direction; it does not establish chemistry
performance. Scientific improvement must be evaluated on our own tasks.

The project should invest primarily in the scientific environment. External
agents should be able to use it directly without entering a second internal LLM
controller. A built-in chemistry assistant can remain another client of the
same capabilities.

## 3. Architecture and ownership

```mermaid
flowchart TD
    User[Chemist and scientific task] --> Agent[Agent runtime]
    Agent <--> Workspace[Investigation workspace]
    Agent <--> Interface[Scientific interfaces: MCP, Python API, CLI]
    Interface --> Taxonomy[reactive_taxonomy]
    Interface --> Registry[condition_registry]
    Interface --> Recommender[condition_recommender]
    Interface --> Retro[core_retrosynthesis]
    Interface --> Sources[External evidence and model adapters]
    Recommender --> Taxonomy
    Recommender --> Registry
```

- `reactive_taxonomy` owns molecular and reaction observations, interpretation,
  correspondence, edits, and structural checks.
- `condition_registry` owns substance identities, contextual roles, canonical
  recipes, and its existing protocol-draft contracts.
- `condition_recommender` owns dataset conversion and admission, indexing,
  compatibility, retrieval, ranking, aggregation, and recommendation evidence.
- `core_retrosynthesis` owns its existing domain planning and verification
  operations; the interface should reuse these rather than copy them.
- Application adapters expose capabilities and manage investigation artifacts.
  They do not introduce chemistry rules or an alternative recommendation path.
- External literature and predictive-model adapters remain explicit boundaries;
  deterministic domain libraries remain free of network calls.

MCP is an access mechanism. It is not the reasoning engine, the scientific
contract, or a reason to deploy every capability as a separate service. Begin
with a thin interface over existing APIs and separate deployment only when
operational requirements justify it.

## 4. Meaningful atomic operations

An atomic capability answers one scientific question and returns a result that
can be inspected and reused. It need not correspond to one Python function or
one low-level graph manipulation.

Expose both composable operations and convenient complete workflows. They must
share domain implementations. Do not require an agent to reconstruct routine
normalization or bypass mandatory compatibility checks by manually assembling
low-level calls.

The names below are proposed interface concepts, not claims about current APIs.

| Area | Proposed operations | Result emphasis |
| --- | --- | --- |
| Reaction analysis | `observe_reaction`, `inspect_reaction_center`, `enumerate_interpretations`, `check_transformation` | Structures, atom references, edits, competing interpretations, checks and limitations. |
| Dataset access | `search_precedents`, `get_precedent`, `get_procedure`, `compare_precedents`, `inspect_dataset_coverage` | Complete records, structural query criteria, source provenance, missingness and independent support. |
| Conditions | `resolve_recipe`, `recommend_conditions`, `assess_recipe`, `compare_recipes`, `inspect_component_roles` | Canonical recipes, contextual roles, compatibility, observed variants, and unsupported assumptions. |
| Retrosynthesis | `propose_disconnections`, `generate_precursors`, `evaluate_step`, `find_step_precedents`, `assemble_route`, `revise_route` | Alternative intermediates and steps, structural checks, evidence, dependencies, and unresolved concerns. |
| Evidence access | `get_evidence`, `get_artifact`, `render_structure` | Inspectable source details, stable references, and visual structures tied to machine objects. |

Dataset access should support exploratory searches that do not require a named
reaction family. Exploratory results must retain admission status and must not
silently become trusted recommendation precedents. Known incompatible recipes
may be inspected as counterevidence but remain excluded from eligible rankings.

## 5. An agent-accessible working environment

Useful tools require useful context and feedback. Provide:

- Searchable documentation with preconditions, examples, limitations, and
  distinctions between observations and predictions.
- Typed input and output schemas, stable object references, and explicit units.
- Compact summaries with pagination and access to full records and artifacts.
- A Python API for custom comparisons and batch analysis using the same domain
  operations as the MCP and CLI interfaces.
- Persistent files or records for candidate routes, recipes, evidence,
  hypotheses, and decisions.
- Errors that identify what failed, what partial result remains usable, and
  which inputs or capabilities are missing.
- Execution limits, cancellation, and resumable job references for expensive
  searches or model calls.

An agent may write a derived analysis script when no predefined operation
answers a question. Preserve its code, input references, dependencies, and
outputs. Such an analysis is a derived artifact, not an automatic change to
production chemistry definitions or trusted datasets.

The agent's conversation is not the sole scientific record. Domain objects and
evidence must survive context compaction, runtime changes, and resumed sessions.

## 6. Evidence, hypotheses, and checks

Deterministic code establishes the result of its implemented checks under its
documented assumptions. Passing a graph check does not establish experimental
feasibility, and missing rule coverage does not establish chemical impossibility.

Every assessment should distinguish at least:

| Check outcome | Meaning |
| --- | --- |
| `pass` | A specified check passed within its stated scope. |
| `fail` | A specified check found a violation or contradiction. |
| `unknown` | Evidence is insufficient or ambiguous. |
| `unsupported` | The capability does not implement the required assessment. |
| `error` | The operation failed to execute successfully. |

Preserve existing domain statuses and map them explicitly at the interface;
these categories are not a replacement chemistry schema.

Scientific claims should reference supporting and conflicting evidence and
state assumptions. Record provenance separately from review status: an observed
source report, deterministic inference, model prediction, and agent hypothesis
are different origins, even when a chemist has reviewed them.

Allow agents to propose unfamiliar transformations and condition adaptations in
the investigation workspace. These proposals do not create verified signatures,
override structural contradictions, enter trusted indices, or become observed
experimental records. Promotion requires the existing admission and review
rules. Unknown cases remain available for investigation without being presented
as validated recommendations.

A useful result can say:

> Atom accounting passed. Selectivity is outside the current validator's scope.
> Two source precedents provide partial support, with different substrate
> environments. The proposed step remains a hypothesis.

## 7. Persistent investigation contracts

Reuse existing chemistry and recommendation contracts through references rather
than duplicating their contents in a universal reaction dictionary.

Proposed application-owned objects:

- **Investigation:** objective, user constraints, artifact references, open
  questions, budget, and completion or stopping status.
- **Hypothesis:** a proposed transformation, route choice, component role, or
  condition adaptation, with assumptions, evidence links, and review status.
- **Evidence reference:** source record or artifact ID, relevant fields or
  excerpt location, provenance, versions, and applicable limitations.
- **Assessment:** subject reference, named checks and outcomes, contradictions,
  unassessed properties, and links to authoritative domain results.
- **Decision:** chosen action or alternative, concise rationale, evidence used,
  rejected alternatives, and unresolved questions.

Record tool arguments, result references, schema and definition versions, data
snapshot identity, and model/provider metadata where applicable. Preserve an
inspectable decision record; private model reasoning is not required.

Candidate objects should be immutable or explicitly revisioned. Agent state may
reference a selected version without overwriting the source observation.
Deterministic computations should be replayable from recorded inputs. Replaying
an agent investigation may require recorded model outputs; a fresh model run is
not guaranteed to follow the same path.

## 8. Retrosynthesis as an evolving investigation

The agent should work with a route graph containing competing intermediates,
proposed steps, evidence, constraints, and unresolved questions. Keep its control
of investigation separate from the deterministic operations that evaluate or
expand that graph.

An illustrative workflow is:

1. Analyze a target and propose several strategic decompositions.
2. Ask available engines for precursors under selected structural constraints.
3. Retain an additional agent-proposed candidate as a hypothesis when useful.
4. Check structures, correspondence, and reconstruction where supported.
5. Retrieve step precedents and inspect conditions and procedure details.
6. Investigate downstream selectivity, intermediate stability, or availability.
7. Revise an earlier disconnection or step order when evidence warrants it.
8. Present alternative routes with evidence and explicit unresolved limitations.

Structural consistency, experimental support, selectivity, operating conditions,
and route practicality are separate assessment dimensions. Avoid a single
unqualified `valid_route` flag or a score that implies experimentally established
success.

The expected advantage over fixed search is choosing what to investigate and
when to reconsider strategy. This is a hypothesis to test, not a performance
claim supplied by the architecture.

## 9. Condition reasoning beyond recipe retrieval

The canonical recommender remains responsible for structure-backed retrieval,
hard compatibility, recipe aggregation, and ranking. The agent adds an
investigation around its evidence and limitations.

Useful condition reasoning includes:

- Examining a component's possible roles in this specific reaction context.
- Comparing whole recipes and operating procedures, preserving co-occurrence
  rather than mixing frequent components into an unsupported recipe.
- Identifying substrate differences that may limit transfer from precedents.
- Searching for contradictory examples and reported failure outcomes.
- Proposing an adaptation with explicit assumptions and unresolved questions.
- Identifying observations or experiments that would distinguish hypotheses.

Keep an observed recipe, a normalized representation of that recipe, and a
proposed adaptation distinct. Never invent missing temperature, quantities,
order of addition, concentration, or time. A role may be multiple or uncertain.
Precedent association does not establish mechanism or causality, and an
unreported outcome is not an experimental failure.

An agent's explanation should identify the evidence behind its claims.
Recommendation scores, evidence strength, and calibrated probabilities are
different quantities and must not share an ambiguous `confidence` field.

## 10. Integration with the current codebase

The existing `chem_coworker/assistance/` implementation already contains
capability wrappers, evidence projections, action contracts, policy, and
evaluation support. Reuse useful contracts and domain calls.

Its closed action vocabulary is a starting point for bounded assistance, but
should not define the maximum reasoning space of every external agent. Expose
the underlying scientific operations independently. Add proposed hypotheses
and targeted assessments where they fill a demonstrated scientific gap.

Keep `chem_coworker` a thin application shell. Do not restore `chemtools`, revive
ConditionCore/family routing as a parallel recommender, or move chemistry rules
into prompts. Agent instructions can explain how to use tools and interpret
results; executable scientific behavior stays with its owning package.

## 11. Staged implementation

### Stage A: Inventory and interface specification

Map existing public operations to the capability table. Document current
contracts, missing capabilities, ownership, and required evidence. Select two
bounded investigation tasks before adding new abstractions.

### Stage B: Direct access to existing capabilities

Implement a thin interface that an external agent can use for reaction analysis,
precedent inspection, condition recommendation, and existing retro operations.
Provide examples, structured failures, full-result access, and version metadata.
Verify parity with direct package calls.

### Stage C: Investigation artifacts and hypothesis assessment

Add persistent route alternatives, hypotheses, evidence references, and
assessments. Enable an agent to submit a proposed step or recipe for checks
without modifying trusted observations or production definitions.

### Stage D: Evaluate agent-assisted investigations

Use two initial demonstrations:

1. Revise a multistep route after identifying a problematic downstream step.
2. Investigate a condition suggestion whose nearest precedents contain an
   important substrate mismatch.

Compare the deterministic workflow, the workflow plus LLM review, and an
advanced agent using the same scientific capabilities. Use the same underlying
model where possible. Report tool/data access and resource budgets, and repeat
agent runs to measure variation. Use matched-budget comparisons to separate
runtime effects from additional computation.

### Stage E: Consolidate validated interfaces

Retain the interfaces that improve task outcomes and remove temporary duplicate
adapters after parity. Keep runtime-specific integration replaceable. Product
embedding can follow once direct agent use has demonstrated value.

These stages do not waive the roadmap's frozen baseline, blind chemist review,
adjudication, untouched evaluation, corpus validation, or release requirements.
Experimental interfaces must identify the validation status of the artifacts
they expose. Do not tune development against the untouched evaluation set.

## 12. Evaluation and acceptance criteria

Blind chemistry review should assess:

- Step and route plausibility, including dependencies between steps.
- Correctly identified condition-transfer limitations and incompatibilities.
- Whether proposed alternatives address the actual problem.
- Whether cited evidence supports the associated claim, rather than merely
  existing as a valid reference.
- Unsupported assertions, missed contradictions, and appropriate abstention.
- Quality of uncertainty communication and remaining questions.

Operational evaluation should assess successful task completion, tool failures,
recoverability, repeatability of recorded calculations, latency, tokens, and
cost. Test unfamiliar, ambiguous, conflicting, and unsupported cases as well as
straightforward successes.

An integration succeeds when an agent can investigate a problem more effectively
using the exposed capabilities while preserving provenance and domain
boundaries. More tool calls, longer explanations, and successful schema
validation alone are not evidence of better chemistry.
