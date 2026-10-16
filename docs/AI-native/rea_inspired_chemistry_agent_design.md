# REA-Inspired Architecture for Agentic Chemistry and Synthesis Planning

**Date:** 2026-10-16  
**Status:** Proposed system design  
**Scope:** Synthesis planning, reaction precedent investigation, condition recommendation, route assessment, and evidence management.

> The agent directs the investigation; a chemistry service layer makes tools consistent and preserves evidence; specialized engines perform the underlying calculations, searches, and predictions.

This design borrows REA’s architecture and evidence-handling discipline rather than its reverse-engineering code wholesale. REA’s actual implementation is more elaborate than a simple three-layer diagram: it separates interfaces, shared application workflows, provider adapters, domain logic, and investigation records. That separation is the useful architectural reference. [REA architecture][rea-architecture]

The chemistry components and tool names below are proposals, not functionality already provided by REA. The design is intended to fit a synthesis-planning system combining a single-step retrosynthesis model, local reaction data, precedent retrieval, and ConditionCore.

## 1. Translating the Three Layers into Chemistry

| Layer | Main responsibilities | Application in the chemistry system |
|---|---|---|
| **1. Chemistry agent** | Interprets the objective, develops strategic hypotheses, chooses investigations, compares alternatives, and explains conclusions. | Decides whether to retrieve a known scaffold synthesis, examine a particular disconnection, change the order of steps, or investigate a weak step. |
| **2. Chemistry service layer** | Provides stable tool interfaces, validates inputs and outputs, coordinates engines, tracks evidence, manages sessions, and enforces constraints. | A chemistry equivalent of REA: molecule and reaction records, precursor generation tools, precedent tools, condition tools, and route assessment. |
| **3. Chemistry engines and data sources** | Perform structure processing, prediction, graph search, retrieval, rule evaluation, and inventory lookup. | RDKit, single-step models, a route-search framework, local reaction databases, ConditionCore, and reagent inventory. |

The middle layer should contain real domain logic—not merely forward requests. It should know, for example, what constitutes the same molecule, how a reaction prediction differs from an experimental record, and which evidence supports a particular route step.

Expose the same application functions through both **MCP for agents** and **Python/CLI interfaces for testing and conventional workflows**. REA similarly shares application workflows and evidence contracts between its MCP server and CLI. [REA CLI][rea-cli]

### 1.1. Separate strategic planning from route-tree search

Give the agent responsibility for questions such as:

> Should we construct this fused core first, or assemble it through a late annulation?

Give the search engine responsibility for:

> Expand precursor candidates, track stock termination, remove duplicate states, and search within this computational budget.

Syntheseus is a relevant architectural reference: it separates reaction models from search algorithms and provides interfaces for substituting either component. It is a framework, not a single retrosynthesis model. [Syntheseus][syntheseus]

**Recommended division:** agent-level strategic decisions, engine-level combinatorial search. Calling an LLM for every individual search-tree expansion should not be the starting design.

## 2. Six REA Design Ideas to Borrow

### 2.1. Provider-independent interfaces without hiding provider differences

REA uses provider adapters and caller-visible contracts to separate its workflows from particular analysis engines. [REA architecture][rea-architecture]

The chemistry equivalent would be a stable operation such as:

```text
propose_precursors(product_id, constraints, providers, limit)
```

Behind that interface, the system could use an internal model, another single-step model, or a rule/template engine.

All providers should return a common representation of the proposed transformation, but **normalization must not erase important differences**.

| Required information | Why it matters |
|---|---|
| Exact product and precursor structures | Prevents differences in stereochemistry, charge, or salt treatment from disappearing. |
| Provider and model revision | Identifies which engine produced the proposal. |
| Original score and its definition | Prevents a ranking score from being interpreted as an experimental success probability. |
| Reaction-center information, when available | Supports comparison by disconnection rather than only by precursor strings. |
| Supported and unsupported features | Makes clear whether conditions, stereochemistry, or selectivity were actually assessed. |
| Original provider output | Allows debugging when normalization produces a questionable result. |

Require **explicit fallback reporting**. When one engine fails and another is used, that change must appear in the result.

A uniform interface should make providers interchangeable operationally—not pretend their predictions have identical meaning.

### 2.2. Make every important result an evidence-bearing record

REA returns normalized results with evidence identities, preserves provider representations when available, and retains records for subsequent inspection. [REA MCP contracts][rea-contracts]

For chemistry, introduce a **reaction-step dossier**.

Instead of returning:

> This coupling is feasible; use condition set A.

return a structured record containing:

| Dossier component | Proposed contents |
|---|---|
| **Transformation** | Exact reactants, product, intended bond changes, and reaction class. |
| **Origin of the proposal** | Model prediction, retrieved experiment, rule-generated candidate, or agent hypothesis. |
| **Experimental support** | Exact or analogous reaction records, with document/example locations. |
| **Transferability assessment** | Relevant similarities and differences between the proposed substrates and the precedents. |
| **Conditions** | Candidate condition sets, their origins, and known missing details. |
| **Checks** | Structure consistency, stereochemical consistency, forward-model results, and applicability checks. |
| **Open questions** | Unresolved selectivity, compatibility, operational, or evidence gaps. |

**Evidence should attach to individual claims, not merely to the route as a whole.**

For example, a paper might support the ring-forming transformation but not the proposed substrate substitution pattern. A different source might support a condition choice. Those are separate claims and should retain separate links.

For experimentally reported reactions, map suitable fields to the ORD schema rather than inventing all experimental vocabulary. ORD already represents inputs, conditions, observations, workups, outcomes, and provenance. Keep predicted steps and planning claims in separate records linked to those experimental records. [ORD schema][ord-schema]

### 2.3. Treat unknown, not searched, and not found as different states

REA explicitly distinguishes unavailable observations from observed empty results. An unsupported or incomplete analysis is not treated as proof of absence. [REA MCP contracts][rea-contracts]

The same distinction is essential in chemistry:

| Result | Correct interpretation |
|---|---|
| No exact precedent found in the local database | No matching record was found **within that database and search scope**. |
| Literature search was not performed | External precedent support remains unknown. |
| A model does not assess stereochemistry | Stereochemical feasibility was not assessed. |
| No compatible conditions returned | The current method found no candidates; this does not establish impossibility. |
| A relevant experiment failed | Negative evidence for the recorded substrates and conditions—not necessarily for every implementation of the transformation. |

Store important uncertainties as objects that the agent can revisit:

```text
Question: Is the required site-selectivity supported?
Affected step: step_04
Current status: unresolved
Existing evidence: analogous reaction, but a different competing site
Next useful investigation: retrieve a closer substrate example
```

This creates a better planning loop than repeatedly asking the model whether it is confident.

It also preserves room for genuinely new chemistry: **an unsupported hypothesis can remain a candidate, provided it is not mislabeled as established.**

### 2.4. Return summaries first, then allow detailed inspection

REA can retain complete analysis results while exposing summaries and selected views for follow-up inspection. [REA MCP contracts][rea-contracts]

Avoid sending the agent dozens of complete procedures and hundreds of precursor candidates at once.

A precedent search should initially return a compact comparison:

```text
Precedent ID
Reaction-center match
Relevant substrate differences
Reported outcome
Condition completeness
Source-verification status
```

The agent can request the full record for the few examples that matter.

Likewise, route search should initially expose **distinct strategies and their bottlenecks**, not every node in the search graph.

The design objective is to support:

> Show me the strongest evidence for the weakest step in route B.

That is a more useful interface than repeatedly returning an entire route database.

### 2.5. Preserve reproducible sessions and versioned results

REA’s snapshot system keys reusable results to the target, operation, parameters, provider, and settings. [REA CLI][rea-cli]

Apply the same principle to chemistry, with additional domain-specific identity rules.

A planning session should preserve:

- The original target, normalization policy, and stereochemical requirements.
- Model revisions, dataset snapshots, search settings, and constraints.
- Condition-rule versions and material-availability checks.

**Do not use only a target SMILES as the cache key.** A result obtained with one stock database, one stereochemical policy, and one model revision is not equivalent to a result obtained with different settings.

Preserve the original structure alongside standardized representations. Standardization should be an explicit, inspectable transformation—not an invisible preprocessing step that changes the planning objective.

For inventory, retain the observation time and quantity basis. For cached literature results, retain search coverage and source versions.

### 2.6. Separate workflow instructions from executable chemistry logic

REA’s guided prompts provide investigation strategies; they do not themselves invoke tools. They also instruct agents to distinguish observations, inferences, and unknowns. [REA MCP prompts][rea-prompts]

The chemistry equivalent could be a small set of workflows:

| Workflow | Purpose |
|---|---|
| **Precedent-first planning** | Retrieve exact or related scaffold syntheses before broad disconnection search. |
| **Bottleneck investigation** | Identify the weakest transformation and gather the evidence needed to assess it. |
| **Route comparison** | Compare strategic differences, experimental support, material access, and unresolved risks. |
| **Condition assessment** | Evaluate substrate features against precedents and ConditionCore rules. |

Critical requirements should be enforced in code.

A prompt can say, “Preserve target stereochemistry.” The service should independently check that the proposed final product matches the requested structure.

A prompt can say, “Cite experimental evidence.” The service should verify that the cited evidence identifier exists and refers to the claimed record.

**Instructions guide behavior; contracts and validators enforce invariants.**

## 3. Chemistry-Specific Verification: Not a Single Yes/No Test

This is the most important departure from a simplistic reverse-engineering analogy.

RDKit’s sanitization checks molecular representation and computes properties needed for molecular processing. Passing those checks is not a demonstration that a proposed synthesis will work. [RDKit Book][rdkit-book]

Similarly, ASKCOS distinguishes retrosynthesis, condition recommendation, and reaction-outcome prediction as separate functions; its description cautions that displayed forward-prediction probabilities should be interpreted qualitatively. [ASKCOS description][askcos-paper]

Do **not** expose a single field such as:

```text
reaction_valid: true
```

Separate the assessments instead:

| Assessment | What the proposed system should report |
|---|---|
| **Structural validity** | Whether the submitted molecular representations pass the defined checks. |
| **Transformation consistency** | Whether the proposed reactants and product encode the intended transformation. |
| **Stereochemical consistency** | Whether specified stereochemistry is preserved or its formation is explicitly addressed. |
| **Model support** | Which models support the proposal, under which inputs and settings. |
| **Experimental support** | Whether exact or analogous experiments have been located and source-checked. |
| **Condition applicability** | Which substrate features are supported, extrapolated, or unresolved for each condition set. |
| **Operational readiness** | Whether materials, procedure details, equipment requirements, and review are sufficient for the intended experiment. |

These are separate dimensions, not a ladder where passing one automatically proves the next.

For transformation consistency, define the scope carefully: checking mapped product-atom accounting is different from requiring a fully balanced reaction record containing every reagent and byproduct.

### 3.1. Keep conditions inside the planning loop

Integrate ConditionCore **during assessment of strategic steps**, not only after a route has been selected.

```text
Proposed transformation
        ↓
Candidate conditions + applicability assessment
        ↓
New compatibility or selectivity concerns
        ↓
Revise conditions, step order, precursors, or disconnection
```

Otherwise, the planner can select an attractive abstract route whose critical steps have no adequately supported implementation.

Treat condition-related warnings as reasons to investigate or rerank—not automatically reject every unfamiliar substrate.

## 4. Proposed Initial Tool Interface

Start with a relatively small agent-facing surface. The following names are **proposed interfaces**, not existing REA commands.

| Tool | Purpose |
|---|---|
| `inspect_molecule` | Validate the target, preserve identity, identify relevant structural features, and expose possible reactive sites. |
| `search_precedents` | Search exact structures, related cores, transformations, and reaction-center environments with explicit coverage. |
| `inspect_precedent` | Retrieve a particular reaction record, procedure, source location, and extraction uncertainties. |
| `propose_precursors` | Obtain normalized single-step candidates from selected providers. |
| `search_routes` | Run bounded multi-step search under defined constraints and return route candidates plus termination reasons. |
| `recommend_conditions` | Return condition sets with precedent/rule origins and applicability information. |
| `assess_step` | Run selected structural, model-based, evidence, and condition checks; return a dossier. |
| `check_materials` | Distinguish local stock, catalog listing, and checked availability for the required material form and quantity. |

Session management, evidence inspection, route comparison, and report generation can sit alongside these tools.

### 4.1. Separate hypothesis generation from assessment

A model or agent should be allowed to propose an unusual disconnection. That proposal then enters the same assessment pipeline as every other candidate.

### 4.2. Choose useful tool granularity

Avoid both extremes. The agent should not need to call every low-level RDKit operation separately, but it should also not be trapped behind one opaque `plan_synthesis()` call.

A composite `assess_step` tool can coordinate several checks while exposing each result and its origin.

### 4.3. Role of the A–B/reactivity taxonomy

The A–B/reactivity taxonomy can organize reaction classes, reactive-site candidates, and rule applicability. A taxonomy label should not substitute for the exact molecular graph or atom-level reaction-center description.

## 5. End-to-End Synthesis-Planning Workflow

For a complex heterocyclic target, use the following workflow.

### 5.1. Establish the target and constraints

The service registers the exact target structure, stereochemical requirements, allowed starting materials, and relevant operational constraints. It returns a structural summary for the agent to inspect.

### 5.2. Retrieve precedents and develop strategic hypotheses

The agent first investigates exact or related scaffold syntheses, following a precedent-first approach.

Retain an **independent disconnection branch** so that precedent retrieval does not dictate every route. That branch could use the internal single-step model, another provider, or an agent-proposed core-construction strategy.

### 5.3. Expand promising strategies with engines

The route-search engine explores the selected branches within explicit budgets.

It should report why exploration stopped: stock termination, depth limit, budget exhaustion, unsupported structures, or another reason.

Label stock-terminated routes as such—not simply as “validated routes.”

### 5.4. Investigate the bottleneck, not every step equally

Suppose two routes contain routine finishing steps but differ in their core-forming transformation.

The agent should concentrate evidence retrieval and condition assessment on those core-forming steps. It can request closer precedents, inspect original procedures, or propose a different construction when the evidence is weak.

### 5.5. Feed findings back into the route

A condition assessment might reveal that the proposed sequence leaves an unresolved compatibility issue. The agent can change the step order, choose different precursors, or ask the search engine to explore another branch.

This is the central benefit of the architecture: **tool results can change the strategy**, rather than merely decorate a route already chosen by the LLM.

### 5.6. Produce one consistent output package

Generate the route diagram, step table, source references, condition proposals, and unresolved-question list from the same saved route graph and evidence records.

For SVG/Markdown workflow comparisons, the drawn structures and written route should not be independently regenerated and accidentally disagree.

Include an alternative route only when it represents a meaningful strategic alternative with enough supporting analysis—not merely different conditions on the same weak transformation.

## 6. Implementation and Evaluation

### 6.1. Build the evidence layer before adding more models

The first implementation should wrap the existing components:

```text
RDKit-based structure handling
        +
Single-step retrosynthesis model
        +
Local precedent retrieval
        +
ConditionCore
        +
Material lookup
```

The first deliverable should be a reliable reaction-step dossier, stable molecular identities, and reproducible planning sessions. A straightforward database and structured route records are enough initially; a large multi-agent system or specialized graph database is not a prerequisite.

Next, add the agent’s iterative workflow and a bounded route-search engine. Additional models should come later, when the system can measure what unique, useful candidates they contribute.

### 6.2. Test adapters separately from chemistry quality

REA separates integration/conformance tests from real-provider verification and agent-investigation evaluation. Its documentation warns that tests using injected providers do not establish real-provider behavior. [REA testing][rea-testing]

Make the same distinction:

| Evaluation level | Question |
|---|---|
| **Software correctness** | Does the adapter preserve structures, scores, conditions, errors, and provenance correctly? |
| **Planning quality** | Are the resulting routes and claims chemically useful and appropriately supported? |
| **Experimental performance** | Do selected proposals work under the recorded experimental conditions? |

For workflow comparison, hold the target set, backend models, data access, and computational budgets as constant as practical. Compare the current workflow with the same tools accessed through the proposed service layer and guided agent.

Important outcomes include:

- Fewer structure and stereochemistry errors, and fewer unsupported “validated” claims.
- Better-supported critical steps and clearer unresolved bottlenecks.
- More meaningful strategic diversity and less chemist review effort.

Do not use route count or stock termination alone as the success metric.

### 6.3. Keep planning and execution behind different permissions

Keep the first version of synthesis planning read-only. An approved route can later become an experimental proposal for ChemRoboX, but instrument execution should require separate procedure checks and authorization.

The three logical layers need not be three machines. Private reaction data, inventory, and the evidence store can remain local or on an internal server; heavier prediction services can run on a separate compute server.

Any information sent to a cloud agent should be an explicit data-sharing decision.

## 7. Core Recommendation

**Borrow REA’s service boundaries, structured evidence, explicit uncertainty, progressive inspection, and reproducible workflows. Do not merely copy its MCP wrapper pattern.**

For the chemistry system, the most valuable reusable object is the **reaction-step dossier**, and the most valuable agent behavior is **investigating and resolving the weakest parts of a route**.

The resulting system should do more than generate plausible routes. It should make clear:

> What was calculated, what was predicted, what was actually reported, what remains uncertain, and what investigation would most improve the plan.

## References

These are the source links cited in the preceding design discussion. The chemistry architecture, proposed interfaces, and implementation recommendations are design proposals rather than descriptions of existing REA chemistry functionality.

1. [REA — architecture][rea-architecture]
2. [REA — CLI and reusable analysis workflows][rea-cli]
3. [Syntheseus — synthesis-planning framework][syntheseus]
4. [REA — MCP contracts][rea-contracts]
5. [Open Reaction Database — schema documentation][ord-schema]
6. [REA — MCP prompts and investigation workflows][rea-prompts]
7. [RDKit Book — molecular processing and sanitization][rdkit-book]
8. [ASKCOS description — arXiv:2501.01835v1][askcos-paper]
9. [REA — testing documentation][rea-testing]

[rea-architecture]: https://raw.githubusercontent.com/morluto/rea/main/docs/architecture.mermaid
[rea-cli]: https://raw.githubusercontent.com/morluto/rea/main/docs/cli.md
[syntheseus]: https://github.com/microsoft/syntheseus
[rea-contracts]: https://raw.githubusercontent.com/morluto/rea/main/docs/mcp-contracts.md
[ord-schema]: https://docs.open-reaction-database.org/en/latest/schema.html
[rea-prompts]: https://raw.githubusercontent.com/morluto/rea/main/docs/mcp-prompts.md
[rdkit-book]: https://www.rdkit.org/docs/RDKit_Book.html
[askcos-paper]: https://arxiv.org/html/2501.01835v1
[rea-testing]: https://raw.githubusercontent.com/morluto/rea/main/docs/testing.md
