# ChemCoworker Chemistry System Implementation Plan

## Scope

This note captures a chemistry-core improvement analysis for the ChemCoworker system (classification, workflow routing, tool execution contracts, chemistry validation, confidence, and knowledge intake integrity), not just CLI/UI changes.

## Summary

The current system is already strong in architecture (tool registry, workflow registry, event bus, critic pass, validators), but several high-impact chemistry correctness gaps remain:

- Tool dependency contracts (`requires`, chemistry prerequisites) are mostly advisory.
- Reaction typing is not propagated as a canonical taxonomy object across the pipeline.
- Confidence is not evidence-calibrated.
- Validator coverage is sparse (mainly `recommend_conditions`).
- Workflow routing likely has unreachable chemistry paths (`forward_synthesis`).
- Intake/note filing is not consistently taxonomy-validated before writing chemistry knowledge.

These issues are systematic and should be addressed at the contract/runtime layer, not only in prompts or UI.

## Key Findings

### 1. Tool chemistry dependencies are not enforced at runtime

- `ToolPlugin.requires` is informational (`chem_coworker/tools/_base.py`)
- Native loop executes whatever tool calls the LLM emits (`chem_coworker/agent.py`)
- This can allow chemistry-invalid execution order (e.g., conditions recommendation before reliable reaction typing / evidence)

Impact:

- Reduced robustness
- More hallucination risk under partial evidence
- Harder to guarantee chemistry-first behavior

### 2. Reaction type handling is not fully taxonomy-canonical end-to-end

- `detect_reaction_type` returns a `reaction_type` string and a display label derived from string formatting
- It does not consistently attach canonical taxonomy metadata (name/category/aliases/constraints)

Impact:

- Downstream tools and outputs may drift from taxonomy terminology
- Harder to reason about ambiguity and family-level constraints systematically

### 3. Single-label reaction typing loses ambiguity information

The pipeline mostly carries a single `reaction_type` forward. In practice, many reactions are ambiguous or partially specified.

Impact:

- Premature commitment to one family
- Weak fallback behavior for conditions/troubleshooting/forward planning
- Confidence cannot be calibrated well

### 4. Workflow routing mismatch / unreachable chemistry workflow path

- `WorkflowRegistry` defines `forward_synthesis`
- `TaskClassifier` does not appear to emit a corresponding task type

Impact:

- Intended chemistry workflow specialization may not be used
- Features in specialized workflows can silently remain inactive

### 5. Task classification is keyword-heavy instead of chemistry-signal-first

Current classifier is practical, but relies heavily on keywords/regex.

Impact:

- Misrouting for chemistry prompts that are structurally obvious but lexically ambiguous
- Reduced generality for new chemistry tasks phrased differently

### 6. Validator coverage is too narrow

Validator pattern is a strong design, but only a small subset of tools uses it (notably `recommend_conditions`).

Impact:

- Chemistry caveats are inconsistently surfaced
- Critical correctness checks remain implicit in prompts or post-hoc critique

### 7. `ChemResponse.confidence` is not evidence-calibrated

Confidence currently behaves more like a placeholder/default than a chemistry evidence score.

Impact:

- UI/API consumers cannot distinguish strong vs weak chemistry conclusions
- Harder to automate downstream decisions safely

### 8. Structured output extraction is hardcoded by tool name

`_extract_structured` in the agent is tool-name-specific and manual.

Impact:

- Adding chemistry tools requires agent edits
- Limits extensibility and consistency

### 9. Intake / notes writing is not fully taxonomy-validated before KB writes

The notes intake pipeline writes chemistry knowledge to files, but reaction type hints/detections are not consistently canonicalized against the taxonomy catalog before filing.

Impact:

- Knowledge base contamination risk
- Duplicate/alias note fragmentation

### 10. No mismatch policy between hinted and detected chemistry labels in intake

If a user hint conflicts with extracted/detected reaction types, current behavior may still write notes.

Impact:

- Incorrect note placement
- Lower trust in accumulated domain knowledge

## Recommended Systematic Improvements (Priority Order)

## Phase 1 (Highest ROI, chemistry correctness)

Status: Implemented (initial version)

Delivered:

- Runtime contract partitioning/defer logic for native tool calls (enforces `prerequisites`/`requires` before execution)
- Contract-violation synthetic tool results fed back into the LLM tool loop (instead of silently executing invalid chemistry order)
- Taxonomy-canonical reaction type metadata returned by `detect_reaction_type`
- Structured response extraction updated to prefer canonical `reaction_type_id` and include taxonomy metadata
- Targeted tests for contract enforcement and canonical metadata propagation

### A. Enforce tool prerequisites and data contracts at runtime

Implement runtime checks in the native tool loop / executor:

- Reject/defer tool calls whose `requires` keys are not satisfied
- Optionally auto-insert prerequisite tools (`prerequisites`) when missing
- Emit explicit chemistry warnings when the LLM tries to skip anchor tools

Why first:

- Converts chemistry constraints from documentation into execution guarantees

### B. Canonicalize reaction type to taxonomy object across the pipeline

Standardize on taxonomy-resolved reaction identity:

- canonical `reaction_type_id`
- display name
- category
- aliases
- constraints (when relevant)
- evidence/confidence

Use `chemtools.taxonomy.reaction_catalog` as the single source of truth.

Why first:

- Aligns all downstream chemistry logic and UI outputs with taxonomy-driven design

## Phase 2 (Reliability + generality)

Status: Implemented (baseline version)

Delivered:

- Added chemistry validators for `detect_reaction_type` (unknown/low-confidence caveats)
- Added chemistry validators for `analyze_bond_changes` (missing formed bonds / low mapping confidence / unclear key bond type)
- Expanded tool metadata on core chemistry tools (`provides`) to support stronger contract + context propagation
- Added deterministic evidence-based confidence aggregation in `ChemCoworker` (reaction typing, bond changes, conditions support, warnings, critic severity)
- Added targeted tests for validators and confidence aggregation heuristics

### C. Propagate top-N reaction-type hypotheses (not only one label)

Carry ambiguity through:

- `reaction_type_candidates: [{id, score, evidence}]`
- choose final answer strategy based on confidence spread
- allow downstream tools to use candidate-aware logic

### D. Expand validator coverage to core chemistry tools

Add validators for:

- reaction typing ambiguity / low confidence
- bond-change mapping confidence / inconsistency
- forward synthesis product plausibility
- retrosynthesis disconnection plausibility and support gaps

### E. Evidence-based confidence aggregation

Compute confidence from chemistry evidence, e.g.:

- reaction detector confidence
- atom mapping / bond-change confidence
- HTE support size and match quality
- tool agreement/disagreement
- critic findings severity

## Phase 3 (Routing + extensibility)

Status: Implemented (baseline version)

Delivered:

- Fixed practical `forward_synthesis` workflow reachability by adding `forward_synthesis` task classification output
- Added chemistry-signal-first classifier heuristics for forward product-prediction queries
  - distinguishes multi-reactant product prediction from retrosynthesis target planning
  - avoids routing full reaction SMILES (`>>`) to the forward-synthesis workflow
- Added classifier routing tests and workflow integration assertions
- Added tool-declared structured output projections (`ToolPlugin.structured_projection`) and made agent structured aggregation projection-first (with small legacy fallback)
- Added tests covering built-in and custom tool structured projections

### F. Fix workflow routing / task typing gaps

- Add explicit task-type detection for `forward_synthesis`
- Validate workflow reachability with tests
- Ensure classifier/workflow registry stay synchronized

### G. Upgrade classifier to chemistry-signal-first routing

Use chemistry structure cues in addition to keywords:

- reaction SMILES vs target-only molecule
- explicit reagent/condition slot requests
- taxonomy-family names and motif labels
- troubleshooting signatures (yield, selectivity, decomposition)

### H. Replace hardcoded `_extract_structured` with schema/projection contracts

Let each tool declare what structured outputs it contributes.

Benefits:

- easier tool addition
- less agent churn
- more consistent machine-readable chemistry API outputs

## Phase 4 (Knowledge integrity / intake)

Status: Implemented (baseline version)

Delivered:

- Taxonomy-validated intake filing for hint and detected reaction labels (canonicalized before write)
- Unknown reaction label handling policy (`general`, `quarantine`, `reject`)
- Intake mismatch policies (`warn`, `confirm`, `reject`, `force`)
- `intake --dry-run` support (compute target note files without writing)
- Machine-readable intake output (`intake --output-format json`)
- Tests for canonicalization, mismatch handling, unknown-label rejection, dry-run, and CLI JSON intake mode

### I. Taxonomy-validate all intake reaction labels before writing notes

- Canonicalize user hint
- Canonicalize detected labels
- Store canonical IDs
- quarantine unknowns or route to `general`

### J. Add intake mismatch policies

Support:

- `warn` (default)
- `confirm`
- `reject`
- `force`

This protects chemistry knowledge quality without blocking advanced users.

### K. Add `intake --dry-run` and machine-readable intake output

Before writing:

- inspect detected/canonical types
- review warnings
- programmatic QA pipelines for note ingestion

## Proposed Implementation Milestones

### Milestone 1: Contract Enforcement + Canonical Taxonomy Reaction Type

Status: Completed (baseline implementation)

Deliverables:

- Runtime `requires` enforcement in executor/native loop
- Canonical taxonomy reaction type metadata returned by detection-related tools
- Tests for dependency enforcement and taxonomy canonicalization

Expected impact:

- Large improvement in chemistry correctness and reproducibility

### Milestone 2: Validators + Confidence

Status: Completed (baseline implementation)

Deliverables:

- Validators for reaction typing + bond changes + forward/retro core outputs
- Confidence aggregation method in `ChemResponse`
- Warning severity structure (optional extension)

Expected impact:

- Better safety and interpretability

### Milestone 3: Routing + Structured Output Contracts

Status: Completed (baseline)

Deliverables:

- `forward_synthesis` routing fix
- chemistry-signal-first classifier upgrades
- tool-declared structured output projections

Expected impact:

- Generality and maintainability

### Milestone 4: Intake Knowledge Integrity

Status: Completed (baseline)

Deliverables:

- taxonomy-validated intake filing
- mismatch policy
- dry-run + JSON output for intake

Expected impact:

- Cleaner, more trustworthy chemistry knowledge base

## Suggested Tests (Chemistry-Focused)

- Tool call requiring `reaction_type` fails/deferred if detector not run
- Reaction type aliases resolve to canonical taxonomy ID
- Ambiguous detection returns multiple candidates and lowers confidence
- Low-confidence bond mapping emits validator warnings
- `forward_synthesis` prompts route to the correct workflow
- Intake with conflicting hint vs detection follows mismatch policy
- Unknown reaction type in intake is quarantined or rejected (per policy)

## Notes

- The existing architecture (ToolPlugin metadata, validators, workflow registry, critic pass) is already a strong foundation for these changes.
- Prioritize moving chemistry constraints into deterministic runtime checks and taxonomy-backed contracts rather than relying on prompt instructions.
