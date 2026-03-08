# ChemCoworker Skills System Proposal

## Purpose

This document turns the OpenClaw skills/runtime research into a concrete design for `chem_coworker`.

The goal is not to copy OpenClaw mechanically. The goal is to add the parts that improve capability isolation, environment-aware loading, and prompt efficiency while preserving this project's core rule:

- chemistry logic stays in `chemtools/taxonomy`, deterministic tools, and validators
- skills orchestrate tool usage and evidence standards
- skills do not become a second source of chemistry truth

## Executive Summary

`chem_coworker` already has a strong runtime core:

- executable tools are typed and self-describing via `ToolPlugin`
- tools are centrally registered in `chem_coworker/tools/__init__.py`
- LLM-visible tool subsets are already curated by workflow in `chem_coworker/workflow.py`
- the native loop already enforces some contract gating in `chem_coworker/agent.py`

What is missing is an externalized capability layer.

OpenClaw's main architectural advantage is not "more tools". It is that it separates:

- tool implementation
- capability instructions
- environment eligibility
- per-workspace policy

ChemCoworker should adopt that same separation.

## What OpenClaw Gets Right

Based on current OpenClaw docs, the useful patterns are:

1. Skills are filesystem packages
   A skill is a directory with a `SKILL.md` manifest plus instructions.

2. Skills are loaded with precedence
   OpenClaw supports bundled, user, and workspace-level skill directories, with later sources able to override earlier ones.

3. Skills are filtered by eligibility
   Skills can declare requirements, and only usable skills are exposed.

4. Skills are loaded economically
   The agent first sees a compact catalog, not the full body of every skill.

5. Skills are part of policy, not just prompts
   Per-workspace or per-agent skill scope lets one runtime expose different capabilities safely.

6. Runtime loop protections matter
   Repetition detection and execution guardrails are treated as runtime responsibilities, not prompt wishes.

For ChemCoworker, all six patterns are relevant.

## Current ChemCoworker Seams

The new skill layer should attach to existing abstractions instead of replacing them.

### Existing abstractions worth preserving

- `ToolPlugin` already describes executable tools, contracts, validators, and LLM exposure.
- `ToolRegistry` already filters tool names and exposes a public/private distinction.
- `WorkflowDefinition.llm_visible_tools` already acts like a coarse tool policy profile.
- `ChemCoworker._build_native_tools()` already composes the tool surface seen by the model.
- `ChemCoworker._run_native_tool_loop()` is already the right place for on-demand skill hydration.

### Existing gaps

- workflow behavior is mostly defined in Python constants and system prompts
- there is no external capability packaging for domain playbooks
- there is no environment-aware gating above raw tools
- there is no workspace-local specialization layer for chemistry teams or project contexts
- there is no compact "available skills" catalog given to the model

## Proposed Architecture

Add a new package:

```text
chem_coworker/
  skills/
    __init__.py
    manifest.py
    loader.py
    eligibility.py
    registry.py
    catalog.py
```

Add a new repo folder for bundled skills:

```text
chem_coworker_skillpacks/
  bundled/
    reaction-analysis/
      SKILL.md
    condition-recommendation/
      SKILL.md
    retrosynthesis-route-planning/
      SKILL.md
    literature-curation/
      SKILL.md
```

Add optional local override locations:

```text
./skills/
~/.chem_coworker/skills/
```

Load precedence:

1. bundled skills
2. user skills
3. workspace skills

Later layers override earlier skills by `id`.

## Design Principle: Tools Execute, Skills Teach

The boundary should stay hard:

- `ToolPlugin` remains the only executable primitive
- `WorkflowDefinition` remains the main task router
- `SkillManifest` adds usage policy, instructions, and eligibility metadata

Skills should answer:

- when should this capability be used
- which tool subset is relevant
- what evidence is required before answering
- what failure modes should be surfaced
- what chemistry caveats must appear in the answer

Skills should not answer:

- how Suzuki coupling is typed
- how functional groups are detected
- how a route is validated
- what deterministic chemistry rule is true

Those stay in `chemtools`.

## Proposed Manifest Schema

Each skill directory contains `SKILL.md` with YAML frontmatter followed by markdown instructions.

Example:

```md
---
id: condition_recommendation
name: Condition Recommendation
version: 1
summary: Use taxonomy-backed reaction typing and HTE evidence to recommend reaction conditions.
category: chemistry
tool_policy: conditions_specialist
workflow_targets:
  - forward_chemistry
  - forward_synthesis
triggers:
  keywords:
    - recommend condition
    - suggest condition
    - what catalyst
    - what solvent
  requires_any:
    - reaction_smiles
    - reaction_context
tool_allowlist:
  - analyze_reaction
  - recommend_reaction_conditions
  - reagent_assistant
tool_preferred_order:
  - analyze_reaction
  - recommend_reaction_conditions
provides_context:
  - reaction_type
  - recommendations
requires_context:
  - reaction_type
eligibility:
  rdkit: true
  python_modules:
    - rdkit
  data_files:
    - data/HTE_db
  env_vars: []
  taxonomy_ids: []
prompting:
  inject_mode: on_demand
  max_tokens_hint: 700
answer_contract:
  require_tool_evidence: true
  require_taxonomy_alignment: true
  must_surface_warnings: true
  forbid_knowledge_only_conditions: true
priority: 80
---

## When to use

Use this skill for requests that ask for catalysts, ligands, bases, solvents, temperature,
or screening strategy for a specific transformation.

## Required behavior

1. Establish reaction identity from taxonomy-backed tool evidence before recommending conditions.
2. Prefer HTE-backed recommendations over unaudited general knowledge.
3. If evidence is weak, say so explicitly and downgrade confidence.
4. Do not invent exact conditions when no successful condition tool evidence exists.

## Failure modes

- Ambiguous reaction class
- Missing HTE coverage
- Reactant roles unresolved
- Product-target mismatch
```

### Parsed model

Recommended dataclass shape:

```python
@dataclass
class SkillManifest:
    id: str
    name: str
    version: int
    summary: str
    category: str
    tool_policy: str | None = None
    workflow_targets: list[str] = field(default_factory=list)
    triggers: SkillTriggers = field(default_factory=SkillTriggers)
    tool_allowlist: list[str] = field(default_factory=list)
    tool_preferred_order: list[str] = field(default_factory=list)
    provides_context: list[str] = field(default_factory=list)
    requires_context: list[str] = field(default_factory=list)
    eligibility: SkillEligibility = field(default_factory=SkillEligibility)
    prompting: SkillPrompting = field(default_factory=SkillPrompting)
    answer_contract: SkillAnswerContract = field(default_factory=SkillAnswerContract)
    priority: int = 50
    instructions_md: str = ""
    source_path: str = ""
```

### Chemistry-specific eligibility fields

OpenClaw uses general runtime requirements. ChemCoworker should extend that idea with chemistry-aware checks:

- `rdkit`
- `python_modules`
- `data_files`
- `env_vars`
- `taxonomy_ids`
- `reaction_families`
- `provider_models`

These let the runtime hide skills that cannot actually run on the current machine.

## New Runtime Components

### 1. `SkillLoader`

Responsibilities:

- discover skill directories
- parse `SKILL.md`
- validate frontmatter
- merge by precedence
- attach source metadata

Recommended API:

```python
class SkillLoader:
    def load_all(self, workspace_root: Path | None = None) -> list[SkillManifest]:
        ...
```

### 2. `SkillEligibilityChecker`

Responsibilities:

- validate each manifest's runtime requirements
- report why a skill is ineligible
- produce both eligible and suppressed skill lists

Recommended output:

```python
@dataclass
class SkillEligibilityResult:
    skill_id: str
    eligible: bool
    reasons: list[str] = field(default_factory=list)
```

### 3. `SkillRegistry`

Responsibilities:

- store eligible skills
- expose a compact skill catalog
- filter by workflow, query triggers, or tool policy
- resolve a specific skill on demand

Recommended API:

```python
class SkillRegistry:
    def catalog_for_workflow(self, workflow_name: str) -> list[SkillManifest]:
        ...

    def match_for_query(self, query: str, task_type: str, smiles_present: bool) -> list[SkillManifest]:
        ...

    def get(self, skill_id: str) -> SkillManifest | None:
        ...
```

### 4. `SkillCatalogFormatter`

Responsibilities:

- generate the compact metadata-only block seen by the model
- avoid full instruction injection unless selected

Example output:

```text
Available skills:
- condition_recommendation: HTE-backed condition selection for known transformations
  Workflows: forward_chemistry, forward_synthesis
  Tools: analyze_reaction, recommend_reaction_conditions, reagent_assistant
- retrosynthesis_route_planning: route planning and proposal validation
  Workflows: retrosynthesis
  Tools: retrosynthesis_step, plan_route, validate_synthesis_proposal
```

## Integration Points In Existing Code

### `chem_coworker/workflow.py`

Add skill policy hooks to `WorkflowDefinition`:

```python
@dataclass
class WorkflowDefinition:
    ...
    llm_visible_tools: Optional[List[str]] = None
    default_skill_ids: Optional[List[str]] = None
    tool_policy: Optional[str] = None
```

Recommended initial mapping:

- `retrosynthesis`
  - `tool_policy="retro_specialist"`
  - `default_skill_ids=["retrosynthesis_route_planning", "route_validation"]`
- `forward_synthesis`
  - `tool_policy="forward_specialist"`
  - `default_skill_ids=["forward_prediction", "condition_recommendation"]`
- `forward_chemistry`
  - `tool_policy="general_chemistry"`
  - `default_skill_ids=["reaction_analysis", "condition_recommendation", "literature_curation"]`

### `chem_coworker/agent.py`

Add skill support in three places.

#### A. Agent initialization

Initialize a skill registry once:

```python
self.skill_registry = build_default_skill_registry(...)
```

#### B. Native prompt assembly

Before `llm.bind_tools(native_tools)`, generate a compact skill catalog for the current workflow and add it to the system message.

Do not inject full skill instructions here.

#### C. On-demand skill hydration

Inside `_run_native_tool_loop()`, when the model appears to need a skill or when workflow defaults require one:

- load the matching skill body
- append it as an additional `SystemMessage` or `HumanMessage` instruction block
- continue the loop with the same tool surface

This mirrors OpenClaw's "catalog first, details on demand" approach without requiring a new model capability.

### `chem_coworker/tools/_base.py`

Keep `ToolPlugin` as-is for execution, but consider adding optional metadata that improves skill matching:

```python
skill_tags: List[str] = field(default_factory=list)
tool_policies: List[str] = field(default_factory=list)
```

This lets skills and workflows refer to stable policy names instead of hardcoding tool names everywhere.

### `chem_coworker/tools/__init__.py`

Keep registry ownership here, but add helper methods:

- `filtered_names_for_policy(policy_name)`
- `describe_tools_for_policy(policy_name)`

This gives skills a stable intermediate layer between markdown and concrete tool IDs.

## Prompting Model

ChemCoworker should use a two-stage skill prompting strategy.

### Stage 1: Always include compact skill catalog

This catalog contains only:

- skill id
- summary
- workflow targets
- core tool subset
- high-level usage note

### Stage 2: Inject full skill only if needed

Inject the full markdown instructions when:

- the workflow declares it as default
- the query strongly matches the skill trigger
- the model repeatedly makes poor tool choices in a domain the skill covers

This avoids token bloat while keeping domain instructions reachable.

## Tool Policy Profiles

OpenClaw's separation of capabilities and agent scope is worth copying directly.

Add named tool policies:

- `general_chemistry`
- `conditions_specialist`
- `retro_specialist`
- `literature_curator`
- `internal_debug`

Each policy resolves to a curated tool set using the existing registry.

This is cleaner than putting large tool lists directly in every workflow and every skill.

## Safety And Chemistry Guardrails

The main risk in a markdown skill system is logic drift. The mitigation is to make skills policy-only.

Required rules:

1. Skills may reference taxonomy IDs and reaction families, but may not redefine them.
2. Skills may require evidence, but may not override deterministic tool outputs.
3. Skills may tighten answer behavior, but may not bypass validators.
4. Skills may hide unavailable capabilities, but may not force-enable hidden tools.
5. Any chemistry claims inside skills should be treated as instructions for tool use, not authoritative chemistry data.

## Loop Protection Improvements

OpenClaw documents repetitive tool-call detection. ChemCoworker should add the same class of guardrail in `_run_native_tool_loop()`.

Add a loop detector keyed by:

- tool name
- normalized args
- success/failure state
- whether new context keys were produced

If the same failing or non-progressing call pattern repeats multiple times:

- stop re-executing it
- inject a runtime warning into the conversation
- ask the model to conclude with current evidence

This is independent of skills, but should be implemented alongside them.

## CLI Additions

Add:

```text
python -m chem_coworker._cli.app skills list
python -m chem_coworker._cli.app skills show condition_recommendation
python -m chem_coworker._cli.app skills doctor
```

`skills doctor` should report:

- eligible skills
- suppressed skills
- missing requirements
- missing datasets
- missing API keys

This is especially useful for chemistry deployments where runtime assets matter.

## Initial Skill Set

Start small. The first bundled skills should correspond to existing major workflows:

1. `reaction_analysis`
2. `condition_recommendation`
3. `retrosynthesis_route_planning`
4. `forward_prediction`
5. `literature_curation`

Do not start with dozens of skills.

## Suggested Rollout

### Phase 1: Read-only skill infrastructure

- implement manifest parsing
- implement eligibility checks
- implement catalog rendering
- expose CLI inspection commands
- no prompt injection yet

### Phase 2: Prompt-side catalog

- add compact skill catalog to workflow prompt assembly
- no full on-demand hydration yet

### Phase 3: On-demand skill hydration

- allow workflow-default and query-matched skills to inject full instructions
- add tests for prompt-size control

### Phase 4: Tool policies

- add named tool policy profiles
- refactor workflow tool lists to use policies

### Phase 5: Loop protections

- add repeated-call detection
- add non-progress bailout behavior

## Testing Strategy

Add tests for:

- manifest parsing and validation
- precedence and override behavior
- ineligible skill suppression
- workflow-target filtering
- query-trigger matching
- compact catalog formatting
- on-demand hydration decisions
- loop protection behavior

Representative test files:

```text
tests/test_skill_manifest.py
tests/test_skill_loader.py
tests/test_skill_eligibility.py
tests/test_skill_registry.py
tests/test_skill_prompting.py
tests/test_native_loop_repetition_guard.py
```

## Why This Is Worth Doing

This design gives ChemCoworker four practical gains:

1. Cleaner specialization
   Chemistry playbooks move out of Python constants and giant prompts.

2. Better environment awareness
   The model sees only capabilities that can actually run in the local chemistry environment.

3. Better workspace customization
   Different teams or projects can add local skill overlays without editing the main codebase.

4. Better prompt efficiency
   Detailed instructions are loaded only when needed.

## What Not To Copy From OpenClaw

Do not import these patterns blindly:

- marketplace-first thinking
- user-authored skills with unrestricted execution semantics
- putting chemistry logic directly into markdown
- treating skills as equivalent to tools

For this repository, skills are a policy and prompting layer above deterministic chemistry tools.

## Recommended First Implementation Cut

The smallest high-value version is:

1. `SkillManifest`
2. `SkillLoader`
3. `SkillEligibilityChecker`
4. `SkillRegistry`
5. compact skill catalog injection

That is enough to get most of the OpenClaw benefit without destabilizing the current runtime.

## Sources

OpenClaw docs consulted on March 8, 2026:

- Skills: <https://docs.openclaw.ai/tools/skills>
- Tools overview: <https://docs.openclaw.ai/tools>
- Agent runtime: <https://docs.openclaw.ai/concepts/agent>
- Multi-agent sandbox tools: <https://docs.openclaw.ai/tools/multi-agent-sandbox-tools>
- MCP tools: <https://docs.openclaw.ai/tools/mcp-tools>

Relevant local integration points:

- `chem_coworker/tools/_base.py`
- `chem_coworker/tools/__init__.py`
- `chem_coworker/workflow.py`
- `chem_coworker/agent.py`
