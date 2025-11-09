"""
Rule Builder Specification – Phase 0 Outcomes
=============================================

This document captures the schema inventory and workflow requirements gathered
from the current rule databases under `data/rule_db_v2/`. It serves as the
reference for implementing the deterministic builder (Phase 1) and later
assistant integrations.
"""

# Overview

- **Schema Version:** All existing files declare `"schema_version": "2.0"` and
  `"source_type": "rule"`. The builder must default to (and validate against)
  schema 2.0 until a change is introduced centrally.
- **Top-level fields:** `schema_version`, `source_type`, `metadata`, `reaction`,
  `applies_if`, `default_rule`, `base_rules`, `modifiers`.
- **Families covered:** Suzuki, Sonogashira, SNAr, Pd/Cu C–N coupling, C–O
  coupling, amide formation, reductive amination, RCM. Any builder output must
  slot directly into the unified index consumed by
  `chemtools.recommend.unified.UnifiedRecommender`.

# Field Inventory

| Field | Required | Notes |
|-------|----------|-------|
| `metadata.id` | ✓ | Unique slug, reused elsewhere (e.g., unified index). |
| `metadata.name` | ✓ | Friendly name displayed in UIs. |
| `metadata.version` | ✓ | Semantic string; builder should enforce non-empty. |
| `metadata.created_date` | ✓ | `YYYY-MM-DD`. |
| `metadata.status` | ✓ | Typically `active`/`draft`. |
| `metadata.tags` | ✓ | Non-empty list of keywords. |
| `reaction.family` | ✓ | Must match taxonomy names (e.g., `Suzuki_Miyaura`). |
| `reaction.reference_reactions` | ✓ | ≥3 reaction SMILES strings observed across files. |
| `reaction.scope.scope_type` | optional | usually `broad`, `targeted`, or `htc`. |
| `reaction.scope.compatible_functional_groups` | optional | list of feature slugs. |
| `reaction.scope.incompatible_functional_groups` | optional | list of feature slugs. |
| `reaction.notes` | optional | Long-form rationale/protocol summary. |
| `applies_if` | ✓ | `all`/`any` keys referencing feature flags (see `chemtools.util.functional_groups`). |
| `default_rule.id` | optional | Present in most files; builder will auto-assign if absent. |
| `default_rule.description` | optional | Free text. |
| `default_rule.conditions` | ✓ | Dict; keys vary per family (e.g., `pd_source`, `ligand`, `solvent`). |
| `base_rules[]` | ✓ | Non-empty list; each entry needs `name`, `id`, `description`, `reactant_features`, `conditions`. |
| `base_rules[].reactant_features` | ✓ | Uses `all`/`any` lists referencing feature slugs. |
| `modifiers[]` | optional | Each entry: `id`, `when` (feature/symptom list), `suggest` (string). |

## Condition Field Variants

- **Cross-couplings (Suzuki, Sonogashira, C–N/C–O):** expect `pd_source`,
  `precatalyst`, `ligand`, `base`, `base_equiv`, `solvent`, `temperature_C`,
  `time_h`, optional `additive`, `atmosphere`.
- **Condensation/amide-like:** common keys include `coupling_system`,
  `stoichiometry_equiv`, `solvent`, `solvent_volume_vol_per_g`,
  `temperature_C`, `time_h`, `order_of_addition`, `notes`.
- **Reductive amination:** adds `reducing_agent`, `acid`, `drying_agent`.

The builder must treat condition blocks as opaque dictionaries validated for
non-empty string values while allowing family-specific keys.

# Interview & Data Capture Workflow

To keep new rule entries consistent, the builder will ask for data in the
following order:

1. **Identification**
   - Reaction family (`reaction.family`)
   - Unique ID slug (`metadata.id`)
   - Friendly name + version + tags

2. **Evidence Pack**
   - ≥3 reference reaction SMILES (`reaction.reference_reactions`)
   - Protocol/notes block (free text) summarizing experimental trends
   - Optional scope annotations (compatible/incompatible functional groups)

3. **Applicability**
   - Feature gates (`applies_if.all`, `applies_if.any`)
   - Default rule description and baseline conditions

4. **Rule Set**
   - Number of base rules to seed
   - For each: human label, ID, triggering features (`reactant_features`), and
     a structured condition table

5. **Modifiers**
   - Symptom or feature triggers (prefixed with `symptom:` when needed)
   - Suggested adjustments/rationale

6. **Validation & Review**
   - Run deterministic validation (schema, feature IDs, required fields)
   - Present diff/summary before writing to disk

# Validation Requirements

Deterministic checks (Phase 1 deliverable):

- Ensure top-level required keys exist.
- Validate metadata strings (non-empty, correct formats).
- Confirm at least one reference reaction and that all entries resemble
  reaction SMILES (`reactants>>products` with `.` separators).
- Ensure `applies_if` uses known logical keys (`all`, `any`) and references
  non-empty feature IDs.
- For each rule:
  - `conditions` dict must be non-empty.
  - `reactant_features` uses `all`/`any`/`and` semantics consistent with
    `RuleSpec`.
  - `id` uniqueness within the file.
- For modifiers:
  - Each `when` entry is non-empty.
  - `suggest` is descriptive text.

The builder module will expose `validate(strict=True)` to raise descriptive
errors and `list_issues()` to return warnings for deferred cleanup.

# Interactive Session & CLI

- `chem_assistant/rule_builder_session.py` implements the wizard flow (metadata →
  references → scope → applicability → default rule → base rules → modifiers)
  on top of `RuleBuilder`. Prompt handling uses injectable `prompt_fn` /
  `output_fn` hooks so CLIs, tests, and LLM agents can share the same logic.
- `python -m chem_assistant.rule_builder_cli` launches the wizard from the
  terminal.
  - `--load PATH`: edit an existing rule DB (use `--output` to write elsewhere).
  - `--family NAME --output PATH`: create a new database for the given family.
- After the wizard completes, the CLI runs deterministic validation and prints
  a diff (`RuleBuilder.diff()`). Blocking validation errors prevent saving until
  resolved.

# LLM-Assisted Drafting (`rule_builder_autofill_tool`)

- Tool path: `chem_assistant/chemtools_wrapper.py` (`rule_builder_autofill_tool`).
- Input schema: `RuleBuilderAutoInput` collects metadata (`metadata_id`,
  `metadata_name`, `metadata_version`, optional date/status/tags), reference
  reactions, `protocol_text`, optional `notes`, `desired_focus`,
  `applies_if_hints`, `modifier_hints`, and `max_base_rules`.
- Workflow:
  1. Deterministic metadata + reference reactions populate a new `RuleBuilder`.
  2. Protocol text, hints, and focus instructions feed the LLM template
     (`llmtools.prompts.RULE_BUILDER_EXTRACTION`) to produce structured JSON.
  3. The structured payload maps through the deterministic builder (scope,
     applies_if, default rule, base rules, modifiers).
  4. Validation runs automatically and returns both the serialized database and
     any issues that need manual cleanup.
- Output contract:
  ```json
  {
    "success": true,
    "rule_database": {...},
    "issues": [{"field": "...", "message": "...", "severity": "warning"}],
    "message": "Draft created with warnings"
  }
  ```
- Usage ideas:
  - Generate first drafts from automation notebooks or protocol dumps, then
    refine interactively via the CLI wizard.
  - Provide `applies_if_hints` to force required gating features and
    `modifier_hints` for high-priority lab symptoms.

## CLI & Automation Helpers

- `python scripts/rule_builder_autofill_demo.py` provides a small CLI wrapper
  around the autofill tool. Supply metadata, references (inline or via existing
  JSON), protocol text (string or file), and optional hints; the script will:
  1. Run the LLM-assisted draft.
  2. Print validation findings.
  3. Write the resulting JSON to `--output` if provided.
- Combine this script with the `builder` CLI command for a two-step workflow:
  autofill to generate a seed, then launch the wizard to review/curate before
  committing to `data/rule_db_v2/`.
