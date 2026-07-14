# Repository Guidelines

## Mission

Build a clean, deterministic, type-agnostic reaction-condition recommendation
system. The new system is founded on `reactive_taxonomy`; named reaction families
are optional chemistry annotations, not mandatory routing keys.

The molecular graph is the source of truth. Reaction observations must come from
parsed structures, reactive sites, atom correspondence, bond edits, and local
environments. Source labels and reaction names may contribute evidence, but they
must never override contradictory structural evidence.

The implementation plan in
`docs/new/type_agnostic_reaction_recommendation_implementation.md` is the primary
design reference for this work.

## New-System Architecture

The clean system consists of three standalone packages:

- `reactive_taxonomy/`: molecular and reaction chemistry, including parsing,
  functional groups, reactive sites, reaction grammars, graph operators, bond
  edits, environments, spectators, family evidence, and reaction signatures.
- `condition_registry/`: condition-substance identity, aliases, families,
  contextual roles, provenance, validation, and canonical condition recipes.
- `condition_recommender/`: dataset conversion, admission, indexing, chemistry
  compatibility, retrieval, scoring, recipe aggregation, and explanations.

The dependency direction is:

```text
reactive_taxonomy       condition_registry
          \                 /
           condition_recommender
                    |
               app / CLI / API
```

Package ownership is strict:

- `reactive_taxonomy` must not import `condition_registry`,
  `condition_recommender`, or legacy `chemtools`.
- `condition_registry` must not import `condition_recommender` or legacy
  `chemtools`. It may consume a narrow, typed reaction context without owning
  reaction chemistry.
- `condition_recommender` may import `reactive_taxonomy` and
  `condition_registry`.
- Application layers may compose all three packages but must not contain core
  chemistry or recommendation rules.

## Legacy Boundary

`chemtools/` and its existing application paths are legacy code. They may be read
to understand prior behavior or to migrate deterministic chemistry ideas, but the
new standalone packages must not depend on them.

- Do not add imports from `chemtools` to any new-system package.
- Do not add adapters that make legacy models part of the new public contracts.
- Do not preserve duplicate recommendation or conversion paths permanently.
- When parity is established, remove migrated legacy paths instead of maintaining
  backward compatibility, unless compatibility is explicitly requested.
- New features belong in the standalone package that owns the responsibility,
  not in `chemtools`.

## Core System Principles

### Chemistry first

- Validate graph transformations, valence, reactive handles, functional-group
  constraints, and condition compatibility before similarity scoring.
- Use hard chemistry filters before ranking or aggregation.
- Do not replace deterministic chemistry with reaction-name matching, string
  heuristics, or opaque embeddings.

### Type agnostic, not chemistry agnostic

- Every usable reaction should receive a generic `ReactionSignature` when exact
  reconstruction or valid mapped edits provide sufficient evidence.
- `named_family` is optional and must carry confidence and evidence.
- Missing or ambiguous family identity must not by itself reject a chemically
  verified record.
- High-confidence family evidence may narrow the first retrieval tier; lower
  confidence may affect scoring but must not block generic fallback.

### Separate system layers

Keep these stages explicit:

1. Observation: components, functional groups, sites, atom correspondence, edits,
   environments, and spectators.
2. Interpretation: transformation class, partner roles, product transformation,
   family candidates, confidence, and evidence.
3. Recommendation: admission, compatibility filtering, retrieval, similarity,
   recipe aggregation, ranking, and explanation.

An interpretation failure must not erase valid observations.

### Evidence and uncertainty

Use evidence in this order:

1. validated supplied atom mapping and observed bond edits;
2. exact product reconstruction from a registered taxonomy operator;
3. uniquely supported reactive-site grammar;
4. unresolved or conflicting candidates retained for review.

Never invent atom correspondence, force a named family, or present a predicted
product as observed fact. Preserve ambiguity, conflicts, warnings, confidence,
and provenance in typed outputs.

### Declarative definitions

- Put handles, grammars, rendering rules, feature vocabularies, compatibility
  rules, and scoring weights in validated, versioned definitions.
- Keep graph edits, descriptor calculations, validation, and other executable
  behavior in explicit Python registries.
- Never dynamically import executable code named by arbitrary JSON.
- Include schema and definition versions in serialized analyses and records.

## Current Implementation Priority

Implement the plan in phases. Do not skip ahead when a later phase depends on an
unfinished contract.

The current priority is Phase A: the generic reaction-signature foundation.

1. Add `ReactionAtomReference`, `ReactionEdit`, `ReactionPartner`,
   `ProductTransformation`, and `ReactionSignature`.
2. Normalize selected candidates, mapped bond changes, family environments,
   product connections, and spectators into the generic schema.
3. Generate deterministic, versioned L0-L4 signature keys and `signature_id`.
4. Attach `reaction_signature` to `ReactionAnalysis` without changing current
   family results.
5. Support a valid mapped unknown-family reaction with `named_family=None`.
6. Serialize the nested signature into recommendation records and review output
   without changing admission or retrieval behavior yet.
7. Establish parity tests before starting the family registry or unified
   converter phases.

Temporary compatibility fields are allowed only when they have clear removal
criteria and regression tests.

## Reaction-Signature Requirements

- Signature identity must use normalized chemistry fields, schema versions, and
  definition versions—not display labels or source reaction names.
- Signatures must be invariant to irrelevant reactant ordering and serialization
  details.
- Atom references must retain side, component index, atom index, atom-map number,
  element, charge, aromaticity, hybridization, and local-environment identity.
- Support formed, broken, order-changed, and schema-level hydrogen changes.
- `ProductConnection` may remain temporarily as a single-bond compatibility view;
  `ProductTransformation` is the general contract.
- Partner roles are optional interpretations. The base schema must not assume
  every reaction has an electrophile, nucleophile, or transfer partner.
- Mapped evidence and operator reconstruction must be reconciled. Contradictions
  must produce review evidence rather than silent precedence.

## Condition and Recommendation Requirements

- Normalize conditions through `condition_registry`; do not infer identities or
  roles from source column names alone.
- Preserve raw condition identifiers, resolved substance IDs, contextual roles,
  role confidence, and provenance.
- Use canonical nested JSON or Parquet as the primary converted artifact. CSV is
  a review/export view.
- Run compatibility rules before similarity.
- Retrieval must follow an explicit family-to-generic fallback ladder and report
  the level used.
- Aggregate precedents by canonical resolved recipe, not raw condition strings.
- Recommendations must cite precedent IDs and explain matching edits, handles,
  environments, mismatches, cautions, fallback level, and uncertainty.

## Project Structure

- `reactive_taxonomy/`: foundation of the new chemistry system.
- `condition_registry/`: standalone condition identity and recipe system.
- `condition_recommender/`: standalone conversion and recommendation system.
- `tests/reactive_taxonomy/`, `tests/condition_registry/`, and
  `tests/condition_recommender/`: new-system tests.
- `docs/new/`: current architecture and implementation documents.
- `app/`: application integration; keep domain logic out of this layer.
- `data/` and `datasets/`: curated or source datasets; avoid committing large or
  proprietary artifacts.
- `results/`: local generated analysis; do not commit large outputs.
- `chemtools/`: legacy system, outside the dependency graph of new packages.

## Coding Standards

- Python 3.10 or newer, PEP 8, four-space indentation.
- Type hints on all public functions and dataclass fields.
- Module docstrings and docstrings for public or chemically complex behavior.
- Use `snake_case` for functions and variables, `PascalCase` for classes, and
  `UPPER_SNAKE_CASE` for constants.
- Prefer immutable typed dataclasses for public contracts.
- Keep deterministic libraries free of mutable global state and network calls.
- Keep functions focused; put shared behavior in the owning package rather than
  copying it between family-specific modules.
- Remove obsolete paths during refactors after parity tests pass.

## SMARTS Compilation

All new-system SMARTS compilation must use the standalone centralized cache:

```python
from reactive_taxonomy.chemistry.smarts_cache import compile_smarts

pattern = compile_smarts("[CX4][Cl,Br,I]", validate=False)
if pattern is not None and molecule.HasSubstructMatch(pattern):
    pass
```

Rules:

- Never call `Chem.MolFromSmarts()` directly outside the cache implementation.
- Store SMARTS strings in versioned definitions or module-level constants.
- Compile lazily through `compile_smarts()`.
- Use `validate=True` when validating new or critical definitions.
- Use `compile_smarts_batch()` when validating or warming a collection.
- RDKit is required for molecular operations; do not implement chemistry-free
  fallbacks unless explicitly requested.

## Testing and Validation

- Run the complete suite with `pytest -q` before handing off a change.
- Keep test module names import-safe across subdirectories; use package markers or
  importlib mode so the full suite collects in one invocation.
- Add fast deterministic unit tests for every changed chemistry contract.
- For reaction changes, include positive, negative, ambiguous, and conflicting
  evidence cases.
- Required signature regressions include Suzuki, C-N, C-O, C-S, partner-order
  invariance, mapped unknown-family reactions, invalid maps, and deterministic
  IDs.
- Treat dataset snapshot changes as chemistry changes. Explain changes in
  coverage, admission tiers, labels, or rejection reasons.
- Verify definition loaders and schemas whenever JSON definitions change.
- Do not accept updated snapshots merely because the code changed.

Useful commands:

```powershell
pytest -q
pytest -q tests/reactive_taxonomy
pytest -q tests/condition_registry
pytest -q tests/condition_recommender
```

## Commit and Review Guidelines

- Use Conventional Commits: `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, or
  `chore:` with an optional package scope.
- Keep subjects imperative and at most 72 characters.
- PRs should state the chemistry motivation, contract/schema impact, definition
  version changes, migration/removal impact, and test results.
- If API routes change, verify the application boots and document the API impact.
- Do not commit secrets. Use environment variables for configuration.

## Definition of Clean-System Progress

A change advances the new system only when it:

- strengthens the standalone package contracts;
- derives behavior from molecular evidence and versioned chemistry definitions;
- preserves uncertainty and provenance;
- adds or improves chemistry regression coverage;
- avoids new dependencies on legacy `chemtools`; and
- moves toward one canonical conversion and recommendation path rather than
  maintaining parallel family-specific systems.
