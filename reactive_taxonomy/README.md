# Reactive Taxonomy

`reactive_taxonomy` is the standalone molecular and reaction-chemistry package
for the type-agnostic condition-recommendation system. Molecular graphs are the
source of truth. Reaction names and source labels are optional annotations and
never determine structural evidence.

## Reaction workflow

```text
Reaction SMILES
    ↓
Parsed molecular graph facts
    ↓
Cross-side atom correspondence
    ↓
Normalized edits or explicit alternative edit hypotheses
    ↓
Structural observation
    ├── minimum reaction core
    ├── graph topology and ring changes
    └── deterministic L0–L4 reaction signature
    ↓
Optional molecular annotations and reaction-pattern interpretation
    ↓
One concise and one detailed chemist-facing label
```

The minimum core and reaction signature are built only from graph structure,
atom correspondence, normalized edits, and local graph environments. Molecular
motifs, reactive-site hypotheses, transformation patterns, synthesis patterns,
named families, and display labels cannot create or modify edits.

When correspondence is ambiguous, the package retains typed edit hypotheses
with provenance and abstains from producing a canonical signature. An optional
external atom mapper may resolve that ambiguity, but cannot silently override a
conflicting internal observation.

## Public layers

- `observe_molecular_structure(smiles)` returns element, charge, bond, ring, and
  component facts.
- `interpret_molecular_reactivity(structure)` adds optional motifs, reactive-site
  hypotheses, contexts, and reactivity profiles.
- `analyze_molecule(smiles)` composes both molecular layers.
- `observe_reaction(reaction_smiles)` returns the interpretation-free reaction
  observation, minimum core, topology, and evidence alternatives.
- `interpret_reaction(observation)` adds optional pattern and family annotations.
- `render_reaction(observation, interpretation)` produces the unified concise and
  detailed display label.
- `featurize_reaction(reaction_smiles)` composes the canonical public analysis.

## Definitions

Identity-affecting definitions are intentionally narrow. The taxonomy manifest
marks only structural signature definitions as identity inputs. Site patterns,
molecular motifs, reaction patterns, source-label mappings, and rendering rules
are annotation definitions.

Reaction patterns are split into two optional tiers:

- transformation patterns describe common net edit shapes such as substitution,
  addition, elimination, bond-order change, coupling, and ring closure;
- synthesis patterns add stronger chemist-facing interpretations such as
  Suzuki-like coupling, amide formation, reductive-amination-like change, or
  cycloaddition.

Pattern definitions contain no graph operators, predicted products, structural
slots, or reconstruction instructions.

## Validation

```powershell
python -m reactive_taxonomy.cli validate
pytest -q tests/reactive_taxonomy
```

All SMARTS compilation goes through
`reactive_taxonomy.chemistry.smarts_cache.compile_smarts`.

The concise architectural walkthrough is in
[`docs/reaction_featurization_workflow.md`](../docs/reaction_featurization_workflow.md).
