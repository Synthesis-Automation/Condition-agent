# Reactive Taxonomy

Standalone, deterministic compound and reaction featurization built directly
on RDKit. This package must not import `chemtools`.

Public entry points:

- `featurize_molecule(smiles)`
- `featurize_reaction(reaction_smiles)`
- `validate_taxonomy()`

Definitions own handle detection, functional groups, rendering, reaction
grammars, and descriptor weights. Python implements graph interpretation,
candidate resolution, reaction operators, and typed result contracts.

Dependency direction:

```text
applications and recommenders -> reactive_taxonomy -> RDKit
```

The recommender and dataset conversion layers belong outside this package.
