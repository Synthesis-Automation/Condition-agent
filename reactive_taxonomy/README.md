# Reactive Taxonomy

Standalone, deterministic compound and reaction featurization built directly
on RDKit. This package must not import `chemtools`.

Public entry points:

- `featurize_molecule(smiles)`
- `featurize_reaction(reaction_smiles)`
- `resolve_source_label(source_label)`
- `validate_taxonomy()`
- `validate_source_label_mappings()`

Integrated CLI tester:

```powershell
python -m reactive_taxonomy.cli validate
python -m reactive_taxonomy.cli self-test
python -m reactive_taxonomy.cli molecule "Brc1ccc(N)cc1C#N"
python -m reactive_taxonomy.cli reaction "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/molecule_features.jsonl
python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/sample_compounds_featurized.csv
python -m reactive_taxonomy.cli batch examples/sample_reactions.csv --mode reaction --output results/reaction_features.jsonl
```

Use concise mode for a chemist-readable view containing only the primary
labels and interpretation:

```powershell
python -m reactive_taxonomy.cli molecule "Brc1ccc(N)cc1C#N" --concise
python -m reactive_taxonomy.cli reaction "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" --concise
```

For reactive sites bearing alkyl groups, the readable output separates the
substitution at the reactive center from attachment-carbon sterics. For
example, tert-butylamine remains a primary amine (`R–NH2`) but its attached
alkyl group is reported as `tertiary` and alpha-branched.

Alkyl reactive centers retain `Alkyl` as their broad machine context while
recording synthetically important attachment subtypes. Benzylic, allylic, and
propargylic leaving groups are rendered as labels such as `Benzyl–Cl`,
`Allyl–Br`, and `Propargyl–OMs`.

Add `--format json` to the single-record and validation commands for the full
typed result. Batch mode prints a coverage summary and optionally writes
source-traceable JSONL or CSV. The format is inferred from the output extension
or selected with `--output-format`. CSV contains readable summary columns;
JSONL retains the complete typed analysis. Use `--column NAME` when a CSV does
not use a standard `smiles` or `reaction_smiles` header.

Definitions own handle detection, functional groups, rendering, reaction
grammars, and descriptor weights. Python implements graph interpretation,
candidate resolution, reaction operators, and typed result contracts.

Dependency direction:

```text
applications and recommenders -> reactive_taxonomy -> RDKit
```

The recommender and dataset conversion layers belong outside this package.

## Source-label normalization

Legacy dataset labels are resolved through the versioned
`definitions/source_label_crosswalk.v1.json` file. Resolution separates a
stable machine label from a chemist-facing display and optional environment
constraints. For example, `RNH2 a-branch` resolves to machine label `R-NH2`,
display label `R–NH₂ (α-C branched)`, and an `alpha_branched=true` qualifier.
Unsupported labels are returned unchanged with `mapping_status=unresolved`.
