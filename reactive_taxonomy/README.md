# Reactive Taxonomy

Standalone, deterministic compound and reaction featurization built directly
on RDKit. This package must not import `chemtools`.

Public entry points:

- `featurize_molecule(smiles)`
- `featurize_reaction(reaction_smiles)`
- `resolve_source_label(source_label)`
- `validate_taxonomy()`
- `validate_source_label_mappings()`

Every reaction with normalized edit evidence also receives a structured
`display_label`. Exact grammar labels remain available as overlays, while
mapped unknown-family reactions fall back to deterministic labels such as
`C–N bond formation` or `C=C → C–C; 2 × H gain at C`. Display style is
explicitly excluded from reaction-signature identity.

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

Phase 1 molecular-feature evaluation, including machine metrics and a blind
chemist-review packet:

```powershell
python -m reactive_taxonomy.molecular_feature_evaluation_cli `
  results/molecular_feature_evaluation
```

The versioned answer key is
`benchmarks/molecular_features/benchmark_manifest.v1.json`. Generated artifacts
include `machine_report.json`, detailed case results, `chemist_review.csv`,
highlighted `review_structures.html`, and `disagreements.csv`. The machine gate
does not complete the separate human chemist gate.

Phase 2 reaction-edit evaluation, including mapped edits, exact reconstruction,
hydrogen changes, evidence conflicts, and a blind reaction-center review packet:

```powershell
python -m reactive_taxonomy.reaction_edit_evaluation_cli `
  results/reaction_edit_evaluation
```

Its versioned answer key is
`benchmarks/reaction_edits/benchmark_manifest.v1.json`, and it writes the same
five evaluation artifacts under `results/reaction_edit_evaluation/`.

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

Unsaturated-bond labels expose endpoint substitution rather than collapsing all
sites to `C=C` or `C≡C`. Examples include `H2C=CH2`, `H2C=CHR1`,
`R1R2C=CR3R4`, `R1–C≡C–H`, and `R1–C≡C–R2`; defined alkene E/Z
stereochemistry is retained. Canonical `PI|Alkene` and `PI|Alkyne` signatures
remain display-independent.

Organic nitriles are also bond-localized π handles with typed carbon and
nitrogen endpoints (`PI|Nitrile`, `R–C≡N`). Organic azides are separate
multi-atom dipolar handles with attachment, proximal-, central-, and
terminal-nitrogen roles (`DG|Azide|Organic`, `R–N3`). Their declared addition,
reduction, or cycloaddition modes describe possible chemistry only; a reaction
interpretation still requires product/edit evidence. Cyanide salts, isocyanides,
and inorganic azide are not promoted to these organic-handle contracts.

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
