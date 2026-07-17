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
`display_label`. An exact, edit-consistent grammar label is preferred when one
is available. Otherwise, versioned generic edit patterns render transformations
such as `C–N substitution`, `C=C hydrogenation`, or `N=N reductive cleavage`.
Any unmatched edit set still receives a deterministic literal summary such as
`C–N bond formation`. The structured result retains the grammar or pattern ID,
the literal structural label, evidence, and confidence. Display style and label
definitions are explicitly excluded from reaction-signature identity.

For a single validated bond-order transformation, the preferred display adds
retained-neighbor and endpoint-hydrogen context. For example, the generic
interpretation `C=C hydrogenation` is retained in `transformation_label`, while
`reaction_label` becomes `Ar–CH=CH2 → Ar–CH2–CH3`. The nested display also
retains `structural_label`, `reactant_context_label`, and
`product_context_label`. Multi-edit cases that cannot be rendered safely keep
their literal edit summary.

Reductive amination is represented as one declarative two-partner grammar plus
one registered graph operator, rather than a reaction-name special case. The
grammar accepts available aldehyde, ketone, or formaldehyde carbonyl sites and
free primary or secondary amines. The operator removes the carbonyl oxygen,
forms C–N, adds H at carbon, and consumes N–H. Exact product reconstruction can
therefore render labels such as
`HeteroAr–CH=O + Ar–NH2 → HeteroAr–CH2–NH–Ar`; product mismatches remain
unverified, and multiple indistinguishable assignments remain unselected.

The high-ROI common-reaction layer also includes explicit alkyl C–N, C–O, and
C–S substitution plus terminal-alkene Heck coupling. Alkyl substitution uses
the shared leaving-group and X–H site contracts and the existing
handle-replacement operator. It intentionally leaves `named_family` unset
because the graph alone does not distinguish SN1, SN2, or protection chemistry.
The Heck grammar uses alkene endpoint hydrogen counts to select a terminal
attachment site, removes Ar–X, forms Ar–C, and records alkene C–H loss. It does
not invent E/Z stereochemistry: a stereospecified product is exact only when
that stereochemistry is supported by the input/operator result.

Friedel–Crafts acylation is represented by an activated acyl electrophile and
an available aromatic C–H site. It reuses the generic handle-replacement
operator to remove the acyl leaving group, form the acyl-C–aryl bond, and record
aromatic C–H loss. Carboxylic acids are excluded because their acyl sites are
latent rather than activated. Regioisomers are selected only by exact product
reconstruction; unresolved sites and mapping conflicts remain visible.

Intramolecularity is a shared reaction-topology dimension rather than a second
set of family grammars. Grammar roles declare `same`, `different`, or
`same_or_different` component relationships. The same graph operators therefore
handle intermolecular joining and same-component ring closure. Every normalized
edit signature now carries `ReactionTopology`, including reaction scope, role
component membership, formed-bond scope, reactant tether distance, formed ring
size, and graph cycle-rank delta. Topology contributes to L0–L2 identity while
L3 remains a topology-agnostic bond-edit fallback. Mapped unknown reactions get
the same topology analysis without requiring a named family.

For example, `NCCc1ccccc1Br>>c1ccc2c(c1)CCN2` is rendered as
`intramolecular (5-membered ring) Ar–Br / R–NH2 → Ar–NH–R` using the same
`sp2_c_n_substitution` grammar as the corresponding intermolecular reaction.

When an unmapped reaction has exactly one conserved heavy-atom scaffold and one
product, a conservative correspondence fallback may supply edit evidence after
registered grammar reconstruction has failed. It accepts only mappings whose
best alternatives imply the same normalized chemistry, requires every product
heavy atom to be accounted for, and reports chemically distinct alternatives as
`ambiguous_atom_correspondence`. Multi-substrate assembly remains grammar- or
mapping-dependent rather than being guessed by this fallback.

Integrated CLI tester:

```powershell
python -m reactive_taxonomy.cli validate
python -m reactive_taxonomy.cli self-test
python -m reactive_taxonomy.cli molecule "Brc1ccc(N)cc1C#N"
python -m reactive_taxonomy.cli reaction "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/molecule_features.jsonl
python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/sample_compounds_featurized.csv
python -m reactive_taxonomy.cli batch examples/sample_reactions.csv --mode reaction --output results/sample_reaction_featurized.csv

python -m reactive_taxonomy.cli batch examples/dataset_300/C_N_Coupling.csv --mode reaction --output results/data_300.csv


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

Organic heteroatom-pair bonds share a bond-localized handle family with two
typed attachments: `HB|Azo` (`R1–N=N–R2`), `HB|Disulfide`
(`R1–S–S–R2`), and `HB|Peroxide` (`R1–O–O–R2`). Hydroperoxides,
hydrogen peroxide, thiols, ordinary ethers, and unattached diazene are excluded.
Mapped reductions can therefore report the observed central-bond cleavage and
the corresponding hydrogen gains without requiring a named reaction family.

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
