# Reaction-Label Condition Recommender

This CLI recommends reaction conditions from the cleaned HTE label dataset.
It uses `reactive_taxonomy` to identify the reaction grammar and reactive
functional-group signatures, then returns the five highest-ranked condition
recipes from compatible precedents.

## Basic usage

Run the command from the repository root:

```powershell
python -m condition_recommender.recommend_cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
```

The input must be a complete reaction SMILES containing reactants and a
product:

```text
reactant1.reactant2>>product
```

The CLI returns JSON containing the detected reaction grammar, query FG
signatures, candidate count, and ranked condition recommendations.

## Options

Return a different number of recommendations:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" --top-k 10
```

Use another cleaned dataset:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" `
  --records datasets/reaction_label/v2.1_cleaned.csv
```

Show all options:

```powershell
python -m condition_recommender.recommend_cli --help
```

## Cleaned CSV schema

The recommender uses one flat CSV. Its main columns are ordered for convenient
filtering and manual review:

- `source_reaction_type`
- `reactive_site_1_normalized_label` and `reactive_site_2_normalized_label`
- `reactive_site_1_display_label` and `reactive_site_2_display_label`
- `yield_pct` and `z_score`
- catalyst, ligand, base, solvent, additive, and procedure columns

Additional `reactive_site_1_*` and `reactive_site_2_*` columns preserve the
original source labels, taxonomy signatures, qualifiers, and mapping status.
An unresolved label is retained verbatim and has `mapping_status=unresolved`.

Regenerate the cleaned CSV from the repository root with:

```powershell
python -m scripts.clean_reaction_label_csv `
  datasets/reaction_label/v2.1.csv `
  datasets/reaction_label/v2.1_cleaned.csv
```

## Ranking logic

The recommender:

1. Verifies the reaction and assigns a grammar with `reactive_taxonomy`.
2. Filters the dataset to compatible source reaction families.
3. Matches `reactive_site_1` and `reactive_site_2` as an unordered pair.
4. Ranks primarily by FG-signature similarity.
5. Uses attachment and branching qualifiers to refine close matches.
6. Uses yield, z-score, and precedent support as lower-weight signals.
7. Aggregates identical condition recipes and returns distinct results.

The ranking weights are configured in
`condition_recommender/definitions/label_retrieval.v1.json`.

## Supported reaction grammars

The cleaned dataset currently supports:

- Suzuki–Miyaura coupling
- C–N substitution and coupling
- Amide formation
- C–O substitution and coupling
- C–S substitution and coupling
- Sonogashira coupling
- Negishi coupling

Unsupported or structurally ambiguous reactions return an error instead of
falling back to an unrelated reaction family.

## Output

Each recommendation includes:

- rank and overall score;
- signature and qualifier similarity;
- expected yield and mean z-score;
- number of supporting precedents;
- source row numbers;
- base, catalyst, solvents, ligand, additives, and procedure text;
- a short explanation of the match.

## Important limitations

The source CSV contains reaction labels rather than reaction structures for
its precedent rows. Recommendations are therefore weak-label precedents, not
structure-verified literature matches. Reactive-site order in the source data
is not reliable and is treated as unordered. Condition names have also not yet
been normalized to condition-registry identities.
