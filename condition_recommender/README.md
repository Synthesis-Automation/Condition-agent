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

## Ranking logic

The recommender:

1. Verifies the reaction and assigns a grammar with `reactive_taxonomy`.
2. Filters the dataset to compatible source reaction families.
3. Matches FG A and FG B as an unordered pair.
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
- base, catalyst, solvent, ligand, additives, and time/temperature text;
- a short explanation of the match.

## Important limitations

The source CSV contains reaction labels rather than reaction structures for
its precedent rows. Recommendations are therefore weak-label precedents, not
structure-verified literature matches. FG A/B order is not reliable and is
treated as unordered. Condition names have also not yet been normalized to
condition-registry identities.

