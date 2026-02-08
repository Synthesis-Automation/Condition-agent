# Chan-Lam Model (`/ml`)

This folder contains a trainable Chan-Lam condition-selection pipeline using:

- motif/spectator descriptors from the converted HTE table
- molecule-level motif + functional-group tokens
- steric/electronic summaries from `chemtools.featurizers`
- condition variables (catalyst/base/solvent/reagents)

## Train

```powershell
python -m ml.train_chanlam_model `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-dir "results/ml/chanlam"
```

## Outputs

- `results/ml/chanlam/chanlam_rf_model.joblib`
- `results/ml/chanlam/metrics.json`
- `results/ml/chanlam/training_set_predictions.csv`
- `results/ml/chanlam/feature_importance_top200.csv`

## Baseline Comparison

```powershell
python -m ml.evaluate_chanlam `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-dir "results/ml/chanlam_eval"
```

## Tuning + External Holdout

```powershell
python -m ml.tune_and_holdout_chanlam `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-dir "results/ml/chanlam_tune_holdout"
```

Artifacts:

- `results/ml/chanlam_tune_holdout/tuning_and_holdout_report.json`
- `results/ml/chanlam_tune_holdout/external_holdout_sulfonamides.txt`

## Simplified Model (Core Features Only)

This variant uses only:

- Reactant descriptors: motif/electronic/steric (+ spectator/formed motif tokens)
- Conditions: `catalyst`, `base`, `solvent`

Evaluate:

```powershell
python -m ml.simple_chanlam_model evaluate `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-json "results/ml/chanlam_simple/evaluation_report.json"
```

Recommend conditions for a new reaction:

```powershell
python -m ml.simple_chanlam_model recommend `
  --reaction-smiles "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1" `
  --output-csv "results/ml/chanlam_simple/new_reaction_top_conditions.csv" `
  --top-k 10
```

Notes:

- `--model-type`: `rf_reg` or `lgbm_rank` (default is `rf_reg`)
- `--feature-profile`: `core_full` (default) or `base_motif_spectator`
- Condition-property features are enabled by default; disable with `--without-condition-props`
- In current Chan-Lam runs, `rf_reg` outperformed `lgbm_rank` on Spearman and APYR.

Base-only condition profile (reactant motifs/spectators + base):

```powershell
python -m ml.simple_chanlam_model evaluate `
  --feature-profile base_motif_spectator `
  --without-condition-props `
  --output-json "results/ml/chanlam_simple/base_motif_spectator_rf.json"
```

```powershell
python -m ml.simple_chanlam_model recommend `
  --reaction-smiles "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1" `
  --feature-profile base_motif_spectator `
  --without-condition-props `
  --output-csv "results/ml/chanlam_simple/new_reaction_top_conditions_base_only.csv" `
  --top-k 10
```

## Metrics

- Repeated random-split CV: Spearman and APYR summary
- Leave-one-sulfonamide-out CV (LOGO): global Spearman + APYR

## Hybrid CSV+ML PoC

Offline backtest of a hybrid recommender that blends:
- ML prediction
- empirical condition priors from CSV (shrinkage)
- nearest-reaction retrieval score

```powershell
python -m ml.hybrid_recommender_poc `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-json "results/ml/chanlam_hybrid_poc/poc_report_full_logo.json" `
  --max-folds 0
```

Artifacts:
- `results/ml/chanlam_hybrid_poc/poc_report_full_logo.json`
- `results/ml/chanlam_hybrid_poc/poc_report_full_logo.predictions.csv`
- `results/ml/chanlam_hybrid_poc/poc_report_full_logo.per_fold.csv`
