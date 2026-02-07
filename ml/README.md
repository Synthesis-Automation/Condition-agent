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

## Metrics

- Repeated random-split CV: Spearman and APYR summary
- Leave-one-sulfonamide-out CV (LOGO): global Spearman + APYR
