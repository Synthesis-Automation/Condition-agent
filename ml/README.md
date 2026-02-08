# Chan-Lam ML Rebuild (`/ml`)

Clean rebuilt system for condition recommendation.

## Architecture

- `ml/contracts.py`: Config contracts
- `ml/data.py`: Dataset loading + LOGO splits + leakage checks
- `ml/chemistry.py`: Descriptor and reaction featurization primitives
- `ml/condition_features.py`: Engineered condition-property features
- `ml/features.py`: Feature profile definitions
- `ml/models.py`: Stage-1 prior + Stage-2 ML models
- `ml/recommender.py`: Two-stage recommender
- `ml/eval.py`: LOGO backtest harness
- `ml/metrics.py`: APYR/Spearman metrics
- `ml/cli.py`: Unified CLI

## CLI

### Train

```powershell
python -m ml.cli train `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --artifact "results/ml/chanlam_rebuild/recommender.joblib"
```

### Evaluate (LOGO)

```powershell
python -m ml.cli evaluate `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --output-json "results/ml/chanlam_rebuild/evaluation_report.json" `
  --max-folds 0
```

### Tune Production Defaults (Phase 3)

Runs LOGO backtest, then grid-searches `shortlist_k` and blend weights
(`w_stage1`, `w_stage2 = 1 - w_stage1`) from fold-level predictions.

```powershell
python -m ml.cli tune-defaults `
  --input-csv "data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv" `
  --eval-output-json "results/ml/chanlam_rebuild/evaluation_for_default_tuning.json" `
  --output-json "results/ml/chanlam_rebuild/default_tuning_report.json" `
  --max-folds 0
```

Current production defaults are set from this tuning step in `ml/contracts.py`.

Composite objective tuning (recommended):

- default objective mode is `composite`
- balances APYR + top1@70 + Spearman + top3@70
- includes small preference for non-zero stage2 weight and an explicit
  `--min-stage2-weight` guard (default `0.10`) to avoid stage1-only collapse

Use `--objective-mode apyr_lexicographic` to recover the previous APYR-first behavior.

### Recommend (new reaction)

```powershell
python -m ml.cli recommend `
  --artifact "results/ml/chanlam_rebuild/recommender.joblib" `
  --reaction-smiles "NS(=O)(=O)c1cccc(F)n1.COc1ccc(B(O)O)cc1>>COc1ccc(NS(=O)(=O)c2cccc(F)n2)cc1" `
  --output-csv "results/ml/chanlam_rebuild/new_reaction_top_conditions.csv" `
  --top-k 10
```

## Feature Profiles

- `core_full`: catalyst/base/solvent + descriptor features
- `base_motif_spectator`: base + motif/spectator minimal profile

Use with `--feature-profile`.

## Model Types

- `rf_reg` (default)
- `lgbm_rank`

Use with `--model-type`.
