# Chan-Lam ML Rebuild: Phase 4 Benchmark

## Setup

- Dataset: `data/HTE_db/literature/ChanLam_dataset_converted (2)_canonical.csv` (3892 rows)
- Evaluation protocol: LOGO (44 folds), artifact `results/ml/chanlam_rebuild/evaluation_phase4_defaults.json`
- Active config: `core_full`, condition props enabled, `shortlist_k=30`, `rf_reg` (`n_estimators=220`), blend `w_stage1=0.8`, `w_stage2=0.2`

## Strict Side-by-Side Metrics

| Metric | Stage1 | Stage2 | Final (0.8/0.2) |
|---|---:|---:|---:|
| APYR (mean) | 84.486 | 81.696 | 83.418 |
| Spearman (mean) | 0.568 | 0.588 | 0.465 |
| Top1 hit @ 50% | 75.00% | 67.05% | 75.00% |
| Top1 hit @ 70% | 43.18% | 43.18% | 43.18% |
| Top3 hit @ 50% | 81.82% | 87.50% | 77.27% |
| Top3 hit @ 70% | 53.41% | 57.95% | 53.41% |
| Top5 hit @ 50% | 81.82% | 92.05% | 88.64% |
| Top5 hit @ 70% | 56.82% | 61.36% | 61.36% |

## Key Deltas

- Final vs Stage2: APYR `+1.722`, Spearman `-0.123`, Top1@70 `+0.00 pp`, Top3@70 `-4.55 pp`, Top5@70 `+0.00 pp`.
- Final vs Stage1: APYR `-1.068`, Spearman `-0.103`, Top1@70 `+0.00 pp`, Top3@70 `+0.00 pp`, Top5@70 `+4.55 pp`.
- Winner (APYR): Stage1
- Winner (Spearman): Stage2
- Winner (Top1@70): tie
- Winner (Top3@70): Stage2
- Winner (Top5@70): tie (Stage2/Final)

## Files

- Full evaluation JSON: `results/ml/chanlam_rebuild/evaluation_phase4_defaults.json`
- Per-fold metrics: `results/ml/chanlam_rebuild/evaluation_phase4_defaults.per_fold.csv`
- Prediction rows: `results/ml/chanlam_rebuild/evaluation_phase4_defaults.predictions.csv`
- Compact benchmark table (CSV): `results/ml/chanlam_rebuild/benchmark_phase4_summary.csv`
