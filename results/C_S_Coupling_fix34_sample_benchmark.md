# Fix 3/4 Sample Benchmark (C_S_Coupling)

- Generated: 2026-02-15 14:58:29
- Source: `data\HTE_db\literature\C_S_Coupling_canonical.csv`
- Sample: 800 rows (seed=20260215)
- Sample input: `results\C_S_Coupling_fix34_sample_input.csv`
- Sample output: `results\C_S_Coupling_fix34_sample_output.csv`
- Audit CSV: `results\C_S_Coupling_fix34_sample_audit.csv`

## Problematic Metrics (Before vs After)

| Metric | Before | After |
|---|---:|---:|
| `total_rows` | 800 | 797 (-3 vs before) |
| `unknown_detected_reaction_type` | 10 | 0 (-10 vs before) |
| `unresolved_detected_reaction_type` | 2 | 0 (-2 vs before) |
| `missing_reaction_key` | 0 | 0 (0 vs before) |
| `missing_reaction_events` | 0 | 0 (0 vs before) |
| `reaction_key_empty_or_degenerate` | 0 | 0 (0 vs before) |
| `reaction_smiles_no_arrow` | 0 | 0 (0 vs before) |
| `reaction_smiles_missing_side` | 0 | 0 (0 vs before) |

## Routing / Scope Audit on Sample

- `routing_excluded`: 3
- `out_of_scope`: 0
- `cleanup_applied`: 16
- coordination components removed: 0
- counterion components removed: 16

## Manual Review Sample: Out-of-Scope Rows

- No out-of-scope rows in this sample.
