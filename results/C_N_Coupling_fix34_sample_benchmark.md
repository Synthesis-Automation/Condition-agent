# Fix 3/4 Sample Benchmark (C_N_Coupling)

- Generated: 2026-02-15 14:53:45
- Source: `data\HTE_db\literature\C_N_Coupling_canonical.csv`
- Sample: 800 rows (seed=20260215)
- Sample input: `results\C_N_Coupling_fix34_sample_input.csv`
- Sample output: `results\C_N_Coupling_fix34_sample_output.csv`
- Audit CSV: `results\C_N_Coupling_fix34_sample_audit.csv`

## Problematic Metrics (Before vs After)

| Metric | Before | After |
|---|---:|---:|
| `total_rows` | 800 | 796 (-4 vs before) |
| `unknown_detected_reaction_type` | 4 | 4 (0 vs before) |
| `unresolved_detected_reaction_type` | 0 | 0 (0 vs before) |
| `missing_reaction_key` | 0 | 0 (0 vs before) |
| `missing_reaction_events` | 0 | 0 (0 vs before) |
| `reaction_key_empty_or_degenerate` | 0 | 0 (0 vs before) |
| `reaction_smiles_no_arrow` | 0 | 0 (0 vs before) |
| `reaction_smiles_missing_side` | 0 | 0 (0 vs before) |

## Routing / Scope Audit on Sample

- `routing_excluded`: 2
- `out_of_scope`: 2
- `cleanup_applied`: 35
- coordination components removed: 0
- counterion components removed: 41

## Manual Review Sample: Out-of-Scope Rows

| sample_row | reason | formed_motifs | reaction_smiles (preview) |
|---:|---|---|---|
| 238 | `formed_motif_conflict_ar_ar` | `Ar-Ar|HeteroAr-H` | `CC1(C)CC(=O)C=C(Nc2ccc([N+](=O)[O-])cc2C#Cc2ccccc2Br)C1>>CC1(C)CC(=O)C=C(n2c(-c3ccccc3Br)cc3cc([N+](=O)[O-])ccc32)C1` |
| 293 | `formed_motif_conflict_ar_ar` | `Ar-Ar|HeteroAr-H` | `C/C=C/N(C(=O)C#Cc1ccccc1)c1cccc2cccc(I)c12>>C/C=C\n1c(=O)c2c(-c3ccccc3)nnn2c2cccc3cccc1c32` |
