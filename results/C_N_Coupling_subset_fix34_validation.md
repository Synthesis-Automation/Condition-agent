# C_N_Coupling Subset Validation for Fixes 3 and 4

- Generated: 2026-02-15 14:48:51
- Input subset: `results/C_N_Coupling_subset_fix34_input.csv`
- Intent: verify coordination/counterion cleanup and out-of-scope C-N routing on representative rows

| Original Index | Decision | Reason | Cleaned | Coord Removed | Counterion Removed |
|---:|---|---|---|---:|---:|
| 1 | kept | - | no | 0 | 0 |
| 2 | kept | - | no | 0 | 0 |
| 15 | kept | - | no | 0 | 0 |
| 17 | kept | - | no | 0 | 0 |
| 18 | kept | - | no | 0 | 0 |
| 27 | kept | - | no | 0 | 0 |
| 29 | kept | - | no | 0 | 0 |
| 30 | kept | - | no | 0 | 0 |
| 164 | kept | - | no | 0 | 0 |
| 1065 | kept | - | yes | 1 | 0 |
| 1066 | kept | - | yes | 1 | 2 |
| 1069 | kept | - | yes | 1 | 1 |
| 1072 | kept | - | yes | 1 | 0 |
| 1073 | kept | - | yes | 1 | 0 |
| 1074 | kept | - | yes | 1 | 1 |
| 4120 | out_of_scope | formed_motif_conflict_ar_ar | no | 0 | 0 |
| 4121 | out_of_scope | formed_motif_conflict_ar_ar | no | 0 | 0 |
| 8116 | kept | - | no | 0 | 0 |
| 23042 | routing_excluded | exclude_organometallic_or_coordination_complex | no | 0 | 0 |

## Summary

- Kept: 16
- Routing excluded: 1
- Out of scope: 2
