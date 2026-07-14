# Mixed reaction dataset audit

- Files: 118
- Rows: 238214
- Canonicalizable reactions: 236160
- Unique canonical reactions: 186703
- Repeated reaction groups: 19785
- Multi-recipe groups: 11642
- Cross-dataset groups: 7400

## Chemistry sample

- Sampled rows: 1180
- Valid reactions: 1158
- Signatures: 130
- Signature rate: 11.02%

## Field coverage

| Field | Non-empty | Rate |
|---|---:|---:|
| `reaction_id` | 238214 | 100.00% |
| `source_declared_family` | 238214 | 100.00% |
| `reaction_smiles` | 238213 | 100.00% |
| `yield_pct` | 235077 | 98.68% |
| `temperature_c` | 212 | 0.09% |
| `time_h` | 1102 | 0.46% |
| `reference` | 234322 | 98.37% |
| `reactant_cas` | 233878 | 98.18% |
| `product_cas` | 233878 | 98.18% |
| `reagent_cas` | 199295 | 83.66% |
| `catalyst_cas` | 145699 | 61.16% |
| `solvent_cas` | 223618 | 93.87% |
| `experimental_procedure` | 19193 | 8.06% |
| `stages` | 238214 | 100.00% |
| `steps` | 238214 | 100.00% |
| `notes` | 164220 | 68.94% |

Chemistry coverage is measured on the deterministic per-file sample; metadata, identity, registry, and duplicate metrics cover every row.
