# TEMPLATE_EVAL_REPORT

## Build integrity

- fragments: 10
- templates: 3
- expansions: 16
- generated patterns: 16
- SMARTS compile successes: 16
- SMARTS compile failures: 0

## Consistency with atomic features

- generated tokens found in atomic features: 16/16 (100.0%)
- generated SMARTS found in atomic features: 16/16 (100.0%)
- duplicate generated tokens: 0
- duplicate generated SMARTS (under multiple tokens): 0

## Functional match tests

- passed: 6/6

| test | smiles | expected tokens | status | notes |
|---|---|---|---|---|
| bromobenzene electrophile | `Brc1ccccc1` | aromatic_bromide_present | PASS |  |
| chlorobenzene electrophile | `Clc1ccccc1` | aromatic_chloride_present | PASS |  |
| iodobenzene electrophile | `Ic1ccccc1` | aromatic_iodide_present | PASS |  |
| phenylboronic acid | `OB(O)c1ccccc1` | aromatic_boronic_acid_present | PASS |  |
| phenyl Bpin (one common SMILES) | `c1ccccc1B1OC(C)(C)C(C)(C)O1` | aromatic_bpin_present | PASS | unexpected=['aromatic_boronic_acid_present'] |
| vinyl bromide (optional) | `C=CBr` | vinyl_bromide_present | PASS |  |

