# Reactive Taxonomy Reaction Benchmark

- Source rows: 610
- Valid reactions: 607/610
- Exact product reconstructions: 200
- Reactant-grammar-only candidates: 54
- Unresolved/unsupported: 356
- Runtime: 3.586 s (170.1 reactions/s)

## Focus-family results

| Declared type | Rows | Exact | Reactant-only | Unresolved |
| --- | ---: | ---: | ---: | ---: |
| Suzuki | 46 | 31 | 13 | 2 |
| Buchwald-Hartwig amination | 125 | 105 | 18 | 2 |
| Buchwald-Hartwig | 2 | 2 | 0 | 0 |
| SNAr | 5 | 4 | 0 | 1 |
| Ullmann CO coupling | 5 | 5 | 0 | 0 |
| CS coupling | 5 | 5 | 0 | 0 |
| Sonogashira | 9 | 6 | 3 | 0 |
| Negishi | 5 | 5 | 0 | 0 |
| Kumada | 5 | 5 | 0 | 0 |
| Stille | 3 | 3 | 0 | 0 |
| ChanLam coupling | 3 | 1 | 2 | 0 |
| Amidation | 2 | 2 | 0 | 0 |
| Esterification | 2 | 2 | 0 | 0 |

## Interpretation

- Exact reconstruction is confirmation by graph equality, not a named-mechanism guess.
- Named families remain compatible-family metadata because most rows do not include catalyst or condition evidence.
- Remaining Suzuki misses include unsupported MIDA/irregular boron encodings and multi-event examples.
- Remaining C–N misses are mostly bis-couplings, product/reactant size inconsistencies, or invalid product fixtures.
- Three sample reactions contain invalid product SMILES and are returned as invalid.
- Multi-event synthesis, site-guided MCS fallback, and tautomer-aware verification remain later extensions.
