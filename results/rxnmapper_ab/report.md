# RXNMapper A/B evaluation

Evaluation date: 2026-08-01

## Scope

- Fresh matched conversion sample: 190 rows, comprising the first 10 rows from
  each of the 19 `examples/small_300` source files.
- Both sample variants used the same checkout, definitions, source rows,
  one-worker conversion, shard size, review export, and index builder.
- The only conversion variable was `use_rxnmapper`.
- Leakage-safe retrieval evaluation used grouped random holdout, seed 17,
  test fraction 0.2, and top-k 5.
- Query-time evaluation used the same current 3,059-row mapper-off index for
  both variants, so only query-time mapping changed.

## Conversion A/B

| Metric | RXNMapper off | RXNMapper on |
| --- | ---: | ---: |
| Records | 190 | 190 |
| Wall time | 41.066 s | 75.147 s |
| Eligible index rows | 84 | 84 |
| Index membership shared | 84 | 84 |

RXNMapper-on was 1.83x as slow (+83.0%) on this CPU-only environment. It
created 35 signatures and 114 cores absent from the off conversion, but none
entered the trusted index.

## Leakage-safe held-out index evaluation

Both variants produced the same 69-row train and 15-row test split and exactly
the same results:

| Metric | Both variants |
| --- | ---: |
| Coverage | 8/15 (53.33%) |
| Abstention | 7/15 (46.67%) |
| Hard-incompatible recommendations | 0 |
| Explanation completeness | 100% |
| Yield MAE | 18.75 percentage points |

The sample is too small to interpret recipe-recovery or yield metrics as model
quality estimates. This test supports index parity, not broad calibration.

## Query-time A/B on one fixed index

| Metric | RXNMapper off | RXNMapper on |
| --- | ---: | ---: |
| Queries | 190 | 190 |
| Wall time | 46.867 s | 58.596 s |
| Valid recommendations | 114 (60.0%) | 133 (70.0%) |
| Abstentions/errors | 76 | 57 |

- RXNMapper preserved all 114 queries that already succeeded internally.
- It rescued 19 additional queries, all with `external_mapping_only` status.
- Among the 114 queries where both variants succeeded, top-1 recipe agreement
  was 88.60% and exact top-k agreement was 87.72%.
- Retrieval level agreed on 81.58% of all queries.
- Query-time RXNMapper increased wall time by 25.0%.

## Interpretation

RXNMapper is not beneficial during trusted dataset conversion for this matched
sample: it increased cost without changing index membership or leakage-safe
retrieval behavior. At query time it provides useful optional coverage for
ambiguous reactions, but those rescued results retain external-mapping review
provenance and should not be treated as verified observations.

Recommended operating policy:

1. Build canonical datasets and trusted indexes with RXNMapper disabled.
2. Run internal query analysis first.
3. Invoke RXNMapper only after the internal path would otherwise abstain.
4. Preserve external evidence status, cautions, and review requirements.

## Full-corpus note

A fresh mapper-off build of all 5,510 rows produced 3,059 eligible index rows.
The existing 3,032-row production mapper-on artifact predates the current
checkout, so its 27-row difference is confounded by intervening code changes
and is not attributed to RXNMapper in this report.
