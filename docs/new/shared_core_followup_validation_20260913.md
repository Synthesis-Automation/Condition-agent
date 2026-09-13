# Shared-core follow-up: disagreement review and larger evaluation

This follow-up preserves the user-selected v2 default and addresses the next
validation priorities. Production chemistry definitions and projection keys are
unchanged. No Full/Compact rebuild is required for the evaluation tooling changes.

## Implemented foundation

`reactive_taxonomy.shared_core_diagnostics` compares existing protected graphs
and tests explicitly non-authorizing counterfactuals: unchanged ring membership,
unchanged center stereochemistry, hydrogen substitution with the same hydrogen
change, and their combination. Connectivity, element, charge, aromaticity,
hybridization, isotope, radical state and edit direction remain protected.
The diagnostic cannot turn a rejected precedent into an eligible recommendation.

The evaluation contract is now **1.6**. Exact map-independent core equivalence
joins evidence groups even when observations cite different publications.
Previous behavior joined that key only for unreferenced rows. Publication-based
independent-support counting in production remains unchanged; evaluation leakage
control and evidence-support counting have different responsibilities.

The evaluator also supports deterministic query caps applied after full group
assignment. Unscored held-out rows never enter training. It records query latency,
candidate-budget warnings, shared-core traces and recommended precedent IDs.
Projection workers and progress callbacks support larger reproducible runs.

`condition_recommender.evaluation_panel.freeze_evaluation_panel` freezes selected
source positions, all train/test/query IDs, protocol, code/definition hashes and
the panel checksum before outcomes are calculated. It refuses to overwrite
frozen inputs. Prior-evidence exclusion is explicit direct token exclusion;
it does not claim transitive exclusion through unsampled source records.

## Disagreement review

The existing consumed-panel audit contained 82 top-precedent comparisons:

- 13 have a match under at least one tested static-context counterfactual;
- 68 remain different under those counterfactuals;
- one has no usable query projection.

Counterfactual counts overlap: seven hydrogen-substitution matches, four
unchanged-ring matches, and 13 matches under the combined profile. No case is
explained solely by unchanged center stereochemistry. These are diagnostic
categories, not chemist compatibility labels, and the two splits reuse one panel.

The packet is `results/shared_core_followup/disagreement_review/review_packet.html`.
It displays graph-participating reactants and observed products, with complete
supplied reactions, original edits and the reported recipe in expandable details.
Automated diagnoses are visible, so it is explicitly **not blind**.
`analyst_triage.md` prioritizes cases DIFF-051, DIFF-055, DIFF-065 and DIFF-070.
The separate chemist form remains unfilled; assistant analysis is not recorded as
independent adjudication. The HTML page accepts decisions and notes and downloads
a review CSV; unreviewed cases stay blank.

The review demonstrates why a blanket removal of ring or hydrogen labels is
premature: some counterfactual matches have different retained functional context,
including imine/oxime products and differently activated carbon partners. Any
future relaxation needs a bounded graph rule and negative controls covering
those distinctions. The 68 unresolved comparisons are not automatically judged
chemically unsuitable.

## Frozen larger experiment

The panel contains **20,000 observations** sampled from the current 571,157-row
Full index. Prior development and consumed panels were excluded by observation,
publication, canonical reaction and map-independent exact core tokens.
Seed: `2026091305`; top-k: five; independent evidence target: two.

| Split | Training observations | All held-out observations | Evaluated queries per engine |
| --- | ---: | ---: | ---: |
| Publication/reaction grouped | 15,944 | 4,056 | 500 |
| Scaffold disjoint | 18,597 | 1,403 | 500 |

Both engines use exactly the same query IDs within each split. Entire test groups
remain excluded when the query cap is applied. Shared projections are built from
the training index only. This is about ten times the previous training size,
but remains a sampled-corpus evaluation, not a Full-training performance estimate.
The two splits are complementary views of one panel, not independent replications.

All artifacts and reproducible execution scripts are in
`results/shared_core_followup/`. The original 2,000-row panel and its reports are
preserved. No chemistry or ranking rules are tuned against the new outcomes.

## Results

| Split | Earlier-engine coverage | V2 coverage | Earlier top-five recipe recovery | V2 top-five recipe recovery |
| --- | ---: | ---: | ---: | ---: |
| Publication/reaction grouped | 90.4% | 80.2% | 13.8% | 8.8% |
| Scaffold disjoint | 88.6% | 75.6% | 21.2% | 11.6% |

The original numerical limits are still missed: coverage loses 10.2 and 13.0
percentage points; top-five recipe recovery loses 5.0 and 9.6 points. The larger
training pool increases absolute coverage relative to the earlier small-panel
experiment, but these are different panels and protocols, not a controlled
training-size learning curve. The new measurements do not erase the gap.

Publication, canonical-reaction and mapping-equivalence overlap is zero in both
splits. Scaffold overlap is zero in the scaffold-disjoint split; shared scaffolds
are intentionally permitted in the grouped split. Both engines report zero
violations detected by the existing hard compatibility rules. This is not an
independent measure of chemical precision or experimental success.

Paired coverage: grouped has 396 queries covered by both engines, 56 only by the
earlier engine, five only by v2 and 43 by neither. Scaffold-disjoint has 373, 70,
five and 52 respectively. V2 reaches its candidate budget on two of 500 queries
per split, so the reported cap is not a sufficient explanation for most misses.
Recorded query timings come from concurrent evaluation runs and are diagnostic
wall times, not a controlled production latency benchmark.

A post-hoc inspection of the top earlier-engine precedent for the 126 new
baseline-only coverage cases found:

| Diagnostic | Grouped | Scaffold disjoint |
| --- | ---: | ---: |
| Protected graph or other context still differs | 44 | 51 |
| Graph evidence unavailable on at least one side | 8 | 7 |
| Static-context counterfactual warrants review | 4 | 12 |
| Already qualifies under the v2 comparator | 0 | 0 |

This audits only the top earlier-engine precedent, not every possible candidate.
It identifies no confirmed lookup omission among those comparisons and does not
establish that every rejected analogue is chemically unsuitable. These new cases
are now consumed diagnostic evidence. They were not used to change chemistry,
ranking or thresholds; the frozen code/definition hashes still match execution.

## Decision and remaining chemistry work

Keep the user-selected v2 operational default and retain the explicit baseline
for controlled comparisons. Do not broaden matching solely to recover the older
coverage number. The next chemistry changes should separate retained functional
context from permissible substitution differences, starting with the four
development cases in the packet. Each accepted case needs a graph-defined rule,
counterexamples and validation of actual source/condition requirements.

Independent chemist adjudication remains open. The interactive review packet
makes that decision concrete; no user approval or experimental evidence has been
invented. After a chemistry rule changes, freeze a further untouched panel rather
than reusing this one as validation. Rebuild projections only if the projection
contract changes; reconvert source records only if their observations change.

## Verification

- Complete `pytest -q`: **1,465 passed in 419.12 seconds**.
- Four frozen runs completed successfully: 500 evaluated queries per engine per
  split, with exact agreement against frozen query IDs and split counts.
- Review-browser check passed: 82 editable cases, priority ordering, correct CSV
  quoting for commas/quotes/newlines, preservation of unreviewed blanks, and zero
  page JavaScript errors. Test decisions were not persisted as reviewer evidence.
- The rendered review packet was visually inspected. Ruff Python name checks and
  `git diff --check` pass.

Machine-readable results: `results/shared_core_followup/validation_summary.json`.
Post-hoc comparisons: `results/shared_core_followup/posthoc_miss_diagnostics.json`.
The prior default-promotion changes in the working tree are preserved.
