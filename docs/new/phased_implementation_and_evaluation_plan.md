# Phased Implementation and Evaluation Plan

## Purpose

Implement the type-agnostic reaction-condition system in independently
verifiable phases. Each phase must produce both machine-testable artifacts and a
chemist-readable review packet. Do not advance until both evaluation gates pass.

## Phase 1 implementation status

**Machine implementation: complete. Human review gate: accepted for progression
on 2026-07-16.**

The project-owner review found the current outputs chemically reasonable for
the 23 simple curated compounds and reported no systematic disagreement. This
acceptance permits Phase 2 work; it does not remove the stated need for broader,
independent chemist-authored cases before freezing the feature vocabulary for
production.

### Major achievements

- Added the versioned benchmark
  `benchmarks/molecular_features/benchmark_manifest.v1.json`.
- Curated 23 cases spanning primary and secondary amines, attachment-carbon
  branching, aryl and heteroaryl leaving groups, boron transfer handles, acyl
  and carbonyl centers, O-H and S-H sites, sulfonamide ownership, epoxide,
  alkene, alkyne, steric/electronic context, a true negative, and invalid input.
- Assigned every case exactly once to development, validation, or untouched test
  partitions (7/8/8 cases in v1).
- Added exact functional-group and reactive-site count evaluation.
- Added exact atom-localization checks for 21 curated reactive sites.
- Added targeted environment checks for nitrogen substitution, attached-group
  sterics, alpha branching, electronic class, nearby groups, and ortho
  hindrance.
- Added repeat-run determinism and equivalent-SMILES invariance checks.
- Added taxonomy-definition validation to the same release gate.
- Added a blind `chemist_review.csv`, atom-indexed highlighted-structure HTML,
  detailed machine case results, and a disagreement-resolution template.
- Added deterministic regression tests for the benchmark and every generated
  artifact.

### Current machine results

The current `molecular_features.v1` run reports:

| Metric | Result |
| --- | ---: |
| Curated cases | 23 |
| Critical case pass rate | 100% |
| Functional-group precision | 100% |
| Functional-group recall | 100% |
| Reactive-site precision | 100% |
| Reactive-site recall | 100% |
| Reactive-atom localization | 21/21 |
| Targeted environment checks | 6/6 |
| Repeat-run determinism | 100% |
| Equivalent-SMILES invariance | 3/3 |
| Taxonomy validation errors | 0 |

These results establish internal consistency against the curated v1 answer key;
they are not a substitute for independent chemical review or a larger external
validation set.

### Run the Phase 1 evaluation

From the repository root:

```powershell
python -m reactive_taxonomy.molecular_feature_evaluation_cli `
  results/molecular_feature_evaluation
```

The command writes:

```text
results/molecular_feature_evaluation/
  machine_report.json
  case_results.jsonl
  chemist_review.csv
  review_structures.html
  disagreements.csv
```

`machine_report.json` is the automated gate. The machine portion passes only
when every configured threshold passes. `case_results.jsonl` retains expected
and observed values for debugging but must not be given to blind reviewers.

### Human chemist evaluation

1. Open `review_structures.html` and inspect the atom-indexed structures. Orange
   atoms are detected reactive sites; blue atoms belong to other detected
   functional groups.
2. Complete the blank fields in `chemist_review.csv` without consulting
   `case_results.jsonl`.
3. Judge reactive-site atoms, functional groups, substitution class,
   attachment sterics, electronic context, missing features, and false
   detections.
4. Copy every disagreement into `disagreements.csv` and classify it as an
   implementation defect, definition defect, ambiguous chemistry, benchmark
   defect, or descriptor limitation.
5. Resolve confirmed defects through a definition/code change and a named
   regression test. Do not merely update the expected snapshot.
6. Mark the human gate complete only after systematic disagreements are
   resolved and any remaining ambiguity is explicitly documented.

### Phase 1 limitations

- The benchmark is deliberately small and curated; it does not estimate
  performance over broad chemical space.
- Steric descriptors are deterministic 2D graph-local features, not conformer,
  Sterimol, or buried-volume calculations.
- Electronic context is the interpretable qualitative `tag_distance_v1`
  descriptor, not a partial-charge or quantum calculation.
- Reaction-relative spectators and retained-group evidence belong to Phase 2,
  because they cannot be determined from a molecule alone.
- More chemist-authored external cases are required before freezing the feature
  vocabulary as a broad production contract.

## Phase 2 implementation status

**Machine implementation: complete. Human chemist gate: pending review.**

Phase 2 adds schema-level hydrogen-change extraction to the existing typed
formed, broken, and order-changed edits. Mapped evidence and exact operator
reconstruction are reconciled; contradictions remain explicit review evidence.

The versioned `reaction_edits.v1` benchmark contains 12 cases split equally
between development, validation, and untouched test partitions. It covers
single and multi-edit reactions, hydrogen gain and loss, exact reconstruction,
mapped/unmapped parity, mapping defects, a no-edit negative, invalid product
valence, and conflicting mapping/reconstruction evidence.

### Current Phase 2 machine results

| Metric | Result |
| --- | ---: |
| Curated cases | 12 |
| Critical case pass rate | 100% |
| Edit precision | 100% |
| Edit recall | 100% |
| Edit detail accuracy | 100% |
| Evidence accuracy | 100% |
| Edit-resolution accuracy | 100% |
| Product-reconstruction accuracy | 100% |
| Invalid-valence rejection | 100% |
| Mapped/unmapped parity | 100% |
| Conflict retention | 100% |
| Repeat-run determinism | 100% |

These metrics establish consistency against the curated v1 answer key, not
broad reaction-space accuracy.

### Run the Phase 2 evaluation

```powershell
python -m reactive_taxonomy.reaction_edit_evaluation_cli `
  results/reaction_edit_evaluation
```

Review `review_structures.html` without consulting `case_results.jsonl`.
Reactant edit atoms are orange and corresponding mapped product atoms are
green. Complete the blank fields in `chemist_review.csv`, record every
disagreement in `disagreements.csv`, and convert confirmed defects into code or
definition changes with named regression tests.

### Structured display-label foundation

The Phase 2 output now exposes a versioned, observation-first display-label
contract for Phase 3. Normalized edits produce deterministic clauses for bond
formation, cleavage, order changes, and hydrogen gain or loss. Exact grammar
labels remain optional higher-specificity overlays; unknown-family mapped
reactions receive generic labels without inventing a family. Conflicting edit
evidence is stated in the label rather than hidden. Display labels are serialized
into generic conversion records but never contribute to signature identity.

This presentation foundation does not complete the Phase 3 signature gate;
partner, spectator, transformation-class, L0-L4, and optional-family review is
still required after the Phase 2 human gate.

## Implementation phases

| Phase | Implementation | Automated evaluation | Human chemist evaluation |
| --- | --- | --- | --- |
| 1. Molecular features | Functional groups, reactive sites, amine substitution, attachment sterics, electronic context, and nearby groups | Curated positive/negative tests, atom-order invariance, SMARTS validation, and deterministic IDs | Blind review of detected sites and local environments on representative molecules |
| 2. Reaction edits | Formed, broken, order-changed, hydrogen, and multi-bond edits from mapping or exact reconstruction | Product reconstruction, valence checks, mapped/unmapped parity, and conflicting-evidence tests | Review reaction-center highlighting and edit descriptions against supplied products |
| 3. Reaction signatures | Stable partners, spectators, transformation classes, L0-L4 keys, and optional family evidence | Partner-order invariance, deterministic signatures, unknown-family cases, and serialization round trips | Review concise transformations, partner roles, spectators, evidence, and confidence |
| 4. Reaction coverage | Add transformation plugins in waves: acyl, redox, addition/elimination, then cyclization | Per-family held-out coverage, reconstruction accuracy, false-positive rate, and full regression suite | Review random accepted, rejected, and ambiguous cases for every new plugin |
| 5. Conditions and recipes | Normalize identities, contextual roles, canonical `RCR1` recipes, and label-dataset condition names | Identity coverage, role-confidence calibration, recipe deduplication, and provenance tests | Review catalysts, ligands, bases, solvents, and ambiguous multi-role substances |
| 6. Recommendation | Compatibility filtering, structural retrieval, weak-label fallback, and recipe ranking | Leakage-safe top-k recovery, yield MAE, coverage, fallback levels, and hard-violation count | Blind ranking of recommended recipes, explanations, cautions, and unsuitable conditions |

## Required evaluation artifacts

Every phase must write:

- `machine_report.json`: coverage, accuracy, failures, confidence, and definition
  versions.
- `chemist_review.csv`: structures, detected features, evidence, predictions, and
  blank reviewer fields.
- `review_structures.html` or rendered images: highlighted reactive atoms, edits,
  and spectators.
- `disagreements.csv`: machine/chemist disagreements and their resolution status.
- `benchmark_manifest.v1.json`: immutable case IDs and train, validation, and test
  assignments.

## Evaluation rules

- Maintain separate development, validation, and untouched test sets.
- Include positive, negative, ambiguous, and conflicting cases.
- Sample accepted, rejected, and low-confidence records for human review.
- Never accept snapshot changes automatically.
- Record inter-reviewer agreement when multiple chemists participate.
- Convert confirmed disagreements into regression tests or definition updates.
- Report results by reaction class and evidence quality, not only as global
  averages.

## Immediate sequence

1. Build the molecular-feature benchmark and chemist review packet.
2. Validate reaction-edit extraction independently of reaction names.
3. Freeze `ReactionSignature` only after phases 1 and 2 pass.
4. Expand reaction coverage one transformation-plugin wave at a time.
5. Normalize structure-rich and label-only condition data.
6. Calibrate recommendation only after chemistry coverage is adequate.

## Release gate

Correct molecular features and bond edits come first. Reaction names, retrieval
weights, and recommendation scores cannot compensate for errors in the molecular
or transformation evidence.
