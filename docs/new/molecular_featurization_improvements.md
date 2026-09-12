# Molecular featurization review fixes

This change addresses the September 2026 molecular-featurization review. The
pre-change revision and molecular benchmark are frozen locally under
`results/molecular_features_review/baseline_summary.json`. The preceding full
suite passed 1,287 tests; the existing molecular benchmark passed 23 cases, but
contained only six environment checks and three equivalent-SMILES checks.

## Chemistry and uncertainty

- Homogeneous C=C/C≡C environments measure motif distance from both endpoints.
  Endpoint tie-breaking uses graph symmetry ranks without atom-map numbers.
  Equivalent SMILES and changed atom numbering must preserve feature tokens
  and continuous scores.
- Allylic/propargylic assignments require an adjacent carbon-carbon pi bond.
  Carbonyl and nitrile activation no longer receives these labels.
- Positively charged aromatic nitrogen has a separate cationic classification;
  hydrogen count or substitution degree cannot turn pyridinium into a pyrrole
  donor. Explicit isotopic hydrogens count consistently in descriptor fields.
- Unsupported contexts retain observed atom facts, but steric/electronic
  scores are null and classes are `unknown`. They contribute no similarity
  credit, and rendering explicitly reports missing descriptors.
- A validated substituent matrix in `reactivity_descriptor_rules.v1.json`
  separates inductive and resonance contributions for directly attached
  alkoxy, hydroxy, amino, nitro, cyano, CF3, carbonyl and halogen groups on the
  same six-membered aromatic ring. Resonance is limited to ortho/para positions;
  groups separated by saturated linkers do not receive direct resonance credit.
  Magnitudes remain **uncalibrated chemistry priors**, not measured electronics.
  Other ring sizes and remote fused-ring pathways are outside this matrix.

The existing `Brc1c(C)cccc1C#N` regression now reports `electron_poor` instead
of `slightly_poor`: the ortho cyano group supplies both inductive and resonance
withdrawal. The test checks those explicit contributions as well as the class.

## Recommendation feature projection

`condition_recommender.molecular_features` owns a typed, independently serialized
`ReactionMolecularFeatures` projection. It selects existing molecular hypotheses
only where their loci intersect **observed reactant edits**. Each retained site
cites component/atom indices, edit indices, hypothesis identity, profile status
and definition hashes. Ambiguous edit hypotheses cannot create observed features.

Conversion persists `molecular_features` separately from `reaction_signature`.
Both in-memory and SQLite indexes retain it; environment lookup, live queries,
evaluation and chemist-review queries pass the projection to scoring. Reaction
signature construction, identity, admission and hard compatibility remain
independent of molecular annotations.

Site comparison uses a globally optimal one-to-one assignment, with compatible
roles, context kinds, available site types and center elements. This prevents
greedy matching from making similarity depend on partner order. Unsupported
profiles are excluded. Where both structural substituent and molecular profiles
are available, they split the existing environment component equally according
to the versioned similarity definition; where only one is available, it supplies
that component. The total environment weight remains 0.11.

Retrosynthesis template context also omits unsupported numeric profiles while
retaining the molecular observations. Streamed SQLite metadata records the
versions actually present in admitted rows, including mixed compatible versions;
its artifact identity matches the in-memory builder.

The development regression for para-OMe versus para-CF3 aryl amination verifies
that structural signatures remain equal while molecular scoring resolves their
electronic difference. This is a controlled wiring check, **not evidence of
improved condition recovery on a real corpus**. The independent chemist review
and untouched evaluation gates in the primary roadmap still apply.

## API, storage and migration

| Contract | Current version / behavior |
| --- | --- |
| Composed molecule analysis | 3.1 |
| Site reactivity profile | 1.1; nullable unsupported scores |
| Reactivity descriptor definition | 1.1; algorithm `molecular_profiles.v1.1` |
| Recommendation molecular features | 1.0 |
| Recommendation record / converter | 10.2 / `generic_conversion.v10.2` |
| Persisted generic index | 6.4 |
| Generic similarity definition | 1.1; schema remains 1.0 |
| Fallback descriptor definition | 2.1; incorporates molecular definition hashes |
| Reaction signature identity | Unchanged |

`detect_reactive_site_hypotheses(rdkit_molecule)` now preserves atom and bond
indices in the caller's molecule, including renumbered/disconnected molecules.
Component indices identify `Chem.GetMolFrags(molecule)` components; returned
atom indices remain in the original whole-molecule coordinate system.

The duplicated `site.context_features["environment"]` view is removed. Molecular
consumers obtain profiles from `analysis.reactive_site_environments`, joining on
`hypothesis_id`; CLI consumers have migrated. `include_context_features=False`
skips profile construction and returns no environments. An empty `site_types`
selection now means no sites. The RDKit detection helper returns hypotheses and
their context records without performing full descriptor construction.

**Regenerate converted records and rebuild persisted indexes.** Old index
schemas and stale molecular/fallback definition hashes are rejected. Record
schemas 10.0/10.1 remain readable when their other chemistry contracts are current,
but absent molecular projections provide no new feature evidence. The existing
rich-signature profile input remains a temporary reader compatibility behavior,
covered by the profile/retrieval regression tests. Remove that input together
with the temporary `ReactionPartner` annotation fields at roadmap gate 6 after
current corpora have been reconverted and all callers use the separate projection.

## Validation

New regressions cover randomized SMILES, atom renumbering, disconnected
components, charged heterocycles, explicit isotopic hydrogen, unsupported
profiles, positional/linker negative controls, definition validation, ambiguous
and conflicting edit provenance, SQLite persistence, scoring and identity
separation. Existing molecular benchmark expectations are not retuned.

Local validation and timing artifacts are written beneath
`results/molecular_features_review/`. Timings use warm molecular analysis of four
illustrative development molecules and must not be interpreted as full reaction
conversion throughput. Notation templates are cached immutably; profiles are no
longer eagerly serialized into every hypothesis.

The warm benchmark measured a median **9.16 ms to 7.83 ms per molecule**,
or **14.4% less time**, across five batches of 120 analyses. The unchanged
molecular benchmark passes **23/23 cases**. The controlled amination example
has equal similarity without molecular features (0.610000 each); adding the
features produces 0.720000 for matching para-OMe and 0.708578 for para-CF3.
Detailed contributions are in the local `development_ablation.json` artifact.

Final validation: `python -m pytest -q` passed **1,323 tests in 356.71 seconds**,
including web API/featurization, retrosynthesis and streamed-index parity.
This adds 36 regression cases to the frozen 1,287-test baseline. Ruff's `F`
checks on changed Python files and `git diff --check` also pass. These checks
do not replace the pending blind chemist review and untouched corpus evaluation.
