# Strategy-led planning and molecular analysis update

Implemented after the stepwise pyridone investigation, 2026-10-11.

## Behavior

The existing general molecular analysis, focused inspection, SSR suppliers and
recorded scientific workspace remain the public workflow. No new autonomous
planner or unvalidated route-ranking score is introduced.

1. Analyze target identity and graph observations. The audit recognizes 2-pyridone
   and aromatic cyclic C(=O)-N without assigning ordinary amide reactivity tags.
   Direct carbonyl–nitrogen observations survive missing motif interpretation.
2. Keep distinct construction strategies in recorded decision notes, including
   entry materials, preparation burden, evidence references, weakest step, next
   decisive question and a stop/revisit reason. Finishing variants sharing a
   bottleneck are not counted as independent strategy discoveries.
3. Use literature and local precedents to establish accessible intermediates and
   construction evidence. Focus SSR on concrete precursor questions; allow bounded
   open exploration where it could reveal an overlooked strategy.
4. Inspect precursor chemistry as well as product similarity. Both step and
   condition precedent inspection now return `transfer_review`, distinguishing
   exact structures, the same product from different precursors, different products,
   and unresolved input. Graph-derived site inventory differences flag donor and
   leaving-group changes, including boron versus tin. They do not identify the
   actual reaction center, exclude spectators or certify condition transfer.
5. Spend the next investigation on the weakest consequential step, inspect recipe
   rule coverage, and retain source procedure/quantity gaps. `no_known_conflict`
   is never promoted to experimental feasibility.

These decision checks are included directly in the Normal and Tools + guidance
workspace prompt. The optional retrosynthesis adviser contains further details.
The ablation modes that intentionally omit guidance retain that behavior.

## Molecular contract

`audit_target` / workspace `analyze_molecule` now return `target_audit.v2`:

- map-free canonical identity, formula, molecular weight, charge and stereo;
- perceived rings, motifs, sites and independent carbonyl–nitrogen observations;
- `canonical_atom_to_input_atom`, input atom-map numbers, and per-component atom
  and bond projections;
- explicit scope, limited pattern coverage, definition hashes and RDKit version.

Stereo, ring and carbonyl observation indices use the whole canonical molecule.
Motif/site atom and bond indices remain component-local, now consistently derived
from canonical SMILES. Project through `components` for multi-component structures.
Input atom references must not be reused as focused disconnection IDs. The same
canonicalization helper is used by focused inspection/disconnection. Tautomers,
chemical forms, salts and specified stereo remain distinct.

The compact summary exposes coverage and omitted counts; full observations remain
in the saved artifact. Existing `inspect_reactive_sites` and `compare_molecules`
provide detailed follow-up. Custom code should answer a named gap and be recorded
with its inputs, source, versions, outputs and evidence references. Recurrent checks
should be promoted into the owning domain package after positive/negative tests.
Existing workspace evidence and SMARTS caches are reused; no new result cache.

Motif definition schema remains `1.0`, with definition version `1.1` and changed
content hash. The new motifs have empty reactivity tags; the ordinary amide rule
is unchanged. This is not exhaustive coverage of every pyridone tautomer/subtype.
Reaction-signature identity definitions are unchanged; motif definitions belong
to the annotation layer. Saved annotation/descriptor provenance still changes.

## Registry corrections and provenance

Four existing records changed; no substance IDs were added:

| CAS | Correction |
| --- | --- |
| 72287-26-4 | Correct Pd(dppf)Cl2 composition; add common formula aliases. |
| 95464-05-4 | Replace erroneous bromoaniline graph with Pd(dppf)Cl2·CH2Cl2 composition; add catalyst capability and explicit adduct aliases. Remove incorrect `Pd(dppf)2Cl2` alias. |
| 127-08-2 | Add potassium acetate alias. |
| 7778-53-2 | Add potassium phosphate tribasic alias. |

Supplier identity/formula evidence:
[unsolvated Pd(dppf)Cl2](https://www.sigmaaldrich.com/DE/en/product/aldrich/697230)
and [dichloromethane adduct](https://www.sigmaaldrich.com/US/en/product/aldrich/379670),
accessed 2026-10-11. Formula checks require C34H28Cl2FeP2Pd and
C35H30Cl4FeP2Pd respectively. The registry uses charge-balanced disconnected
coordination-component representations; these encode composition, not solution
speciation or metal coordination geometry. Formula agreement alone is not a
general identity-validation method. The base aliases refer to the unchanged
existing potassium salts.

The existing content-addressed registry version changes automatically. Unsolvated
and adduct identities are not merged. Old name-only `Pd(dppf)2Cl2` records now remain
unresolved rather than being silently assigned the adduct. CAS-based identity is
retained. Corpus-wide consequences of reconversion have not been measured.

## Capability and migration

`w.capabilities()['operation_readiness']` reports missing required artifacts from
the operation catalog. The [completion fix](workspace_completion_fix_20261011.md)
also checks SQLite condition-index metadata at startup and reports incompatibility;
corpus content and other prerequisites remain validated on call.
Configure a new investigation with the intended processed release and required
indexes/catalogs instead of changing a pinned investigation or repeatedly calling
a tool with missing inputs.

Workspace contracts: analyze_molecule 2, inspect_step_precedents 5,
inspect_condition_precedents 3. Condition evidence schema is now
`condition_evidence_comparison.v2`; historical v1 evidence remains readable without
inventing a transfer review. New typed structural review schema is
`precedent_transfer_review.v1`.

Restart the application and create a new investigation for changed tools/definitions.
Historical artifacts are not rewritten. Existing converted corpora/indexes contain
old definition identity; do not relabel them as current or weaken their compatibility
checks. A versioned reconversion/rebuild and normal validation gates are required
before claiming that a corpus incorporates these chemistry/registry corrections.

## Validation boundary

Regression tests cover the investigation target, N-H/N-methyl and Kekule pyridones,
ordinary amides/lactams, negative heterocycles, separate tautomers, graph projections,
stereo conflicts, donor changes, unknown input and registry formula/charge checks.
Workspace integration tests cover saved evidence and compact review output.

A local recorded smoke check is available at
[`planning_update_20261011/review.html`](../../results/ai_native/planning_update_20261011/review.html).
It includes the numbered target, full v2 audit, canonical atoms 5–6 inspection,
the boron/tin precursor contrast, registry lookups, scripts and recorded artifacts.
These generated files are local and are not required for another checkout.

This implements the reusable engineering changes. It does not establish faster
route discovery, experimental success or improved multi-target chemist acceptance.
A controlled multi-target comparison, independent review and untouched evaluation
remain follow-up work. The current release gates are unchanged.

Validation result: `pytest -q` — **2980 passed, 4 skipped** in 803.54 seconds.
The initial run exposed a summary-size regression; the projection was corrected
to keep atom tables in the full artifact, followed by a successful complete run.
No chemistry snapshots were regenerated. `git diff --check` also passed.
