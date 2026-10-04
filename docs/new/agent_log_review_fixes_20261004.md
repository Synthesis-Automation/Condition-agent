# Scientific workspace log review fixes — 2026-10-04

The reviewed route investigation completed but left its methylation and oxidation
graphs unresolved. Source acquisition and scientific assessment were useful;
source quantity reconciliation and repeated nested-input repairs were gaps.

## Source quantity evidence

`condition_registry.quantity_audit.audit_source_quantities` compares complete,
adjacent parenthetical mass/amount pairs with the molecular weight of an explicitly
supplied graph. `quantity_consistency.v1@1.0` declares units and a 5% discrepancy
tolerance; the definition hash is included in registry baselines and each check.
Both reports, literal text, passage offsets and graph-conditional arithmetic are
retained. No value is selected, corrected or converted into an approved recipe.

Literature preparation v2 (operation contract 3) returns `quantity_checks` and
copies conflict cautions into its immutable `literature_reaction.v2` block.
Answer preflight reports `source_quantity_conflict`. Empty checks mean no eligible
pair was found; nonadjacent values, volumes, solution annotations and unsupported
units are not guessed. Invalid graphs and ambiguous organic mixtures remain
`not_assessed`. This does not verify source identity, compound assignment or purity.

For example, the [patent's Reference Example 4](https://patents.google.com/patent/US20130143877A1/en)
lists methyl iodide as 162 g and 2.60 mol. Conditional on `CI`, these quantities
disagree: 2.60 mol corresponds to about 369 g. The original reports remain intact.

## Reaction inputs and compatibility

`reaction_recipe_assessment.v3` adds the existing typed taxonomy
`ReactionCompletenessAssessment` to direct recipe results. This exposes product
element deficits without recomputing chemistry or changing admission, ranking,
mapping, operator selection or compatibility rules. The proposed-recipe operation
is contract 2. Summaries and answer preflight retain the deficits.

An input containing one `CI` while gaining two product carbons has a carbon
deficit. An oxidation graph gaining oxygen without a supplied oxygen-bearing
contributor has an oxygen deficit. Recipe quantities and solvent identities do
not supply graph atoms. Inspect actual source evidence and represent justified
contributors/multiplicity explicitly; never invent donors or mappings to pass.
Complete elemental accounting alone still does not prove atom correspondence,
selectivity, compatibility or feasibility.

## Input help and identity resolution

`w.help('fetch_source')` now documents the existing helper. On-demand help for
recipe assessment and literature preparation includes nested schemas derived
from their actual typed contracts. Provenance must be a JSON object. Invalid
provenance now reports the component/stage path before resolution. Literature
claims may omit empty `limitations`; required attribution and literal quotations
remain enforced. Answer-helper help exposes claims without an `id` field and
instructs callers to repair failed preflight before finalization.

Two explicit common-name aliases were added to existing substances:

- Acetic acid → `cas:64-19-7`, retaining canonical `AcOH`, structure and roles.
  Identity source: [NIST Acetic acid](https://webbook.nist.gov/cgi/cbook.cgi?ID=C64197).
- Methyl iodide → `cas:74-88-4`, retaining canonical Iodomethane and its unassigned
  role. Identity source: [NIST Iodomethane](https://webbook.nist.gov/cgi/cbook.cgi?ID=C74884).

No records were added or merged, and no contextual roles were inferred from the
new aliases. Registry hashes change; source corpus, indexes, labels and admission
tiers were not rebuilt or relabeled. Saved investigation evidence is unchanged.
New recipes using these names now resolve to their existing substance IDs, so
their canonical recipe IDs can change from the previous raw-name identities.
Future corpus conversion may improve identity coverage for these exact aliases;
no corpus-wide coverage or admission improvement has been measured here.
Restart the server and start a new conversation to use the new scientific baseline.

## Validation and release scope

The focused pre-change baseline passed 68 tests. Regressions cover consistent and
conflicting quantities, supported units, invalid/ambiguous inputs, salts, literal
provenance, deterministic replay, incomplete reaction graphs, partner ordering,
closed claims, help discovery and registry validation.

Final validation: `python -m pytest -q` passed **2,705 tests**, with **1 skipped**
(820.32 seconds). The focused integration check passed 178 tests; the final
quantity-boundary and recipe checks passed 33 tests. `git diff --check` passed.
Local logs and recorded evidence live under
`results/ai_native/log_review_fixes_20261004/` and are not committed.

A fresh-baseline verification using exact saved browser bytes and source arguments
found one methyl-iodide mass/amount conflict, product-element deficits of C:1 and
O:1 for the first and third steps, and no unresolved identities in the three
proposed recipes. The Fischer recipe now reports no known conflict; applicable
reaction-capability coverage still remains incomplete. The original saved answer
remains readable and unchanged, without fabricated historical checks.

Engineering regression checks do not satisfy independent chemistry review or
untouched evaluation gates. These changes expose uncertainty and improve inputs;
they do not validate the reviewed route or its proposed operating conditions.
