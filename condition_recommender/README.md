# Condition Recommender

`condition_recommender` contains several condition-recommendation approaches at
different stages of maturity. They share reaction observations from
`reactive_taxonomy` and canonical condition identities and recipes from
`condition_registry`, but they are not interchangeable.

There is currently no automatic fallback from one recommendation path to
another. Callers must choose the path explicitly and preserve its provenance,
warnings, and uncertainty.

## Recommendation paths and status

| Path | Public entry point | Current status | Intended use |
| --- | --- | --- | --- |
| Expert structural rules | `recommend_rule_conditions()` | Phase I C-N engine and definitions implemented; one production rule and multiple review-only drafts | Auditable starter protocols selected from explicit molecular facts |
| Generic structure-backed retrieval | `recommend_generic_conditions()` | Functional pilot; coverage and calibration depend on converted structure-rich records | Type-agnostic precedent retrieval and canonical recipe aggregation |
| Weak-label retrieval | `recommend_conditions_from_labels()` | Functional transitional path with important data limitations | Datasets that contain reaction labels but not precedent structures |

Supporting audit, conversion, indexing, and evaluation commands prepare data for
these paths; they do not recommend conditions by themselves.

## Which path should I use?

- Use expert rules when the reaction is inside a reviewed rule scope and an
  explicit starter protocol is desired.
- Use generic structure-backed retrieval when a compatible converted record set
  or persisted generic index is available.
- Use weak-label retrieval only when structure-rich precedents are unavailable,
  and retain its uncertainty warnings.
- Abstain when the selected path cannot support the query. Do not silently call
  a weaker path and present the result as equivalent evidence.

## Expert rule-based recommendation

### Design

Expert rules live under `definitions/rule_sets/`. They:

- match allowlisted structural facts from `reactive_taxonomy`;
- never route using source reaction names, display labels, or legacy rule IDs;
- reference explicit recipe templates owned by `condition_registry`;
- apply chemistry scope and compatibility checks before returning a recipe;
- report complete match evidence and failed predicates;
- distinguish reviewed production rules from review-only drafts.

Legacy rulebook identifiers may appear only in provenance. They do not determine
whether a rule matches.

### Phase I C-N scope

The first rule set proposes palladium conditions for verified, single-event,
intermolecular `sp2_c_n_substitution` reactions.

The current hierarchy is:

1. combined hindered-Ar-Cl override;
2. Ar-Cl, hindered-Ar-Br, or alpha-branched-amine override;
3. nucleophile-class default;
4. bounded free-amine Ar-Br fallback;
5. typed abstention when no supported rule matches.

Defaults currently distinguish:

- primary alkyl amines;
- secondary alkyl amines, including cyclic secondary amines;
- primary aryl and heteroaryl amines;
- secondary aryl-containing amines;
- aromatic N-H partners;
- primary carboxamides.

Structural overrides and defaults are selected deterministically. Within the
first matching tier, only the highest-priority rule is selected. A matching
draft override blocks a lower-priority production default, preventing unsafe
fallback around a recognized special case.

The Phase I rules abstain from unsupported handles and chemistry, including
Ar-F, uncurated pseudohalides, intramolecular reactions, multi-event reactions,
and unsupported nitrogen classes. An electron-poor heteroaryl halide may also
support SNAr; proposing palladium conditions does not assert that the product
proves a Buchwald-Hartwig mechanism.

### Production status

The primary-amide/aryl-chloride rule is currently the only active rule. It uses
a located primary-literature protocol and returns a fully specified
tBuBrettPhos Pd G3/K3PO4/t-BuOH recipe.

```python
from condition_recommender import recommend_rule_conditions

result = recommend_rule_conditions(
    "Clc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1"
)
```

Production mode excludes every draft rule and draft recipe variant.

### Review the draft C-N rules

The draft C-N defaults and structural overrides are standalone clean-system
rules driven by structural taxonomy facts. Their screening recipes specify
identities, loadings, equivalents, partner stoichiometry, temperature, time,
concentration, and atmosphere, but remain drafts until checked against located
primary procedures.

For chemistry review, export the concise rulebook CSV:

```powershell
python -m condition_recommender.rule_review_cli export `
  results/rule_review/pd_sp2_cn_rulebook.csv
```

The export has one row per rule and explicit recipe variant, with 20 columns
chosen for chemistry review. It includes a readable structural-match summary,
condition names and quantities, operating parameters, rationale, cautions,
source locator, and editable review fields. Alternative condition sets are
separate rows; slash-separated combinations are not generated.

Internal schema fields, individual structural predicates, registry IDs,
definition versions, and provenance JSON are deliberately not flattened into
this sheet. Those belong to the validated canonical rule and recipe-template
definitions, not the chemist review view.

The final columns are blank review fields:

- `review_decision`;
- `reviewer`;
- `review_date`;
- `review_notes`;
- `proposed_change`.

The CSV is a generated review view, not a runtime definition. Editing it does
not change recommendations. Record chemistry feedback in those columns, then
curate and validate the canonical rule JSON and registry recipe-template JSON.
Regenerate the CSV after canonical changes. The `rule_id` and
`recipe_variant_id` link each review row back to those canonical definitions.

The exporter includes both active and draft definitions because its purpose is
review. It does not need `--include-draft`. That flag belongs to reaction
recommendation, where it explicitly opts into returning unreviewed recipes.

Enable them only for review:

```powershell
# Primary-alkyl-amine Ar-Br default
python -m condition_recommender.rule_recommend_cli `
  "Brc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft

# Ar-Cl structural override
python -m condition_recommender.rule_recommend_cli `
  "Clc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft

# Ortho-hindered Ar-Br structural override
python -m condition_recommender.rule_recommend_cli `
  "Cc1cccc(C)c1Br.CN>>Cc1cccc(C)c1NC" `
  --include-draft

# Alpha-branched-amine structural override
python -m condition_recommender.rule_recommend_cli `
  "Brc1ccccc1.CC(C)(C)N>>CC(C)(C)Nc1ccccc1" `
  --include-draft
```

The result reports:

- the selected rule and whether it is a default, structural override, or
  fallback;
- canonical recipe components with quantities;
- temperature, time, concentration, atmosphere, and partner stoichiometry;
- compatibility score, rationale, cautions, and review flags.

This concise review view is the default. Add `--json` to emit the complete
audit payload, including structural match evidence, failures for nonmatching
rules, excluded variants, provenance, and definition versions:

```powershell
python -m condition_recommender.rule_recommend_cli `
  "Clc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft `
  --json
```

Without `--include-draft`, structural draft matches remain visible in the
complete JSON `match_traces`, but their recipes are not returned.

Load and validate the rule definitions directly:

```python
from condition_recommender.rules import load_condition_rule_set

rule_set = load_condition_rule_set()
```

Loading performs strict schema, taxonomy-vocabulary, environment, role,
substance-identity, template-reference, status, and provenance validation.

### Rule path limitations and remaining work

The following work is still required before the C-N drafts can be activated:

- verify every draft against a located primary procedure;
- complete an independent C-N applicability and recipe benchmark;
- use the generated rulebook in blind chemist review and track disagreements;
- calibrate which coordination and functional-group risks should be hard
  exclusions rather than cautions;
- add separately reviewed regimes for currently unsupported handles,
  intramolecular C-N formation, and additional deactivated nitrogen classes.

## Generic structure-backed recommendation

This is the intended type-agnostic precedent path. It operates on canonical
converted records rather than reaction-name partitions.

### Build reference-safe pilot samples

Develop and calibrate against deterministic samples before converting the full
corpus:

```powershell
python -m condition_recommender.sample_cli `
  data-processor/reaction_dataset `
  results/generic_sampling/v1
```

The sampler writes smoke, development, validation, and untouched-test CSVs plus
a versioned manifest. Rows connected by the same normalized publication
reference or exact normalized reaction text remain in one primary partition.
The smoke sample is a subset of development. Original source file and row
provenance survive conversion of the sampled CSVs.

Convert the smoke sample while developing:

```powershell
python -m condition_recommender.generic_conversion_cli `
  results/generic_sampling/v1/smoke.csv `
  results/generic_conversion/smoke_v1
```

### Audit and convert a structure-rich dataset

Audit source data without modifying it:

```powershell
python -m condition_recommender.audit_cli `
  data-processor/reaction_dataset `
  results/reaction_dataset_audit `
  --chemistry-sample-per-file 100
```

Convert a CSV file or directory:

```powershell
python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset `
  results/generic_conversion
```

The converter writes canonical nested `records.jsonl`, tiered CSV review views,
and JSON/Markdown coverage reports. Exact reconstructed signatures and valid
mapped signatures may be admitted even when `named_family` is absent. Source
family labels are provenance only. Chemistry, condition, outcome, and index
eligibility are assessed independently. The legacy verified/review/rejected
tier remains a compatibility summary, not the indexing contract.

### Build an index and recommend

```powershell
python -m condition_recommender.generic_index_cli `
  results/generic_conversion/records.jsonl `
  results/generic_conversion/generic_index.json

python -m condition_recommender.generic_recommend_cli `
  "<reaction_smiles>" `
  --records results/generic_conversion/generic_index.json
```

Retrieval follows a chemistry hierarchy: exact signature, relaxed handle
signature, high-confidence family, generic transformation, then compatible bond
edits narrowed by local-environment neighbors, followed by a broad compatible
bond-edit fallback. Hard bond-edit and recipe-compatibility gates run before
similarity and aggregation. Minimum support is counted by independent
publication, or by canonical reaction when no reference is available; repeated
scope examples from one paper cannot satisfy the threshold by themselves.

Results aggregate by canonical `RCORE1` recipe core and expose the observed
`RCR1` operating-condition variants. Support distinguishes raw observations,
canonical reactions, reference-local condition series, independent references,
and datasets. Repeated rows from one publication contribute one evidence unit
to similarity and yield aggregation. Expected yield is omitted when no usable
outcomes support the recipe; it is never invented from missing values. In that
case the unavailable yield term receives zero applied weight and the available
ranking terms are renormalized.

Generic configuration is split by responsibility:

- `definitions/generic_retrieval.v1.json` selects the narrowest adequately
  supported candidate pool;
- `definitions/generic_similarity.v1.json` weights interpretable molecular
  comparison features;
- `definitions/generic_ranking.v1.json` combines similarity, compatibility,
  condition certainty, outcome evidence, and reference-aware support.

Ranking and retrieval parameters are selected on development splits and
promoted only after a separate validation gate. Similarity weights remain an
explicit chemistry prior rather than a statistical fit.

Every result includes a typed retrieval trace. Each attempted tier reports raw
rows, independent support units, compatibility-filtered rows, exclusions, the
configured minimum support, and why the tier was selected or skipped.
Each recommendation also includes a typed score trace with similarity feature
values and contributions, ranking components and contributions, applied
weights, definition versions, and evidence/outcome counts.

For an application or batch process, load the validated index once:

```python
from condition_recommender import GenericConditionRecommender

recommender = GenericConditionRecommender.from_path(
    "results/generic_conversion/generic_index.json"
)
first = recommender.recommend("<reaction_smiles>")
second = recommender.recommend("<another_reaction_smiles>")
```

`recommend_generic_conditions()` remains the convenient one-shot entry point.

Reaction topology participates in the more specific signature tiers. A
topology-agnostic fallback is allowed only with an explicit
`REACTION_TOPOLOGY_FALLBACK_USED` warning. Multi-event queries and precedents
are compared using their complete event and net-edit sets.

### Held-out evaluation

```powershell
python -m condition_recommender.evaluation_cli `
  results/generic_conversion/generic_index.json `
  results/generic_conversion_evaluation `
  --test-fraction 0.2 --seed 17 --top-k 5
```

Rows connected by a normalized publication reference or canonical reaction
remain wholly in train or test. The report includes coverage, fallback levels,
top-1/top-k recipe recovery, conditional recovery when the recipe was observed
in training, yield MAE, compatibility exclusions, and the count of
hard-incompatible recommendations.

Use `--split-mode scaffold_disjoint`, `source_disjoint`, or `forward_time` for
the stricter leakage diagnostics. Compare chemistry-gated baselines with:

```powershell
python -m condition_recommender.baseline_cli <index> <output>
```

Calibrate only from development and validation indices:

```powershell
python -m condition_recommender.calibration_cli `
  <development-index> <validation-index> <output>
```

Generate the blind review gate before opening the untouched test:

```powershell
python -m condition_recommender.chemist_review_cli <index> <output>
```

An independent chemist reviews the self-contained `review_packet.html` and
records decisions in `review_form.csv` without opening `answer_key.jsonl`.
The JSONL packet and highlighted SVGs remain available for programmatic or
individual-case inspection. Every candidate decision must be one of
`compatible`, `compatible_with_caution`, `incompatible`, or `uncertain`. After
the form is complete, unblind and bind the signed review to the exact packet
artifacts:

```powershell
python -m condition_recommender.chemist_review_adjudication_cli `
  <review-output> <review-output>/review_summary.json `
  --reviewer "<name>" --independent-reviewer --sign-off
```

Use one `--unresolved-defect "<description>"` argument per unresolved
systematic issue and omit `--sign-off` until those issues have been addressed.
Adjudication refuses to read the answer key while any review decision is blank
or invalid.

### Restartable full conversion

After chemist review and untouched-test gates pass:

```powershell
python -m condition_recommender.sharded_conversion_cli `
  data-processor/reaction_dataset results/generic_conversion/v2/full `
  --shard-size 100 --workers 6

python -m condition_recommender.conversion_integrity_cli `
  results/generic_conversion/v2/full/shard_manifest.json

python -m condition_recommender.generic_index_cli `
  results/generic_conversion/v2/full/records.jsonl.gz `
  results/generic_conversion/v2/full/generic_index.json

python -m condition_recommender.generic_index_integrity_cli `
  results/generic_conversion/v2/full/generic_index.json
```

Canonical shards and catalogs are deterministic compressed nested JSONL. The
manifest binds each shard to source and output checksums plus every
participating schema and chemistry/registry definition.

Compose the final machine and independent-human gates with
`condition_recommender.release_validation_cli`. A review summary is accepted
only when it has the adjudication schema, independent reviewer sign-off, no
unresolved systematic defect, and hashes for the packet, form, report, and
answer key.

### Generic path limitations

- Performance depends on the quality and coverage of the supplied converted
  records.
- Current evaluation datasets are pilots and do not establish broad production
  validity.
- Registry coverage and role normalization are still incomplete for some
  source datasets.
- Calibration is sample-level and is not broad production-accuracy evidence.
- Production release remains gated by independent chemist review and the
  untouched-test report.

## Weak-label recommendation

The weak-label path uses a cleaned flat CSV whose precedent rows contain source
reaction labels and reactive-site labels rather than reaction structures.

```powershell
python -m condition_recommender.recommend_cli `
  "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
```

Options:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" --top-k 10

python -m condition_recommender.recommend_cli "<reaction_smiles>" `
  --records datasets/reaction_label/v2.1_cleaned.csv
```

The query reaction is structure-verified with `reactive_taxonomy`, but precedent
rows are not. The implementation:

1. assigns a query reaction grammar;
2. restricts source reaction types through a declarative crosswalk;
3. matches two reactive participants as an unordered pair;
4. ranks by functional-group signature and structural qualifiers;
5. uses yield, z-score, and support as lower-weight signals;
6. aggregates canonical `condition_registry` recipes.

Weights live in `definitions/label_retrieval.v1.json`.

Weak handle signatures may be family-level prefixes when the source label lacks
structural detail. Known incompatible handle families or `XH` center classes
are hard mismatches; an unknown/blank participant may still contribute a
partial match through the other participant. This allows `XH|Csp3` precedents
to support an exact activated-carbon query without allowing N–H or O–H rows to
leak into the candidate recipes.

For the cleaned HTE data, `CH-Activation` is routed to the intermolecular
`Ar–X + Ar–H → Ar–Ar` grammar. The same source reaction type also contains
`Ar–H + alkyne` records, so the generic `alkyne` label is normalized to
`PI|Alkyne`; those records then fail participant compatibility instead of
entering direct-arylation recommendations. Friedel–Crafts remains a separate
structure grammar and is not backed by this CSV because the dataset has no
corresponding acylation source type.

Source reaction context disambiguates the same `alkyne` label in Sonogashira
rows as the terminal reactive handle `XH|Csp|H1|Alkynyl`. Two-partner
retrieval requires both precedent signatures, and transfer handles from
different element families (for example Zn, Mg, Sn, and B) are hard
incompatible even when the generic transformation grammar is shared.

Important limitations:

- precedent reactions are not structure-verified;
- source reactive-site order is unreliable and treated as unordered;
- unresolved condition names remain explicit with identity uncertainty;
- reaction topology cannot be recovered from the precedent rows;
- only structurally verified intermolecular queries are accepted;
- recommendations are weak-label precedents, not verified literature matches.

The path returns explicit warnings including
`WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED` and
`UNRESOLVED_CONDITION_IDENTITIES_RETAINED_WITH_UNCERTAINTY`.

The cleaned CSV uses a component-oriented condition contract instead of fixed
`catalyst`/`base`/`solvent` slots:

- `condition_recipe_id` is the canonical `RCR1` aggregation key;
- `condition_display` is a compact human-review label;
- `temperature_c`, `time_h`, `concentration_m`, and `atmosphere` are convenient
  single-stage review fields;
- `procedure_text` preserves source text and never contributes opaque text to
  recipe identity.

The paired `<csv-stem>.condition_recipes.jsonl.gz` catalog stores each unique
full recipe exactly once, including resolved identities, contextual roles,
provenance, stages, declared absences, and uncertainty. It is losslessly
compressed and read transparently by the recommender. This keeps the CSV small
without making the recommendation path parse display strings.

Multi-stage procedures remain ordered in the recipe catalog; their
top-level temperature/time cells stay blank rather than collapsing multiple
stages. Source columns are role hints only and cannot override contradictory
registry chemistry.

Regenerate the cleaned CSV with:

```powershell
python -m scripts.clean_reaction_label_csv `
  datasets/reaction_label/v2.1.csv `
  datasets/reaction_label/v2.1_cleaned.csv
```

## Shared compatibility behavior

`definitions/compatibility.v1.json` compares unchanged query spectator groups
with contextual recipe roles and condition families.

Hard conflicts currently include explicit oxidation, reduction, acid, moisture,
and atmosphere conflicts. Hydrolysis, acid/base, metal coordination, and
catalyst-poisoning risks generally remain eligible with auditable penalties and
cautions. Within one penalty category, only the strongest matching penalty is
applied. The definition loader validates query tags against
`reactive_taxonomy`, and recipe buckets and family IDs against
`condition_registry`, so a future taxonomy extension fails explicitly if the
matching compatibility vocabulary is stale.

Compatibility evidence does not repair an unsupported transformation or make a
draft recipe production-ready.

## Input contract

Recommendation APIs require a complete product-specified reaction SMILES:

```text
reactant1.reactant2>>product
```

Invalid, unsupported, ambiguous, or chemically conflicting inputs return a
typed error or abstention rather than an unrelated recommendation.

## Development status

The package is moving toward one canonical structure-backed recommendation
workflow. Until parity and evaluation gates are complete:

- expert C-N rules remain a separate explicit path;
- generic retrieval remains the preferred cross-family direction;
- weak-label retrieval remains transitional;
- application code must not merge results from these paths without preserving
  their evidence level and provenance.
