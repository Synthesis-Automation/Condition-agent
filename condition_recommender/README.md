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
| Suzuki-specific structural retrieval | `recommend_conditions()` | Narrow earlier pilot; not the intended long-term general path | Existing verified Suzuki pilot datasets |
| Weak-label retrieval | `recommend_conditions_from_labels()` | Functional transitional path with important data limitations | Datasets that contain reaction labels but not precedent structures |

Supporting audit, conversion, indexing, and evaluation commands prepare data for
these paths; they do not recommend conditions by themselves.

## Which path should I use?

- Use expert rules when the reaction is inside a reviewed rule scope and an
  explicit starter protocol is desired.
- Use generic structure-backed retrieval when a compatible converted record set
  or persisted generic index is available.
- Use the Suzuki-specific path only for the existing Suzuki pilot.
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

The remaining C-N defaults and structural overrides contain explicit screening
recipes derived from the legacy Buchwald-Hartwig rulebook plus expert
refinement. Their identities, loadings, equivalents, partner stoichiometry,
temperature, time, concentration, and atmosphere are complete, but they remain
drafts until checked against located primary procedures.

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

- `rule_kind`: `default`, `structural_override`, or `fallback`;
- selected tier and priority;
- structural match evidence;
- failures for every nonmatching rule;
- explicit canonical recipe variants;
- compatibility scores, penalties, and hard exclusions;
- draft and fallback warnings;
- rule, taxonomy, template, and recipe definition versions.

Without `--include-draft`, structural draft matches remain visible in
`match_traces`, but their recipes are not returned.

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
- add blind chemist review artifacts and disagreement tracking;
- calibrate which coordination and functional-group risks should be hard
  exclusions rather than cautions;
- add separately reviewed regimes for currently unsupported handles,
  intramolecular C-N formation, and additional deactivated nitrogen classes.

## Generic structure-backed recommendation

This is the intended type-agnostic precedent path. It operates on canonical
converted records rather than reaction-name partitions.

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
family labels are provenance only.

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
edits and environments. Hard bond-edit and recipe-compatibility gates run before
similarity and aggregation.

Results aggregate by canonical `RCR1` recipe and report retrieval level,
support, dataset diversity, expected yield, precedent IDs, compatibility
evidence, explanations, and cautions.

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

Canonical reaction groups remain wholly in train or test. The report includes
coverage, fallback levels, top-1/top-k recipe recovery, conditional recovery
when the recipe was observed in training, yield MAE, compatibility exclusions,
and the count of hard-incompatible recommendations.

### Generic path limitations

- Performance depends on the quality and coverage of the supplied converted
  records.
- Current evaluation datasets are pilots and do not establish broad production
  validity.
- Registry coverage and role normalization are still incomplete for some
  source datasets.
- Retrieval calibration across all supported transformation classes is not
  finished.

## Suzuki-specific structural pilot

`recommend_conditions()` is an earlier, Suzuki-only structural precedent path:

```python
from condition_recommender import recommend_conditions

result = recommend_conditions(
    "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
    records_path="results/suzuki_recommendation_pilot/verified.csv",
)
```

It requires an exactly verified Suzuki reaction and a verified pilot CSV. It
does not support general reaction classes, and temperature and time are not
modeled. New cross-family work should target the generic structure-backed path
instead of expanding this specialized implementation.

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
6. aggregates identical raw condition recipes.

Weights live in `definitions/label_retrieval.v1.json`.

Important limitations:

- precedent reactions are not structure-verified;
- source reactive-site order is unreliable and treated as unordered;
- condition names are not yet normalized through `condition_registry`;
- reaction topology cannot be recovered from the precedent rows;
- only structurally verified intermolecular queries are accepted;
- recommendations are weak-label precedents, not verified literature matches.

The path returns explicit warnings including
`WEAK_LABEL_PRECEDENTS_NOT_STRUCTURE_VERIFIED` and
`CONDITION_NAMES_NOT_REGISTRY_NORMALIZED`.

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
applied.

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
- the Suzuki-specific pilot remains narrow;
- weak-label retrieval remains transitional;
- application code must not merge results from these paths without preserving
  their evidence level and provenance.
