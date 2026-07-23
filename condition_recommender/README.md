# Reaction-Label Condition Recommender

## Declarative expert rules

The clean expert-rule definitions live under `definitions/rule_sets/`. They
match allowlisted structural facts from `reactive_taxonomy`, reference recipe
templates owned by `condition_registry`, and do not match source reaction names
or display labels. The first definition set covers palladium condition proposals
for verified `sp2_c_n_substitution` events.

```python
from condition_recommender.rules import load_condition_rule_set

rule_set = load_condition_rule_set()
```

Loading performs strict shape, taxonomy-vocabulary, role, environment,
substance-identity, and cross-package template validation. The Phase I C-N
rules cover bounded intermolecular Ar/HeteroAr substitution and use a
deterministic hierarchy:

1. combined hindered-Ar-Cl structural override;
2. Ar-Cl, hindered-Ar-Br, or alpha-branched-amine structural override;
3. primary/secondary alkyl amine, primary/secondary aryl amine, aromatic N-H,
   or primary-amide class default;
4. bounded free-amine Ar-Br fallback.

Within the first matching tier, only the highest-priority rule is selected.
A matching draft override blocks a lower-priority production default, avoiding
unsafe fallback around a known special case. Ar-F, uncurated leaving groups,
intramolecular reactions, and unsupported nitrogen classes abstain. The rule
engine proposes a palladium condition regime from structural evidence; it does
not claim that a generic `sp2_c_n_substitution` product proves a
Buchwald-Hartwig mechanism.

The narrow primary-amide/aryl-chloride rule remains the only active rule and
has a complete primary-literature protocol. The other defaults and overrides
are complete, explicitly materializable screening drafts derived from the
legacy Buchwald-Hartwig rulebook plus expert refinement. They require
primary-procedure review before activation.

The production API excludes drafts by default:

```python
from condition_recommender import recommend_rule_conditions

result = recommend_rule_conditions(
    "Clc1ccccc1.CC(N)=O>>CC(=O)Nc1ccccc1"
)
```

It returns the reviewed tBuBrettPhos Pd G3/K3PO4/t-BuOH protocol with catalyst
loading, equivalents, partner stoichiometry, temperature, time, concentration,
atmosphere, compatibility evidence, and source provenance. The active rule is
restricted to primary amides and aryl chlorides; it does not extrapolate to
bromides or secondary amides.

Structurally evaluate the draft rules during review with an explicit opt-in:

```powershell
python -m condition_recommender.rule_recommend_cli `
  "Brc1ccccc1.CN>>CNc1ccccc1" `
  --include-draft
```

The review result contains the selected specificity tier and rule kind
(`default`, `structural_override`, or `fallback`), structural match evidence,
failed predicates for every nonmatching rule, the referenced recipe template,
canonical recipes for its explicit variants, compatibility scores and evidence,
hard-excluded variants, and a draft warning. Without `--include-draft`,
structural matches are reported in the trace but draft recipes are excluded from
recommendations. Template alternatives that are not part of an explicit variant
are never materialized implicitly.

## Mixed dataset audit

Audit a CSV file or directory before conversion without modifying the source
data:

```powershell
python -m condition_recommender.audit_cli `
  data-processor/reaction_dataset `
  results/reaction_dataset_audit `
  --chemistry-sample-per-file 100
```

Metadata, canonical reaction identities, duplicate groups, yields, and condition
registry coverage are computed over every row. Full reaction featurization uses
a deterministic per-file sample so large corpora remain practical. The command
writes `audit_report.json` and `audit_report.md`.

## Generic mixed-dataset conversion

Convert a single CSV or an entire directory without selecting a named reaction
family:

```powershell
python -m condition_recommender.generic_conversion_cli `
  data-processor/reaction_dataset `
  results/generic_conversion
```

For a bounded smoke test, add `--max-rows 1000`. This limit reads source rows in
file order and is not a statistically representative sample. The converter
writes canonical nested records to `records.jsonl`, tiered CSV review views, and
JSON/Markdown coverage reports. Exact reconstructed signatures and valid mapped
signatures can be verified even when `named_family` is absent. Source family
labels are retained only as provenance.

## Type-agnostic recommendation

Recommend canonical resolved recipes directly from converted JSONL records:

```powershell
python -m condition_recommender.generic_recommend_cli `
  "c1ccc2[nH]cnc2c1.COc1ccc(B(O)O)cc1>>COc1ccc(-n2cnc3ccccc32)cc1" `
  --records results/generic_conversion_chan_lam_pilot/records.jsonl
```

Retrieval uses the first signature tier with adequate support: exact signature,
relaxed handle signature, high-confidence named family, generic transformation,
or compatible bond edits and environments. A hard bond-edit gate prevents
fallback to precedents with incompatible net chemistry. Results aggregate by
canonical `RCR1` recipe and report support, dataset diversity, expected yield,
precedent IDs, explanations, and condition-identity cautions.

Reaction topology is part of the L0-L2 signature tiers and of fallback
similarity. Intramolecular precedents are preferred for intramolecular queries;
ring-size proximity refines cyclization matches. L3 remains a deliberately
topology-agnostic bond-edit fallback. If that fallback crosses reaction scope,
the result reports `REACTION_TOPOLOGY_FALLBACK_USED` and identifies the mismatch
in each affected recommendation.

Multi-event reactions are indexed by their complete deterministic event
multiset at L0-L2. The L3 hard gate continues to compare the complete net bond
edits, so a precedent covering only one event cannot silently stand in for a
mixed transformation. Canonical JSONL records retain nested events and event
relations; CSV review output exposes `reaction_event_count` and
`reaction_event_scope`.

Before ranking, `compatibility.v1.json` compares unchanged spectator-group tags
from `reactive_taxonomy` with contextual role buckets and curated family IDs in
the resolved recipe. Explicit oxidation, reduction, acid, and moisture conflicts
are excluded. Hydrolysis, acid/base, metal-coordination, and catalyst-poisoning
risks remain eligible with auditable penalties and cautions. Overlapping risks
within one category use only the strongest penalty. The same definition checks
oxidative or hydrogen atmospheres, strengthens hydrolysis penalties for hot
aqueous recipes, and validates mandatory catalyst roles for named metal-coupling
regimes. Unresolved condition identities turn a missing-role decision into a
penalty rather than an unsupported hard rejection.

## Persisted index and held-out evaluation

Build a deterministic, integrity-checked index once and reuse it for queries:

```powershell
python -m condition_recommender.generic_index_cli `
  results/generic_conversion_chan_lam_pilot/records.jsonl `
  results/generic_conversion_chan_lam_pilot/generic_index.json
```

The index stores only admitted signature/recipe precedents and has a content-based
`GRI1` identity. Both `records.jsonl` and the persisted index are accepted by the
generic recommendation CLI.

Generic index schema 1.1 records the reaction-signature schema, taxonomy
definition versions, recommendation-record schema, and converter version. Stale
or mixed converted records are rejected instead of silently falling through to
less-specific keys. Re-run generic conversion and rebuild the index after any
signature schema or identity-definition change.

Run a canonical-reaction-group held-out benchmark:

```powershell
python -m condition_recommender.evaluation_cli `
  results/generic_conversion_chan_lam_pilot/generic_index.json `
  results/generic_conversion_chan_lam_evaluation `
  --test-fraction 0.2 --seed 17 --top-k 5
```

Every `CRX1` group is assigned wholly to train or test, preventing duplicate
observations and alternate recipes for the same reaction from leaking across the
boundary. The report includes retrieval coverage, fallback levels, top-1/top-k
recipe recovery, recipe-seen conditional recovery, yield MAE, compatibility
exclusions, and a zero-tolerance count of hard-incompatible recommendations.

This CLI recommends reaction conditions from the cleaned HTE label dataset.
It uses `reactive_taxonomy` to identify the reaction grammar and reactive
functional-group signatures, then returns the five highest-ranked condition
recipes from compatible precedents.

## Basic usage

Run the command from the repository root:

```powershell
python -m condition_recommender.recommend_cli "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
```

The input must be a complete reaction SMILES containing reactants and a
product:

```text
reactant1.reactant2>>product
```

The CLI returns JSON containing the detected reaction grammar, query FG
signatures, candidate count, and ranked condition recommendations.

## Options

Return a different number of recommendations:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" --top-k 10
```

Use another cleaned dataset:

```powershell
python -m condition_recommender.recommend_cli "<reaction_smiles>" `
  --records datasets/reaction_label/v2.1_cleaned.csv
```

Show all options:

```powershell
python -m condition_recommender.recommend_cli --help
```

## Cleaned CSV schema

The recommender uses one flat CSV. Its main columns are ordered for convenient
filtering and manual review:

- `source_reaction_type`
- `reactive_site_1_normalized_label` and `reactive_site_2_normalized_label`
- `reactive_site_1_display_label` and `reactive_site_2_display_label`
- `yield_pct` and `z_score`
- catalyst, ligand, base, solvent, additive, and procedure columns

Additional `reactive_site_1_*` and `reactive_site_2_*` columns preserve the
original source labels, taxonomy signatures, qualifiers, and mapping status.
An unresolved label is retained verbatim and has `mapping_status=unresolved`.

Regenerate the cleaned CSV from the repository root with:

```powershell
python -m scripts.clean_reaction_label_csv `
  datasets/reaction_label/v2.1.csv `
  datasets/reaction_label/v2.1_cleaned.csv
```

## Ranking logic

The recommender:

1. Verifies the reaction and assigns a grammar with `reactive_taxonomy`.
2. Filters the dataset to compatible source reaction families.
3. Matches `reactive_site_1` and `reactive_site_2` as an unordered pair.
4. Ranks primarily by FG-signature similarity.
5. Uses attachment and branching qualifiers to refine close matches.
6. Uses yield, z-score, and precedent support as lower-weight signals.
7. Aggregates identical condition recipes and returns distinct results.

The ranking weights are configured in
`condition_recommender/definitions/label_retrieval.v1.json`.

## Supported reaction grammars

The cleaned dataset currently supports:

- Suzuki–Miyaura coupling
- C–N substitution and coupling
- Amide formation
- C–O substitution and coupling
- C–S substitution and coupling
- Sonogashira coupling
- Negishi coupling

Unsupported or structurally ambiguous reactions return an error instead of
falling back to an unrelated reaction family.

## Output

Each recommendation includes:

- rank and overall score;
- signature and qualifier similarity;
- expected yield and mean z-score;
- number of supporting precedents;
- source row numbers;
- base, catalyst, solvents, ligand, additives, and procedure text;
- a short explanation of the match.

## Important limitations

The source CSV contains reaction labels rather than reaction structures for
its precedent rows. Recommendations are therefore weak-label precedents, not
structure-verified literature matches. Reactive-site order in the source data
is not reliable and is treated as unordered. Condition names have also not yet
been normalized to condition-registry identities. Because reaction scope and
ring formation cannot be reconstructed from these rows, the weak-label
recommender accepts only structurally verified intermolecular queries and returns
`QUERY_TOPOLOGY_NOT_SUPPORTED_BY_LABEL_DATASET` for intramolecular, mixed,
unimolecular, or unresolved topology. Use the generic structure-backed path for
those reactions.
