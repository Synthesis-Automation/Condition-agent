# Higher-Level Route Dataset Audit and POC

## Decision

Use the released compressed route JSONL as the source for route learning. Do
not use the two reaction CSV exports as the canonical route source. The route
JSONL preserves patent and route membership, the original mapped three-part
reaction SMILES, reagents, and the paper's higher-level abstractions. The CSV
writer deliberately converts reactions to `reactants>>products`, losing the
middle reagent field and all route context.

The raw Figshare snapshot remains unchanged. A deterministic curator in
`core_retrosynthesis.route_curation` produces a separate, normalized POC
artifact.

Sources:

- paper: <https://doi.org/10.1021/acscentsci.5c02014>
- author repository: <https://github.com/jihye-roh/higherlev_retro>
- released data: <https://doi.org/10.6084/m9.figshare.28306673>

## What the release actually contains

The local Figshare version 2 archive contains:

| File | Records | Intended use |
| --- | ---: | --- |
| `routes/uspto.higher-level.routes.jsonl.gz` | 976,597 | Route membership, original steps, reagents, and higher-level subtrees |
| `reactions/uspto_original.csv` | 1,414,143 | Deduplicated source reaction IDs for single-step training |
| `reactions/uspto_higher-level.csv` | 1,549,008 | Algorithmically abstracted reactions for single-step training |
| `route_testset/uspto190_canon_reactions.smi` | 640 | Reactions associated with the USPTO-190 route benchmark |

All 976,597 route records have unique route IDs and valid top-level count
relationships. They cover 103,332 patent IDs. There are 1,761,206 source-step
occurrences across original trees but only 1,414,143 unique source reaction
IDs, so route occurrences overlap materially.

Most records are not useful multistep examples:

| Original steps | Routes |
| ---: | ---: |
| 1 | 613,452 |
| 2 | 176,183 |
| 3 | 84,869 |
| 4 | 44,065 |
| 5 | 25,109 |
| 6 | 14,257 |
| More than 6 | 18,662 |

Only 363,145 records are multistep. The large number of single-step records is
expected from the authors' extraction script, which explicitly adds reactions
not included in an extracted tree as one-step routes.

## Important data semantics

### Routes are inferred, not reported chronology

The source pipeline groups reactions by patent and connects products and
reactants by molecular identity to infer trees. It does not provide a manually
verified experimental chronology or a direct record of a chemist's planning
intent. A connected tree is therefore route evidence, not proof that the
patent executed exactly that sequence.

### Reaction array order is not route order

The author's abstraction code collects subtree reactions from a Python `set`
before serializing them. Consequently, `subtrees[].reactions[]` has no reliable
synthesis or retrosynthesis order. Any consumer that treats array position as
time order will learn incorrect sequences.

The local curator ignores source array order and reconstructs the chain from
canonical atom-map-free product/precursor identity. It writes explicit
`retrosynthetic_position`, internal intermediate, and terminal precursor
fields.

### Higher-level labels are weak supervision

`abstracted_reaction_smiles` is produced by the paper's abstraction algorithm.
It is not a human annotation saying that a chemist intentionally selected a
protecting group, tandem operation, or strategic disconnection. These fields
can be useful as auxiliary weak labels, but must not be presented as observed
chemist intent or treated as unquestionable ground truth.

### One-step CSVs have heavy exact duplication

The original CSV has 815,456 unique exact mapped reaction strings among
1,414,143 rows; 598,687 rows repeat an exact reaction string. The higher-level
CSV has 871,732 unique exact strings among 1,549,008 rows; 677,276 rows repeat
an exact string. IDs remain unique because they represent source occurrences.

For a single-step model, exact chemistry should be deduplicated or weighted so
patent enumeration does not dominate. For a route model, occurrences may be
retained as evidence, but training/evaluation grouping must remain at patent or
near-duplicate route-cluster level.

### Mapping is useful but not uniformly source-complete

Among one deterministically selected 3-6-step reduced candidate per eligible
patent, 2,907 of 18,892 routes failed the strict requirement that every mapped
product atom have a mapped reactant source. These reactions are allowed by the
authors' preprocessing when a small number of product atoms have no source,
but they are unsafe for our graph-edit and signature contracts. The POC rejects
them rather than inventing correspondence.

## POC curation contract

The first corpus deliberately favors clean implementation and strong evidence:

1. Read only `uspto.higher-level.routes.jsonl.gz`.
2. Require 3-6 original steps.
3. Require `depth == reaction_count` and one subtree, giving a linear first
   contract while the general branched schema remains separate work.
4. Require at least one step reduction in the paper's higher-level route. This
   ensures the POC contains an actual weak strategic signal.
5. Parse reactants, reagents, products, and non-empty abstractions with RDKit.
6. Require complete and unique atom maps and require product maps to be a
   subset of reactant maps.
7. Reconstruct a connected chain with exactly one target and no repeated
   intermediate.
8. Exclude exact USPTO-190 reactions and routes whose target occurs as a
   product in the released USPTO-190 reaction list.
9. Choose at most one route per patent using a stable content-independent hash.
10. Assign train, validation, and test splits by patent hash.

The step/depth/linearity restrictions are POC sampling choices, not chemistry
quality claims. Branched, longer, unreduced, and multi-subtree records remain in
the raw snapshot for later phases.

## Generated POC

The generated artifact is:

```text
datasets/external/higher_level_retrosynthesis/figshare_v2/curated/
  routes.poc.v1.jsonl.gz
  routes.poc.v1.jsonl.gz.manifest.json
```

It contains 5,000 routes from 5,000 distinct patents:

| Split | Routes |
| --- | ---: |
| Train | 4,003 |
| Validation | 525 |
| Test | 472 |

| Original steps | Routes |
| ---: | ---: |
| 3 | 2,636 |
| 4 | 1,389 |
| 5 | 667 |
| 6 | 308 |

The output SHA-256 is
`0bb7da925f0b261f3164cb95aa61afcb78be47f81e4a4016ff86bae12c099c97`.
The manifest records source/test hashes, RDKit version, policy version, every
rejection count, split counts, and abstraction-reduction distribution.

Each route record contains source identity, patent split, target identity,
route-level counts, validation status, and ordered steps. Each step retains:

- the original mapped `reactants>reagents>product` string;
- the algorithm-generated abstracted reaction, when present;
- separated reactant, reagent, and mapped-product strings;
- canonical product and precursor identities;
- explicit internal and terminal precursors; and
- retrosynthetic position.

A reproducible 50-route chemist review is generated at
`results/core_retrosynthesis/route_reviews/higher_level_routes_random50_chemist_sequence.html`.
It uses simple random sampling with seed `20260814`. Each route starts with a
forward synthetic sequence of structures, added reactants, recorded conditions,
intermediates, and final product; the toolbar can reverse it into retrosynthetic
order. Original and non-empty higher-level reactions remain below each sequence.
The report also supports locally persisted review status, notes, filtering, and
JSON export.

The curated records have also been converted to the shared nested route-tree
schema. The 50-route tree artifact matches the HTML sample exactly, and all
5,000 POC routes converted with zero rejections. Design, schema, invariants,
commands, and artifact hashes are documented in
`docs/new/route_tree_contract_and_conversion.md`.

All route steps have additionally been analyzed through the generic reaction
signature and reaction-core pipeline. The resulting route-core sidecars preserve
cross-step atom-lineage ambiguity and define inspectable two-/three-event motifs.
The minimized 50-route chemistry review is generated at
`results/core_retrosynthesis/route_reviews/higher_level_routes_random50_route_core.html`.
Design, audit results, commands, and hashes are documented in
`docs/new/route_core_projection_design_and_status.md`.

## How to use it first

The highest-ROI first experiment is route-action ranking, not a new end-to-end
route generator:

1. At each observed target or intermediate, ask the existing validated
   single-step engine for candidate disconnections.
2. Mark the recorded next step as the positive only when its graph edit and
   precursors can be reconstructed by our contracts.
3. Use other structurally valid candidates at the same state as hard negatives.
4. Train or tune a route-context reranker using structure-derived state,
   disconnection, complexity, compatibility, and terminal-material evidence.
5. Use the higher-level reduction flag only as an auxiliary weak target.
6. Evaluate by patent-disjoint split and report both next-action recovery and
   whole-route outcomes.

This leaves candidate generation and chemistry validation deterministic while
using data to learn which valid action is strategically preferable.

This first experiment is now implemented as the versioned route-action replay
benchmark. On the random 50-route review sample, 18 of 188 observed steps are
strictly eligible positives; all 18 receive validated candidates. Exact precursor
recall is 33.3% at top 1 and 66.7% at top 10/25, while SITE1 recall is 44.4% at
top 1 and 72.2% at top 10/25. The low 9.6% positive eligibility is the immediate
data/chemistry bottleneck and must be improved without weakening validation.
Design and audit details live in
`docs/new/route_action_replay_benchmark_design_and_status.md`.

## Required cleanup after the POC

Before production use:

- add a general occurrence-preserving branched route graph importer;
- reconcile every original step with `ReactionSignature` and retain mapping
  conflicts for review rather than silently dropping or remapping them;
- normalize middle-field conditions through `condition_registry` while
  retaining the raw reagent string and provenance;
- cluster patent families and near-duplicate molecules/routes before final
  benchmarking, because patent-level splitting alone cannot prevent chemistry
  analog leakage across related patents;
- quantify route-extraction confidence and identify disconnected enumeration,
  alternative preparations, and duplicate intermediates within patents;
- deduplicate or reweight exact single-step chemistry separately from route
  occurrence evidence;
- preserve tandem, multi-site, and one-pot reactions as atomic observed steps
  unless an external source provides a defensible internal sequence; and
- obtain additional human-curated route or total-synthesis data if the target
  is chemist intent, protecting-group strategy, or route selection among known
  alternatives. This release alone does not directly annotate those decisions.
