# Route-Core Projection: Design and Implementation Status

**Status date:** 2026-08-14

## Decision

`RouteCoreProjection` is the canonical minimized chemistry view of one concrete
observed or planned route. It composes the generic single-step reaction core and
signature across the route tree, retains cross-step atom-lineage ambiguity, and
provides deterministic route and motif identities.

It is deliberately not an executable multistep template. The full
`ReactionRouteTree` remains source truth; route-core data is a derived projection
for chemistry review, evaluation, retrieval, motif mining, and learned search
priors.

Package ownership remains directional:

```text
reactive_taxonomy reaction observations and per-step cores
                         ↓
core_retrosynthesis route topology and route-core composition
```

## Projection schema

Schema `1.0` contains:

- source tree, route, patent, split, target, depth, and provenance;
- one `RouteCoreStep` per observed reaction occurrence;
- complete serialized `ReactionSignature` and `ReactionCoreProjection`;
- minimum and rendering reaction SMILES from the display-only projection;
- chemistry validity, evidence, quality, warnings, and definition versions;
- internal and terminal precursor occurrence links;
- observed remaining-step labels;
- explicit producer-to-consumer atom-lineage candidates;
- exact, typed, and shape route-core keys; and
- inspectable two- and three-event typed and shape motifs.

Exact route-core identity includes terminal molecular identity. Typed and shape
identity retain transformation topology while omitting exact leaf identity.
Multiple source routes may therefore share one content-derived route-core ID.

## Cross-step atom lineage

Atom-map numbers are local to each patent reaction and are never assumed to be
stable across steps. For every observed intermediate, the builder:

1. finds the corresponding producer product and consumer reactant components by
   canonical molecular identity;
2. performs a chirality-aware full-graph isomorphism;
3. translates producer-local atom maps to consumer-local atom maps;
4. stores every bounded distinct correspondence candidate; and
5. reports unique, symmetry-ambiguous, bounded-ambiguous, unresolved, or invalid
   evidence.

Symmetry ambiguity is valid connectivity evidence, not an error. The projection
does not select an arbitrary global atom lineage.

## Minimized chemistry display

Each step reuses the standalone reaction display projection. Active atoms,
reaction interfaces, connecting ring paths, and versioned local context remain;
remote unchanged structures are represented by labeled attachment points.

The display is a view only. It is not used as reaction identity, admission
evidence, or an executable reaction program. Original mapped reactions remain in
the route tree and step record.

## Full-corpus audit

The 5,000-route POC produced:

| Measure | Result |
| --- | ---: |
| Routes | 5,000 |
| Reaction occurrences | 18,647 |
| Signatures, cores, and minimized drawings | 18,644 |
| Passing core quality | 1,233 |
| Review core quality | 17,085 |
| Blocked core quality | 326 |
| Unavailable cores | 3 |
| Fully chemistry-resolved routes | 4,681 |
| Cross-step lineage links | 13,647 |
| Unique lineage | 3,674 |
| Symmetry-ambiguous lineage | 9,064 |
| Bounded-ambiguous lineage | 909 |
| Fully lineage-connected routes | 4,457 |
| Unique exact route-core identities | 4,523 |
| Unique typed route-core identities | 4,463 |
| Unique shape route-core identities | 4,423 |
| Unique typed 2/3-event motifs | 17,991 |
| Unique shape 2/3-event motifs | 16,566 |
| Route conversion rejections | 0 |

The three unavailable reactions contain atom-map element contradictions and
unaccounted product atoms. They remain in their route projections with explicit
warnings. Blocked steps are likewise retained; no route topology or observation
is erased by chemistry interpretation quality.

## Generated artifacts

```text
datasets/external/higher_level_retrosynthesis/figshare_v2/curated/
  routes.poc.random50.core.v1.jsonl.gz
  routes.poc.random50.core.v1.jsonl.gz.manifest.json
  routes.poc.core.v1.jsonl.gz
  routes.poc.core.v1.jsonl.gz.manifest.json

results/core_retrosynthesis/route_reviews/
  higher_level_routes_random50_route_core.html
```

Hashes:

- 50-route route-core JSONL:
  `2b62f37b9cf4eea7f5c89df5d2a8a8a818e75d9066bc149f96ae9936439264e5`
- 5,000-route route-core JSONL:
  `7396640f621c27666b8ac4dd065b3d51f0c6c7ddaa709224bdf646cd72ee7f51`

The 50-route HTML contains 188 minimized core-event drawings, 138 cross-step
lineage connectors, 50 full-target anchors, and no drawing-error placeholders.
It supports forward and retrosynthetic direction, chemistry and lineage filters,
review notes, and JSON review export.

## Commands

```powershell
python -m core_retrosynthesis build-route-cores `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.tree.v2.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.core.v1.jsonl.gz `
  --workers 8

python -m core_retrosynthesis render-route-core-review `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.core.v1.jsonl.gz `
  results/core_retrosynthesis/route_reviews/higher_level_routes_random50_route_core.html `
  --sample-size 50 --seed 20260814
```

## Deliberate next work

1. Review the blocked, unavailable, and bounded-lineage examples in the HTML and
   targeted audit exports.
2. Expand the implemented route-action replay POC after improving strict
   single-step eligibility coverage; see
   `route_action_replay_benchmark_design_and_status.md`.
3. Compute patent-disjoint empirical support for typed and shape motifs.
4. Add motif compatibility as optional ranking evidence; do not make it a hard
   routing rule.
5. Consider an executable route operator only after recurrent motifs survive
   held-out evaluation and every constituent one-step action remains independently
   applicable and forward-validatable.
