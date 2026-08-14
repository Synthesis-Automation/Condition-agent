# Shared Route-Tree Contract and Observed-Route Conversion

## Decision

`core_retrosynthesis.route_contract` owns the canonical representation of one
concrete retrosynthesis route. The contract is evidence-neutral: observed
patent reactions and predicted search actions use the same alternating tree,
while typed evidence records how each reaction and route was obtained.

The contract represents one chosen route, not the complete search space. A
future set of alternative disconnections belongs in a separate AND/OR search
graph whose realized paths produce this tree contract.

## Schema

Schema version `2.0` is an occurrence-preserving, nested, bipartite tree:

```text
ReactionRouteTree
└── root: MoleculeOccurrenceNode
    └── reaction: RouteReactionNode
        ├── child: MoleculeOccurrenceNode
        │   └── reaction: RouteReactionNode | null
        └── child: MoleculeOccurrenceNode
            └── reaction: RouteReactionNode | null
```

A molecule occurrence is not molecular identity. Two nodes may have the same
canonical SMILES and still have different occurrence IDs because they occupy
different branches or stoichiometric positions.

### `ReactionRouteTree`

The tree records:

- stable content-derived `tree_id`;
- `route_kind`: `observed` or `planned`;
- canonical target and nested root occurrence;
- reaction count and maximum retrosynthetic depth;
- deterministic fingerprint tokens;
- dataset, source route, patent, and split provenance;
- source record and schema versions;
- route-connectivity method and uncalibrated confidence;
- higher-level abstraction counts; and
- route-level warnings.

### `MoleculeOccurrenceNode`

Each occurrence records:

- stable `occurrence_id`;
- canonical molecular identity;
- depth from the root target;
- terminal state and evidence; and
- the optional reaction that disconnects this occurrence.

The former `molecule_node_id` remains as a temporary compatibility alias and
serialized compatibility field.

### `RouteReactionNode`

Each reaction occurrence records:

- stable reaction-node and step IDs;
- retrosynthetic depth;
- complete reaction SMILES;
- typed `RouteStepEvidence`;
- optional prediction-only `PlannedRouteAction`; and
- all precursor molecule occurrences.

The former `proposed_reaction_smiles`, `operator_id`, and
`disconnection_site_key` views remain available during the v2 migration.
Observed reactions leave prediction-only annotations empty.

### Evidence

`RouteStepEvidence.evidence_kind` is one of:

- `observed`: a source reaction record exists;
- `inferred`: a relationship is reconstructed rather than directly recorded;
- `predicted`: the route search proposed the action.

Evidence kind does not hide route-level inference. For the higher-level USPTO
release, reaction records are `observed`, while
`connectivity_method=reconstructed_from_canonical_molecule_identity` states
that their route connections were inferred within a patent. Confidence remains
`null` because the source does not provide calibrated probabilities.

Algorithm-generated higher-level reactions are stored as weak annotations with
`abstraction_status`; they are not presented as observed chemist intent.

## Structural invariants

`validate_route_tree()` verifies:

- exactly one root at depth zero;
- alternating molecule and reaction layers;
- reaction depth equals parent molecule depth plus one;
- precursor molecule depth equals reaction depth;
- unique molecule-occurrence and reaction-node IDs;
- no terminal molecule has a child reaction;
- every reaction has at least one precursor;
- declared reaction count and maximum depth match traversal; and
- fingerprint tokens are deterministically ordered.

`ReactionRouteTree.from_dict()` reconstructs the immutable typed tree and runs
the validator, so persisted rows cannot silently bypass the public contract.

Chemistry validation remains separate and earlier in the pipeline. The route
curator validates mapping, parsing, molecular identity, and connectivity before
the tree converter is called.

## Observed conversion

`core_retrosynthesis.route_conversion` converts curated rows without trusting
their step-array order:

1. index each reaction by canonical product identity;
2. start at the declared target;
3. recursively connect canonical precursor identities to recorded products;
4. create a unique molecule occurrence for every precursor position;
5. retain the original mapped reaction and reagent fields;
6. attach explicit observed/inferred/algorithmic evidence;
7. reject cycles, reused reaction occurrences, unreachable steps, duplicate
   products, and count conflicts;
8. validate the constructed tree; and
9. serialize deterministic gzip JSONL with a checksum manifest.

The converter supports convergent branches now. The current 5,000-route POC is
linear because linearity was a curation policy, not a schema limitation.
One-pot and tandem source reactions remain atomic reaction nodes unless a
future source provides defensible internal steps.

## Generated artifacts

```text
datasets/external/higher_level_retrosynthesis/figshare_v2/curated/
  routes.poc.random50.tree.v2.jsonl.gz
  routes.poc.random50.tree.v2.jsonl.gz.manifest.json
  routes.poc.tree.v2.jsonl.gz
  routes.poc.tree.v2.jsonl.gz.manifest.json
```

The 50-tree artifact uses seed `20260814` and matches the HTML review sample.
Its SHA-256 is
`3934ccb87eb6aa6a46a2c6abe0431d6cdec760d23d62e324674f3e14b9f128bc`.

The full artifact contains:

- 5,000 routes from 5,000 patents;
- 18,647 reaction occurrences;
- 35,668 molecule occurrences;
- 5,000 unique tree IDs;
- zero conversion rejections; and
- zero structural-validation failures.

Its SHA-256 is
`5e9d45a45f58678cadb6f4b37cae0a89ab4817191af1fc10790da505810479a8`.

## Commands

Convert the full curated corpus:

```powershell
python -m core_retrosynthesis convert-route-trees `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.v1.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.tree.v2.jsonl.gz
```

Generate the sample matching the visual review:

```powershell
python -m core_retrosynthesis convert-route-trees `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.v1.jsonl.gz `
  datasets/external/higher_level_retrosynthesis/figshare_v2/curated/routes.poc.random50.tree.v2.jsonl.gz `
  --sample-size 50 --seed 20260814
```

Load validated typed trees:

```python
from core_retrosynthesis import iter_route_trees

for tree in iter_route_trees("routes.poc.tree.v2.jsonl.gz"):
    assert tree.route_kind == "observed"
```

## Deliberate next work

The contract is ready for route-learning experiments, but two chemistry layers
remain intentionally empty:

- `reaction_signature_id` will be populated only after each step reconciles
  mapped edits with `reactive_taxonomy`; and
- `condition_recipe_id` will be populated only after the reagent field is
  normalized through `condition_registry`.

Neither field is inferred from names or raw strings during route conversion.
