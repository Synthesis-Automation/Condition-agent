# Lightweight ML Design for Reaction Mechanism Hints

## Purpose

Develop an optional, lightweight machine-learning feature that estimates broad
mechanistic regimes from graph-derived reaction evidence. The feature is intended
to provide useful hints such as:

- polar or ionic;
- radical;
- pericyclic-compatible;
- transition-metal-mediated;
- photochemical or excited-state; and
- mixed or unresolved.

The feature must not claim to determine the experimentally established mechanism.
Many overall transformations are compatible with several pathways, and many
published reactions do not have a uniquely established mechanism. The output is
therefore a versioned, evidence-bearing, uncertainty-preserving interpretation.

Mechanism hints are not part of reaction identity. They must not change
`ReactionSignature`, its L0-L4 keys, `signature_id`, atom correspondence, observed
bond edits, admission, or chemistry compatibility decisions.

## Design principles

1. Molecular graph observations remain the source of truth.
2. The model predicts broad, non-exclusive mechanism hints rather than a detailed
   sequence of elementary steps.
3. Missing labels mean unknown, not negative.
4. Model inputs must use chemistry-bearing fields, not hashes, source labels, or
   reaction names.
5. Structure-only and condition-aware inference are separate tasks.
6. Model predictions cannot override contradictory structural evidence.
7. Low-support and out-of-domain cases must abstain or remain ambiguous.
8. Training and inference artifacts must be deterministic, versioned, and
   auditable.
9. Initial deployment is shadow-only and must not alter recommendation behavior.

## Why a hint is feasible, but an actual mechanism is not

Reaction SMILES and reaction graphs contain substantial information about reaction
classes and transformation patterns. Transformer reaction fingerprints can learn
reaction-class structure from reaction SMILES, while condensed reaction graphs and
atom-aligned graph models make bond changes and reaction centers explicit.

These representations nevertheless describe overall reactants and products, not
transition states, kinetic isotope effects, intermediate lifetimes, crossover
behavior, or complete potential-energy surfaces. Distinct polar, radical,
concerted, and metal-mediated pathways can produce the same net graph edits.

Mechanism-specific datasets address this limitation by supplying extra labels:

- PMechDB stores curated elementary polar steps with atom mapping, source/sink
  annotations, and arrow-pushing information.
- RMechDB separately curates elementary radical steps.
- MechFinder combines extracted reaction templates with expert-authored
  mechanistic templates. Its published mech-USPTO-31K dataset includes polar and
  pericyclic pathways but excludes radical and organometallic reactions.
- Kinetic mechanism classification uses time-dependent experimental measurements,
  demonstrating that mechanistic information can exist outside endpoint
  structures.

The appropriate goal for this system is therefore a broad mechanistic
compatibility assessment with explicit uncertainty.

## Mechanism-hint ontology

The first version should use independent labels rather than a mutually exclusive
softmax class. A reaction can support several labels at the same time.

### Primary regime hints

| Hint | Meaning |
| --- | --- |
| `polar_ionic` | A closed-shell donor/acceptor or ionic pathway is structurally compatible. |
| `radical` | A single-electron or radical pathway is supported or compatible. |
| `pericyclic_compatible` | The connected edit topology is compatible with a concerted pericyclic transformation. |
| `transition_metal_mediated` | A transition-metal-mediated pathway is structurally or conditionally supported. |
| `photochemical_excited_state` | Photochemical or excited-state participation is supported. |
| `mixed` | More than one mechanistic regime has material independent support. |
| `unresolved` | Evidence is insufficient, conflicting, or out of domain. |

`mixed` and `unresolved` should normally be derived assessment states rather than
independently trained labels.

### Optional later refinements

Once the broad labels are validated, optional subordinate hints may include:

- polar substitution, addition, elimination, or rearrangement;
- radical addition, fragmentation, substitution, or rearrangement;
- cycloaddition-, electrocyclization-, or sigmatropic-compatible;
- palladium-, nickel-, copper-, iron-, or other metal-mediated; and
- photoredox, direct photochemical, or electrochemical activation.

These refinements must not be introduced until the parent label has adequate
coverage and calibration.

## Public output contract

A possible immutable contract is:

```python
@dataclass(frozen=True)
class MechanismHint:
    hint_id: str
    support_score: float
    status: Literal[
        "supported",
        "possible",
        "contradicted",
        "insufficient_evidence",
        "out_of_domain",
    ]
    evidence_strength: Literal["direct", "contextual", "weak", "none"]
    matched_edit_indices: Tuple[int, ...]
    matched_core_event_ids: Tuple[str, ...]
    supporting_features: Tuple[str, ...]
    opposing_features: Tuple[str, ...]
    missing_evidence: Tuple[str, ...]
    nearest_labeled_precedent_ids: Tuple[str, ...]


@dataclass(frozen=True)
class MechanismHintAssessment:
    evidence_scope: Literal["structure_only", "structure_and_conditions"]
    hints: Tuple[MechanismHint, ...]
    primary_hint_id: Optional[str]
    ambiguity: Literal["low", "moderate", "high", "unresolved"]
    out_of_domain: bool
    model_id: str
    model_version: str
    feature_definition_version: str
    training_snapshot_id: str
    warnings: Tuple[str, ...]
    schema_version: str = "1.0"
```

`support_score` is initially a model support value. It must not be presented as a
physical probability or as the probability that a proposed mechanism is true.
Probability language is allowed only after label-specific calibration has been
validated.

Scores do not need to sum to one because labels may overlap. `primary_hint_id`
should be populated only when a validation-derived threshold and separation margin
are satisfied.

## Input representations

### Required structure-only inputs

The model should consume a structured feature bundle assembled from:

```text
ReactionSignature
+ ReactionCoreProjection
+ molecular reactive-site environments
+ reaction pattern matches
```

A signature-only model should be retained as an ablation baseline, but the core and
local contexts are expected to contain important discriminating evidence.

### Excluded inputs

The following must not be used as model features:

- `signature_id` and L0-L4 digest values;
- `reaction_id`, reference IDs, or source row numbers;
- source-declared reaction type;
- source reaction name or rendered reaction label;
- `named_family` or family confidence;
- dataset name when it can reveal the target label;
- condition identities in the structure-only model; and
- product or mechanism labels derived from the evaluation split.

Digest values are stable identifiers but contain no transferable chemical
similarity. Named families and source labels create target leakage and violate the
graph-first evidence policy.

## Feature extraction

Feature extraction should be deterministic and generic. It should serialize
existing graph facts into sparse categorical/count features without encoding
mechanism-specific rules.

### Reaction edits

- edit type: formed, broken, order changed, or hydrogen changed;
- endpoint element pair;
- old and new bond order;
- endpoint aromaticity and hybridization;
- endpoint formal charge;
- intra- versus intercomponent edit;
- local-environment identities and normalized environment tokens;
- stereo change type; and
- edit evidence quality.

Example tokens:

```text
edit:formed:C-N:SINGLE
edit:broken:C-Br:SINGLE
edit:order_changed:C-C:DOUBLE>SINGLE
edit:hydrogen_change:C:+1
endpoint:C:aromatic:SP2
stereo:atom:inverted
```

### Reaction-core features

- exact edit tokens and participant tokens;
- active-atom and event counts;
- atom-state changes in formal charge, radical electrons, aromaticity,
  hybridization, isotope, and stereochemistry;
- core event membership;
- shortest paths between edit centers;
- event relations such as shared atom, same component, or independent;
- remote-subgraph class and continuity;
- attachment-port chemistry;
- substituent branching, cyclicity, and benzylic, allylic, or propargylic flags;
- mapping coverage, evidence status, and core quality.

Example tokens:

```text
state_change:formal_charge:N:-1>0
state_change:radical:C:0>1
core_event_count:1
event_relation:shared_active_atom
remote_class:heteroaryl
substituent:benzylic
core_quality:pass
```

### Reaction topology

- reaction scope;
- formed-bond scopes;
- tether distances;
- ring-count delta;
- formed-ring sizes;
- ring element sequence;
- ring bond-order sequence; and
- number and connectivity of changing bonds in each event.

Pericyclic-compatible classification may benefit from a future generic
`ReactionEditPathDescriptor` that represents continuous or cyclic arrangements of
bond changes. This descriptor must be an observation of edit topology, not a
pericyclic label or executable reaction template.

### Molecular environments

- reactive-center element, charge, radical count, hybridization, hydrogen count,
  aromaticity, and ring membership;
- context kind: aromatic, alkyl, alkenyl, alkynyl, acyl, sulfonyl, phosphoryl,
  heteroatom, or other;
- electronic activation axis, class, score, and contribution types;
- steric accessibility and approach-burden class;
- conjugation and lone-pair availability;
- leaving-group, transfer-carrier, redox, elimination, strain, and coordination
  modifiers;
- molecular motif IDs and tags; and
- distance of nearby motifs from active atoms and attachment ports.

### Numeric features

Useful numeric fields include:

- counts of each edit and state-change type;
- active-atom, event, component, and spectator counts;
- ring sizes and ring-count delta;
- path lengths between edits;
- electronic activation and steric scores;
- branching counts;
- mapping coverage and checked-edit fraction; and
- evidence and completeness confidence.

Continuous features must use training-set normalization parameters stored in the
model artifact. Missing numeric values require explicit missingness features rather
than silent zero substitution.

## Lightweight model

### Primary model

Train one sparse logistic-regression model per mechanism hint. This supports
partially labeled multi-label data because each classifier can be fitted only on
records explicitly labeled positive or negative for that hint.

Recommended starting configuration:

- scikit-learn `LogisticRegression` or `SGDClassifier(loss="log_loss")`;
- one binary classifier per hint;
- L1 or elastic-net regularization for sparse, inspectable coefficients;
- label-confidence sample weights;
- class balancing assessed per label;
- fixed random seed; and
- validation-derived decision and abstention thresholds.

A tree ensemble can be evaluated as a secondary baseline, but sparse linear models
are preferred initially because they are easier to inspect, export, reproduce, and
calibrate.

### Nearest-neighbour companion

Use similarity-weighted labeled neighbours as an explanation and out-of-domain
guard. The existing generic signature/core similarity can retrieve the nearest
labeled precedents without becoming a hard mechanism rule.

The neighbour layer should report:

- precedent ID;
- structural similarity;
- known mechanism labels and label provenance;
- shared edit/core/context features; and
- important mismatches.

Low maximum neighbour similarity or poor feature coverage should set
`out_of_domain=True` and cap or suppress the primary hint.

## Training-label schema

Training labels must distinguish positive, negative, and unknown states:

```json
{
  "reaction_id": "example-1",
  "labels": {
    "polar_ionic": {
      "status": "positive",
      "confidence": 0.9,
      "source": "curated_mechanism_dataset",
      "evidence": ["curated elementary polar pathway"]
    },
    "radical": {
      "status": "unknown",
      "confidence": 0.0,
      "source": "not_assessed",
      "evidence": []
    },
    "pericyclic_compatible": {
      "status": "negative",
      "confidence": 0.8,
      "source": "chemist_review",
      "evidence": ["reviewed non-concerted pathway"]
    }
  },
  "label_schema_version": "1.0"
}
```

An absent label must never be converted automatically into a negative. Each binary
classifier should exclude unknown records. Label confidence becomes the training
sample weight.

### Candidate label sources

1. Curated mechanistic datasets with explicit provenance.
2. Chemist-reviewed broad mechanism labels.
3. Published mechanism annotations linked to the exact reaction record.
4. Condition-derived weak labels for training only, with reduced confidence.
5. Existing named-family or source labels as weak evidence only, never as
   structure truth.

Condition-derived labels are acceptable for training a structure-only student
model if condition fields are excluded from its inputs. Their provenance and weak
status must be preserved.

### External-data cautions

- Elementary-step datasets and overall synthetic reactions are different domains.
- PMechDB polar positives do not automatically provide radical or pericyclic
  negatives.
- mech-USPTO-31K intentionally excludes radical and organometallic reactions;
  missing labels for those regimes are unknown.
- Reaction-class labels from patent datasets are not equivalent to mechanism
  labels.
- Licensing must be reviewed before copying or redistributing external data.

External datasets should seed training and evaluation, followed by calibration on
complete reactions representative of this project's intended use.

## Structure-only and condition-aware models

### Structure-only model

Inputs:

- reaction signature;
- reaction core;
- molecular contexts; and
- structure-derived reaction patterns.

Use:

- query-side mechanism hints;
- optional retrieval/ranking features after validation;
- review output; and
- out-of-domain detection.

This model must not consume recommended or observed condition identities.

### Condition-aware model

Inputs may additionally include:

- resolved catalyst and ligand substance/family IDs;
- oxidant, reductant, acid, base, additive, and solvent families;
- contextual role confidence;
- temperature and staged process information;
- atmosphere;
- declared absences; and
- future irradiation and electrochemical metadata.

Use:

- classifying historical precedents;
- explaining why a recorded or proposed recipe supports a mechanistic regime; and
- post-recommendation consistency assessment.

The condition-aware model must not be used to infer a query mechanism before
conditions are recommended. Doing so would leak the prediction target into the
query representation.

## Model artifact and inference

Train with scikit-learn, but do not use an opaque pickle as the canonical deployed
artifact. Export a validated JSON and optional numeric-array artifact containing:

```text
model ID and schema version
mechanism-label vocabulary
feature vocabulary
feature normalization parameters
per-label coefficients and intercepts
calibration parameters
decision and abstention thresholds
training snapshot checksum
training-source summaries
signature/core/descriptor definition versions
training library versions
validation metrics and evaluation split ID
```

Inference then requires only deterministic feature extraction, sparse dot
products, and sigmoid/calibration functions. The deployed runtime does not need to
load arbitrary executable objects.

Artifact validation must reject:

- unknown labels;
- duplicated or unsorted feature definitions;
- incompatible schema or definition versions;
- non-finite coefficients;
- missing thresholds or calibration metadata; and
- a feature vocabulary inconsistent with the coefficient dimensions.

## Explanation generation

For a sparse linear model, explanation can report the largest positive and negative
feature contributions:

```text
coefficient(feature) * observed_feature_value
```

An explanation should include:

- the top supporting features;
- the top opposing features;
- missing high-value feature groups;
- matched edit indices and core event IDs;
- nearest labeled precedents;
- model and definition versions;
- out-of-domain status; and
- a warning that the result is a mechanistic hint, not an established pathway.

Feature explanations must use chemistry-readable projections. Internal categorical
tokens may be included for audit but should not be the only explanation shown to a
chemist.

## Confidence, ambiguity, and abstention

### Calibration

Evaluate calibration independently for every label. Candidate approaches include
sigmoid calibration or isotonic calibration when sufficient validation examples
exist. Store calibration parameters in the exported artifact.

### Thresholds

Select thresholds from held-out validation data to meet a defined precision or
false-positive target. Do not choose a universal arbitrary threshold for all
labels.

### Ambiguity

Report more than one hint when multiple labels pass their thresholds. Derive
ambiguity from:

- number of supported labels;
- separation between the strongest scores;
- calibration uncertainty;
- conflicting high-contribution features;
- nearest-neighbour agreement; and
- structural evidence quality.

### Out-of-domain handling

Mark a prediction out of domain when one or more of the following applies:

- no sufficiently similar labeled core exists;
- too many extracted features are unseen by the model;
- reaction-core quality is blocked or under review;
- atom correspondence or edit hypotheses remain materially ambiguous;
- the reaction is a multi-event combination absent from training; or
- the feature/definition version is incompatible with the model.

Out-of-domain status should suppress `primary_hint_id` unless direct structural
evidence independently supports a limited statement such as an explicitly encoded
radical state change.

## Evaluation design

Random reaction-level splits are inadequate because closely related cores and
scaffolds can occur on both sides. Evaluation should include:

1. group splits by reaction-core typed key;
2. scaffold or remote-subgraph group splits;
3. source-dataset holdout;
4. unknown-family reactions;
5. mapped versus inferred-correspondence subsets;
6. single-event versus multi-event subsets; and
7. complete versus incomplete condition records for the conditioned model.

Report per label:

- positive and negative support counts;
- precision, recall, and F1;
- precision-recall AUC;
- Brier score and calibration error;
- precision at the deployed threshold;
- coverage at the abstention threshold;
- out-of-domain performance; and
- nearest-neighbour agreement.

Required chemistry review sets should include:

- polar and radical pathways with similar net bond changes;
- cycloaddition versus stepwise ring formation;
- electrocyclization and sigmatropic examples;
- explicitly encoded radicals and charge changes;
- transition-metal/radical mixed systems;
- tandem or multi-event reactions;
- disputed mechanisms with multiple positive labels;
- invalid or ambiguous atom mappings; and
- reactions whose source family contradicts structural evidence.

## Baselines and ablations

Evaluate at least these baselines:

1. class-prevalence baseline;
2. similarity-weighted nearest neighbours;
3. signature-only sparse logistic regression;
4. signature plus reaction core;
5. signature, core, and molecular contexts; and
6. condition-aware model for historical records only.

The incremental value of family names and source labels may be measured in a
separate diagnostic experiment, but those fields must remain excluded from the
deployed graph-first model.

## Package ownership

### `reactive_taxonomy`

Owns:

- structure-only mechanism-hint contracts;
- generic feature extraction from reaction observations, signatures, cores, and
  molecular contexts;
- model-artifact validation and deterministic structure-only inference; and
- attachment of optional hints after the observation and signature are complete.

It must not import `condition_registry`, `condition_recommender`, or legacy
`chemtools`.

### `condition_registry`

Owns:

- condition substance identity;
- contextual roles and role confidence;
- condition family identity; and
- explicit mechanistically relevant substance/family properties when curated.

It must not infer a reaction mechanism.

### `condition_recommender`

Owns:

- labeled training-record conversion;
- reconciliation of structure-only hints with resolved condition recipes;
- condition-aware mechanism assessment;
- optional use of validated hints in similarity or ranking; and
- precedent-based explanations.

Application layers only render or expose the resulting assessment.

## Serialization

Add `mechanism_hint_assessment` as an optional nested field in reaction analyses and
recommendation records. CSV and review exports may expose concise fields such as:

```text
mechanism_primary_hint
mechanism_supported_hints
mechanism_ambiguity
mechanism_out_of_domain
mechanism_model_version
mechanism_hint_json
```

The full nested JSON remains canonical. Do not add mechanism fields to the
reaction-signature identity payload.

## Rollout plan

### Stage 1: contracts and labeled data

- Define the broad multi-label ontology.
- Add typed label and assessment contracts.
- Define the training-record schema and provenance requirements.
- Import or reference licensed seed datasets without treating missing labels as
  negatives.
- Create a small chemist-reviewed target-domain validation set.

### Stage 2: deterministic feature extraction

- Build sparse features from signatures, cores, and contexts.
- Add feature-vocabulary validation and snapshot tests.
- Confirm reactant-order and serialization invariance.
- Explicitly test that excluded leakage fields do not enter the feature vector.

### Stage 3: lightweight baselines

- Implement nearest-neighbour label aggregation.
- Train one sparse logistic model per hint.
- Export coefficients and vocabulary to a validated artifact.
- Compare signature-only, core, and context ablations.

### Stage 4: shadow inference

- Attach structure-only mechanism hints after signature construction.
- Serialize hints into review records.
- Leave admission, retrieval, compatibility, and recommendation behavior
  unchanged.
- Perform chemist review of confident, ambiguous, contradictory, and
  out-of-domain cases.

### Stage 5: conditioned assessment

- Add condition-derived features in `condition_recommender`.
- Train and evaluate a separate condition-aware model.
- Compare pre-condition and post-condition assessments for historical records.
- Do not expose condition-aware scores to query-side recommendation.

### Stage 6: controlled recommendation integration

- Add mechanism-hint similarity as a low-weight, soft scoring component only if
  evaluation demonstrates value.
- Do not use hints as hard chemistry filters.
- Report retrieval/ranking changes and evaluate them on held-out recommendation
  outcomes.

## Initial success criteria

The first version is successful when it:

- predicts broad, multi-label mechanism hints without source-label leakage;
- preserves unknown and disputed mechanisms;
- produces deterministic results for the same versioned input and artifact;
- explains predictions using graph-derived features and labeled precedents;
- abstains reliably on low-quality and out-of-domain reactions;
- shows improved precision over similarity-only and prevalence baselines on
  core-grouped validation splits;
- leaves reaction signatures and recommendation behavior unchanged during shadow
  deployment; and
- adds no dependency on legacy `chemtools`.

## References

1. Schwaller, P. et al. *Mapping the space of chemical reactions using
   attention-based neural networks*. Nature Machine Intelligence 3, 144-152
   (2021). https://doi.org/10.1038/s42256-020-00284-w
2. Schwaller, P. et al. *Extraction of organic chemistry grammar from
   unsupervised learning of chemical reactions*. Science Advances 7, eabe4166
   (2021). https://doi.org/10.1126/sciadv.abe4166
3. Heid, E. and Green, W. H. *Machine learning of reaction properties via learned
   representations of the condensed graph of reaction*. Journal of Chemical
   Information and Modeling 62, 2101-2110 (2022).
   https://doi.org/10.1021/acs.jcim.1c00975
4. Zeng, K. et al. *A general-purpose framework for chemical reaction
   representation with atomic correspondence and flexible condition adaptation*.
   Journal of Cheminformatics 18, 96 (2026).
   https://doi.org/10.1186/s13321-026-01201-w
5. Zhang, Q. et al. *A large-scale reaction dataset of mechanistic pathways of
   organic reactions*. Scientific Data 11, 829 (2024).
   https://doi.org/10.1038/s41597-024-03709-y
6. Miller, R. J. et al. *PMechDB: A Public Database of Elementary Polar Reaction
   Steps*. Journal of Chemical Information and Modeling 64, 2637-2646 (2024).
   https://doi.org/10.1021/acs.jcim.3c01810
7. Tavakoli, M. et al. *RMechDB: A Public Database of Elementary Radical Reaction
   Steps*. Journal of Chemical Information and Modeling 63, 2464-2470 (2023).
   https://doi.org/10.1021/acs.jcim.2c01359
8. Burés, J. and Larrosa, I. *Organic reaction mechanism classification using
   machine learning*. Nature 613, 689-695 (2023).
   https://doi.org/10.1038/s41586-022-05639-4

