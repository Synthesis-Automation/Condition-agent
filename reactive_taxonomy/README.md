# Reactive Taxonomy

Standalone, deterministic compound and reaction featurization built directly
on RDKit. This package must not import `chemtools`.

Public entry points:

- `featurize_molecule(smiles)`
- `featurize_reaction(reaction_smiles)`
- `resolve_source_label(source_label)`
- `validate_taxonomy()`
- `validate_source_label_mappings()`

Every reaction with normalized edit evidence also receives a structured
`display_label`. An exact, edit-consistent grammar label is preferred when one
is available. Otherwise, versioned generic edit patterns render transformations
such as `C–N substitution`, `C=C hydrogenation`, or `N=N reductive cleavage`.
Any unmatched edit set still receives a deterministic literal summary such as
`C–N bond formation`. The structured result retains the grammar or pattern ID,
the literal structural label, evidence, and confidence. Display style and label
definitions are explicitly excluded from reaction-signature identity.

Mapped observations also receive a template-free `ReactionCoreProjection`.
V2 keeps every edit-participating atom, selects primary centers only for
explanation, removes each remote connected component once, and records all cut
connections as typed attachment ports. Remote classes such as aryl,
heteroaryl, alkyl, alkenyl, alkynyl, acyl, ring-aliphatic,
heteroatom-linked, or generic R are derived from the removed molecular graph.
Exact fragment SMILES and functional-group annotations remain available even
when the concise label uses `Ar`, `HetAr`, or `R`.

Four versioned keys separate different uses of the minimized observation:
`RCX2` retains exact mapped-edit and fragment identity, `RCT2` retains typed
remote-subgraph identity, `RSH2` is the mapping-robust retrieval shape, and
`RCS2` is a diagnostic primary-center transition. `RSH2` also includes
participant handles/sites, retained remote shape, and event count so reactions
with similar centers but different partner chemistry do not collapse. `RCS2`
is never used for retrieval. Construction does not load reaction grammars or
the reaction-template registry. It abstains when no edited atom is observed on
both sides and does not change `RS3`, admission, or named-family results.
The reaction CLI and desktop featurizer show the minimized label, evidence,
shape and center keys, core size, remote subgraphs, attachment counts, and
warnings. Batch reaction CSV exports the same audit fields; JSON retains the
complete nested projection.

The desktop featurizer also presents mapped cores graphically. The renderer in
`visualization` draws active atoms and replaces only unchanged retained remote
subgraphs with deterministic `Ar`, `HetAr`, `R`, or indexed placeholders. It
always exposes the exact fragment-SMILES legend and the core evidence status;
external mapping is visibly review-only. This rendering does not affect core
identity, signatures, admission, or retrieval.

Mapped, single-event references can also be compiled into the versioned
`reaction_templates.v1.json` registry. This authoring registry stores a
map-number- and reactant-order-invariant edit fingerprint, optional family
annotations, provenance, and the mapped reference; it does not store or copy a
final reaction signature. Its semantic layer is intentionally small: each
template has a reaction label, a product label, and automatically derived roles
that link reference atom maps to canonical taxonomy site types. Query signatures
remain structure-derived, while matching active templates contribute
interpretation candidates. Use
`python -m reactive_taxonomy.reaction_template_cli --help` for import, list,
show, validation, and query matching commands. The PyQt6 wrapper is available
as `python -m app.reaction_template_registry_gui`.

When an incomplete query has no normalized edit set, active templates may
first compile minimal edit-bearing reactant roles from the mapped reference.
Repeated reference participants become explicit multiplicity requirements.
The bounded executor may reuse one supplied query species only for such a
repeatable role, applies the stored bond and schema-H edits, validates valence,
and requires exact canonical reconstruction of the reported main product.
This emits `exact_template_reconstruction_with_inferred_multiplicity` evidence
without inventing atom correspondence or an `RS3` signature.

Successful execution also emits a typed `TemplateReactionInterpretation`.
It retains the actual query component/atom/site bindings and renders a concise
label such as `R–CH=O + 2 × R–OH → acetal`. Steric accessibility, electronic
class, nearby groups, and context flags are read from the existing
`SiteReactivityProfile` for each bound site; templates do not contain descriptor
rules or SMARTS. Repeated supplied species remain separate bindings, while a
demonstrated inferred repeat is reported once with its multiplicity and warning.

Role display labels and map-level element alternatives are the only curated
generalization controls. For example, the registered Darzens reference labels
its derived `activated_sp3_carbon` role as `α-halo ester` and explicitly permits
`Br`, `Cl`, or `I` at the reference leaving-group map. The role also requires
the existing `alpha_to:ester` taxonomy token. The executor therefore requires
the query atom to be the ester-activated C–H site and must reconstruct the
reported epoxide exactly. Consequently, an unactivated alkyl halide, an
α-halo nitrile, or a ring-opened product does not match. These controls are
available during import as `--role-label ROLE=LABEL`,
`--role-tokens ROLE=alpha_to:ester`, and `--atom-elements MAP=Cl,Br,I`.

If exact reconstruction fails, a lower conservative fallback compares the
edited-centre before/after bond-state multisets. Neighbor hydrogen count and
external heavy-atom degree remain in this comparison, allowing an acetal
centre to remain distinct from a hemiacetal. This fallback emits
`template_center_transition_hypothesis` evidence.

Every valid parsed reaction also receives a
`ReactionCompletenessAssessment`. It records heavy-atom and element counts,
product-atom coverage, mapping coverage, suspected missing reactants or
insufficient partner multiplicity, and typed warnings. Exact reconstruction,
complete product mapping, conservative scaffold correspondence, or bounded
global multi-reactant correspondence may verify product-atom provenance.
Definite product-heavy-atom excess blocks signature
generation; unresolved provenance remains explicit for downstream review.
Reactant-side excess is retained without rejection because main-product
reaction records commonly omit byproducts. Missing reactants or stoichiometric
partner copies are never synthesized.

When one conserved scaffold uniquely exchanges one connected branch, an
incomplete record may retain an observation-only
`PartialProductTransformation`. Product-only atoms are grouped into a rooted
`ProductOriginGap` with atom references, internal bonds, attachment topology,
and a deterministic `PFG1` key. A product-heavy-atom provenance ledger links
conserved atoms to reactant atoms and leaves gap source atoms null. This
represents omitted multi-atom sources such as cyanide or azide without
inventing a reagent. Structurally matching
middle-side agents are recorded only as source support; ambiguous or absent
sources remain explicit, and the record is not promoted to a verified
`ReactionSignature`.

Exact grammar labels use one reaction-wide namespace for repeated generic
fragments. Distinct graph-derived fragments are numbered by semantic role and
retain the same alias across the arrow, for example
`Ar1–Br + Ar2–NH2 → Ar1–NH–Ar2`. Existing site-local notation is reconciled
with that namespace, so `R–Br + R1R2–NH` becomes
`R1–Br + R2R3–NH → R1–NR2R3`. A symbol remains unnumbered when only one
distinct fragment of that kind participates, and standalone molecule labels
remain unchanged.

For a single validated bond-order transformation, the preferred display adds
retained-neighbor and endpoint-hydrogen context. For example, the generic
interpretation `C=C hydrogenation` is retained in `transformation_label`, while
`reaction_label` becomes `Ar–CH=CH2 → Ar–CH2–CH3`. The nested display also
retains `structural_label`, `reactant_context_label`, and
`product_context_label`. Multi-edit cases that cannot be rendered safely keep
their literal edit summary.

Reductive amination is represented as one declarative two-partner grammar plus
one registered graph operator, rather than a reaction-name special case. The
grammar accepts available aldehyde, ketone, or formaldehyde carbonyl sites and
free primary or secondary amines. The operator removes the carbonyl oxygen,
forms C–N, adds H at carbon, and consumes N–H. Exact product reconstruction can
therefore render labels such as
`HeteroAr–CH=O + Ar–NH2 → HeteroAr–CH2–NH–Ar`; product mismatches remain
unverified, and multiple indistinguishable assignments remain unselected.

The high-ROI common-reaction layer also includes explicit alkyl C–N, C–O, and
C–S substitution plus terminal-alkene Heck coupling. Alkyl substitution uses
  the shared leaving-group and X–H site contracts and the registered
  `center_replacement` operator. It intentionally leaves `named_family` unset
because the graph alone does not distinguish SN1, SN2, or protection chemistry.
The Heck grammar uses alkene endpoint hydrogen counts to select a terminal
attachment site, removes Ar–X, forms Ar–C, and records alkene C–H loss. It does
not invent E/Z stereochemistry: a stereospecified product is exact only when
  that stereochemistry is supported by the input/operator result.

Alcohol C–O displacement is represented by one shared declarative rewrite for
O–H, N–H, and S–H partners. The operator removes the alcohol oxygen from the
reported main-product projection, forms the new C–X bond, and consumes partner
X–H. Defined tetrahedral centers enumerate inversion and retention as separate
structural outcomes; the observed product selects the supported outcome and the
signature records it explicitly. No Mitsunobu or other named family is forced
from this graph pattern alone because omitted reagents may be needed to support
that optional annotation.

The v2 grammar registry declares compatible molecular roles, not executable
operators or reaction archetypes. Every grammar selects one bounded
connectivity rewrite. `substitution`, `addition`, and `elimination` are derived
from the emitted or observed bond/H changes and describe graph topology, not
mechanism. Mapped reactions receive the same edit-derived archetype even when
no grammar or named family is available.

The internal Phase 1 connectivity observation preserves stronger evidence than
the compatibility edit view. `BondTransition` distinguishes definite
`bond`/`no_bond` states from `endpoint_absent` and `unknown`; every transition
retains observed-product, main-product-projection, exact-reconstruction,
correspondence-inference, or unresolved scope. Aggregated `HydrogenDelta` and
formal-charge `AtomStateTransition` observations are included in an internal
`ConnectivityEditGraph` with a deterministic `CEG1` shadow key. This graph is
dual-written inside edit normalization for evaluation only and does not yet
participate in `ReactionSignature`, serialized analyses, retrieval, or
recommendation behavior.

V2 uses one separately versioned, bounded connectivity-rewrite definition and
generic executor for all grammars. Suzuki C-C coupling, C-N/O/S
release-and-connect, Br2/ICl/Si-H/Si-B addition, beta elimination, and simple
bond-order changes use localized bond-state, schema-H, charge,
endpoint-permutation, product-seed, and authorized-projection instructions.
There is no shadow executor or legacy fallback.

Molecular featurization emits immutable `ReactiveLinkSite`,
`BondCapacitySite`, and `ConnectionEndpointSite` records in
`connectivity_sites`. Existing leaving-group, transfer-group, X-H, aromatic
C-H, addition-donor, and unsaturated-bond detectors remain annotation
provenance. Declarative adapters expose
real or virtual endpoints, carrier provenance, bounded bond capacity,
connection requirements, context, availability, and component-qualified IDs.
The shadow rewrite DSL now consumes these normalized interfaces, allowing
N/O/S-H, Si-H, B-H, and explicit A-B alkene addition to share the same
connectivity program. The adapters are not serialized into public analyses and
the legacy operators remain authoritative for that Phase 3 boundary.

Phase 4 extends the normalized interfaces to activated acyl/sulfonyl centers,
anionic partners, eliminable pairs, carbonyls, and alcohols. Versioned rewrite
definitions now cover reusable release-and-connect, two-anchor coupling,
alkene/alkyne addition, beta elimination, and simple bond-order changes.
Authority is declared per grammar after corpus-wide exact parity; 28 grammars
now execute through connectivity rewrites. Direct aromatic C-H substitution,
Chan-Lam coupling, reductive amination, Heck coupling, and multi-event
sequences remain legacy-backed where projection semantics or corpus parity are
not complete.

`pair_addition` combines an unsaturated-bond acceptor with either an existing
N/O/S–H site or a curated `addition_donor`. The latter represents both explicit
A–B bonds such as Br–Br and implicit A–H bonds such as Si–H and B–H. The
operator enumerates both constitutional orientations for A≠B, collapses
symmetry-equivalent products, reduces C=C or C≡C by one bond order, and emits
formed, broken, order-changed, and schema-level hydrogen edits. Broad addition
grammars require product verification and do not infer Markovnikov selectivity,
syn/anti selectivity, or a named family.

`pair_elimination` is the inverse net topology. The initial conservative
`eliminable_pair` detector exposes a carbon bearing Cl/Br/I and each adjacent
carbon bearing at least one hydrogen. Each beta site is a separate outcome, so
the reported product resolves regioselectivity. Omitted HX is not invented as
an observed byproduct.

Friedel–Crafts acylation is represented by an activated acyl electrophile and
an available aromatic C–H site. It reuses the generic handle-replacement
operator to remove the acyl leaving group, form the acyl-C–aryl bond, and record
aromatic C–H loss. Carboxylic acids are excluded because their acyl sites are
latent rather than activated. Regioisomers are selected only by exact product
reconstruction; unresolved sites and mapping conflicts remain visible.

Intramolecularity is a shared reaction-topology dimension rather than a second
set of family grammars. Grammar roles declare `same`, `different`, or
`same_or_different` component relationships. The same graph operators therefore
handle intermolecular joining and same-component ring closure. Every normalized
edit signature now carries `ReactionTopology`, including reaction scope, role
component membership, formed-bond scope, reactant tether distance, formed ring
size, and graph cycle-rank delta. Topology contributes to L0–L2 identity while
L3 remains a topology-agnostic edit fallback. L2 and L3 include schema-level
hydrogen gain or loss as well as formed, broken, and order-changed heavy-atom
bonds, so reactions with the same heavy-atom order change but different
hydrogen balance do not collapse to one identity. Mapped unknown reactions get
the same topology analysis without requiring a named family.

For example, `NCCc1ccccc1Br>>c1ccc2c(c1)CCN2` is rendered as
`intramolecular (5-membered ring) Ar–Br / R–NH2 → Ar–NH–R` using the same
`sp2_c_n_substitution` grammar as the corresponding intermolecular reaction.

Normalized edits are also partitioned into deterministic `ReactionEvent`
objects. Edits that share atom provenance form one event; disconnected edit
groups form a multi-event reaction. This represents mixed chemistry without
forcing one family label. For example, validated C-O and C-S substitutions at
different sites are rendered as `C-O substitution + C-S substitution`, while
two equivalent substitutions are counted as `2 x C-S substitution` in ASCII
style. Each event retains its edits, partners, sites, topology, evidence, and
optional interpretation. The reaction signature stores the ordered event
multiset and cross-event relations such as `shared_component`.

Multi-event inference still follows the normal evidence hierarchy. Valid atom
mapping is preferred. For balanced unmapped substitutions, distinct grammar
assignments may be composed and accepted only when the combined graph operators
exactly reconstruct the observed product. An atom-unbalanced record is not
assigned invented partner copies or atom correspondence merely to explain its
product.

When an unmapped reaction has exactly one conserved heavy-atom scaffold and one
product, a conservative correspondence fallback may supply edit evidence after
registered grammar reconstruction has failed. If that narrower path fails, a
bounded global fallback can combine conserved subgraphs from several reactants
into non-overlapping product assignments. Both paths require every product
heavy atom to be accounted for and accept alternatives only when the
minimum-edit mappings imply the same normalized chemistry. The global path is
limited by component, atom, candidate, and search-combination bounds and rejects
additional substantial products, product element excess, overflow, and
chemically distinct best mappings. Products with explicit stereochemistry are
also rejected because this fallback does not validate stereochemical
correspondence. Accepted global cases keep
`global_atom_correspondence` evidence, review-level chemistry confidence, and no
forced named family. Their labels describe the graph transition, while
recognized generic edit patterns or `generic_graph_transformation` provide a
mechanism-neutral transformation class.

If product atoms cannot be sourced from the reported reactant side, a bounded
single-attachment correspondence may instead retain an observation-only
product-origin gap. The fallback descriptor represents it with a deterministic
`PTS1` partial-transformation key and typed `FSR1` fragment-source
requirements. These are structural observations only: they do not create a
reaction signature, named family, or atom correspondence. The recommender may
later compare the requirements with its separately owned, curated condition
capability registry.

Integrated CLI tester:

```powershell
python -m reactive_taxonomy.cli validate
python -m reactive_taxonomy.cli self-test
python -m reactive_taxonomy.cli molecule "Brc1ccc(N)cc1C#N"
python -m reactive_taxonomy.cli reaction "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"

python -m reactive_taxonomy.cli reaction "Oc1cccc(C(F)(F)F)c1.Oc1cccc(C(F)(F)F)c1.N#Cc1c(F)c(F)c(Oc2cccc(Oc3ccccc3)c2)c(F)c1F>>N#Cc1c(Oc2cccc(C(F)(F)F)c2)c(F)c(Oc2cccc(Oc3ccccc3)c2)c(F)c1Oc1cccc(C(F)(F)F)c1"


python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/molecule_features.jsonl
python -m reactive_taxonomy.cli batch examples/sample_compounds.csv --mode molecule --output results/sample_compounds_featurized.csv
python -m reactive_taxonomy.cli batch examples/sample_reactions.csv --mode reaction --concise --output results/sample_reaction_featurized.csv


python -m reactive_taxonomy.cli batch examples/test_reactions.csv --mode reaction --concise --output results/test_reactions_featurized.csv


python -m reactive_taxonomy.cli batch examples/dataset_300/C_N_Coupling.csv --mode reaction --output results/data_300.csv
python -m reactive_taxonomy.cli batch examples/dataset_300/C_N_Coupling.csv --mode reaction --concise --output results/data_300.csv

python -m reactive_taxonomy.cli batch examples/dataset_300/C_N_Coupling.csv --mode reaction --concise --output results/C_N_Coupling_300.csv

python -m reactive_taxonomy.cli batch data-processor/reaction_dataset/Amide_formation.csv --mode reaction --concise --output results/Amide_formation.csv



```

Reaction-mode CSV output places a semicolon-separated `spectator_groups`
column immediately after `reaction_label`. Values are stable functional-group
IDs such as `nitrile`. Add `--concise` to write only `reaction_smiles`,
`reaction_label`, and `spectator_groups`, in that order.

Phase 1 molecular-feature evaluation, including machine metrics and a blind
chemist-review packet:

```powershell
python -m reactive_taxonomy.molecular_feature_evaluation_cli `
  results/molecular_feature_evaluation
```

The versioned answer key is
`benchmarks/molecular_features/benchmark_manifest.v1.json`. Generated artifacts
include `machine_report.json`, detailed case results, `chemist_review.csv`,
highlighted `review_structures.html`, and `disagreements.csv`. The machine gate
does not complete the separate human chemist gate.

Phase 2 reaction-edit evaluation, including mapped edits, exact reconstruction,
hydrogen changes, evidence conflicts, and a blind reaction-center review packet:

```powershell
python -m reactive_taxonomy.reaction_edit_evaluation_cli `
  results/reaction_edit_evaluation
```

Its versioned answer key is
`benchmarks/reaction_edits/benchmark_manifest.v1.json`. In addition to the five
standard evaluation artifacts, it writes
`connectivity_shadow_report.json` with per-case shadow keys, observation-scope
counts, compatibility parity, unsupported bond domains, and canonicalization
overflow.

Use concise mode for a chemist-readable view containing the primary labels and
the shared context-aware reactivity profile:

```powershell
python -m reactive_taxonomy.cli molecule "Brc1ccc(N)cc1C#N" --concise
python -m reactive_taxonomy.cli molecule "BrBr" --concise

python -m reactive_taxonomy.cli reaction "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" --concise

python -m reactive_taxonomy.cli reaction "BrBr.C=CCCCCCCCCCCC>>CCCCCCCCCCCC(Br)CBr" --concise

```

For reactive sites bearing alkyl groups, the readable output separates the
substitution at the reactive center from attachment-carbon sterics. For
example, tert-butylamine remains a primary amine (`R–NH2`) but its attached
alkyl group is reported as `tertiary` and alpha-branched.

Reaction partner contexts retain generic `Ar` and `HeteroAr` handle context,
while the nested profile records the actual aromatic family, ring sizes,
fusion, anchor-relative heteroatoms, ortho occupancy and burden, and
context-specific electronic-demand evidence. Alkyl profiles report
methyl/primary/secondary/tertiary identity, branching, cyclicity,
benzylic/allylic/propargylic activation, and beta-H status. Heteroatom profiles
separate N/O/S/P identity from attached-group burden and lone-pair
availability. GUI, CLI, concise review, condition rules, signature tokens, and
retrieval all consume this same typed record.

Compact output omits atom indices, raw graph-shell counts, methods, and scores.
Expanded rendering retains scores, contributors, confidence, and provenance.
Reaction signatures use stable categorical profile tokens under schema `3.0`
and the `RS3` namespace; raw scores and display labels do not enter identity.

Unsaturated-bond labels expose endpoint substitution rather than collapsing all
sites to `C=C` or `C≡C`. Examples include `H2C=CH2`, `H2C=CHR1`,
`R1R2C=CR3R4`, `R1–C≡C–H`, and `R1–C≡C–R2`; defined alkene E/Z
stereochemistry is retained. Canonical `PI|Alkene` and `PI|Alkyne` signatures
remain display-independent.

Organic nitriles are also bond-localized π handles with typed carbon and
nitrogen endpoints (`PI|Nitrile`, `R–C≡N`). Organic azides are separate
multi-atom dipolar handles with attachment, proximal-, central-, and
terminal-nitrogen roles (`DG|Azide|Organic`, `R–N3`). Their declared addition,
reduction, or cycloaddition modes describe possible chemistry only; a reaction
interpretation still requires product/edit evidence. Cyanide salts, isocyanides,
and inorganic azide are not promoted to these organic-handle contracts.

Organic heteroatom-pair bonds share a bond-localized handle family with two
typed attachments: `HB|Azo` (`R1–N=N–R2`), `HB|Disulfide`
(`R1–S–S–R2`), and `HB|Peroxide` (`R1–O–O–R2`). Hydroperoxides,
hydrogen peroxide, thiols, ordinary ethers, and unattached diazene are excluded.
Mapped reductions can therefore report the observed central-bond cleavage and
the corresponding hydrogen gains without requiring a named reaction family.

Alkyl reactive centers retain `Alkyl` as their broad machine context while
recording synthetically important attachment subtypes. Benzylic, allylic, and
propargylic leaving groups are rendered as labels such as `Benzyl–Cl`,
`Allyl–Br`, and `Propargyl–OMs`.

Add `--format json` to the single-record and validation commands for the full
typed result. Batch mode prints a coverage summary and optionally writes
source-traceable JSONL or CSV. The format is inferred from the output extension
or selected with `--output-format`. CSV contains readable summary columns;
JSONL retains the complete typed analysis. Use `--column NAME` when a CSV does
not use a standard `smiles` or `reaction_smiles` header.

Definitions own handle detection, functional groups, rendering, reaction
grammars, and descriptor weights. Python implements graph interpretation,
candidate resolution, reaction operators, and typed result contracts.

Dependency direction:

```text
applications and recommenders -> reactive_taxonomy -> RDKit
```

The recommender and dataset conversion layers belong outside this package.

## Source-label normalization

Legacy dataset labels are resolved through the versioned
`definitions/source_label_crosswalk.v1.json` file. Resolution separates a
stable machine label from a chemist-facing display and optional environment
constraints. For example, `RNH2 a-branch` resolves to machine label `R-NH2`,
display label `R–NH₂ (α-C branched)`, and an `alpha_branched=true` qualifier.
Unsupported labels are returned unchanged with `mapping_status=unresolved`.

Qualified mappings may deliberately stop at a validated handle-family prefix.
`Alkyl-H acidic`, for example, resolves to `XH|Csp3` because the source label
does not identify the activating group or hydrogen count. Structure-derived
sites remain exact, such as `XH|Csp3|H3|alpha_to:nitrile`. The taxonomy emits
these sites only for explicit alpha-to-aldehyde, ketone, ester, amide, nitrile,
sulfone, or nitro rules; ordinary alkane C–H sites remain excluded.

The source aliases `ArH` and `Ar-H` both resolve canonically to the
atom-localized aromatic C–H signature `CH|ArH` and display label `Ar–H`.
This handle is shared by structurally distinct grammars, including
intermolecular direct arylation and Friedel–Crafts acylation; the handle alone
does not imply either reaction.
