# Organic Chemistry Taxonomy System Overview (Living Doc)

Last updated: 2026-01-10
Owner: Chemistry data platform team

Chemistry introduction

We define a controlled vocabulary of organic groups and
compound motifs (e.g., Ar-Br, Ar-B(OH)2, Ar-NH2), then use those motifs to
describe reaction families (e.g., Suzuki, Buchwald-Hartwig) and the roles of
reagents (catalyst, ligand, base, solvent). The goal is to capture organic
reactivity in a standard, machine-readable form so that condition extraction,
recommendation, and dataset analytics reflect chemical meaning rather than
black-box labels.

Purpose and scope

- Document how taxonomy data and code work together across typing, extraction, and recommendation.
- Provide a single reference for extension and alignment decisions.
- Cover reaction taxonomy, organic groups and compounds, reagent taxonomy, and downstream usage.
- Clarify where the taxonomy ends and where external heuristics or ML models begin.
- Provide guardrails for changes that can impact training data or downstream analytics.

---

## System boundaries and non-goals

- This taxonomy defines identity, roles, and reaction typing; it does not define kinetics, yields, or optimal conditions.
- Free text or vendor labels are treated as inputs, not authoritative ground truth.
- Model training or ranking logic lives outside the taxonomy; this system provides stable, interpretable features.
- Inventory, purchasing, or safety compliance is out of scope except for basic reagent roles.

## Design philosophy

- **Atoms over Labels**: We prefer structural evidence (SMARTS) over brittle string matching.
- **Hierarchical Priority**: Chemical scaffolds are ranked (e.g., Aromatic > Aliphatic) to resolve ambiguity in polyfunctional molecules. This ensures a Carbon attached to both a Phenyl and an Ethyl group is correctly typed as "Aryl-X" rather than "Alkyl-X".
- **Separation of Concerns**: SMARTS patterns define _identity_, Logic sets define _role_, and Templates define _connectivity_.
- **Substituent-Level Deduplication**: The default filtering mode collapses matches at the substituent atom level (`site_filter="substituent"`). This prevents "label explosion" (e.g., getting both `Ar-NR2` and `RCH2-NR2` for the same nitrogen).
- **Transformation-Awareness**: Reaction keys are represented as net-chemical changes: `[Reacted] -> [Formed] || [Spectators]`. This allows the system to distinguish between a pyridine that reacts (N-arylation) and a pyridine used as a base (spectator).
- **Stateless & Deterministic**: The system relies on fixed rules and cached SMARTS compilation for 100% reproducibility.
- **Machine-Ready Chemistry**: These tags make it easier for LLMs to access structured chemistry context and reason from well-defined signals instead of raw data such as smiles strings.

### SMARTS compilation policy

- Use `chemtools/util/smarts_cache.py` for all SMARTS compilation.
- Prefer lazy compilation during detection rather than module import time.
- Use `compile_smarts(..., validate=True)` only when introducing or auditing new patterns.
- Avoid direct calls to RDKit `Chem.MolFromSmarts()` in taxonomy code.

---

## Chemical Logic & Motif Resolution

### The A-B Model (Scaffold-Link-Substituent)

Most organic molecules are modeled as `A-B` pairs:

- **A (Scaffold)**: The core framework (Ar, Alkyl, Alkenyl, etc.). These are intentionally "Variety-Heavy" to capture local steric and electronic environments (e.g., distinguishing primary `RCH2-` from tertiary `R3C-`).
- **B (Substituent)**: The functional group (Br, NH2, B(OH)2, etc.). These are "Interface-Heavy," defining the actual chemical transformation.
- **Link**: The bond connecting them (single, double, via oxygen, etc.).

#### Naming Convention: R vs Alkyl

To ensure chemical clarity, we distinguish between generic wildcards and specific alkyl groups:

- **Alkyl**: Specifically refers to an $sp^3$ carbon scaffold (e.g., `Alkyl-Br`).
- **R**: A generic wildcard representing any non-hydrogen group or atom.
- **Hybrid Cases**: In names like `RCH2` (Primary Alkyl), `R2CH` (Secondary), and `R3C` (Tertiary), the `R` represents "any non-H substituent" attached to the central alkyl carbon. These names are preserved as they accurately describe the degree of substitution.

#### Design Consideration: Asymmetry vs Complexity (A-B vs A-L-B)

There is an inherent asymmetry in this system: **Scaffolds** (A) have dozens of nuanced variations, while **Substituents** (B) are relatively monolithic. This reflects the reality of synthetic chemistry—we care deeply about the _environment_ of a functional group (is it hindered? is it benzylic?).

- **Why not make everything A-L-B?**: Introducing a universal Linker layer (e.g., treating `RCH2-` as `R-link(CH2)-`) would make the system more "symmetrical" but significantly more complex. It would turn a 2-layer hierarchy into a 3-layer one, requiring more complex templates and making SMARTS detection slower.
- **The "Gestalt" Principle**: We treat `RCH2-` as a first-class Scaffold because "Primary Alkyl" is a holistic chemical concept. Breaking it down into components loses the semantic "Gestalt" that a researcher or a recommender engine expects (e.g., "This reaction works for primary alkyls").
- **When to use Links**: Templates like `via_oxygen` are reserved for future use (not now). E.g., in places where the linking atom (O, S, NH) fundamentally changes the chemistry class (e.g., Sulfonate leaving groups) or where the "B" part is actually an entire fragment (e.g., `Ar-O-R`).

This modularity allows the system to generate thousands of motifs (e.g., `Ar-NH2`, `R-NH2`, `Acyl-NH2`) from a small set of primitive definitions.

### Recommended feature layers (conditions recommendation, minimal maintenance)

The goal is a small, stable feature stack that keeps chemistry signal high and
cost low. Start with the core, then add overlays only if they improve model
quality.

- **Core (required)**:
  - Motif delta (`reacted`, `formed`, `spectator`) with confidence.
  - Reagent roles (`base`, `ligand`, `catalyst`, `solvent`) from v2 taxonomy.
- **Lightweight overlays (recommended)**:
  - Steric tags (e.g., `ortho`, `benzylic`, `hindered`).
  - Electronic tags (e.g., `EWG`, `EDG`, `heteroaryl`).
  - Scaffold spectators (e.g., `pyridine`, `dmf`) to prevent role confusion.
- **Optional, low-maintenance extras**:
  - Simple topology tags (ring count, heteroatom count, substitution degree).
  - Basic physchem descriptors (MW, HBD/HBA, TPSA, logP) when RDKit is available.

Avoid high-dimensional fingerprints or full graph descriptors unless accuracy
gains justify the maintenance cost.

### Priority Hierarchy

When multiple scaffolds claim the same functional group atom, we apply a priority-based resolution:

1. **Priority 5 (Aromatic/Acyl)**: `Ar`, `AromN`, `Acyl`, `Alkenyl`.
2. **Priority 3 (Aliphatic)**: `RCH2`, `R2CH`, `R3C`.
3. **Priority 1 (Specific/Generic)**: `Alkyl`.
4. **Priority 0 (Wildcards)**: `R` (non-H), `Any_Scaffold` (any atom).

This hierarchy ensures that chemical "intent" is preserved—aromaticity and unsaturation always take precedence over simple alkyl chains during classification.

### Stoichiometry & Transformation Calculus

To determine what happened in a reaction, we use raw motif counts from reactants and products:

1. **Count Reactant Motifs**: `Ar-Br: 1, Ar-NH2: 1, Pyridine: 1`
2. **Count Product Motifs**: `Ar-NR2: 1, Pyridine: 1`
3. **Net Change (Counter Math)**:
   - `-1 Ar-Br`, `-1 Ar-NH2` (Reacted)
   - `+1 Ar-NR2` (Formed)
   - `Pyridine` (Spectator - count unchanged)
4. **Persistent Spectators**: Common heterocycles (Pyridine, Quinoline) are mathematically protected. If they are present in the reactant but missing from the product (due to extraction noise), they are automatically re-assigned to the spectator slot.

---

## Key data assets (taxonomy/data)

- `organic_groups.v1.3.json`
  - Canonical organic groups (scaffold/substituent) with attachpoints in SMARTS.
  - Examples (scaffold): `Ar` (aryl), `Alkyl` ($sp^3$ carbon), `R` (non-H wildcard), `Alkenyl`.
  - Examples (substituent): `Br`, `Cl`, `NH2`, `OH`, `B(OH)2`.
- `organic_compounds.v1.3.json`
  - Compound motifs built from group templates or direct SMARTS.
  - Examples: `Ar-Br` (Aryl Bromide), `Ar-AromN` (N-Aryl Heterocycle), `Alkyl-OH` (Alcohol).
- `scaffold_motifs.v1.3.json`
  - Direct SMARTS motifs for whole-molecule spectators or complex cores that do not fit simple A-B templates.
  - Examples: `Pyridine`, `Quinoline`, `THF`, `DMF`.
- `group_logic.json`
  - Group sets and priorities (e.g., `X = {Cl, Br, I, OTf}`, `LeavingGroup`).
- `compound_logic.json`
  - Motif sets for reaction typing.
  - Examples: `@sp2_electrophiles` (includes `Ar-X`, `Alkenyl-X`), `@amine_nucleophiles`.
- `reaction_types.v4.0.json`
  - Canonical reaction types with motif-slot constraints and aliases.
  - High-precision matching: `buchwald_hartwig` requires an electrophile slot AND a nucleophile slot to be occupied.
- `reagent_roles.v2.json`
  - Reagent roles and priorities (v2).

Notes

- Reagent roles are defined in `reagent_roles.v2.json`. Family allowlists are
  handled outside the v2 core (legacy registry and heuristics).

---

## End-to-end data flow (high level)

1. **Input normalization**: Canonicalize reactant and product SMILES, map atoms when available.
2. **Group detection**: Identify scaffold and substituent groups using SMARTS.
3. **Compound synthesis**: Combine groups into motifs via templates and direct SMARTS.
4. **Resolution**: Apply subsumption, priority, and site filtering to avoid duplicates.
5. **Delta calculation**: Compute reacted, formed, and spectator motifs.
6. **Reaction typing**: Match motif delta against reaction type slots and constraints.
7. **Feature bundling**: Assemble condition features (motifs, roles, sterics, electronics).
8. **Outputs**: Provide stable IDs to downstream recommenders, analytics, and LLM tools.

## How typing works (organic groups -> compounds -> reactions)

1. **Primitive Detection**: Organic groups define scaffold/substituent SMARTS with attachpoint maps (`:1` for scaffold atom, `:2` for substituent atom).
   - Source: `chemtools/taxonomy/data/organic_groups.v1.3.json`.
2. **Template Synthesis**: Organic compounds combine groups into motifs using Templates (e.g. `{A}-{B}`).
   - Source: `chemtools/taxonomy/data/organic_compounds.v1.3.json`.
3. **Resolution & Filtering**: Motif detection compiles SMARTS and resolves overlapping matches.
   - **Subsumption**: If motif A is a structural subset of motif B, B wins.
   - **Priority**: Aromatic scaffolds win over alkyl ones at the same attachment site.
   - **Site Filtering**: Controlled by `site_filter`. The "substituent" mode ensures each Nitrogen, Oxygen, or Carbon center gets exactly one classification.
4. **Reaction Assignment**: Matches motif profiles against reaction type slot requirements.
   - **Transformation Match**: Checks the net change in motifs (delta-Counter).
   - **Role Evidence**: Uses reagents (bases, ligands) as supporting evidence if available.
   - Engine: `chemtools/featurizers/reaction_detection.py`.

### Edge cases and safeguards

- **Missing products**: If products are missing or unmapped, fall back to reactant-only motif profiles and lower confidence.
- **Multi-site ambiguity**: Priority rules resolve overlapping scaffold assignments at the same site.
- **Spectator protection**: Heterocycles like pyridine are protected from false-negative drops.
- **RDKit optionality**: Feature extraction must degrade gracefully when RDKit is unavailable.

---

## Core modules and responsibilities

- `chemtools/featurizers/motif_registry.py`
  - Build compound motif registry from groups + compounds + templates.
- `chemtools/featurizers/motif_detect.py`
  - SMARTS detection, subsumption filtering, priority tie-breaks.
- `chemtools/featurizers/reaction_detection.py`
  - Match motif profiles against reaction type slot requirements.
- `chemtools/taxonomy/reaction_catalog.py`
  - Load reaction types and alias map; expand motif sets.
- `chemtools/taxonomy/calculable_spec.py`
  - Generate calculable feature spec from organic compounds + overlays.
- `chemtools/util/functional_groups.py`
  - Functional group detection sourced from calculable spec.
- `chemtools/taxonomy/reagent_v2.py`
  - Reagent role classification (CAS -> SMARTS -> name, if families provided).
- `chemtools/reagent/taxonomy_store.py`
  - Heuristics and registry helpers for reagent classification.
- `chemtools/featurizers/unified.py`
  - Unified reaction bundle: motifs, sterics, electronics, roles.
- `chemtools/reaction_type_detection.py`
  - Unified detect API returning family + confidence + slot evidence.

---

## Condition extraction and recommendation

- Reagent roles for conditions come from the v2 taxonomy.
  - Classifier: `chemtools/taxonomy/reagent_v2.py`.
  - Heuristics: `chemtools/reagent/taxonomy_store.py`.
- HTE recommender uses motif typing for reactant keying.
  - Engine: `chemtools/HTE/recommender.py`.
- Unified recommendation index merges protocols, datasets, and HTE.
  - Builder: `chemtools/recommend/index_builder.py`.
  - Recommender: `chemtools/recommend/unified.py`.

## Spectator ranking system (HTE recommender)

Spectator motifs and groups provide a secondary ranking signal for condition
recommendation. This is separate from reaction typing and is only used for
scoring and ranking matches.

- **Spectator derivation**: Spectators are motifs that appear in both reactants
  and products (net unchanged). When product SMILES is unavailable, spectators
  default to query motifs minus reacted motifs from the transformation key.
- **Spectator group extraction**: Convert spectator motifs to group IDs by
  taking the suffix of motif IDs (e.g., `Ar-NH2` -> `NH2`) and include whole
  scaffold motifs (e.g., `Pyridine`). A stoplist removes generic groups like
  `Ar`, `R`, `Any`, `Alkyl`, `Alkenyl`, `Alkynyl`, and `H`.
- **Transformation match score** (core + spectator):
  - Require reacted motifs to be a subset of query motifs.
  - `match_score = 0.5 + 0.5 * spectator_score`
  - `spectator_score` uses Jaccard similarity for overlap, with fallbacks:
    - both empty: 1.0
    - query has spectators but DB does not: 0.3
    - DB has spectators but query does not: 0.5
- **Spectator group weighting** (when `spectator_groups` are present in HTE data):
  - `match_score *= (0.7 + 0.3 * spectator_similarity)`
  - `spectator_similarity` is Jaccard with empty-set fallbacks:
    - both empty: 1.0
    - query empty: 0.7
    - DB empty: 0.3

Evaluation snapshot (C-S coupling canonical dataset, `results/spectator_groups_eval.json`):

- With spectator groups: hit@10 0.42, MRR 0.50, avg rank 3.43
- Without spectator groups: hit@10 0.18, MRR 0.32, avg rank 5.16
- Rank deltas: 107 improved, 23 worsened, 250 gained, 13 lost

---

## Canonicalization and naming alignment

- Reaction type labels must map to taxonomy IDs.
  - Dataset fields should store the taxonomy `id`; `name` is display-only and resolved via aliases.
  - Aliases resolved in `chemtools/taxonomy/reaction_catalog.py`.
  - Legacy names mapped in `chemtools/featurizers/analysis/reactions.py`.
- Dataset naming guidance is in `docs/NAMING_CONVENTION.md`.

---

## Runtime contracts and schema alignment

- Public outputs should align with `chemtools/contracts.py` and documented feature schemas.
- Breaking changes require a schema version bump and an explicit changelog entry.

## Extension workflow (quick checklist)

### Adding a new Chemical Motif

1. **Define the Group**: If it's a new functional group (e.g., `Sulfinate`), add it to `chemtools/taxonomy/data/organic_groups.v1.3.json`.
2. **Define the Compound**: Combine it with a scaffold (e.g., `Ar-Sulfinate`) in `organic_compounds.v1.3.json`.
3. **Set Priority**: Ensure the group has a priority score relative to existing ones to handle overlap.
4. **Test Resolution**: Run `python chemtools/featurizers/motif_detect.py` on a molecule containing the group to verify it is detected and prioritized correctly.

Add or update organic groups and compounds

1. Add group SMARTS to `chemtools/taxonomy/data/organic_groups.v1.3.json`.
2. Add compound motifs to `chemtools/taxonomy/data/organic_compounds.v1.3.json`.
3. Update set logic in `chemtools/taxonomy/data/group_logic.json` or
   `chemtools/taxonomy/data/compound_logic.json` as needed.
4. Re-run detection tests and spot-check motifs.

Add a new reaction type

1. Add a new entry to `chemtools/taxonomy/data/reaction_types.v4.0.json`.
2. Use motif sets from `compound_logic.json` when possible.
3. Ensure aliases cover common naming variants.

Add or update reagent families

1. Update legacy reagent registry files under `data/reagent_db`.
2. Keep role IDs aligned with `reagent_roles.v2.json`.
3. Confirm name normalization and SMARTS are compatible with classifier rules.

Build or refresh unified recommendation index

1. Standardize datasets (see `data-processor/convert_hte_to_canonical.py`).
2. Run `python data-processor/build_unified_recommendation_index.py`.

---

## Core principles

- Single source of truth: taxonomy JSON files in `chemtools/taxonomy/data`.
- Deterministic, cacheable: SMARTS compilation uses `chemtools/util/smarts_cache.py`.
- Canonical IDs only: detection outputs should resolve to taxonomy IDs.

---

## Feature schema draft (condition recommendation)

Define a stable, versioned feature bundle so recommendation logic can rely on a
consistent contract while taxonomy evolves.

- **Schema ID**: `cond_features.v1`
- **Feature groups**:
  - `reaction_delta`: net motif changes (formed, reacted, spectator).
  - `motif_context`: primary motif per site plus secondary context tags.
  - `role_evidence`: reagent roles with confidence + provenance.
  - `steric_overlay`: steric tags (ortho, hindered, bulky).
  - `electronic_overlay`: EW/ED group tags, aromaticity, heteroaryl.
  - `reaction_family`: canonical taxonomy ID + confidence.
- **Confidence fields**:
  - `reaction_delta.confidence`: penalize missing products, low mapping quality.
  - `role_evidence.confidence`: based on name/SMARTS match strength.
- **Evidence fields**:
  - `reaction_delta.evidence`: motif counts and change rationale.
  - `role_evidence.provenance`: `name_match`, `smarts_match`, `cas_match`.
- **Versioning**:
  - Change log in `chemtools/taxonomy/v2_docs/A_TAXONOMY_SYSTEM_OVERVIEW.md`.
  - Bump schema ID on breaking changes only.

---

## Actionable TODO (near-term)

- **Tie-breaker policy**: Define deterministic ordering within each priority
  class (e.g., Aromatic > Acyl > Alkenyl) and document a secondary-label policy
  for multi-context sites.
- **Secondary context tags**: Allow a site to carry one primary motif plus
  optional secondary tags (e.g., `N-arylated` + `benzylic_adjacent`) to preserve
  condition-relevant context without label explosion.
- **Spectator config**: Move the "persistent spectator" list to a named
  taxonomy data file and document how it is used.
- **Transformation quality score**: Add a scoring heuristic for net-change
  reliability (missing products, mapping gaps, low motif coverage).
- **SMARTS policy**: Make `compile_smarts()` the only allowed compilation
  mechanism in examples and module docs.
- **Dynamic compound typing workflow**: Define candidate -> validation ->
  approval -> version bump steps, with rollback criteria.
- **Coverage dashboards**: Track motif hit rate, family coverage, role coverage,
  and unknown buckets to guide taxonomy expansion.

---

## Future directions (aligning to program goals)

- **Transformation-Aware HTE Indices**:
  - Using net-change logic (`Ar-Br|Ar-NH2 -> Ar-NR2`) to index 20k+ reactions.
  - Filtering "chemical noise" using persistent spectator protection for heterocycles.
- **Dynamic Compound Typing**:
  - Feed undocumented motif hits back into `organic_compounds.v1.3.json`.
  - Expand logic sets to cover common motif groupings.
- **Steric & Electronic Overlays**:
  - Automatically injecting steric tags (e.g., `ortho-substituted`) into taxonomy IDs to create more granular categories for recommendation.
- **Dataset extraction and full reaction picture**
  - Use unified reaction bundles to extract motifs, roles, and reagents.
  - Build dashboards for coverage, gaps, and taxonomy expansion.

---

## Open questions and risks

- Reagent coverage depends on legacy registry content and heuristic patterns; CAS
  data lives in files such as `data/reagent_db/condensation_agent.json`.
- Motif overlap: confirm priority rules are consistent for ambiguous sites.
- Reaction typing coverage: expand motif sets as new families are added.

---
