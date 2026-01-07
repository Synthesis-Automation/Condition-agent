# Organic Chemistry Taxonomy System Overview (Living Doc)

Last updated: 2026-01-07
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

---

## Design philosophy

- SMARTS live at the lowest stable level (groups and direct motifs), not at every
  derived layer. This keeps pattern maintenance centralized and consistent.
- Compound motifs are generated dynamically from group templates (A-B patterns),
  so expanding group coverage automatically expands compound typing.
- Reaction typing is defined by motif-slot constraints and shared motif sets
  (for example `@sp2_electrophiles`), avoiding duplicated SMARTS across families.
- Deterministic outputs are favored over model-only predictions; ML can be added
  later but must resolve to canonical taxonomy IDs.
- All SMARTS compilation goes through the cache helper for speed and reuse.
- The system is layered: groups -> compounds -> reactions, so higher levels are
  compositional and inherited from lower-level definitions.
- It generates human-readable tags and descriptors (e.g., sterically demanding
  aryl chloride) that map to well-defined motif and steric/electronic analysis,
  keeping tag meaning explicit and reproducible.
- These meaning tags make it easier for LLMs to access structured chemistry
  context and reason from well-defined signals instead of raw data such as smiles strings.

---

## Key data assets (taxonomy/data)

- `organic_groups.v1.3.json`
  - Canonical organic groups (scaffold/substituent) with attachpoints in SMARTS.
  - Examples (scaffold): `Ar` (aryl attachment), `R` (alkyl attachment), `Alkenyl`.
  - Examples (substituent): `Br`, `Cl`, `NH2`, `OH`, `B(OH)2`.
- `organic_compounds.v1.3.json`
  - Compound motifs built from group templates or direct SMARTS.
  - Examples: `Ar-Br`, `Ar-B(OH)2`, `Ar-NH2`, `R-OH`.
- `scaffold_motifs.v1.3.json`
  - Direct SMARTS motifs that do not fit simple A-B templates.
- `group_logic.json`
  - Group sets and priorities (e.g., X, LeavingGroup).
- `compound_logic.json`
  - Motif sets for reaction typing (e.g., @sp2_electrophiles).
- `reaction_types.v4.0.json`
  - Canonical reaction types with motif-slot constraints and aliases.
  - Examples: `suzuki_miyaura`, `buchwald_hartwig`, `sonogashira`, `snar_cn`.
- `reagent_roles.v2.json`
  - Reagent roles and priorities (v2).

Notes

- Reagent roles are defined in `reagent_roles.v2.json`. Family allowlists are
  handled outside the v2 core (legacy registry and heuristics).

---

## How typing works (organic groups -> compounds -> reactions)

1. Organic groups define scaffold/substituent SMARTS with attachpoint maps.
   - Source: `chemtools/taxonomy/data/organic_groups.v1.3.json`.
2. Organic compounds combine groups into motifs.
   - Source: `chemtools/taxonomy/data/organic_compounds.v1.3.json`.
   - Templates default to `single_bond` and `via_oxygen`.
3. Motif detection compiles SMARTS and picks the most specific hits per site.
   - Module: `chemtools/featurizers/motif_detect.py`.
4. Reaction typing matches motif hits to reaction slot requirements.
   - Types + slots: `chemtools/taxonomy/data/reaction_types.v4.0.json`.
   - Motif sets: `chemtools/taxonomy/data/compound_logic.json`.
   - Engine: `chemtools/featurizers/reaction_detection.py`.

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

---

## Canonicalization and naming alignment

- Reaction type labels must map to taxonomy IDs.
  - Aliases resolved in `chemtools/taxonomy/reaction_catalog.py`.
  - Legacy names mapped in `chemtools/featurizers/analysis/reactions.py`.
- Dataset naming guidance is in `docs/NAMING_CONVENTION.md`.

---

## Extension workflow (quick checklist)

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

## Future directions (aligning to program goals)

- Canonicalization and standardization
  - Enforce alias resolution to taxonomy IDs across all detection paths.
  - Maintain consistent naming in dataset ingestion and curated DBs.
- Dynamic compound typing
  - Feed undocumented motif hits back into `organic_compounds.v1.3.json`.
  - Expand logic sets to cover common motif groupings.
- HTE database and recommendations
  - Use motif typing and reaction typing to create reaction-level HTE indices.
  - Add similarity search and uncertainty ranking for conditions.
- Dataset extraction and full reaction picture
  - Use unified reaction bundles to extract motifs, roles, and reagents.
  - Build dashboards for coverage, gaps, and taxonomy expansion.

---

## Open questions and risks

- Reagent coverage depends on legacy registry content and heuristic patterns; CAS
  data lives in files such as `data/reagent_db/condensation_agent.json`.
- Motif overlap: confirm priority rules are consistent for ambiguous sites.
- Reaction typing coverage: expand motif sets as new families are added.

---

## Changelog

- YYYY-MM-DD: Initial living doc created.
