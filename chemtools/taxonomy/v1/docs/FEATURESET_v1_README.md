# Organic chemistry SMARTS feature set (v1)

Generated: 2025-12-18

## Files

### Core JSON (v1)

- `chemtools/taxonomy/v1/specs/organic_groups.v1.json` — chemist-style group labels (Ar-, Vinyl-, -CHO, -CN, -NO2, -B(OH)2, -Bpin, -OTf, ...)
- `chemtools/taxonomy/v1/specs/smarts_templates.v1.json` — template rules and explicit expansions (build-time)
- `chemtools/taxonomy/v1/specs/calculable_features.atomic.v1.json` — atomic SMARTS detectors (runtime)
- `chemtools/taxonomy/v1/specs/calculable_features.derived.v1.json` — boolean logic detectors (runtime)
- `chemtools/taxonomy/v1/specs/calculable_features.compiled.v1.json` — combined atomic + derived
- `chemtools/taxonomy/v1/specs/reactant_types.v1.json` — feature-token taxonomy (cross-coupling + amide formation focus)
- `chemtools/taxonomy/v1/specs/reaction_types.v1.json` — coarse reaction typing: Suzuki, Buchwald–Hartwig amination, amide formation

### Reactivity computed layers (POC)

These are RDKit-graph-derived descriptors anchored on **templated Ar-* features** (Ar ipso atom is atom-map `:1` in templated SMARTS like `[c:1][Br]`), intended as reusable “reactivity feature layers”.

- `chemtools/taxonomy/v1/specs/reactivity_features.computed.v1.json` — sterics POC v1 spec (ortho substitution counts only)
- `chemtools/taxonomy/v1/docs/REACTIVITY_STERICS_POC_v1_README.md` — sterics POC v1 notes + usage
- `chemtools/taxonomy/v1/specs/reactivity_features.computed.v2.json` — sterics v2 spec (ortho substitution + topological bulk → 0–10 steric score)
- `chemtools/taxonomy/v1/docs/REACTIVITY_STERICS_POC_v2_README.md` — sterics v2 notes + usage
- `chemtools/taxonomy/v1/specs/reactivity_features.electronics.computed.v1.json` — electronics POC v1 spec (Gasteiger delta charge → 0–10 electron-poor score)
- `chemtools/taxonomy/v1/docs/REACTIVITY_ELECTRONICS_POC_v1_README.md` — electronics POC v1 notes + usage

### POC CLI entrypoints

- Sterics v1/v2: `reactivity_sterics_poc_v1.py`, `reactivity_sterics_poc_v2.py`
- Electronics v1: `reactivity_electronics_poc_v1.py`

## Scope emphasis

- Cross-coupling: sp2 electrophiles (aryl/vinyl halides and pseudohalides) and organoboron partners.
- Amide formation: carboxylic acids or activated acyl electrophiles (acid chlorides/anhydrides) with amine N–H.

All files are UTF-8 JSON.
