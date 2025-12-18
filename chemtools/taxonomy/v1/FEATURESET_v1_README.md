# Organic chemistry SMARTS feature set (v1)

Generated: 2025-12-18

## Files

- `organic_groups.v1.json` — chemist-style group labels (Ar-, Vinyl-, -CHO, -CN, -NO2, -B(OH)2, -Bpin, -OTf, ...)
- `smarts_templates.v1.json` — template rules and explicit expansions (build-time)
- `calculable_features.atomic.v1.json` — atomic SMARTS detectors (runtime)
- `calculable_features.derived.v1.json` — boolean logic detectors (runtime)
- `calculable_features.compiled.v1.json` — combined atomic + derived
- `reactant_types.v1.json` — feature-token taxonomy (cross-coupling + amide formation focus)
- `reaction_types.v1.json` — coarse reaction typing: Suzuki, Buchwald–Hartwig amination, amide formation

## Scope emphasis

- Cross-coupling: sp2 electrophiles (aryl/vinyl halides and pseudohalides) and organoboron partners.
- Amide formation: carboxylic acids or activated acyl electrophiles (acid chlorides/anhydrides) with amine N–H.

All files are UTF-8 JSON.
