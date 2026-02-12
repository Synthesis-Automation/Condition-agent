# Oxidation Taxonomy Progress

Date: 2026-02-12

## Scope

Refactor oxidation reaction typing to use specific, chemistry-defined reactant/product classes (in the style of `Reduction_carbonyl_to_alcohol`), and expand coverage for common oxidation families.

## Completed

- Replaced the broad legacy alcohol-oxidation type with specific reaction types:
  - `Oxidation_primary_alcohol_to_aldehyde`
  - `Oxidation_secondary_alcohol_to_ketone`
  - `Oxidation_primary_alcohol_to_carboxylic_acid`
  - `Benzylic_CH_oxidation_to_carbonyl`
  - `Benzylic_CH_oxidation_to_carboxylic_acid`
  - `Oxidation_aldehyde_to_carboxylic_acid`
- Added common oxidation families:
  - `Baeyer_villiger_oxidation`
  - `Alkene_oxidative_cleavage_to_carbonyl`
  - `Alkene_oxidative_cleavage_to_carboxylic_acid`
  - `Oxidation_sulfide_to_sulfone`
- Split broad legacy Wacker typing into:
  - `Wacker_oxidation_terminal_alkene_to_ketone`
  - `Wacker_oxidation_internal_alkene_to_ketone`
- Tightened `Epoxidation` substrate scope to explicit alkene motif classes.
- Synced motif-set usage metadata in `compound_logic.json` for alcohol, carbonyl, carboxylic acid, and ester sets.
- Synced transformation behavior in `transformation_patterns.json` for all new oxidation IDs.

## Files Updated

- `chemtools/taxonomy/data/reaction_types.v4.0.json`
- `chemtools/taxonomy/data/compound_logic.json`
- `chemtools/taxonomy/data/transformation_patterns.json`

## Validation

- JSON parse check passed for all modified taxonomy files.
- `reaction_catalog.load_reaction_catalog()` loads successfully with all new oxidation IDs.

## Compatibility Note

- Legacy coarse oxidation IDs were removed from the reaction type registry.

## Next Steps

- Remap historical datasets using removed coarse oxidation IDs to the new specific IDs.
- Update any reporting dashboards or downstream scripts that assume legacy oxidation IDs.
- Add focused regression tests for:
  - primary alcohol -> aldehyde vs acid separation
  - benzylic C-H oxidation carbonyl vs acid separation
  - terminal vs internal Wacker split
  - oxidative alkene cleavage classification
