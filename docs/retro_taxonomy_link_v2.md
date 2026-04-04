# Retro Taxonomy Link v2 (Draft)

## Goal

Separate three concepts that are currently conflated in retro JSONs:

- `taxonomy family` (canonical reaction family in `reaction_types.v4.0.json`)
- `retro transform` (retrosynthetic disconnection logic identifier)
- `template/retron record name` (entry key in retro libraries)

## Problem in v1

`taxonomy_id` in the retro data files
`chemtools/retro/data/retron_patterns.json` /
`chemtools/retro/data/hte_templates.json` currently mixes:

- canonical taxonomy IDs (`Suzuki_miyaura`)
- taxonomy aliases (`Click_cuaac`)
- subtype labels that are broader/narrower than a family (`SNAr`)

This causes drift and makes overlap audits noisy.

## v2 Field Model (Additive / Backward-Compatible)

### `retron_patterns.json` entries

- `name`: stable retron record ID (existing)
- `taxonomy_id`: legacy link field (keep during migration)
- `taxonomy_family_id`: canonical reaction family ID (new)
- `reaction_name`: legacy retro transform label (keep during migration)
- `retro_transform_id`: explicit retro transform identifier (new, copied from `reaction_name`)

### `hte_templates.json` entries

- `name`: stable template record ID (existing)
- `taxonomy_id`: legacy link field (keep during migration)
- `taxonomy_family_id`: canonical reaction family ID (new)

## Migration Plan

1. Add new fields (`taxonomy_family_id`, `retro_transform_id`) without removing v1 fields.
2. Update runtime consumers to prefer `taxonomy_family_id` when present.
   Status: implemented for retro registry/loaders and key forward/retrosynthesis consumers.
3. Audit for ambiguous mappings (e.g. `SNAr`) and resolve explicitly.
4. Once consumers no longer depend on v1 semantics, deprecate `taxonomy_id` and `reaction_name`.

## Tooling

Use:

`python scripts/migrate_retro_json_links_v2.py`

By default it writes preview files (`*.v2.preview.json`) and does not modify source JSONs.
