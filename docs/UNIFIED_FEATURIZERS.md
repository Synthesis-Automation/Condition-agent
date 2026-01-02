# Unified Featurizers

This document describes the unified feature bundles used for both molecule and
reaction analysis. The unified featurizers live in `chemtools/featurizers/unified.py`.

## Entry points

- `featurize_molecule(smiles, registry_paths=None, options=None)`
- `featurize_reaction(reaction_smiles, registry_paths=None, options=None)`

## Molecule bundle (v1)

Shape (top-level):

```json
{
  "schema_version": "v1",
  "kind": "molecule",
  "molecule": { ... },
  "meta": { "rdkit_available": true, "errors": [] }
}
```

Key `molecule` fields:

- `smiles`: input SMILES
- `motifs`: motif hits from organic compound registry
- `steric`, `electronics`, `nearby`, `analyses`: motif-based analysis payloads
- `functional_group_map`, `functional_groups`: taxonomy functional group detection
- `rdkit_props`: core RDKit descriptors (MW, logP, TPSA, HBA/HBD, rings, etc.)

## Reaction bundle (v1)

Shape (top-level):

```json
{
  "schema_version": "v1",
  "kind": "reaction",
  "reaction": { ... },
  "meta": { "rdkit_available": true, "errors": [] }
}
```

Key `reaction` fields:

- `reaction_smiles`: input reaction SMILES
- `normalized`: normalized reaction payload (reactants, agents, products)
- `detection`: reaction-type detection matches
- `reaction_type`: best match summary with confidence
- `reactants`: list of per-reactant molecule bundles
- `aggregates`: reaction-level aggregates (max sterics, avg electronics, motifs, etc.)
- `roles`: reactant roles and expectedness summary from reaction context
- `agent_roles`: reagent-role classification for agents (catalyst, solvent, etc.)

### Reaction feature spec (roles + catalysts + solvent flags)

`roles` is reactant-focused (electrophile/nucleophile/coupling partner) and is
computed from the reaction-type context.

`agent_roles` is agent-focused and uses the reagent taxonomy v2 to classify
agents by SMARTS when possible. It includes:

- `agents`: entries with `smiles`, optional `role_id`, `family_id`, `match_kind`
- `role_counts`, `family_counts`: frequency by role and family
- `role_flags`: per-role booleans (ligand, metal_catalyst, base, solvent, etc.)
- `flags`: convenience booleans for condition roles, including:
  - `has_catalyst`, `has_metal_catalyst`, `has_organo_catalyst`
  - `has_ligand`, `has_base`, `has_acid`
  - `has_solvent`, `has_additive`
  - `has_oxidant`, `has_reductant`, `has_condensation_agent`

## Options

`options` is a dict with feature toggles:

- `include_rdkit_props` (default: True)
- `include_functional_groups` (default: True)
- `include_roles` (default: True)
- `include_agent_roles` (default: True)
- `max_hits_per_compound` (optional cap for motif detection)

## Notes

- RDKit is required for SMARTS-based detection and descriptor calculation.
- The unified bundles are designed to be backward-friendly with the existing
  motif and reaction detection layers.
