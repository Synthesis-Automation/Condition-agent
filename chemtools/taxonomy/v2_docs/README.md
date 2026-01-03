# Taxonomy v2 (Design and Usage)

This package hosts the new taxonomy system. Python modules live at
`chemtools/taxonomy/`, while data assets are separated into:

- `chemtools/taxonomy/data/`: JSON data used by v2 modules.
- `chemtools/taxonomy/v2_docs/`: Spec and design notes.
- `chemtools/archive/taxonomy/`: Legacy registry and data (reference only; not
  part of the v2 core).

The v2 system is intentionally small and deterministic, and includes:

- **Reagent taxonomy v2**: Minimal roles and families with allowlists and SMARTS
  detection (`reagent_roles.v2.json`, `reagent_families.v2_cas.json`, and the
  merged CAS pack).
- **Reaction catalog**: Canonical reaction IDs and aliases from
  `reaction_types.v4.0.json`.
- **Organic groups and compounds**: Template-based motif definitions from
  `organic_groups.v1.2.json` and `organic_compounds.v1.2.json` (plus
  `smarts_templates.v1.json`).
- **Motif registry + detection**: Compiles the organic group/compound templates
  into motif SMARTS and uses them for reaction classification.
- **Steric and electronic factor analysis**: Motif-driven steric and electronic
  descriptors for aryl and alkyl systems.

RDKit is optional at runtime; SMARTS and motif detection degrade gracefully when
it is unavailable.

## Core Modules

### Reagent taxonomy v2

- `chemtools/taxonomy/reagent_v2.py`
  - Loads v2 roles/families and merges the CAS pack into allowlists.
  - Implements the v2 classification algorithm (CAS > SMARTS > name/keywords).

Example:

```python
from chemtools.taxonomy.reagent_v2 import ReagentTaxonomyV2

taxonomy = ReagentTaxonomyV2.from_path()  # defaults to taxonomy/data
match = taxonomy.classify({
    "name": "triethylamine",
    "cas": "121-44-8",
    "smiles": "CCN(CC)CC",
})
print(match)  # -> family_id/role_id + match metadata
```

### Reaction catalog

- `chemtools/taxonomy/reaction_catalog.py`
  - Loads canonical reaction types and alias map.

Example:

```python
from chemtools.taxonomy.reaction_catalog import resolve_reaction_type

reaction_id = resolve_reaction_type("Suzuki coupling")
```

### Reaction detection (motif-based)

- `chemtools/featurizers/reaction_detection.py`
  - Uses motif detection to map reactant SMILES to reaction types.

Example:

```python
from chemtools.featurizers.reaction_detection import detect_reaction_types

result = detect_reaction_types("Brc1ccccc1.O=B(O)O[Na]>>")
print(result.to_dict())
```

### Organic groups, compounds, and motif registry

- `chemtools/featurizers/motif_registry.py`
  - Loads organic groups + compounds + SMARTS templates, then compiles the
    motif registry used across detection and analysis.

### Steric and electronic factor analysis

- `chemtools/featurizers/molecule.py`
  - Computes motif hits and steric/electronic descriptors for aryl and alkyl
    systems using the organic group/compound registry.

Example:

```python
from chemtools.featurizers.molecule import featurize_molecule

analysis = featurize_molecule("c1ccccc1Br")
```

### Rule DB lookup

- `chemtools/taxonomy/rule_db.py`
  - Resolves a reaction family label to a rule DB identifier using v2 metadata.

Example:

```python
from chemtools.taxonomy.rule_db import resolve_rule_db_v2

db_id = resolve_rule_db_v2("Suzuki")
```

## Data Files (taxonomy/data)

- `reagent_roles.v2.json`: Role definitions (id, priority, default family).
- `reagent_families.v2_cas.json`: Families with allowlists and optional SMARTS.
- `reagent_cas_pack_solvent_base.v2.json`: CAS pack merged into allowlists.
- `reaction_types.v4.0.json`: Canonical taxonomy for reaction types.
- `organic_groups.v1.2.json`, `organic_compounds.v1.2.json`,
  `smarts_templates.v1.json`: Motif registry sources.

## Legacy Archive

Legacy registry modules and data are preserved under
`chemtools/archive/taxonomy/`. New code should use v2 modules instead of the
archived registry. For legacy tooling, import from `chemtools.archive.taxonomy`
explicitly.
