# Condition Registry

Standalone substance identity, role, and family registry for condition dataset
conversion and recommendation. It does not import `chemtools`.

Data ownership:

- `definitions/substances.v1.csv`: migrated flat substance registry
- `definitions/pending_substances.csv`: unresolved additions awaiting curation
- `definitions/roles_families.v1.json`: role and family taxonomy
- `definitions/role_resolution.v1.json`: contextual role and recipe-bucket rules

## Contextual recipes

Build role-aware recipes without collapsing multi-role substances:

```python
from condition_registry import build_resolved_recipe

recipe = build_resolved_recipe(
    {
        "catalyst_cas": ("14221-01-3",),
        "reagent_cas": ("584-08-7",),
        "solvent_cas": ("7732-18-5",),
    },
    transformation_class="c_c_transfer_coupling",
)
```

The resulting `RCR1` identity is based on resolved substances, contextual roles,
operating conditions, and definition versions. Unresolved identities remain in
`other_components` with provenance and uncertainty warnings.

## CLI tester

Run the built-in smoke checks:

```powershell
python -m condition_registry.cli self-test
```

Resolve an exact CAS number, canonical name, or curated alias:

```powershell
python -m condition_registry.cli resolve --cas "14221-01-3"
python -m condition_registry.cli resolve --name "Tetrakis(triphenylphosphine)palladium"
python -m condition_registry.cli resolve --name "HCl" --format json
```

Resolution is deliberately conservative. CAS lookups require a valid CAS
checksum, and name lookups only use exact normalized canonical names and
aliases from `substances.v1.csv`. The command exits with status `0` when a
substance resolves and `1` for invalid, unresolved, or ambiguous queries.

Validate substance definitions against the package-owned role/family taxonomy:

```powershell
python -m condition_registry.cli validate
python -m condition_registry.cli validate --format json
```

Validation exits with status `1` when any row has an issue, which makes the
command suitable for CI. The detailed JSON output includes affected rows and
issue codes.

## Migration audit

Write accepted and quarantined rows plus issue reports to an output directory:

```powershell
python -m condition_registry.cli audit results/condition_registry_migration
```

The audit creates `accepted.csv`, `quarantine.csv`, `issues.csv`, and
`migration_report.json`. It is read-only with respect to registry definitions.
