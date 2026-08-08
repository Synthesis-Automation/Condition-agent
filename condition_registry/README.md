# Condition Registry

`condition_registry` owns condition-substance identity, aliases, possible roles,
contextual role resolution, and canonical condition recipes. It is standalone
and does not import `chemtools` or `condition_recommender`.

## Canonical definitions

- `definitions/substances.v2.jsonl`: the single active substance registry
- `definitions/roles.v2.json`: the small allowed-role vocabulary
- `definitions/role_resolution.v2.json`: contextual role and recipe-bucket rules
- `definitions/recipe_templates.v1.json`: expert recipe templates
- `definitions/synthesis_protocol.v1.schema.json`: canonical protocol JSON schema

`definitions/substances.v1.csv` remains only because legacy application paths
outside this package still read it. The clean registry runtime does not.

## Minimal substance format

Each line of `substances.v2.jsonl` is one independent JSON object:

```json
{"id":"cas:64-17-5","name":"Ethanol","cas":"64-17-5","smiles":"CCO","aliases":[{"type":"abbreviation","value":"EtOH"}],"roles":["solvent"]}
```

Only these fields are stored:

| Field | Required | Meaning |
| --- | --- | --- |
| `id` | yes | Stable registry identity; normally `cas:<CAS>` |
| `name` | yes | Canonical display name |
| `cas` | no | Valid CAS number |
| `smiles` | no | Canonical structure when known |
| `aliases` | no | Any number of typed aliases |
| `roles` | no | Possible intrinsic roles from `roles.v2.json` |

An alias contains `type` and `value`, with optional `language` and `shared`.
Set `shared: true` when the normalized alias intentionally belongs to
multiple substances. Generated identifier IDs, normalization profiles,
confidence defaults, and capability IDs are runtime details and are not stored.

Roles are multi-valued because one substance can serve different purposes in
different reactions. They are possibilities, not a claim about a particular
reaction. When no role is known, omit `roles`; do not add an `other` role.
Contextual resolution reports `unassigned`, `ambiguous`, or `conflicting`
instead of forcing a role.

Named reaction families do not belong in this registry. Reaction-family
evidence is owned by `reactive_taxonomy`, while reaction-specific role
preferences are owned by `condition_recommender`.

## Resolution

Resolution is exact and conservative. CAS values require a valid checksum.
Names resolve through canonical names and typed aliases; shared aliases return
an ambiguous result rather than selecting a record by order.

```python
from condition_registry import resolve_identifier

result = resolve_identifier("EtOH", identifier_type="abbreviation")
assert result.substance is not None
assert result.substance.substance_id == "cas:64-17-5"
```

```powershell
python -m condition_registry.cli resolve --cas "64-17-5"
python -m condition_registry.cli resolve --name "EtOH"
python -m condition_registry.cli validate
```

## Curation

Use the curator or the typed Python API. Writes replace the one JSONL file
atomically and verify that the new identity resolves before completing.

```powershell
python app/compound_registry_gui.py
```

```python
from condition_registry import CompoundAdditionRequest, add_compound

result = add_compound(
    CompoundAdditionRequest(
        canonical_name="Ethanol",
        cas="64-17-5",
        smiles="CCO",
    )
)
```

The curation API may validate or return derived physical information, but the
identity registry deliberately does not persist formula, molecular weight,
state, density, boiling point, melting point, curator notes, or generated IDs.

## Recipes

Resolved recipe IDs use the v2 condition contract: `RCR2` includes quantities,
process stages, and operating conditions; `RCORE2` identifies the role-aware
substance set without those variant details. Unresolved identities and roles
remain explicit in serialized recipe components.

Every resolved component includes `cas` when the registry has a CAS number.
The recommender serializes a `synthesis_protocol` beside each converted or
recommended recipe. Its materials contain stable registry IDs, canonical names,
CAS numbers, roles, and quantities; its operating conditions use explicit field
units (`temperature_c`, `time_h`, and `concentration_m`). Observed process stages
become `maintain_conditions` operations.

```json
{
  "protocol_id": "CP1:<sha256>",
  "recipe_id": "RCR2:<sha256>",
  "materials": [
    {
      "material_id": "condition_001",
      "category": "condition",
      "substance_id": "cas:584-08-7",
      "canonical_name": "Potassium carbonate",
      "cas": "584-08-7",
      "role": "base",
      "amount": 2.0,
      "amount_unit": "equiv"
    }
  ],
  "operating_conditions": {
    "temperature_c": 80.0,
    "time_h": 12.0,
    "concentration_m": null,
    "atmosphere": "nitrogen"
  },
  "execution_readiness": "review_required"
}
```

Dataset conditions normally lack ordered additions, vessel configuration,
mixing, quench, and workup. These gaps are listed in
`missing_required_fields`; the registry does not infer them or label such a
draft executable. A future procedure parser or reviewed template can populate
those fields before export to a hardware-specific execution format.

Run the migration audit without changing definitions:

```powershell
python -m condition_registry.cli audit results/condition_registry_migration
```

It writes `accepted.jsonl`, `quarantine.jsonl`, `issues.csv`, and
`migration_report.json`.
