# Condition Registry

Standalone substance identity, role, and family registry for condition dataset
conversion and recommendation. It does not import `chemtools`.

Data ownership:

- `definitions/substances.v1.csv`: migrated flat substance registry
- `definitions/pending_substances.csv`: unresolved additions awaiting curation
- `definitions/roles_families.v1.json`: role and family taxonomy

Run the migration audit with:

```powershell
python -m condition_registry.cli results/condition_registry_migration
```
