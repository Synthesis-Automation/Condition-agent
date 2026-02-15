# Substituent Fragments Editing Guide

This guide describes how to edit `chemtools/taxonomy/data/substituent_fragments.v1.json`
for the composite substituent system.

## Purpose

`substituent_fragments.v1.json` defines composable substituent fragments that generate
groups such as `-CONH2`, `-COOR`, `-SO2NH2` from:

- a `linker` template (`CO`, `SO2`, etc.)
- a `terminal_group` from `organic_groups.v1.3.json`

Generated groups are merged by `chemtools/taxonomy/substituent_composer.py`.

## File Contract

Top-level keys:

- `version`: string
- `notes`: string
- `linkers`: non-empty list
- `groups`: non-empty list

### `linkers[]` schema

Required keys:

- `label`: uppercase linker label (for example `CO`, `SO2`, `OCO`)
- `smarts_template`: SMARTS template that must include `{TAIL}`

Optional keys:

- `priority`: integer
- `description`: string
- `allowed_terminal_groups`: list of allowed terminal group ids
- `blocked_terminal_groups`: list of blocked terminal group ids

### `groups[]` schema

Required keys:

- `id`: substituent id that starts with `-` (for example `-CONH2`)
- `label`: uppercase linker label matching one linker entry
- `terminal_group`: existing group id from `organic_groups.v1.3.json`

Optional keys:

- `aliases`: list of non-empty strings
- `description`: string
- `priority`: integer

## Naming Rules

- Use uppercase linker labels: `CO`, `SO`, `SO2`, `OCO`.
- Prefer canonical ids such as `-COOH`, `-COOR` (legacy names can stay in `aliases`).
- Do not use deprecated `terminal`; use `terminal_group` only.
- Keep ids chemistry-first and stable across scaffolds.
- If both allow/block lists are used, they must not overlap.

## Edit Workflow

1. Edit `chemtools/taxonomy/data/substituent_fragments.v1.json`.
2. Run the checklist:

```bash
python scripts/check_substituent_fragments.py
```

3. If you want full suite coverage:

```bash
python scripts/check_substituent_fragments.py --full-tests
```

## Examples

Add a new linker:

```json
{
  "label": "PO2",
  "priority": 2,
  "smarts_template": "[PX4:2](=O)(O){TAIL}",
  "allowed_terminal_groups": ["-OH", "-OR", "-NH2", "-NHR", "-NR2"],
  "description": "Phosphoryl linker."
}
```

Add a new composed group using an existing terminal:

```json
{
  "id": "-PO2NH2",
  "label": "PO2",
  "terminal_group": "-NH2",
  "aliases": ["-PO2-NH2"],
  "description": "Phosphoramidate-like substituent."
}
```

## Validation Sources

- Schema validation: `validate_substituent_fragments_payload()` in
  `chemtools/taxonomy/substituent_composer.py`
- Taxonomy checks: `python -m chemtools.taxonomy.validate_and_sync --check-only`
- Composer tests: `tests/test_substituent_composer.py`
