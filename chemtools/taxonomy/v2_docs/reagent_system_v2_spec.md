# Reagent taxonomy v2 (minimal spec)

## Goal

Keep only the information required to (1) classify a reagent into a **family** and (2) assign a **role** for condition extraction and rule-matching.

This v2 keeps only the fields needed for deterministic classification. Legacy metadata
may be preserved in the taxonomy file but is not used by the classifier.

---

## Files

### 1) `reagent_roles.v2.json`

Single taxonomy file containing both roles and families:

```json
{
  "roles": [...],
  "families": [...]
}
```

#### `roles`

Array of role objects.

Required fields:

- `id` (string): stable machine id (e.g., `base`, `ligand`)
- `name` (string): display name
- `priority` (int): lower = more specific / preferred (used only for tie-breakers)
Optional:
- `default_family_id` (string): fallback family id when role is known but family is not
- `description` (string)

#### `families`

Array of family objects.

Required fields:

- `id` (string): stable machine id
- `name` (string): display name
- `role_id` (string): must match a role `id`
- `precedence` (int): lower = evaluated earlier / higher priority match
Optional:
- `description` (string)
- `detect` (object): structure-based detection
- `allowlists` (object): name/CAS-based detection (for non-structural materials or when SMILES is unavailable)

#### `detect` (minimal DSL)

Currently supported keys:

- `smarts.any` : list of SMARTS; match if **any** matches
- `smarts.none`: list of SMARTS; reject match if **any** matches

(Extend later only if necessary: `elements.any`, `predicates.any/all`, etc.)

#### `allowlists` (minimal)

Currently supported keys:

- `cas`      : list of CAS strings (exact match)
- `names`    : list of normalized names (exact match)
- `keywords` : list of tokens; match if **all** tokens appear in the normalized reagent name

(Use `cas`/`names` first; use `keywords` only as a fallback because it is fuzzier.)

---

## Classification algorithm (reference)

Input reagent record may contain: `{name?, cas?, smiles?}`.

1) **CAS exact match**

- If `cas` is present, match any family where `cas` is in `allowlists.cas`.

1) **SMILES structural match**

- If `smiles` is present, evaluate `detect.smarts`:
  - If `smarts.any` exists: must match at least one.
  - If `smarts.none` exists: must match none.

1) **Name-based fallback**

- Normalize `name` (lowercase, collapse whitespace, strip punctuation)
- If `allowlists.names` exists: exact match.
- Else if `allowlists.keywords` exists: require all tokens present.

1) **Tie-breaking / scoring**
Return all candidate families, then pick best by:

- Match strength: CAS > SMILES > name exact > keywords
- Then by `precedence` (lower wins)
- Then by role `priority` (lower wins)
- Then by `id` (stable deterministic)

1) **Role assignment**

- Role comes from chosen family’s `role_id`.
- If no family matched, return `role_id=None` (caller may fall back using other heuristics).

---

## Migration checklist (short)

1) Keep role IDs stable.
2) Convert old `include_smarts` / `exclude_smarts` into `detect.smarts.any/none`.
3) Keep `keywords` only where needed; prefer `allowlists.cas` as you build a registry.
4) Move solvent “properties” (polarity/proticity/etc.) into a separate lookup table if you need them later (not part of v2 core).
