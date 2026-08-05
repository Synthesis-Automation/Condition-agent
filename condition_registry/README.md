# Condition Registry

Standalone substance identity, role, and family registry for condition dataset
conversion and recommendation. It does not import `chemtools`.

Data ownership:

- `definitions/substances.v1.csv`: migrated flat substance registry
- `definitions/substance_additions.v1.csv`: explicit, provenance-bearing new
  compounds created by the curation workflow
- `definitions/substance_identifiers.v1.csv`: typed one-to-many aliases and
  external identifiers with provenance
- `definitions/pending_substances.csv`: unresolved additions awaiting curation
- `definitions/roles_families.v1.json`: role and family taxonomy
- `definitions/role_resolution.v1.json`: contextual role and recipe-bucket rules
- `definitions/recipe_templates.v1.json`: typed expert recipe templates

## Substance identifiers and aliases

Each substance can have any number of typed identifiers. During the v1
migration, the loader converts the canonical `name`, `cas`, and `abbreviation`
columns in `substances.v1.csv` into identifier objects and merges additional
rows from `substance_identifiers.v1.csv`. New aliases belong in the identifier
file; do not add numbered alias columns or delimiter-separated aliases to the
substance file.
`Substance.aliases` remains as a computed compatibility view for current
callers; typed `Substance.identifiers` is the source-of-truth contract.

```csv
identifier_id,substance_id,identifier_type,value,language,is_preferred,source,confidence,status,normalization_profile,allow_ambiguous
sid:584-08-7:dipotassium-carbonate,cas:584-08-7,systematic_name,Dipotassium carbonate,en,false,condition_registry_curated,1.0,active,chemical_name_v1,false
```

Supported types are canonical, common, systematic, abbreviation, trade, and
legacy names, plus CAS, InChIKey, and external database identifiers. Alias
resolution remains exact and returns the matched identifier and provenance.
When a normalized alias belongs to multiple substances, resolution returns all
candidate substance IDs rather than choosing one. Such intentional collisions
must set `allow_ambiguous=true` on the supplemental definition row.

Typed resolution is available through the API and CLI:

```python
from condition_registry import resolve_identifier

result = resolve_identifier(
    "Dipotassium carbonate",
    identifier_type="systematic_name",
)
assert result.matched_identifier is not None
assert result.matched_identifier.source == "condition_registry_curated"
```

```powershell
python -m condition_registry.cli resolve --identifier "Dipotassium carbonate" --identifier-type systematic_name
```

## Adding compounds

Use the PyQt curator instead of editing definition CSV files directly:

```powershell
python app/compound_registry_gui.py
```

The form accepts a required CAS number, canonical name, and provenance source;
optional structure, physical properties, two curated role capabilities, and
any number of typed aliases can also be supplied. SMILES is canonicalized with
RDKit, formula and molecular weight are derived when omitted, and supplied
values are checked against the structure. Duplicate CAS numbers, names,
abbreviations, aliases, invalid role/family pairs, and malformed structures are
rejected before anything is written.

Select **Check CAS** to resolve the number against the local registry. If the
compound already exists, its identity, physical properties, roles, aliases,
and available provenance are loaded into the form and the save action changes
to **Validate and Save Changes**. The CAS field is locked while editing because
it determines the stable substance identity. Saving an older legacy record
migrates that record into the provenance-bearing additions definition; an
existing curated addition is updated in place.

Enter a valid CAS number and select **Auto-fill from Web** to query PubChem in
the background, with the NCI/CADD Chemical Identifier Resolver as a fallback.
When available, the form fills the preferred and systematic names,
abbreviation, SMILES, formula, molecular weight, InChIKey, PubChem CID,
physical state, density, boiling point, melting point, synonyms, and source
provenance. Existing condition roles and families are deliberately not inferred
from web names; the curator must review those chemistry annotations. PubChem
experimental physical values are visibly marked for review before saving.

Successful additions are written atomically to
`substance_additions.v1.csv`, with supplemental aliases written to
`substance_identifiers.v1.csv`. Both definitions are restored if either write
or post-write resolution verification fails. Updates atomically coordinate the
legacy, additions, and identifier definitions and restore all three if an edit
cannot be written and resolved successfully. The same workflow is available
programmatically with `add_compound()` and `update_compound()`:

```python
from condition_registry import (
    CompoundAdditionRequest,
    add_compound,
    update_compound,
)

result = add_compound(
    CompoundAdditionRequest(
        canonical_name="Ethanol",
        cas="64-17-5",
        smiles="CCO",
        source="doi:example",
    )
)
assert result.substance.substance_id == "cas:64-17-5"

updated = update_compound(
    result.substance.substance_id,
    CompoundAdditionRequest(
        canonical_name="Ethanol",
        cas="64-17-5",
        smiles="CCO",
        source="curator:reviewed",
    ),
)
```

## Expert recipe templates

Rule-based recommendation definitions reference registry-owned recipe templates
by stable ID. Every template option is an exact `substance_id`; free-text and
unresolved alternatives are rejected by validation. Alternatives are retained
as choices and are never implicitly expanded into a Cartesian product. Only
explicitly authored `variants` can be materialized as recipes.

```python
from condition_registry import get_recipe_template, materialize_recipe_variant

template = get_recipe_template("pd_sp2_cn.primary_alkyl_amine.v1")
assert template is not None
assert template.identity_complete

recipe = materialize_recipe_variant(
    template,
    "tbu_brettphos_k3po4_mecn.v1",
    transformation_class="sp2_c_n_substitution",
    include_draft=True,
)
assert recipe.recipe_id.startswith("RCR1:")
```

Template schema 1.2 stores catalyst loading and reagent equivalents on each
explicit selection, plus partner stoichiometry and operating parameters.
Activation validation requires all of those fields and a located primary
procedure. The primary-amide template is the first active protocol; the other
templates remain `draft`. The Phase I C-N drafts now contain explicit loadings,
equivalents, partner stoichiometry, time, concentration, temperature, and
atmosphere so they can be evaluated as complete screening recipes. They remain
inactive because their legacy-rulebook/expert basis still requires review
against located primary procedures. Draft status is therefore distinct from
both identity resolution and recipe-field completeness.

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

The resulting `RCR1` variant identity is based on resolved substances,
contextual roles, amounts, operating conditions, stages, and definition
versions. Its companion `RCORE1` identity contains only the role-aware
substances, so amounts and operating parameters do not fragment the underlying
condition regime. Unresolved identities remain in the role-hint bucket (or
`other_components` when no useful hint exists) with provenance and uncertainty
warnings. Typed `ConditionComponentInput` values can use exact CAS, name, or
stable substance-ID resolution; `ConditionProcessStage` preserves ordered
multi-stage time/temperature programs.

`CONDITION_RECIPE_COMPONENT_BUCKETS` is the single package-owned component
bucket vocabulary used when flattening, converting, and validating recipes.
`load_condition_vocabulary()` exposes immutable role and family IDs from the
versioned registry definition for downstream validators. Adding a bucket, role,
or family therefore has one registry-owned definition point and produces an
explicit validation failure in stale consumers.
`condition_registry_definition_versions()` exposes content-addressed
definition identities for restartable conversion manifests.

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
