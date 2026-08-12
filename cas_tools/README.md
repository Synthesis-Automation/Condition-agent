# CAS tools

This standalone package provides CAS extraction and web-backed compound lookup.

## Literature CAS/SMILES extractor

Launch the small desktop UI from the repository root:

```powershell
python .\cas_tools\cas_smiles_extractor_gui.py
```

The input defaults to `raw_dataset/literature_reaction_dataset`. The extractor
recursively reads every CSV, examines all JSON-valued columns (including
`reactants_json`, `products_json`, `reagents_json`, `catalysts_json`, and
`solvents_json`), and falls back to matching flat `*_cas`/`*_smiles` columns
when no structured role column exists. It writes a deterministic, deduplicated
CSV with exactly these columns:

```text
cas_no,compound_smiles,reaction_id,citation,source_role
```

The generated file contains exactly one row per CAS number. If a CAS number is
associated with conflicting structures, the extractor selects the SMILES
supported by the most distinct reaction records. Ties are resolved
deterministically. One reaction ID, literature citation, and observed source
role supporting the selected structure are retained. A reactant-role record is
preferred when the selected structure appeared in multiple roles so terminal
starting-material decisions do not silently rely on catalyst, solvent, reagent,
or product observations.

## Canonical molecule index

`molecule_index.py` converts any CSV SMILES column into a reusable SQLite
identity index. Atom maps and serialization differences are normalized to
canonical isomeric SMILES; a full InChIKey provides a second strict lookup key.
Invalid structures are excluded and counted. Selected source columns are
retained as bounded match provenance.

Build the literature index from the repository root:

```powershell
python -m cas_tools.molecule_index_cli `
  cas_tools/literature_cas_smiles_pairs.csv `
  results/literature_molecule_index.sqlite `
  --provenance-column cas_no `
  --provenance-column reaction_id `
  --provenance-column citation `
  --provenance-column source_role
```

The builder writes through a temporary database and atomically replaces the
requested output only after a successful build. Other CSV catalogs can select
a different structure column with `--smiles-column`.

```python
from cas_tools import CanonicalMoleculeIndex

with CanonicalMoleculeIndex("results/literature_molecule_index.sqlite") as index:
    match = index.lookup("C(C)O")
```

## Compound lookup

```python
from cas_tools import lookup_compound_by_cas

result = lookup_compound_by_cas("64-17-5")
if result.found:
    print(result.canonical_name, result.smiles, result.source_ids)
```

The lookup validates the CAS checksum before accessing the network. PubChem
PUG REST supplies compound identity, structure, formula, molecular weight,
InChIKey, CID, and synonyms. PubChem PUG View is queried for optional density,
boiling point, melting point, and physical-description data. Independent
optional requests run concurrently and have explicit timeouts. If PubChem does
not return core identity data, the NCI/CADD Chemical Identifier Resolver is
used as a fallback.

`CompoundLookupResult` is deliberately partial: unavailable endpoints produce
warnings without erasing fields returned by another endpoint. Source IDs and
URLs are retained for curation provenance. Experimental physical values are
representative web values and must be reviewed rather than treated as uniquely
authoritative. The lookup does not infer reaction-condition roles or families.
