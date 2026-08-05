# CAS tools

This standalone package provides CAS extraction and web-backed compound lookup.

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
