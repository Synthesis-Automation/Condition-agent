# Reaction-Key → Reactant Signature Bridge

Purpose: allow recommendation to stay uniform when some datasets only have
reactant motifs (experiments/rules), while others have full `Reaction_Key`.

## Canonical projection

- **Reaction_Key** format (CRK):
  `|Reactants -> Product | bond_formed: ... | bond_broken: ... | spectators: ...`
- **Reactant signature (core)**: tokens from `Reactants` (the CRK summary)
- **Reactant signature (ext)**: tokens from `Reactants + spectators`

Signatures are encoded as a **sorted, deduped** `|`-joined list of motif IDs.
Counts are not encoded yet (multiplicity support is reserved for v2).

Example:

```
Reaction_Key: |Ar-Br|R-NH2 -> Ar-NR2 | bond_formed: C(ar)-N | bond_broken: Br-C(ar) | spectators: Pyridine
Core signature: Ar-Br|R-NH2
Ext signature: Ar-Br|Pyridine|R-NH2
```

## Dataset columns

The loader now populates:

- `Reactant_Signature_Core`
- `Reactant_Signature_Ext`

Rules:

- If `Reaction_Key` exists, signatures are derived from it.
- If `Reaction_Key` is missing, signatures are derived from
  `Reactant_A_Type / Reactant_B_Type / Reactant_C_Type`.

## Matching tiers (recommender)

1) **Exact Reaction_Key (CRK)** match  
2) **Core signature** match (reactant motifs only)  
3) **Ext signature** match (reactants + spectators)  
4) **Reaction type** fallback (if neither key nor signature is available)

This lets experiment/rule datasets participate in transformation-aware
matching without requiring products or full reaction SMILES.
