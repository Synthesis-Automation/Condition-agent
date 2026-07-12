# Suzuki recommendation pilot conversion

- Source rows: 300
- Verified: 181
- Review: 90
- Rejected: 29
- Mean verified yield: 80.42%
- Temperature coverage: 0/300
- Time coverage: 3/300
- Review-label coverage: 90/90

## Admission policy

Verified records require exact product reconstruction, a unique Suzuki–Miyaura
family environment, a yield from 0–100%, and non-empty catalyst, reagent, and
solvent CAS identity groups. Review records retain usable but incomplete or
unverified observations. Invalid, yield-invalid, and unresolved reactions are
rejected.

## Reason counts

- `exact_single_event_suzuki_with_conditions`: 181
- `incomplete_condition_identity`: 6
- `invalid_reaction_smiles`: 2
- `missing_suzuki_family_environment`: 83
- `multiple_products`: 5
- `not_verified_as_suzuki_miyaura`: 83
- `reaction_not_exactly_verified`: 82
- `unresolved_transformation`: 27
