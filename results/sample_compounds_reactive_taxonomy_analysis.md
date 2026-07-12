# Reactive Taxonomy Sample Analysis

- Source: `examples/sample_compounds.csv`
- Rows: 201
- Valid structures: 201/201
- Total sites: 261
- Runtime: 0.443 s (454.1 molecules/s)
- Taxonomy validation errors: 0

## Site families

| Site family | Count |
| --- | --- |
| leaving_group | 116 |
| pronucleophile_XH | 99 |
| electrophilic_center | 28 |
| transfer_group | 18 |

## Availability

| Availability | Count |
| --- | --- |
| available | 116 |
| free | 65 |
| deactivated | 34 |
| transferable | 18 |
| latent | 11 |
| activated | 11 |
| conditional | 6 |

## Coverage by declared sample role

| Role | Rows | Rows with sites | Rows without sites | Total sites |
| --- | --- | --- | --- | --- |
| electrophile | 101 | 101 | 0 | 140 |
| nucleophile | 69 | 67 | 2 | 74 |
| di-nucleophile | 12 | 12 | 0 | 24 |
| reference | 8 | 2 | 6 | 4 |
| bifunctional | 5 | 5 | 0 | 9 |
| di-electrophile | 3 | 3 | 0 | 10 |
| ligand | 2 | 0 | 2 | 0 |
| none | 1 | 0 | 1 | 0 |

## Most common matched patterns

| Pattern | Count |
| --- | --- |
| terminal_carbon_halogen | 109 |
| neutral_nh | 61 |
| neutral_oh | 26 |
| acyl_oxygen_or_halide | 23 |
| carbon_boronic_acid | 10 |
| carbon_sulfonate | 7 |
| terminal_alkyne | 7 |
| neutral_sh | 5 |
| sulfonyl_halide | 5 |
| carbon_zinc | 2 |
| carbon_tin | 2 |
| carbon_boronate | 1 |
| carbon_trifluoroborate | 1 |
| carbon_silicon | 1 |
| carbon_magnesium | 1 |

## Context tokens

| Context | Count |
| --- | --- |
| Ar | 118 |
| Alkyl | 66 |
| HeteroAr | 30 |
| C(O)R | 18 |
| N | 17 |
| Alkenyl | 10 |
| Alkynyl | 8 |
| C(O)OR | 7 |
| C(O)NR | 4 |
| SO2R | 3 |
| O | 1 |

## Rows without detected sites

| Name | Declared role | SMILES |
| --- | --- | --- |
| N,N-Dimethylaniline | nucleophile | c1ccc(N(C)C)cc1 |
| Thioanisole | none | CSc1ccccc1 |
| Phenylmagnesium bromide | nucleophile | [Mg]Brc1ccccc1 |
| Triphenylphosphine | ligand | c1ccc(P(c2ccccc2)c2ccccc2)cc1 |
| 1,8-Bis(dimethylphosphino)naphthalene | ligand | c1ccc2c(c1)c(P(C)C)ccc2P(C)C |
| Ethane | reference | CC |
| Anthracene | reference | c1ccc2cc3ccccc3cc2c1 |
| Coronene | reference | c1ccc2cc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67 |
| Spiro[5.5]undecane | reference | C1CCC2(CC1)CCCC2 |
| (1R,3R,5S)-1,3,5-Trimethylcyclohexane | reference | C[C@H]1CC[C@@H](C)C[C@H]1C |
| (1S,2R,5S)-1,2,5-Trimethylcyclohexane | reference | C[C@@H]1CCC[C@@H](C)[C@@H]1C |

## Audit flags

- Invalid structures: 0
- Sites missing pattern provenance: 0
- Sites with `Other` context: 0
- Rows with warnings: 0
- Multifunctional rows (>1 site): 48

## Representative multifunctional compounds

| Name | Site count | Labels |
| --- | --- | --- |
| 1,3,5-Tribromobenzene | 3 | Ar–Br; Ar–Br; Ar–Br |
| 1,4-Dichlorobenzene | 2 | Ar–Cl; Ar–Cl |
| 3-Bromopyrrole | 2 | HeteroAr–Br; AromN–H |
| 5-Bromoindole | 2 | Ar–Br; AromN–H |
| 2,5-Dichloropyridine | 2 | HeteroAr–Cl; HeteroAr–Cl |
| 3-Aminopyrrole | 2 | AromN–H; HeteroAr–NH2 |
| Trimethylsilylacetylene | 2 | R–C≡C–H; Alkynyl–SiR3 |
| 4-Bromophenylboronic acid | 2 | Ar–Br; Ar–B(OH)2 |
| 4-Bromoaniline | 2 | Ar–Br; Ar–NH2 |
| 4-Bromophenol | 2 | Ar–Br; Ar–OH |
| 4-Bromophenylacetylene | 2 | Ar–Br; R–C≡C–H |
| 1,4-Dibromobenzene | 2 | Ar–Br; Ar–Br |
| 1,4-Phenylenediamine | 2 | Ar–NH2; Ar–NH2 |
| 4-Bromo-N-Boc-aniline | 3 | Ar–Br; RO–C(O)–NHR; N–C(O)OR |
| Methyl 4-bromobenzoate | 2 | Ar–Br; Ar–C(O)OR |
| Pentafluorobromobenzene | 6 | Ar–F; Ar–F; Ar–F; Ar–Br; Ar–F; Ar–F |
| 1,3,5-Tribromobenzene | 3 | Ar–Br; Ar–Br; Ar–Br |
| Pentafluorobromobenzene | 6 | Ar–F; Ar–F; Ar–F; Ar–Br; Ar–F; Ar–F |
| 2,5-Dichloropyridine | 2 | HeteroAr–Cl; HeteroAr–Cl |
| 2,6-Dichloro-iodobenzene | 3 | Ar–I; Ar–Cl; Ar–Cl |
| 4-Iodo-N-Boc-aniline | 3 | Ar–I; RO–C(O)–NHR; N–C(O)OR |
| 4-Bromo-N-Cbz-aniline | 3 | Ar–Br; RO–C(O)–NHR; N–C(O)OR |
| 4-Chloro-N-Fmoc-aniline | 3 | Ar–Cl; RO–C(O)–NHR; N–C(O)OR |
| Phenylalanine | 3 | R–NH2; R–C(O)OH; R–C(O)–OH |
| (S)-Alanine | 3 | R–C(O)OH; R–NH2; R–C(O)–OH |

## Interpretation

- All structures parsed and all emitted sites retained SMARTS provenance.
- No emitted site used the unresolved `Other` context token.
- Rows without sites are mostly intentional v1 negatives: tertiary amines without N-H, ligands, hydrocarbons, and reference scaffolds.
- `Phenylmagnesium bromide` remains a fixture-quality issue because `[Mg]Brc1ccccc1` encodes an Mg-Br-C chain rather than a C-Mg bond.
- The output is site-complete, not reaction-role-selected: multifunctional compounds intentionally return every supported candidate site.
