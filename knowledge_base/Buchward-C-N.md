---
id: Buchwald_Hartwig_CN
title: Buchwald-Hartwig C-N Coupling
tags: C-N, Buchwald, Pd, amination
scope: aryl/heteroaryl halide + amine
source: Internal notes; AbbVie HTE Perspective Table 4; AbbVie HTE SI Tables S1, S4
---

# Reaction Family: Buchwald_Hartwig_CN - Buchwald-Hartwig C-N Coupling

## Summary

- Pd-catalyzed C-N bond formation between aryl halides (or pseudohalides) and amines.

## Applicability

- Substrates: aryl halides (often aryl bromides/chlorides) and amines (aryl or alkyl).
- Functional groups tolerated: depends on catalyst/ligand; tune selectivity with less active ligands when needed.
- Common incompatibilities: air sensitivity; scale-up can be difficult in dry apolar solvents without water.

## Recommended Conditions

| Catalyst | Ligand | Base | Solvent | Temp (C) | Time (h) | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Pd precatalyst (ligand-containing) | Xantphos if selectivity issues | K3PO4 (3 equiv) | apolar solvent + 5-10 equiv water | 80 | ~4 | 8 mol% catalyst, 10 vol (mL/g); avoid Pd2dba3 |

## Alternatives

- Use less active ligands (e.g., Xantphos) to improve selectivity when over-coupling occurs.
- Increase catalyst loading if conversion stalls after ~4 h.

## HTE Insights

- AbbVie HTE SI (Table S1) reports C-N coupling success rate: 67% (313 runs).

### Top conditions overall (AbbVie HTE Perspective Table 4)

| Rank | Catalyst | Ligand | Base | Solvent | Total Top-1 Hits | Total Top-3 Hits |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | RuPhos Pd G3 | RuPhos | LiHMDS | dioxane | 10 | 18 |
| 2 | BrettPhos Pd G3 | BrettPhos | Cs2CO3 | dioxane | 10 | 15 |
| 3 | BrettPhos Pd G3 | BrettPhos | LiHMDS | dioxane | 8 | 15 |
| 4 | BrettPhos Pd G3 | n/a | Cs2CO3 | dioxane | 5 | 15 |
| 5 | Pd2(dba)3 | Me4tBuXPhos | LiHMDS | t-amylOH | 7 | 14 |
| 6 | BrettPhos Pd G3 | RuPhos | K3PO4 | t-amylOH/DMA | 4 | 11 |
| 7 | BrettPhos Pd G3 | BrettPhos | K3PO4 | t-amylOH | 2 | 10 |
| 8 | BINAP Pd G3 | n/a | Cs2CO3 | dioxane | 6 | 9 |
| 9 | RuPhos Pd G3 | n/a | Cs2CO3 | dioxane | 3 | 9 |
| 10 | XPhos Pd G3 | n/a | Cs2CO3 | dioxane | 5 | 8 |

Notes:

- Subscripts restored from the PDF text (e.g., Cs2CO3, K3PO4).
- Pd2(dba)3 and Me4tBuXPhos are reported without subscripts in the PDF text.
- Ligands are listed as shown in the PDF; some rows omit a ligand entry.

- Top conditions for 5-membered heterocyclic halides (Table S4):

| Catalyst | Ligand | Base | Solvent | In Top 10 |
| --- | --- | --- | --- | --- |
| RuPhos Pd G3 | RuPhos | LiHMDS | dioxane | Yes |
| Pd2(dba)3 | Me4tBuXPhos | LiHMDS | t-amyl alcohol | Yes |
| BrettPhos Pd G3 | BrettPhos | Cs2CO3 | dioxane | Yes |
| tBuXPhos Pd G3 | tBuXPhos | P2-Et | DMSO | No |
| TrixiePhos Pd G3 | TrixiePhos | Cs2CO3 | dioxane | No |

## Safety / Handling

- Air-sensitive: degas the reaction mixture before adding base to the catalyst system.

## References

- doi: 10.1021/acs.jmedchem.5c00814
- AbbVie HTE SI: "The Design and Implementation of a High Throughput Experimentation Platform to Accelerate Drug Discovery" (Tables S1, S4).
