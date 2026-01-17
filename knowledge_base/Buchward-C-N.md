---
id: Buchwald_Hartwig_CN
title: Buchwald-Hartwig C-N Coupling
tags: C-N, Buchwald, Pd, amination
scope: aryl/heteroaryl halide + amine
source: Internal notes; AbbVie HTE SI Tables S1, S4
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

- Internal notes (user-provided).
- AbbVie HTE SI: "The Design and Implementation of a High Throughput Experimentation Platform to Accelerate Drug Discovery" (Tables S1, S4).
