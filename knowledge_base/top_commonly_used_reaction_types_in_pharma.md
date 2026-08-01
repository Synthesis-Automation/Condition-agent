# Top Commonly Used Reaction Types in the Pharmaceutical Industry

## Scope

This document provides a practical **top-44 reaction repertoire** for medicinal chemistry and pharmaceutical process development. It is based primarily on the Roche enterprise electronic laboratory notebook (ELN) analysis covering reactions executed from 2010 to 2024, supplemented by major surveys of medicinal chemistry publications, pharmaceutical patents, process-development routes, and published synthetic routes.

The ranking is a **cross-source consensus**, not a pooled numerical league table. The underlying studies measure different things:

- Roche counts executed single-batch ELN records, including incomplete and failed experiments.
- Medicinal chemistry surveys count literature reaction steps or papers containing a reaction.
- Patent studies mine disclosed pharmaceutical transformations.
- Process chemistry surveys emphasize scalable routes to drug candidates.

Consequently, ranks 1–15 are highly stable, ranks 16–31 are broadly common but project-dependent, and ranks 32–44 represent an important industrial long tail.

> **Important:** The Roche study is a 2026 ChemRxiv preprint and had not been peer reviewed at the time of posting.

---

## Principal finding from the Roche ELN study

The Roche study processed 688,659 structurally valid reaction records and classified 560,260 reactions. Four top-level classes account for approximately **71%** of classified reactions:

| Roche top-level class | Share of classified reactions |
| --- | ---: |
| Heteroatom alkylation and arylation | 24.9% |
| Acylation and related processes | 18.8% |
| Deprotection reactions | 15.3% |
| Carbon–carbon bond formation | 12.0% |

The operational repertoire is therefore concentrated around amide formation, N-arylation, deprotection, Suzuki coupling, alkylation, ether formation, reductive amination, and a relatively limited group of enabling functional-group transformations.

---

## Evidence codes

| Code | Source type |
| --- | --- |
| **R** | Roche enterprise ELN study, 2010–2024 |
| **C** | Carey et al., process chemistry survey at GSK, AstraZeneca, and Pfizer |
| **RJ** | Roughley and Jordan medicinal chemistry toolbox analysis |
| **BB** | Brown and Boström analysis of medicinal chemistry methods |
| **P** | Schneider et al. pharmaceutical patent mining |
| **G** | Genheden and Howell analysis of 2.4 million published reactions |
| **HTE** | Ahlbrecht et al. analysis of 66,000 HTE reactions on drug-like molecules |

---

# Tier 1 — Core pharmaceutical reaction repertoire

## 1. Amide bond formation / amidation

- **General transformation:** `RCO2H + R′NH2 → RCONHR′`
- **Typical variants:** HATU-, T3P-, EDC-, CDI-, acid chloride-, or mixed-anhydride-mediated coupling.
- **Why common:** Reliable union of commercially available acid and amine building blocks; especially important for analogue-library synthesis.
- **Roche signal:** N-acylation to amide = **14.0%**; direct carboxylic acid + amine condensation = **9.6%** of all classified reactions.
- **Evidence:** R, C, RJ, BB, P, G

## 2. SNAr amination

- **General transformation:** activated `Ar–F/Cl + HNR2 → Ar–NR2`
- **Typical substrates:** electron-deficient heteroaryl chlorides and fluorides.
- **Why common:** Catalyst-free, operationally simple, and highly effective for heteroaromatic medicinal chemistry.
- **Roche signal:** Part of N-arylation with Ar–X, which accounts for **10.0%**.
- **Evidence:** R, RJ, BB, P

## 3. Amine deprotection

- **General transformation:** protected amine → `N–H`
- **Typical variants:** N-Boc, N-Cbz, N-Fmoc, N-benzyl, and N-SEM deprotection.
- **Why common:** Amine masking is frequently required during multistep synthesis and late-stage diversification.
- **Roche signal:** NH deprotection = **9.5%**; N-Boc deprotection alone = **6.8%**.
- **Evidence:** R, C, RJ, BB

## 4. Suzuki–Miyaura coupling

- **General transformation:** `Ar–X + Ar′–B(OR)2 → Ar–Ar′`
- **Typical variants:** aryl bromide, chloride, or iodide with boronic acid, BPin, BF3K, or related organoboron partner.
- **Why common:** Modular access to biaryl and heterobiaryl structures from widely available building blocks.
- **Roche signal:** Suzuki coupling = **6.2%** of classified reactions; its share of Roche C–C bond formations increased from approximately 42% in 2010 to more than 65% in 2024.
- **Evidence:** R, RJ, BB, P, G, HTE

## 5. N-Alkylation

- **General transformation:** `R2NH + R′–X → R2N–R′`
- **Typical variants:** alkyl halide, mesylate, tosylate, or triflate substitution; heteroaryl N-alkylation; N-methylation.
- **Why common:** Direct installation of side chains on amines and N-heterocycles.
- **Roche signal:** Heteroaryl N-alkylation = **3.4%**; N-substitution with alkyl-X = **1.4%**.
- **Evidence:** R, C, RJ, BB, P

## 6. Reductive amination

- **General transformation:** aldehyde or ketone + amine + reducing agent → substituted amine.
- **Typical reagents:** NaBH(OAc)3, NaBH3CN, NaBH4, borane reagents, or catalytic hydrogenation.
- **Why common:** Rapid formation of C–N bonds with broad building-block availability and limited protecting-group manipulation.
- **Roche signal:** **3.7%**.
- **Evidence:** R, C, RJ, BB, P

## 7. Ester hydrolysis / carboxylate deprotection

- **General transformation:** `RCO2R′ → RCO2H`
- **Typical variants:** methyl, ethyl, tert-butyl, benzyl, and other ester cleavage.
- **Why common:** Generates carboxylic acids for subsequent amidation or salt formation.
- **Roche signal:** RCO2H deprotection = **3.4%**.
- **Evidence:** R, C, RJ, BB

## 8. N-Containing heterocycle formation

- **General transformation:** cyclization or cyclocondensation → N-heterocycle.
- **Typical products:** pyrazoles, imidazoles, triazoles, pyridines, pyrimidines, benzimidazoles, and related ring systems.
- **Why common:** Nitrogen heterocycles are highly prevalent in drug-like compounds.
- **Roche signal:** All heterocycle-forming reactions = **3.9%**; individually classified N-heterocycle formation = **0.7%**, with many ring-forming reactions distributed among smaller or unclassified subclasses.
- **Evidence:** R, C, RJ, BB, P, G

## 9. Williamson ether synthesis / O-alkylation

- **General transformation:** `ROH/ArOH + R′–X → ROR′`
- **Why common:** Simple introduction of alkoxy, linker, and solubilizing substituents.
- **Roche signal:** Williamson ether synthesis = **1.4%**; the broader O-substitution class = **4.8%**.
- **Evidence:** R, C, RJ, BB

## 10. Buchwald–Hartwig amination

- **General transformation:** `Ar–X + HNR2 → Ar–NR2`, Pd-catalyzed.
- **Why common:** Extends aryl amination beyond electron-deficient substrates suitable for SNAr.
- **Roche signal:** Included within the 10.0% N-arylation family; explicitly identified bromo Buchwald–Hartwig amination = **0.5%**, although broader classifier assignments may be distributed among generic N-arylation labels.
- **Evidence:** R, RJ, BB, P, HTE

## 11. SNAr etherification

- **General transformation:** activated `Ar–F/Cl + ROH → Ar–OR`
- **Why common:** Efficient formation of aryl and heteroaryl ethers without palladium catalysis.
- **Roche signal:** **1.4%**.
- **Evidence:** R, BB, P

## 12. Sulfonamide formation

- **General transformation:** `RSO2Cl + HNR′2 → RSO2NR′2`
- **Why common:** Sulfonamides are common pharmacophores and can substantially modify acidity, polarity, and conformation.
- **Roche signal:** N-sulfonylation = **1.5%**.
- **Evidence:** R, C, RJ, BB, P

## 13. Urea formation

- **General transformation:** isocyanate or activated carbamoyl reagent + amine → urea.
- **Why common:** Ureas provide two hydrogen-bond vectors and are useful rigid amide bioisosteres.
- **Roche signal:** N-acylation to urea = **1.0%**.
- **Evidence:** R, C, RJ, BB, P

## 14. Amine protection

- **General transformation:** amine → N-Boc, N-Cbz, N-Fmoc, or another protected amine.
- **Why common:** Controls chemoselectivity in multistep routes.
- **Roche signal:** NH protection = **1.4%**; N-Boc protection = **0.9%**.
- **Evidence:** R, C, RJ, BB

## 15. Alkene hydrogenation

- **General transformation:** `C=C → C–C`
- **Typical methods:** H2 with Pd/C, Pt, Rh, or other homogeneous/heterogeneous catalysts; transfer hydrogenation in selected cases.
- **Why common:** Robust reduction and route simplification; particularly important in development chemistry.
- **Roche signal:** **1.3%**.
- **Evidence:** R, C, RJ, BB, G

---

# Tier 2 — Routine enabling and functional-group transformations

## 16. Nitro reduction to amine

- **General transformation:** `Ar–NO2 → Ar–NH2`
- **Roche signal:** Nitro-to-amine reduction = **1.3%**.
- **Evidence:** R, C, RJ, BB, P

## 17. Aldehyde or ketone reduction

- **General transformation:** `C=O → CH–OH`
- **Typical reagents:** NaBH4, LiAlH4, DIBAL-H, boranes, hydrosilanes, or catalytic hydrogenation.
- **Roche signal:** Ketone-to-alcohol reduction = **1.1%**, plus additional aldehyde reductions.
- **Evidence:** R, C, RJ, BB, G

## 18. Sonogashira coupling

- **General transformation:** `Ar–X + HC≡CR → Ar–C≡CR`
- **Roche signal:** **0.9%**.
- **Evidence:** R, RJ, BB, P, G

## 19. Alcohol deprotection

- **General transformation:** protected alcohol → `OH`
- **Typical variants:** O-TBS, O-benzyl, O-acetyl, and aryl methyl ether cleavage.
- **Roche signal:** ROH deprotection = **1.4%**.
- **Evidence:** R, C, RJ, BB

## 20. Halogenation

- **General transformation:** introduction of Br, I, Cl, or F.
- **Typical roles:** direct pharmacophore installation or preparation of a cross-coupling handle.
- **Roche signal:** Halogenation = **1.8%**, dominated by bromination.
- **Evidence:** R, C, RJ, BB, P

## 21. Alcohol oxidation

- **General transformation:** primary or secondary alcohol → aldehyde or ketone.
- **Typical methods:** Swern, Dess–Martin, TEMPO-based, MnO2, PCC/PDC, or catalytic aerobic oxidation.
- **Roche signal:** Alcohol oxidation family = approximately **1.0%**.
- **Evidence:** R, C, RJ, BB, G

## 22. Esterification

- **General transformation:** `RCO2H + R′OH → RCO2R′`
- **Typical methods:** Fischer esterification, coupling-reagent-mediated esterification, acid chloride, or transesterification.
- **Roche signal:** O-acylation to ester = **0.6%**; explicitly classified esterification = **0.3%**.
- **Evidence:** R, C, RJ, BB

## 23. Ester or carboxylic-acid reduction

- **General transformation:** `RCO2R′/RCO2H → RCH2OH`
- **Roche signal:** Ester-to-alcohol reduction = **0.7%**; carboxylic-acid-to-alcohol reduction is included in other reductions.
- **Evidence:** R, C, RJ, G

## 24. Alcohol protection

- **General transformation:** `ROH →` silyl ether, benzyl ether, ester, acetal, or related protected form.
- **Roche signal:** ROH protection = **0.4%**; O-TBS protection = **0.2%**.
- **Evidence:** R, C, RJ, BB

## 25. Carbamate or carbonate formation

- **General transformation:** amine or alcohol + chloroformate/carbonate reagent → carbamate or carbonate.
- **Roche signal:** Carbamate/carbonate formation = **0.5%**.
- **Evidence:** R, C, RJ, BB, P

## 26. O-Sulfonylation

- **General transformation:** `ROH → ROTs/ROMs/ROTf`
- **Why common:** Converts alcohols into leaving groups or generates aryl/vinyl triflates for coupling.
- **Roche signal:** Hydroxy-to-triflyloxy = **0.3%**; other sulfonates may be distributed among smaller subclasses.
- **Evidence:** R, C, RJ, BB, P

## 27. Alcohol-to-halide conversion

- **General transformation:** `ROH → RCl/RBr/RI`
- **Roche signal:** Alcohol-to-halide = **1.1%**.
- **Evidence:** R, C, RJ, P

## 28. Sulfide oxidation

- **General transformation:** sulfide → sulfoxide or sulfone.
- **Roche signal:** Oxidations at sulfur = **0.6%**; sulfanyl-to-sulfonyl = **0.5%**.
- **Evidence:** R, C, RJ, BB, P

## 29. Carboxylic acid activation to acid chloride

- **General transformation:** `RCO2H → RCOCl`
- **Typical reagents:** SOCl2, oxalyl chloride, PCl3, PCl5, or related chlorinating agents.
- **Roche signal:** **0.4%**.
- **Evidence:** R, C, RJ, P

## 30. Miyaura borylation / halogen-to-boronate conversion

- **General transformation:** `Ar–X → Ar–B(OR)2`
- **Why common:** Generates Suzuki-compatible boronate building blocks.
- **Roche signal:** Bromo Miyaura borylation = **0.5%**; bromo-to-pinacolatoboranyl = **0.2%**.
- **Evidence:** R, BB, P, G

## 31. Chiral resolution or stereoisomer separation

- **General operation:** racemate or diastereomeric mixture → isolated stereoisomer.
- **Typical methods:** chiral chromatography, crystallization of diastereomeric salts, enzymatic kinetic resolution, or preferential crystallization.
- **Roche signal:** Resolution reactions = **1.2%** overall; the proportional use is higher in downstream development than in medicinal chemistry.
- **Evidence:** R, C, RJ

> Chiral resolution is not a covalent transformation, but it is included because it is a frequent route-enabling operation in pharmaceutical synthesis and is treated as a formal class in the Roche ontology.

---

# Tier 3 — Common industrial long-tail reactions

## 32. Grignard or organolithium addition

- **General transformation:** organometallic reagent + carbonyl or another electrophile → C–C bond.
- **Typical use:** Core scaffold construction and alcohol synthesis.
- **Evidence:** C, RJ, BB, P, G

## 33. Wittig or Horner–Wadsworth–Emmons olefination

- **General transformation:** carbonyl compound → alkene.
- **Typical use:** Installation of terminal, internal, or conjugated alkenes.
- **Evidence:** C, RJ, BB, P, G

## 34. Enolate α-alkylation or α-acylation

- **General transformation:** carbonyl compound + electrophile → α-functionalized carbonyl.
- **Roche signal:** Keto α-alkylation and keto α-acylation are each individually identified at approximately **0.2%**.
- **Evidence:** R, C, RJ, BB, P

## 35. Aldol, Claisen, or Knoevenagel condensation

- **General transformation:** enolate or related carbanion + carbonyl partner → new C–C bond.
- **Roche signal:** Knoevenagel condensation = approximately **0.2%**.
- **Evidence:** R, C, RJ, BB, P, G

## 36. Mitsunobu reaction

- **General transformation:** `ROH + H–Nu → R–Nu`, usually with inversion at a stereogenic alcohol center.
- **Typical use:** Ether, ester, N-substitution, or S-substitution formation.
- **Roche signal:** Mitsunobu aryl ether synthesis = **0.4%**.
- **Evidence:** R, RJ, BB, P

## 37. Diazotization / Sandmeyer substitution

- **General transformation:** `ArNH2 → ArN2+ → Ar–X/CN/OH`
- **Typical use:** Aromatic functional-group interconversion.
- **Evidence:** C, RJ, BB, P

## 38. Cyanation / nitrile formation

- **General transformation:** `R/Ar–X → R/Ar–CN`
- **Typical use:** Direct nitrile pharmacophore installation or access to acids, amides, amidines, and amines.
- **Evidence:** C, RJ, BB, P, G

## 39. Friedel–Crafts alkylation or acylation

- **General transformation:** arene + carbon electrophile → substituted arene.
- **Typical use:** Core C–C bond formation, particularly in route-specific or historical processes.
- **Evidence:** C, RJ, BB, P, G

## 40. Heck coupling

- **General transformation:** `Ar–X + alkene → substituted alkene`
- **Typical use:** Aryl–alkene bond formation and annulation.
- **Evidence:** C, RJ, BB, P, G

## 41. Negishi or Kumada coupling

- **General transformation:** `Ar–X + R–Zn/R–Mg → Ar–R`
- **Typical use:** Suzuki alternatives for alkyl, heteroaryl, or otherwise difficult coupling partners.
- **Evidence:** C, RJ, BB, P, G

## 42. Stille coupling

- **General transformation:** `Ar–X + R–SnR3 → Ar–R`
- **Typical use:** Specialized coupling where organoboron chemistry is unsuitable; use is limited by organotin toxicity and removal requirements.
- **Evidence:** C, RJ, BB, P, G

## 43. Ullmann- or Chan–Lam-type C–N/C–O coupling

- **General transformation:** aryl partner + N/O nucleophile → aryl–heteroatom bond.
- **Typical use:** Pd-free thermal coupling or oxidative coupling of boronic acids with amines/alcohols.
- **Evidence:** R, RJ, BB, P, G

## 44. Olefin metathesis, especially ring-closing metathesis

- **General transformation:** alkene redistribution or intramolecular ring closure.
- **Typical use:** Cyclic scaffolds, macrocycles, and constrained linkers.
- **Evidence:** C, RJ, BB, P, G

---

# Additional reactions immediately below the top-44 set

Depending on therapeutic area and route type, the following may also be important:

- Epoxide or aziridine ring opening
- Nitrile hydrolysis
- Curtius, Hofmann, or Lossen rearrangement
- Diels–Alder and other cycloadditions
- Hydroboration–oxidation
- Asymmetric hydrogenation
- Biocatalytic ketone reduction
- Biocatalytic transamination
- C–H functionalization
- Photoredox coupling
- Click reactions, especially azide–alkyne cycloaddition

These are strategically important but generally occur less often than the core modular transformations in large medicinal chemistry datasets.

---

# Discovery chemistry versus process development

The Roche departmental analysis shows that a single global ranking can obscure different operational needs.

| Reaction superclass | Medicinal chemistry | Preclinical chemistry manufacturing | Pharma technical development |
| --- | ---: | ---: | ---: |
| Heteroatom alkylation and arylation | 25.8% | 18.8% | 24.6% |
| Acylation and related processes | 20.5% | 10.8% | 12.0% |
| Deprotection reactions | 16.1% | 12.9% | 10.5% |
| Carbon–carbon bond formation | 12.3% | 9.4% | 12.8% |
| Functional-group interconversion | 9.3% | 15.4% | 14.1% |
| Reductions | 4.2% | 15.3% | 13.0% |
| Resolution reactions | 0.9% | 2.8% | 2.5% |

## Practical interpretation

### Medicinal chemistry and library synthesis

Prioritize:

1. Amidation
2. SNAr and Buchwald–Hartwig amination
3. Suzuki coupling
4. N-alkylation
5. Reductive amination
6. Protection/deprotection
7. Ether formation
8. Sulfonamide and urea formation

These reactions maximize speed, modularity, commercial building-block usage, and analogue throughput.

### Process chemistry and CMC

Increase the priority of:

1. Catalytic and stereoselective reductions
2. Functional-group interconversions
3. Chemoselective oxidation-state adjustment
4. Hydrolysis and activation
5. Crystallization and resolution
6. Robust heterocycle formation
7. Route-convergent C–C bond formation
8. Scalable asymmetric synthesis

Process routes place greater weight on robustness, impurity control, safety, atom economy, cost, isolation, and stereochemical control than on rapid analogue generation.

---

# Reagent-level standardization observed at Roche

## Amidation

The Roche ELN analysis found strong convergence toward a limited operational protocol set:

- HATU appeared in approximately 28,300 N-acylation reactions.
- DIPEA appeared in approximately 38,900 reactions.
- HATU/DIPEA was the dominant coupling-agent/base pairing.
- HATU usage increased from below 30% of annotated amidations in 2010 to above 65% in 2024.

This is an observation of historical Roche practice, not a recommendation that HATU/DIPEA is optimal for every substrate or manufacturing route.

## Suzuki coupling

The corresponding Roche pattern was:

- PdCl2(dppf) was the dominant catalyst, driving approximately 30% of annotated Suzuki reactions overall.
- K2CO3, Na2CO3, and Cs2CO3 collectively covered approximately 80% of annotated base use.
- PdCl2(dppf) increased from approximately 20% of annotated catalyst use in 2010 to roughly 47% in 2024.
- XPhos Pd G3 increased to approximately 10% by 2024.

The later HTE study on 66,000 reactions reinforces an important qualification: historical frequency is not equivalent to universal substrate performance, and less common modern catalysts may outperform standard systems for difficult drug-like substrates.

---

# Yield and failure-data caution

The Roche study retained unfiltered operational records, unlike publication datasets that strongly favor successful reactions. Reported zero-or-unquantified-yield bins were:

- N-acylation to amide: **24%**
- Suzuki coupling: **32%**

However, the paper estimates that explicit true-zero reactions constitute approximately:

- **4%** of the amide subset
- **5%** of the Suzuki subset

The remainder of the nominal zero-yield bins largely reflects missing yield annotations. Therefore:

1. Reaction popularity should not be treated as proof of substrate compatibility.
2. Missing yields must not automatically be interpreted as chemical failure.
3. Condition-recommendation models should distinguish explicit zero yield, missing yield, low conversion, isolation failure, and analytical-only outcomes.

---

# Recommended taxonomy for reaction informatics and automation

For a practical reaction database, separate four layers:

1. **Net transformation** — the bonds broken and formed.
2. **Mechanistic or named reaction family** — Suzuki, SNAr, Buchwald–Hartwig, reductive amination, etc.
3. **Reacting-handle subtype** — Ar–Br versus Ar–Cl, boronic acid versus BPin, primary versus secondary amine, aldehyde versus ketone.
4. **Condition and process context** — catalyst, ligand, base, solvent, temperature, atmosphere, scale, concentration, addition order, and workup.

Examples:

- Keep **Suzuki–Miyaura coupling** as one main type; represent Ar–Br, Ar–Cl, Ar–I, boronic acid, BPin, and BF3K as substrate subtypes.
- Separate **SNAr amination** from **Buchwald–Hartwig amination**, even though both have the net pattern `Ar–X + N–H → Ar–N`.
- Separate **Williamson ether synthesis**, **SNAr etherification**, and **Mitsunobu substitution**.
- Use a parent class such as `amine_deprotection`, with Boc, Cbz, Fmoc, benzyl, and SEM as child labels or modifiers.
- Do not merge reaction frequency with success probability; store attempt count, explicit failure count, and outcome completeness separately.

For an initial pharmaceutical reaction recommendation system, ranks **1–30** should be mandatory primary classes. Ranks **31–44** are appropriate secondary classes for scaffold synthesis, process chemistry, and wider patent coverage.

---

# References

1. Schwarz, K.; Hilleke, M.; Boddy, A. J.; et al. **Charting the Evolving Chemical Synthesis Repertoire at Roche.** ChemRxiv, 2026. [https://doi.org/10.26434/chemrxiv.15002055/v1](https://doi.org/10.26434/chemrxiv.15002055/v1)
2. Carey, J. S.; Laffan, D.; Thomson, C.; Williams, M. T. **Analysis of the Reactions Used for the Preparation of Drug Candidate Molecules.** *Org. Biomol. Chem.* **2006**, *4*, 2337–2347. [https://doi.org/10.1039/B602413K](https://doi.org/10.1039/B602413K)
3. Roughley, S. D.; Jordan, A. M. **The Medicinal Chemist’s Toolbox: An Analysis of Reactions Used in the Pursuit of Drug Candidates.** *J. Med. Chem.* **2011**, *54*, 3451–3479. [https://doi.org/10.1021/jm200187y](https://doi.org/10.1021/jm200187y)
4. Brown, D. G.; Boström, J. **Analysis of Past and Present Synthetic Methodologies on Medicinal Chemistry: Where Have All the New Reactions Gone?** *J. Med. Chem.* **2016**, *59*, 4443–4458. [https://doi.org/10.1021/acs.jmedchem.5b01409](https://doi.org/10.1021/acs.jmedchem.5b01409)
5. Schneider, N.; Lowe, D. M.; Sayle, R. A.; Tarselli, M. A.; Landrum, G. A. **Big Data from Pharmaceutical Patents: A Computational Analysis of Medicinal Chemists’ Bread and Butter.** *J. Med. Chem.* **2016**, *59*, 4385–4402. [https://doi.org/10.1021/acs.jmedchem.6b00153](https://doi.org/10.1021/acs.jmedchem.6b00153)
6. Genheden, S.; Howell, G. P. **An Analysis of Published Synthetic Routes, Route Targets, and Reaction Types (2000–2020).** *Org. Process Res. Dev.* **2024**, *28*, 4225–4239. [https://doi.org/10.1021/acs.oprd.4c00389](https://doi.org/10.1021/acs.oprd.4c00389)
7. Ahlbrecht, J.; Lutz, M. D. R.; Jost, V.; et al. **Which Reaction Conditions Work on Drug-Like Molecules? Lessons from 66,000 High-Throughput Experiments.** *ACS Cent. Sci.* **2026**, *12*, 222–232. [https://doi.org/10.1021/acscentsci.5c02031](https://doi.org/10.1021/acscentsci.5c02031)

---

## Suggested machine-readable class set

```text
amide_formation
snar_amination
amine_deprotection
suzuki_miyaura_coupling
n_alkylation
reductive_amination
ester_hydrolysis
n_heterocycle_formation
williamson_ether_synthesis
buchwald_hartwig_amination
snar_etherification
sulfonamide_formation
urea_formation
amine_protection
alkene_hydrogenation
nitro_reduction
carbonyl_reduction
sonogashira_coupling
alcohol_deprotection
halogenation
alcohol_oxidation
esterification
ester_or_acid_reduction
alcohol_protection
carbamate_or_carbonate_formation
o_sulfonylation
alcohol_to_halide
sulfide_oxidation
acid_chloride_formation
miyaura_borylation
chiral_resolution
grignard_or_organolithium_addition
wittig_or_hwe_olefination
enolate_alpha_functionalization
aldol_claisen_or_knoevenagel
mitsunobu_reaction
diazotization_or_sandmeyer
cyanation
friedel_crafts_reaction
heck_coupling
negishi_or_kumada_coupling
stille_coupling
ullmann_or_chan_lam_coupling
olefin_metathesis
```
