# Role-Aware Featurization – Implementation Plan

## Scope & Goals
- Build a **role-aware** featurizer for amines, alcohols, and aryl halides that outputs a **single consistent vector** per reactant with presence masks.
- Stack three layers: **global descriptors**, **role-specific features**, **reactive-site centered ECFP**.

## Inputs / Outputs
- **Input:** RDKit `Mol` (optionally reaction SMILES with atom-maps).
- **Output:** `{ "vector": np.ndarray, "fields": List[str], "masks": Dict[str, int], "meta": Dict }`.

## Directory Structure
```
chem_feats/
  __init__.py
  registry.yaml          # feature registry & defaults
  smarts.py              # role SMARTS & Finders
  global_feats.py        # role-agnostic RDKit descriptors
  role_feats/
    amine.py
    alcohol.py
    aryl_halide.py
  fingerprints.py        # centered ECFP helpers
  assemble.py            # concat, order, masks, sentinel fill
  tests/
    test_find_roles.py
    test_vectors_shape.py
    test_reference_cases.py
```

## Feature Registry (minimal example)
```yaml
global:
  - {name: global.MW, type: float, default: 0.0}
  - {name: global.logP, type: float, default: 0.0}
amine:
  - {name: amine.present, type: int, default: 0}
  - {name: amine.class_ps3, type: int, default: -1}
alcohol:
  - {name: alcohol.present, type: int, default: 0}
aryl_halide:
  - {name: aryl_halide.present, type: int, default: 0}
fingerprints:
  amine: {bits: 512, radius: 2}
  alcohol: {bits: 512, radius: 2}
  aryl_halide: {bits: 512, radius: 2}
```

## Core Steps
1. **Role detection**
   - Prefer atom-maps; else match SMARTS per role (amine N, alcohol O, aryl ipso C + X).
2. **Global features**
   - MW, logP, TPSA, RB, HBD/HBA, aromatic rings, heteroatom count.
3. **Role-specific features**
   - **Amines:** degree (1°/2°/3°), aniline flag, α-branching, formal charge.
   - **Alcohols:** 1°/2°/3°, benzylic/allylic flags.
   - **Aryl halides:** halide one-hot (F/Cl/Br/I), ortho-block (0–2), ipso degree.
4. **Reactive-site fingerprints**
   - Centered Morgan (radius 2–3) using `fromAtoms=[center_idxs]`.
5. **Assemble vector**
   - Fixed field order from `registry.yaml`; fill missing with sentinels; emit `has_amine/alcohol/aryl_halide` masks.
6. **Validation**
   - Shape consistency across roles; deterministic outputs; speed baseline.

## Public API (Python)
```python
from chem_feats import featurize_mol, featurize_reaction
vec = featurize_mol(mol, roles=["amine","aryl_halide"])
rxn_vecs = featurize_reaction(rxn_smiles)  # returns dict per role + masks
```

## Reference Test Cases
- Aryl halides: p-chloroanisole (ortho_block=0), 2,6-di-tBu-bromobenzene (ortho_block=2).
- Amines: aniline vs tert-butylamine (class + aniline flag).
- Alcohols: benzyl alcohol vs t-BuOH (benzylic/tertiary flags).
- Mixed: `Ar–Cl + aniline` and `Ar–Br + t-BuNH2`.

## Milestones
- **M1 (Day 1–2):** SMARTS finders + global features + registry scaffold.
- **M2 (Day 3–4):** Role-specific features + centered ECFP + assembler.
- **M3 (Day 5):** Tests, SHAP sanity on small model, integration into ConditionCore pipeline.

## Acceptance Criteria
- Fixed-length vectors for any input; masks correctly reflect roles.
- All tests pass; <5 ms per molecule on average; fields documented from registry.


## Extended Role Catalog: Additional Starting‑Material Types

> This section **supplements** the Core Roles with common nucleophiles, electrophiles, and radical precursors frequently encountered as *starting materials*. Use these to improve routing and role‑aware featurization.

### N‑Based Nucleophiles
| Role ID | Short Name | Canonical Examples | Seed SMARTS (indicative) | Key Role‑Specific Features | Typical Reaction Families |
|---|---|---|---|---|---|
| amine | Aliphatic amine (1°/2°/3°) | nBuNH2, i‑Pr2NH, Et3N | `[NX3;H2,H1,H0;!$(N=*)]` | degree(1/2/3); basicity_est; α‑branch; steric_map | Buchwald–Hartwig (aliphatic), SN2, reductive amination |
| aniline | Aromatic amine | aniline, p‑toluidine | `c[NH2,NH,NR]` | aniline_flag; para/meta directors; ortho_block | Buchwald–Hartwig (aryl) |
| sulfonamide | Sulfonamide N–H | TsNH2, NsNH2 | `S(=O)(=O)[NH]` | acidity_est; N‑substitution; hardness | Buchwald–Hartwig (sulfonamide), S_NAr |
| amide_anion | Deprotonated amide/carbamate | NaNHAc, BocNH− | `N-` with adjacent `C(=O)` | resonance_strength; charge; counterion | N‑acylation; metalations; amidations |
| hydrazine | Hydrazines | N2H4, BnNHNH2 | `NN` | degree; bifunctional_flag | Buchwald–Hartwig to hydrazines; Wolff–Kishner precursors |
| hydroxylamine | Hydroxylamine & derivatives | NH2OH, H2NOAc | `NO` with N‑H | oxidation_state; O‑acylated_flag | O/N‑acylations; oxime formation |
| imine_oxime | Imines/oximes | R2C=NH, R2C=NOH | `C=N[*]` | E/Z; electrophilicity_index | Reductive amination, additions |
| azole_NH | NH‑azole donors | imidazole, pyrazole | `n1[c,n]...N` (NH in azole) | pKa_est; ring_basicity | N‑arylation (SNAr/Buchwald), N‑acylation |

### C‑Nucleophiles (Polar)
| Role ID | Short Name | Canonical Examples | Seed SMARTS | Key Features | Typical Families |
|---|---|---|---|---|---|
| 1_3_dicarbonyl | 1,3‑Dicarbonyl donor | malonate, acetoacetate | `C(=O)C[CH2]C(=O)` | pKa_alpha; substituent_class; chelation_sites | Alkylation; Michael; Knoevenagel |
| nitroalkane | Nitroalkane (Henry) | nitromethane, nitroethane | `C[N+](=O)[O-]` | pKa_alpha; β‑branch | Henry; Michael |
| enolate_equiv | Silyl enol ether/enol carbonate | TMS‑enol ether; enol‑OTf/OPiv | `C=C[O;Si]` / `C(=C)OP(=O)` | E/Z; substitution; leaving_group | Mukaiyama aldol; α‑arylation |
| cyanide_source | CN donor | TMSCN, KCN | `C#N` (ionic or silyl) | counterion; safety_flag | Cyanation of carbocations/carbonyls |
| ylides_phosphorane | Phosphonium ylide | Ph3P=CH2 | `P(Ph)3=C` | ylide_basicity; substituent_class | Wittig/HWE olefination |

### Electrophiles Beyond Halides
| Role ID | Short Name | Canonical Examples | Seed SMARTS | Key Features | Typical Families |
|---|---|---|---|---|---|
| alkyl_sulfonate | Alkyl OTs/OMs/ONs | MeOTs, nBuOMs | `C-OS(=O)(=O)[Cl,F,O-]` | LG_class; degree_123 | SN1/SN2; Ni cross‑coupling (OTs) |
| enol_triflate | Enol/vinyl triflate | R2C=CR‑OTf | `C=C-OS(=O)(=O)C(F)(F)F` | E/Z; substitution_degree | Suzuki/Negishi; C–N |
| acyl_imidazole | Acyl‑Imidazole | R‑C(=O)‑N(imid) | `C(=O)N1C=NC=N1` | acyl_electrophilicity; stability | Amide/ester formation, decarboxylative couplings |
| activated_carbonate | p‑NP / PFP carbonates | RO‑C(=O)O‑Ar(NO2) | `C(=O)OAr` (strong EWG) | LG_strength | Carbamate/urethane formation |
| isocyanate | Isocyanate | R‑N=C=O | `N=C=O` | N‑electrophilicity; hazard_flag | Urea/Carbamate formation |
| isothiocyanate | Isothiocyanate | R‑N=C=S | `N=C=S` | S/N selectivity | Thiourea formation |

### Radical / Redox Precursors
| Role ID | Short Name | Canonical Examples | Seed SMARTS | Key Features | Typical Families |
|---|---|---|---|---|---|
| nhp_ester | N‑hydroxyphthalimide ester (RAE) | R‑C(=O)‑ONPhth | `C(=O)ON1C(=O)ccc(=O)c1` | radical_class(alkyl/benzyl/α‑oxy); decarbox_rate | Decarboxylative arylation/amination/borylation |
| pyridinium_salt | Katritzky pyridinium | R‑Py+ | `n1cccc1` with `C(sp3)` on N | alkyl_class; counterion | Alkyl radical generation; C–N/C–B |
| sulfonylate | Sulfinate salts | R‑SO2− M+ | `S(=O)(=O)[O-]` | radical_persistence; EWG_score | Minisci; aryl/alkenyl sulfones |
| diazo | α‑Diazocarbonyl | R‑CO‑CHN2 | `N#N` adjacent to carbonyl | carbene_type (donor/acceptor); safety_flag | Carbene insertions, cyclopropanations |
| oxime_ester | O‑Acyl oxime | R‑C(=NOAc)‑OAc | `C(=NOC(=O))` | N–O BDE; alkyl vs aryl | N‑centred or alkyl radical generation |

### C–H Partners (as Starting Materials)
| Role ID | Short Name | Canonical Examples | Seed SMARTS | Key Features | Typical Families |
|---|---|---|---|---|---|
| aryl_CH | Arene C–H partner | anisole, toluene | `cH` with directing group | dg_class (pyridyl, amide, oxazoline); ortho_sites | Direct arylation (C–H), borylation |
| benzylic_CH | Benzylic C–H | ethylbenzene | `c-CH2-` | BDE_est; radical_stability | Benzylic oxidation/amination |
| allylic_CH | Allylic C–H | cyclohexene | `C=C-CH2-` | BDE_est; substitution_degree | Allylic oxidation/amination |
| aliphatic_CH | Unactivated C–H | cyclohexane | `CH` without activation | BDE_est; site_count | Late‑stage C–H functionalization |

### Normalization & Assignment Hints
- **Amine disambiguation:** treat aniline as a distinct role (`aniline`) from aliphatic `amine`; prefer `aniline` in C(sp2)–N couplings.  
- **Acid derivatives:** map acyl‑imidazoles, NHS/PFP esters, and acid halides under a **superclass `acyl_electrophile`** for shared features (electrophilicity index, hydrolysis risk).  
- **Radical precursors:** ensure **mutual exclusivity** with corresponding acids/alkyl alcohols to avoid double‑counting; emit a `radical_source.present` mask.  
- **C–H roles:** only activate when the router predicts a C–H activation family or when no better electrophile/nucleophile is present.

