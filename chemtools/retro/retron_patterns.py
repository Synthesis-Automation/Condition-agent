"""
Retron SMARTS Pattern Library.

A retron is a substructure in a target molecule that implies a specific
retrosynthetic transform — i.e., a bond whose formation defines the
"last step" in a synthetic route.

Each entry contains:
  name            : unique identifier
  product_smarts  : SMARTS that matches the bond/motif in the TARGET molecule
  reaction_name   : maps to the existing chemtools reaction taxonomy ID
  difficulty      : 0.0 (trivial) → 1.0 (heroic); guides ranking
  description     : human-readable retrosynthetic description
  notes           : chemistry notes and caveats
  precursor_hints : list of precursor type names (for LLM context)

SMARTS retron naming convention:
  [a:1] = aromatic atom with atom-map label
  [C:1] = aliphatic carbon with atom-map label
  The pattern identifies the bond that would be FORMED in the forward reaction.

~50 retrons covering 8 major reaction classes:
  1. C-C bond formation (cross-coupling, additions)
  2. C-N bond formation
  3. C-O bond formation
  4. Carbonyl-based (ester, amide, etc.)
  5. Ring-forming reactions
  6. Oxidation / reduction markers
  7. Protecting group markers
  8. Halogenation
"""
from __future__ import annotations

from typing import List, Dict, Any

# ---------------------------------------------------------------------------
# Retron definitions
# Each entry is a plain dict — intentionally not a dataclass so it can be
# serialized to JSON and used in LLM context directly.
# ---------------------------------------------------------------------------

RETRONS: List[Dict[str, Any]] = [

    # -----------------------------------------------------------------------
    # 1. C–C Bond Formation (Cross-Coupling)
    # -----------------------------------------------------------------------

    {
        "name": "biaryl_suzuki",
        "product_smarts": "[c:1]-[c:2]",
        "reaction_name": "suzuki_miyaura",
        "taxonomy_id": "Suzuki_miyaura",
        "difficulty": 0.15,
        "description": "Biaryl C–C bond ← aryl/heteroaryl halide + arylboronic acid (Suzuki)",
        "precursor_hints": ["aryl_halide", "aryl_boronic_acid"],
        "notes": "Highly reliable. Works for electron-rich and electron-poor arenes. "
                 "Pd/phosphine catalyst (SPhos, XPhos, dppf). Base: K2CO3, Cs2CO3, K3PO4. "
                 "Prefer Br/I over Cl. Watch for protodeboronation of electron-poor boronic acids.",
    },
    {
        "name": "aryl_alkyl_negishi",
        "product_smarts": "[c:1]-[CX4:2]",
        "reaction_name": "negishi_coupling",
        "taxonomy_id": "Negishi",
        "difficulty": 0.35,
        "description": "Aryl–alkyl C–C bond ← aryl halide + alkylzinc (Negishi)",
        "precursor_hints": ["aryl_halide", "alkyl_zinc"],
        "notes": "Pd or Ni catalysis. Organozinc reagents are air-sensitive. "
                 "β-hydride elimination can be a problem with primary alkyl chains.",
    },
    {
        "name": "aryl_alkene_heck",
        "product_smarts": "[c:1]-[CH:2]=[CH:3]",
        "reaction_name": "heck_reaction",
        "taxonomy_id": "Heck",
        "difficulty": 0.30,
        "description": "Aryl vinyl bond ← aryl halide + alkene (Heck/Mizoroki-Heck)",
        "precursor_hints": ["aryl_halide", "alkene"],
        "notes": "Pd catalysis. E-selectivity usually preferred. "
                 "Regioselectivity depends on alkene substitution. "
                 "DMF or DMA solvent, tertiary amine base.",
    },
    {
        "name": "aryl_alkyne_sonogashira",
        "product_smarts": "[c:1]-[C:2]#[C:3]",
        "reaction_name": "sonogashira_coupling",
        "taxonomy_id": "Sonogashira",
        "difficulty": 0.25,
        "description": "Aryl–alkyne bond ← aryl halide + terminal alkyne (Sonogashira)",
        "precursor_hints": ["aryl_halide", "terminal_alkyne"],
        "notes": "Pd/Cu co-catalysis in amine solvent (Et3N, iPr2NH). "
                 "Avoid Cu if alkyne dimerization (Glaser) is a concern. "
                 "TMS-protected alkyne → deprotect after coupling.",
    },
    {
        "name": "beta_hydroxy_carbonyl_aldol",
        "product_smarts": "[CX4:1]([OH1])-[CX4:2]-[C:3](=[O:4])",
        "reaction_name": "aldol_condensation",
        "taxonomy_id": "Aldol_addition",
        "difficulty": 0.40,
        "description": "β-Hydroxy carbonyl ← aldehyde/ketone + enolizable carbonyl (Aldol)",
        "precursor_hints": ["aldehyde", "enolizable_ketone"],
        "notes": "Classic C–C forming reaction. Mukaiyama aldol (TMS enol ether + Lewis acid) "
                 "for stereocontrol. Intramolecular aldol for ring formation. "
                 "Evans oxazolidinone auxiliaries for asymmetric.",
    },
    {
        "name": "alkene_wittig",
        "product_smarts": "[CX3:1]=[CX3:2]",
        "reaction_name": "wittig_reaction",
        "taxonomy_id": "Wittig_reaction",
        "difficulty": 0.30,
        "description": "Alkene ← carbonyl compound + phosphonium ylide (Wittig)",
        "precursor_hints": ["aldehyde_or_ketone", "phosphonium_ylide"],
        "notes": "Non-stabilized ylides give Z-alkene. Stabilized (ester, CN) give E. "
                 "Horner-Wadsworth-Emmons (phosphonate ester) is E-selective and more reactive. "
                 "Byproduct: Ph3P=O (add MgBr2 for HWE).",
    },
    {
        "name": "grignard_alcohol",
        "product_smarts": "[CX4:1]([OH1:2])[CX4:3]",
        "reaction_name": "grignard_reaction",
        "taxonomy_id": "Grignard_addition",
        "difficulty": 0.35,
        "description": "Secondary/tertiary alcohol ← carbonyl + Grignard (RMgX)",
        "precursor_hints": ["carbonyl_compound", "alkyl_or_aryl_halide_for_grignard"],
        "notes": "Anhydrous conditions. THF or Et2O solvent. Sensitive to moisture and air. "
                 "Incompatible with acidic protons, esters (double addition), nitriles. "
                 "Consider organolithium for hindered substrates.",
    },
    {
        "name": "alkyl_alkyl_kumada",
        "product_smarts": "[CX4:1]-[CX4:2]",
        "reaction_name": "kumada_coupling",
        "taxonomy_id": "Kumada",
        "difficulty": 0.45,
        "description": "Alkyl–alkyl C–C bond ← alkyl halide + alkyl Grignard (Kumada)",
        "precursor_hints": ["alkyl_halide", "alkyl_grignard"],
        "notes": "Ni or Pd catalysis. Less common than Suzuki/Negishi. "
                 "Grignard is an air-sensitive, basic reagent. "
                 "Used for unfunctionalized alkyl couplings.",
    },

    # -----------------------------------------------------------------------
    # 2. C–N Bond Formation
    # -----------------------------------------------------------------------

    {
        "name": "aryl_amine_buchwald",
        "product_smarts": "[c:1]-[NX3;H1,H0:2]",
        "reaction_name": "buchwald_hartwig_amination",
        "taxonomy_id": "C_N_Coupling",
        "difficulty": 0.25,
        "description": "Aryl amine ← aryl halide/triflate + amine (Buchwald-Hartwig)",
        "precursor_hints": ["aryl_halide_or_triflate", "primary_or_secondary_amine"],
        "notes": "Pd/bulky phosphine (BINAP, XPhos, DavePhos). Base: Cs2CO3, NaOtBu. "
                 "Works for primary amines (→ secondary arylamines) and secondary amines (→ tertiary). "
                 "Cu-based Ullmann conditions available as Pd-free alternative (higher T).",
    },
    {
        "name": "alpha_amino_reductive_amination",
        "product_smarts": "[CX4:1]([NX3:2])[CX4,CX3:3]",
        "reaction_name": "reductive_amination",
        "taxonomy_id": "Reductive_amination",
        "difficulty": 0.20,
        "description": "Amine ← aldehyde/ketone + amine + reductant (reductive amination)",
        "precursor_hints": ["aldehyde_or_ketone", "primary_or_secondary_amine"],
        "notes": "NaBH3CN or NaBH(OAc)3 in MeOH or DCE. "
                 "For ketones: molecular sieves + NaBH3CN. "
                 "Leuckart-Wallach: formic acid/formamide (no reductant). "
                 "Chiral: asymmetric transfer hydrogenation with Ir/Rh catalyst.",
    },
    {
        "name": "aryl_amine_ullmann",
        "product_smarts": "[c:1]-[NX3:2]-[c:3]",
        "reaction_name": "ullmann_coupling",
        "taxonomy_id": "C_N_Coupling",
        "difficulty": 0.40,
        "description": "Diaryl amine ← two aryl halides + amine/ammonia (Ullmann C–N)",
        "precursor_hints": ["aryl_halide_1", "aryl_halide_2"],
        "notes": "Cu catalysis (CuI, Cu2O). Higher temperature than Buchwald. "
                 "Ligand acceleration: 1,10-phenanthroline, DMEDA. "
                 "Suitable when Pd-free desired.",
    },
    {
        "name": "amide_direct",
        "product_smarts": "[CX3:1](=[O:2])-[NX3:3]",
        "reaction_name": "amide_coupling",
        "taxonomy_id": "Amide_formation",
        "difficulty": 0.20,
        "description": "Amide ← carboxylic acid + amine (amide bond formation)",
        "precursor_hints": ["carboxylic_acid", "amine"],
        "notes": "DCC/DMAP, HATU/HOAt, EDC/HOBt, or T3P coupling reagents. "
                 "Acyl chloride + amine (Schotten-Baumann) for scale-up. "
                 "Watch for epimerization at α-chiral center. "
                 "HATU preferred for hindered or acid-sensitive substrates.",
    },
    {
        "name": "nitrile_reduction_amine",
        "product_smarts": "[CX4:1]-[NH2:2]",
        "reaction_name": "nitrile_reduction",
        "taxonomy_id": "Reduction_nitrile_to_amine",
        "difficulty": 0.30,
        "description": "Primary amine ← nitrile + reductant (nitrile reduction)",
        "precursor_hints": ["nitrile"],
        "notes": "LiAlH4 (THF, 0°C→rt), BH3·THF, or H2/Raney Ni. "
                 "DIBAL-H gives imine intermediate (can be trapped). "
                 "Watch for over-reduction to secondary amine with LAH.",
    },

    # -----------------------------------------------------------------------
    # 3. C–O Bond Formation
    # -----------------------------------------------------------------------

    {
        "name": "aryl_ether_ullmann_o",
        "product_smarts": "[c:1]-[OX2:2]-[CX4,c:3]",
        "reaction_name": "ullmann_ether_synthesis",
        "taxonomy_id": "C_O_Coupling",
        "difficulty": 0.35,
        "description": "Aryl ether ← aryl halide + alcohol/phenol (Ullmann O-arylation)",
        "precursor_hints": ["aryl_halide", "alcohol_or_phenol"],
        "notes": "Cu or Pd catalysis. For phenols: K2CO3 base, DMF. "
                 "Chan-Lam coupling (Cu/O2, boronic acid + alcohol) at room temperature. "
                 "Pd/Xantphos for secondary alcohols.",
    },
    {
        "name": "williamson_ether",
        "product_smarts": "[CX4:1]-[OX2:2]-[CX4:3]",
        "reaction_name": "williamson_synthesis",
        "taxonomy_id": "Alkyl_Nucleophilic_Substitution",
        "difficulty": 0.15,
        "description": "Dialkyl ether ← alkoxide + alkyl halide (Williamson synthesis)",
        "precursor_hints": ["alkyl_halide", "alcohol"],
        "notes": "SN2 mechanism: requires primary or methyl electrophile. "
                 "Base: NaH, K2CO3, Ag2O. Solvent: DMF, THF. "
                 "Elimination competes with secondary/tertiary electrophiles.",
    },
    {
        "name": "mitsunobu_inversion",
        "product_smarts": "[CX4:1]([OH0:2])[CX4:3]",
        "reaction_name": "mitsunobu_reaction",
        "taxonomy_id": "Mitsunobu_reaction",
        "difficulty": 0.30,
        "description": "Ether/ester with inversion ← alcohol + acidic partner + PPh3/DIAD",
        "precursor_hints": ["primary_or_secondary_alcohol", "acidic_pronucleophile"],
        "notes": "Inverts stereochemistry at secondary alcohols. "
                 "DEAD or DIAD + PPh3. Partners: carboxylic acids, phenols, sulfonamides, azides. "
                 "Byproduct: Ph3P=O. Use DEAD-free protocols (cyanomethylene phosphoranes) for scale.",
    },
    {
        "name": "ester_from_acid_alcohol",
        "product_smarts": "[CX3:1](=[O:2])-[OX2:3]-[CX4:4]",
        "reaction_name": "esterification",
        "taxonomy_id": "Esterification",
        "difficulty": 0.15,
        "description": "Ester ← carboxylic acid + alcohol (Fischer esterification or coupling)",
        "precursor_hints": ["carboxylic_acid", "alcohol"],
        "notes": "Fischer esterification: H2SO4 or HCl catalyst, reflux, Dean-Stark. "
                 "DCC/DMAP or coupling reagents for sensitive substrates. "
                 "Acyl chloride + alcohol (Et3N, DMAP) for scale. "
                 "Enzymatic esterification for stereoselective.",
    },
    {
        "name": "epoxide_from_alkene",
        "product_smarts": "[C:1]1[O:2][C:3]1",
        "reaction_name": "epoxidation",
        "taxonomy_id": "Epoxidation",
        "difficulty": 0.25,
        "description": "Epoxide ← alkene (epoxidation)",
        "precursor_hints": ["alkene"],
        "notes": "mCPBA in DCM for unfunctionalized alkenes. "
                 "Sharpless asymmetric epoxidation (Ti/tartrate) for allylic alcohols. "
                 "Jacobsen Mn-salen for non-allylic alkenes. "
                 "Dimethyldioxirane (DMDO) for sensitive substrates.",
    },

    # -----------------------------------------------------------------------
    # 4. Carbonyl-Based Disconnections
    # -----------------------------------------------------------------------

    {
        "name": "ketone_from_friedel_crafts",
        "product_smarts": "[c:1]-[CX3:2](=[O:3])-[CX4,c:4]",
        "reaction_name": "friedel_crafts_acylation",
        "taxonomy_id": "Acylation_friedel_crafts",
        "difficulty": 0.25,
        "description": "Aryl ketone ← aromatic + acyl halide/anhydride (Friedel-Crafts)",
        "precursor_hints": ["aromatic_compound", "acyl_chloride_or_anhydride"],
        "notes": "AlCl3 Lewis acid (1.1+ equiv). "
                 "Works well for electron-rich arenes. Deactivated by EWG. "
                 "Fries rearrangement alternative for phenol esters. "
                 "Haworth synthesis for naphthalene/anthracene derivatives.",
    },
    {
        "name": "aldehyde_from_oxidation",
        "product_smarts": "[CX3H1:1](=[O:2])",
        "reaction_name": "primary_alcohol_oxidation",
        "taxonomy_id": "Oxidation_primary_alcohol_to_aldehyde",
        "difficulty": 0.20,
        "description": "Aldehyde ← primary alcohol (selective oxidation)",
        "precursor_hints": ["primary_alcohol"],
        "notes": "Swern oxidation (oxalyl chloride, DMSO, Et3N, −78°C): mild, no over-oxidation. "
                 "Dess-Martin periodinane (DMP): excellent for sensitive substrates. "
                 "TEMPO/NaOCl for aqueous conditions. "
                 "PCC or PDC in DCM for general use.",
    },
    {
        "name": "ketone_from_oxidation",
        "product_smarts": "[CX3:1](=[O:2])([CX4,c:3])[CX4,c:4]",
        "reaction_name": "secondary_alcohol_oxidation",
        "taxonomy_id": "Oxidation_secondary_alcohol_to_ketone",
        "difficulty": 0.15,
        "description": "Ketone ← secondary alcohol (oxidation)",
        "precursor_hints": ["secondary_alcohol"],
        "notes": "PCC, PDC, Swern, DMP, TEMPO, Oppenauer. "
                 "Highly reliable and straightforward. "
                 "Selectfluor/NaBr in MeCN/H2O also works (see notes).",
    },
    {
        "name": "carboxylic_acid_from_nitrile",
        "product_smarts": "[CX3:1](=[O:2])-[OH1:3]",
        "reaction_name": "nitrile_hydrolysis",
        "taxonomy_id": "Nitrile_hydrolysis",
        "difficulty": 0.20,
        "description": "Carboxylic acid ← nitrile (hydrolysis) or primary alcohol (oxidation)",
        "precursor_hints": ["nitrile", "primary_alcohol"],
        "notes": "Nitrile hydrolysis: H2SO4/H2O reflux or NaOH/H2O2. "
                 "Selective: partial hydrolysis to amide with H2O2/KOH. "
                 "Permanganate or Jones oxidation of primary alcohols for acid.",
    },

    # -----------------------------------------------------------------------
    # 5. Ring-Forming Reactions
    # -----------------------------------------------------------------------

    {
        "name": "cyclohexene_diels_alder",
        "product_smarts": "[C:1]1[C:2][C:3]=[C:4][C:5][C:6]1",
        "reaction_name": "diels_alder",
        "taxonomy_id": "Diels_Alder",
        "difficulty": 0.35,
        "description": "Cyclohexene ← diene + dienophile (Diels-Alder [4+2])",
        "precursor_hints": ["s_cis_diene", "activated_dienophile"],
        "notes": "Concerted [4+2] cycloaddition. endo/exo selectivity by geometry. "
                 "Lewis acid catalysis (BF3, Et2AlCl) accelerates and improves endo selectivity. "
                 "Asymmetric: Evans auxiliaries, chiral Lewis acids. "
                 "Intramolecular DA for polycyclic targets.",
    },
    {
        "name": "cyclopropane_simmons_smith",
        "product_smarts": "[C:1]1[C:2][C:3]1",
        "reaction_name": "simmons_smith",
        "taxonomy_id": "Simmons_Smith_cyclopropanation",
        "difficulty": 0.40,
        "description": "Cyclopropane ← alkene + carbenoid (Simmons-Smith)",
        "precursor_hints": ["alkene", "zinc_carbenoid"],
        "notes": "IZnCH2I (from Et2Zn + CH2I2 or Zn + CH2I2). "
                 "Syn addition, directed by adjacent OH groups. "
                 "Alternative: Corey-Chaykovsky (sulfonium ylide) for terminal epoxides → cyclopropanes.",
    },
    {
        "name": "cyclopentane_pauson_khand",
        "product_smarts": "[C:1]1[C:2][C:3][C:4][C:5]1(=[O:6])",
        "reaction_name": "pauson_khand",
        "taxonomy_id": "Pauson_Khand",
        "difficulty": 0.55,
        "description": "Cyclopentenone ← enyne (Pauson-Khand [2+2+1])",
        "precursor_hints": ["enyne"],
        "notes": "Co2(CO)8 mediated. Intramolecular preferred. "
                 "NMO or TMANO as promoter. "
                 "Asymmetric PKR with chiral Co complexes.",
    },
    {
        "name": "lactam_ring_closure",
        "product_smarts": "[C:1]1(=[O:2])-[NX3:3]-[CX4:4]1",
        "reaction_name": "lactam_formation",
        "taxonomy_id": "Lactam_formation",
        "difficulty": 0.30,
        "description": "Lactam ← amino acid or amino ester (intramolecular amide coupling)",
        "precursor_hints": ["amino_acid_or_amino_ester"],
        "notes": "Intramolecular amide bond formation. "
                 "Ring size determines strategy: 5-membered (γ-lactam) is easiest. "
                 "β-lactam: Staudinger ([2+2]), ketene + imine. "
                 "High-dilution or slow addition conditions for macrolactams.",
    },
    {
        "name": "heterocycle_buchwald_n_arylation",
        "product_smarts": "[n:1]-[c:2]",
        "reaction_name": "n_heterocycle_arylation",
        "taxonomy_id": "C_N_Coupling",
        "difficulty": 0.30,
        "description": "N-aryl heterocycle ← heteroaryl halide + N-H heterocycle (N-arylation)",
        "precursor_hints": ["aryl_halide", "nh_heterocycle"],
        "notes": "Pd/XPhos or Cu/phenanthroline catalysis. "
                 "Pyrazoles, imidazoles, pyrroles, indoles are compatible. "
                 "Chan-Lam (Cu/O2, boronic acid) at room temperature.",
    },
    {
        "name": "macrocycle_rce",
        "product_smarts": "[C:1]=[C:2]",
        "reaction_name": "ring_closing_metathesis",
        "taxonomy_id": "Ring_Closing_Metathesis",
        "difficulty": 0.50,
        "description": "Macrocyclic or medium-ring alkene ← diene (ring-closing metathesis)",
        "precursor_hints": ["diene"],
        "notes": "Grubbs 2nd gen or Hoveyda-Grubbs catalyst. "
                 "High-dilution conditions for macrocycles. Ethylene gas evolution drives reaction. "
                 "E/Z selectivity challenging. Cl2CH2 or toluene solvent.",
    },

    # -----------------------------------------------------------------------
    # 6. Oxidation / Reduction Markers
    # -----------------------------------------------------------------------

    {
        "name": "alcohol_from_reduction",
        "product_smarts": "[CX4:1]([OH1:2])-[CX3:3]=[O:4]",
        "reaction_name": "carbonyl_reduction",
        "taxonomy_id": "Reduction_carbonyl_to_alcohol",
        "difficulty": 0.10,
        "description": "Alcohol ← carbonyl compound (reduction)",
        "precursor_hints": ["aldehyde_or_ketone"],
        "notes": "NaBH4 (MeOH, 0°C→rt) for aldehydes/ketones. "
                 "LiAlH4 for esters and amides. "
                 "DIBAL-H for selective reduction (ester → aldehyde at −78°C). "
                 "Asymmetric: CBS reduction (Corey-Bakshi-Shibata) with BH3·Me2S.",
    },
    {
        "name": "amine_from_reduction",
        "product_smarts": "[CX4:1]-[NX3H2:2]",
        "reaction_name": "imine_reduction",
        "taxonomy_id": "Reduction_imine_to_amine",
        "difficulty": 0.20,
        "description": "Amine ← imine or nitro compound (reduction)",
        "precursor_hints": ["imine_or_nitro_compound"],
        "notes": "NaBH4 or NaBH3CN for imines. "
                 "H2/Pd-C for nitro → amine (transfer hydrogenation also works). "
                 "Zn/AcOH for selective partial reduction of nitro groups. "
                 "Asymmetric imine reduction: Ir(I)/chiral phosphine.",
    },
    {
        "name": "alkene_sharpless",
        "product_smarts": "[CX4:1]([OH1:2])[CX4:3]([OH1:4])",
        "reaction_name": "dihydroxylation",
        "taxonomy_id": "Dihydroxylation",
        "difficulty": 0.25,
        "description": "1,2-Diol ← alkene (dihydroxylation)",
        "precursor_hints": ["alkene"],
        "notes": "OsO4 (cat., NMO as reoxidant, Upjohn). syn addition. "
                 "Sharpless AD (AD-mix-α or β): asymmetric with DHQD or DHQ ligands. "
                 "KMnO4 (basic): syn-diol but over-oxidation risk.",
    },

    # -----------------------------------------------------------------------
    # 7. Protecting Group Markers
    # -----------------------------------------------------------------------

    {
        "name": "boc_protected_amine",
        "product_smarts": "[NX3:1]-[C:2](=[O:3])-[OX2:4]-[C:5]([CH3])([CH3])[CH3]",
        "reaction_name": "boc_protection",
        "taxonomy_id": "Boc_protection",
        "difficulty": 0.10,
        "description": "Boc-amine ← free amine + Boc2O (protection)",
        "precursor_hints": ["primary_or_secondary_amine"],
        "notes": "Boc2O, Et3N or DMAP, DCM. Removed by TFA (DCM) or HCl (dioxane). "
                 "Orthogonal to Cbz (H2/Pd) and Fmoc (piperidine).",
    },
    {
        "name": "cbz_protected_amine",
        "product_smarts": "[NX3:1]-[C:2](=[O:3])-[OX2:4]-[CH2:5]-[c:6]",
        "reaction_name": "cbz_protection",
        "taxonomy_id": "Cbz_protection",
        "difficulty": 0.10,
        "description": "Cbz-amine ← free amine + benzyl chloroformate (Cbz-Cl)",
        "precursor_hints": ["primary_or_secondary_amine"],
        "notes": "CbzCl, Et3N, DCM. Removed by H2/Pd-C (hydrogenolysis). "
                 "Stable to TFA conditions (cf. Boc). Useful for orthogonal protection.",
    },
    {
        "name": "silyl_ether_protection",
        "product_smarts": "[OX2:1]-[Si:2]([C:3])([C:4])[C:5]",
        "reaction_name": "silyl_ether_protection",
        "taxonomy_id": "Silyl_ether_protection",
        "difficulty": 0.10,
        "description": "TMS/TBS/TIPS silyl ether ← alcohol + silyl chloride",
        "precursor_hints": ["alcohol"],
        "notes": "TBSCl/imidazole/DMF (most common). TIPS for sensitive substrates. "
                 "TMS for temporary protection. Removed by TBAF (THF) or HF·pyridine. "
                 "Selectivity: primary > secondary; less hindered > more hindered.",
    },

    # -----------------------------------------------------------------------
    # 8. Halogenation
    # -----------------------------------------------------------------------

    {
        "name": "aryl_fluorine_snarf",
        "product_smarts": "[c:1]-[F:2]",
        "reaction_name": "nucleophilic_aromatic_fluorination",
        "taxonomy_id": "SNAr",
        "difficulty": 0.55,
        "description": "Aryl fluoride ← aryl halide + F– (SNAr) or Balz-Schiemann",
        "precursor_hints": ["aryl_halide_with_ewg", "fluoride_source"],
        "notes": "SNAr requires strong EWG (NO2, CF3). KF or CsF, polar aprotic. "
                 "Balz-Schiemann: diazonium tetrafluoroborate → pyrolysis. "
                 "Pd-catalyzed fluorination (Pd/BrettPhos + Et3N·HF) for wider scope.",
    },
    {
        "name": "aryl_chloride_electrophilic",
        "product_smarts": "[c:1]-[Cl:2]",
        "reaction_name": "electrophilic_aromatic_chlorination",
        "taxonomy_id": "Halogenation_aromatic",
        "difficulty": 0.20,
        "description": "Aryl chloride ← arene + Cl2/electrophilic Cl source (EAS)",
        "precursor_hints": ["aromatic_compound"],
        "notes": "Cl2/AlCl3, NCS (N-chlorosuccinimide)/Lewis acid, or SO2Cl2. "
                 "Electron-rich arenes react well. Regioselectivity determined by existing substituents.",
    },
    {
        "name": "aryl_bromide_electrophilic",
        "product_smarts": "[c:1]-[Br:2]",
        "reaction_name": "electrophilic_aromatic_bromination",
        "taxonomy_id": "Halogenation_aromatic",
        "difficulty": 0.15,
        "description": "Aryl bromide ← arene + Br2 or NBS (EAS)",
        "precursor_hints": ["aromatic_compound"],
        "notes": "Br2/FeBr3 or NBS in DCM/H2O. "
                 "NBS preferred for electron-rich arenes (mild, selective). "
                 "Regioselectivity driven by existing substituents (o/p directors).",
    },
    {
        "name": "allylic_bromination",
        "product_smarts": "[CH2:1]([Br:2])-[CX3:3]=[CX3:4]",
        "reaction_name": "nbs_allylic_bromination",
        "taxonomy_id": "Halogenation_unsaturated",
        "difficulty": 0.25,
        "description": "Allylic bromide ← alkene + NBS (radical allylic bromination)",
        "precursor_hints": ["alkene"],
        "notes": "NBS, benzoyl peroxide or hν, CCl4. "
                 "Wohl-Ziegler conditions. Allylic radical is resonance-stabilized. "
                 "Mixture of regioisomers possible with unsymmetric alkenes.",
    },
    {
        "name": "benzylic_bromination",
        "product_smarts": "[CH2:1]([Br:2])-[c:3]",
        "reaction_name": "benzylic_bromination",
        "taxonomy_id": "Halogenation_aromatic",
        "difficulty": 0.20,
        "description": "Benzyl bromide ← toluene + NBS (radical benzylic bromination)",
        "precursor_hints": ["methylarene"],
        "notes": "NBS, AIBN or hν, CCl4. High selectivity for benzylic position. "
                 "Used to access benzyl electrophiles for further reactions.",
    },

    # -----------------------------------------------------------------------
    # 9. Additional C–N / C–C Specialties
    # -----------------------------------------------------------------------

    {
        "name": "aryl_boronic_acid_miyaura",
        "product_smarts": "[c:1]-[B:2]([OH1:3])[OH1:4]",
        "reaction_name": "miyaura_borylation",
        "taxonomy_id": "Miyaura_borylation",
        "difficulty": 0.30,
        "description": "Arylboronic acid ← aryl halide + B2pin2 (Miyaura borylation)",
        "precursor_hints": ["aryl_halide"],
        "notes": "Pd(dppf)Cl2, B2pin2, KOAc, dioxane, 80°C. "
                 "Alternative to arylboron reagents from ArMgX + B(OEt)3. "
                 "Pin boronate can be hydrolyzed to boronic acid with NaIO4 or H2O2/NaOH.",
    },
    {
        "name": "vinyl_triflate_from_ketone",
        "product_smarts": "[CX3:1]=[CX3:2]-[OX2:3]-S(=O)(=O)-[CF3:4]",
        "reaction_name": "enol_triflate_formation",
        "taxonomy_id": "Enol_triflate_formation",
        "difficulty": 0.30,
        "description": "Vinyl triflate ← ketone + Tf2O or N-phenyltriflimide (via enolate)",
        "precursor_hints": ["ketone"],
        "notes": "LDA → enolate → Tf2O or PhNTf2, THF, −78°C. "
                 "Regioselective: kinetic (LDA, −78°C) vs thermodynamic (Et3N/PhNTf2). "
                 "Used as aryl halide surrogate in cross-coupling.",
    },
    {
        "name": "organotrifluoroborate_from_halide",
        "product_smarts": "[c:1]-[B-:2]([F:3])([F:4])[F:5]",
        "reaction_name": "pot_aryltrifluoroborate",
        "taxonomy_id": "Aryltrifluoroborate_synthesis",
        "difficulty": 0.35,
        "description": "Potassium aryl/vinyl trifluoroborate ← halide + B2pin2 then KHF2",
        "precursor_hints": ["aryl_halide_or_vinyl_halide"],
        "notes": "Pd-catalyzed borylation then KHF2 precipitation. "
                 "BF3K salts are air/moisture-stable alternatives to boronic acids. "
                 "Direct coupling via Suzuki-Miyaura with Pd/base.",
    },
    {
        "name": "c_h_functionalization_direct",
        "product_smarts": "[c:1]-[CX4:2]",
        "reaction_name": "c_h_activation",
        "taxonomy_id": "Arylation_Ar_H",
        "difficulty": 0.65,
        "description": "C–C bond via C–H activation ← arene + coupling partner",
        "precursor_hints": ["aromatic_substrate", "coupling_partner"],
        "notes": "Pd, Rh, Ir, or Ru catalysis. Directing group (DG) required for selectivity. "
                 "CMD mechanism (acetate base, CMD conditions). "
                 "Intramolecular C–H arylation for fused rings. "
                 "Challenging: typically high T, specialized ligands.",
    },

    # -----------------------------------------------------------------------
    # 10. Aliphatic C–N / C–O Specialties
    # -----------------------------------------------------------------------

    {
        "name": "sn2_alkyl_amine",
        "product_smarts": "[CX4:1]-[NX3:2]",
        "reaction_name": "n_alkylation",
        "taxonomy_id": "Alkyl_Nucleophilic_Substitution",
        "difficulty": 0.15,
        "description": "N-alkylated amine ← amine + alkyl halide (SN2 N-alkylation)",
        "precursor_hints": ["amine", "alkyl_halide_primary_or_methyl"],
        "notes": "K2CO3 in DMF or MeCN. Primary alkyl halides only (Br/I preferred). "
                 "Over-alkylation is a major concern: use excess amine or Mitsunobu. "
                 "Gabriel synthesis (phthalimide → phthalimide alkylation → hydrazine deprotection) "
                 "for primary amines without over-alkylation.",
    },
    {
        "name": "lactone_from_hydroxy_acid",
        "product_smarts": "[CX3:1]1(=[O:2])-[OX2:3]-[CX4:4]1",
        "reaction_name": "lactonization",
        "taxonomy_id": "Lactonization",
        "difficulty": 0.25,
        "description": "Lactone ← hydroxy acid (intramolecular esterification)",
        "precursor_hints": ["hydroxy_acid"],
        "notes": "5- and 6-membered lactones form readily (DCC or acid catalyst). "
                 "Macrolactonization: high dilution, Yamaguchi, Corey-Nicolaou, or Shiina conditions. "
                 "Enzymatic (lipase) for asymmetric.",
    },
    {
        "name": "sulfonamide_formation",
        "product_smarts": "[NX3:1]-[SX4:2](=[O:3])(=[O:4])-[CX4,c:5]",
        "reaction_name": "sulfonamide_synthesis",
        "taxonomy_id": "Sulfonamide_synthesis",
        "difficulty": 0.15,
        "description": "Sulfonamide ← amine + sulfonyl chloride",
        "precursor_hints": ["amine", "sulfonyl_chloride"],
        "notes": "Et3N or pyridine in DCM. Very reliable. "
                 "Sulfonyl chlorides can be made from sulfonic acids + PCl5 or SOCl2. "
                 "Sulfonamides as NH-acid: can be deprotonated (pKa ~10) and alkylated.",
    },
]


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def get_retron_by_name(name: str) -> Dict[str, Any]:
    """Return a retron dict by its name, or empty dict if not found."""
    for r in RETRONS:
        if r["name"] == name:
            return r
    return {}


def get_retrons_by_reaction(reaction_name: str) -> List[Dict[str, Any]]:
    """Return all retrons that map to a given reaction taxonomy name."""
    return [r for r in RETRONS if r["reaction_name"] == reaction_name]


def get_retrons_by_difficulty(max_difficulty: float = 1.0) -> List[Dict[str, Any]]:
    """Return all retrons with difficulty <= max_difficulty, sorted ascending."""
    filtered = [r for r in RETRONS if r["difficulty"] <= max_difficulty]
    return sorted(filtered, key=lambda r: r["difficulty"])
