"""
Integrate high-quality extracted retro templates into the production
hte_templates.json and retron_patterns.json files.

This module takes the output of ``extract_templates`` and merges the
best templates into the production taxonomy data files used by the
retrosynthetic analysis engine.

Filters applied:
  - retro_validated == True
  - quality_score >= 0.7
  - count >= 10
  - Pick best template per (cluster_id, primary_source_family)

Mapping:
  - source_files → taxonomy_family_id via a curated lookup table
  - reaction_smarts (forward) → retro_smarts (reversed) + forward_smarts

Output:
  - Updated chemtools/taxonomy/data/hte_templates.json  (backup created)
  - Updated chemtools/taxonomy/data/retron_patterns.json (backup created)

Usage (as CLI):
    # After running extraction (Step 1):
    python -m chemtools.retro.extract_templates --all --limit 0

    # Integrate the extracted templates (Step 2):
    python -m chemtools.retro.integrate_templates

Usage (as library):
    from chemtools.retro.integrate_templates import integrate
    n_templates, n_retrons = integrate()

Adding a NEW reaction family:
    If your new CSV introduces a reaction family not in SOURCE_TO_TAXONOMY,
    add a mapping entry before running integration:

        SOURCE_TO_TAXONOMY["My_new_reaction"] = {
            "taxonomy_family_id": "My_new_reaction",
            "category": "C-C Bond Forming",
            "description": "Description of the reaction",
            "difficulty": 0.3,
        }

    Otherwise the script prints "Skipped (no source mapping): N" and drops
    those templates silently.

Full pipeline (extract + integrate):
    # Step 1: Extract templates from HTE data
    python -m chemtools.retro.extract_templates --all --limit 0

    # Step 2: Integrate high-quality templates into production files
    python -m chemtools.retro.integrate_templates

    # Step 3: Verify
    python -m pytest tests/ -q
"""
from __future__ import annotations

import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_MODULE_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _MODULE_DIR.parent.parent
EXTRACTED_FILE = _PROJECT_ROOT / "results" / "extracted_templates.json"
HTE_FILE = _PROJECT_ROOT / "chemtools" / "taxonomy" / "data" / "hte_templates.json"
RETRON_FILE = _PROJECT_ROOT / "chemtools" / "taxonomy" / "data" / "retron_patterns.json"

# ---------------------------------------------------------------------------
# Source-file stem → (taxonomy_family_id, category, description, difficulty)
# ---------------------------------------------------------------------------
# Maps each extracted source_file stem to the canonical taxonomy_family_id,
# a category label, a brief description, and a default difficulty score.
SOURCE_TO_TAXONOMY: dict[str, dict] = {
    # ---- C-N Bond Forming ----
    "C_N_Coupling": {
        "taxonomy_family_id": "C_N_Coupling",
        "category": "C-N Bond Forming",
        "description": "Aryl-N bond via Buchwald-Hartwig / Pd-catalyzed C-N coupling",
        "difficulty": 0.35,
    },
    "SNAr_amination": {
        "taxonomy_family_id": "C_N_Coupling",
        "category": "C-N Bond Forming",
        "description": "SNAr amination on electron-poor aryl fluorides",
        "difficulty": 0.2,
    },
    "ChanLam_Narylation": {
        "taxonomy_family_id": "Chan_Lam_C_N_Coupling",
        "category": "C-N Bond Forming",
        "description": "Chan-Lam N-arylation (Cu-catalyzed)",
        "difficulty": 0.3,
    },
    "ChanLam_sulfonamide_converted": {
        "taxonomy_family_id": "Chan_Lam_C_N_Coupling",
        "category": "C-N Bond Forming",
        "description": "Chan-Lam sulfonamide coupling",
        "difficulty": 0.3,
    },
    "Amide_formation": {
        "taxonomy_family_id": "Amide_formation",
        "category": "C-N Bond Forming",
        "description": "Amide bond formation (HATU/EDC/DIC coupling)",
        "difficulty": 0.15,
    },
    "CDI-mediated_amidation": {
        "taxonomy_family_id": "Amide_formation",
        "category": "C-N Bond Forming",
        "description": "CDI-mediated amidation",
        "difficulty": 0.2,
    },
    "Anhydride_coupling": {
        "taxonomy_family_id": "Amide_formation",
        "category": "C-N Bond Forming",
        "description": "Anhydride-based amide formation",
        "difficulty": 0.15,
    },
    "Weinreb_amide": {
        "taxonomy_family_id": "Amide_formation",
        "category": "C-N Bond Forming",
        "description": "Weinreb amide formation",
        "difficulty": 0.25,
    },
    "Reductive_amination": {
        "taxonomy_family_id": "Reductive_amination",
        "category": "C-N Bond Forming",
        "description": "Reductive amination (NaBH(OAc)3/NaBH3CN)",
        "difficulty": 0.2,
    },
    "Reduction_NaBHOAc3_NaBH3CN": {
        "taxonomy_family_id": "Reductive_amination",
        "category": "C-N Bond Forming",
        "description": "Reductive amination with NaBH(OAc)3/NaBH3CN",
        "difficulty": 0.2,
    },
    "Urea_and_thiourea_formation": {
        "taxonomy_family_id": "Urea_thiourea_formation",
        "category": "C-N Bond Forming",
        "description": "Urea/thiourea formation from amine + isocyanate",
        "difficulty": 0.15,
    },
    "Sulfonamide_formation": {
        "taxonomy_family_id": "Sulfonamide_synthesis",
        "category": "C-N Bond Forming",
        "description": "Sulfonamide from sulfonyl chloride + amine",
        "difficulty": 0.1,
    },
    "Gabriel_amine_synthesis": {
        "taxonomy_family_id": "Gabriel_amine_synthesis",
        "category": "C-N Bond Forming",
        "description": "Gabriel amine synthesis (phthalimide route)",
        "difficulty": 0.25,
    },
    "Curtius_rearrangement": {
        "taxonomy_family_id": "Curtius_Hofmann_rearrangement",
        "category": "C-N Bond Forming",
        "description": "Curtius rearrangement (acyl azide → isocyanate)",
        "difficulty": 0.35,
    },
    "Hofmann_rearrangement": {
        "taxonomy_family_id": "Curtius_Hofmann_rearrangement",
        "category": "C-N Bond Forming",
        "description": "Hofmann rearrangement (amide → amine)",
        "difficulty": 0.35,
    },
    "CH_amination_and_azidation": {
        "taxonomy_family_id": "CH_amination",
        "category": "C-N Bond Forming",
        "description": "C-H amination / azidation",
        "difficulty": 0.5,
    },
    "EschweilerClarke_methylation": {
        "taxonomy_family_id": "EschweilerClarke_methylation",
        "category": "C-N Bond Forming",
        "description": "Eschweiler-Clarke N-methylation (HCHO/HCO2H)",
        "difficulty": 0.2,
    },
    # ---- C-C Bond Forming ----
    "suzuki_miyaura": {
        "taxonomy_family_id": "Suzuki_miyaura",
        "category": "C-C Bond Forming",
        "description": "Suzuki-Miyaura biaryl coupling",
        "difficulty": 0.15,
    },
    "HeckMizoroki_coupling": {
        "taxonomy_family_id": "Heck",
        "category": "C-C Bond Forming",
        "description": "Heck-Mizoroki vinylation",
        "difficulty": 0.3,
    },
    "Sonogashira_coupling": {
        "taxonomy_family_id": "Sonogashira",
        "category": "C-C Bond Forming",
        "description": "Sonogashira aryl-alkyne coupling",
        "difficulty": 0.25,
    },
    "Stille_coupling": {
        "taxonomy_family_id": "Stille",
        "category": "C-C Bond Forming",
        "description": "Stille coupling with organostannanes",
        "difficulty": 0.35,
    },
    "Negishi_coupling": {
        "taxonomy_family_id": "Negishi",
        "category": "C-C Bond Forming",
        "description": "Negishi coupling with organozinc reagents",
        "difficulty": 0.35,
    },
    "Kumada_coupling": {
        "taxonomy_family_id": "Kumada",
        "category": "C-C Bond Forming",
        "description": "Kumada coupling with Grignard reagents",
        "difficulty": 0.35,
    },
    "Hiyama_coupling": {
        "taxonomy_family_id": "Hiyama",
        "category": "C-C Bond Forming",
        "description": "Hiyama coupling with organosilanes",
        "difficulty": 0.35,
    },
    "LiebeskindSrogl_coupling": {
        "taxonomy_family_id": "Liebeskind_Srogl",
        "category": "C-C Bond Forming",
        "description": "Liebeskind-Srogl thioether C-C coupling",
        "difficulty": 0.45,
    },
    "Decarboxylative_arylation": {
        "taxonomy_family_id": "Decarboxylative_Coupling",
        "category": "C-C Bond Forming",
        "description": "Decarboxylative arylation",
        "difficulty": 0.5,
    },
    "Giese_radical_additions": {
        "taxonomy_family_id": "Giese_radical_addition",
        "category": "C-C Bond Forming",
        "description": "Giese radical addition to alkenes",
        "difficulty": 0.4,
    },
    "Michael_addition": {
        "taxonomy_family_id": "Michael_addition",
        "category": "C-C Bond Forming",
        "description": "Michael 1,4-addition to enones",
        "difficulty": 0.3,
    },
    "Knoevenagel_condensation": {
        "taxonomy_family_id": "Knoevenagel_condensation",
        "category": "C-C Bond Forming",
        "description": "Knoevenagel condensation (active methylene + aldehyde)",
        "difficulty": 0.25,
    },
    "Aldol_classic__Mukaiyama": {
        "taxonomy_family_id": "Aldol_addition",
        "category": "C-C Bond Forming",
        "description": "Aldol condensation (classic & Mukaiyama)",
        "difficulty": 0.4,
    },
    "HornerWadsworthEmmons": {
        "taxonomy_family_id": "Wittig_reaction",
        "category": "C-C Bond Forming",
        "description": "Horner-Wadsworth-Emmons olefination",
        "difficulty": 0.25,
    },
    "NozakiHiyamaKishi": {
        "taxonomy_family_id": "NHK_reaction",
        "category": "C-C Bond Forming",
        "description": "Nozaki-Hiyama-Kishi allylation/vinylation",
        "difficulty": 0.45,
    },
    "Minisci": {
        "taxonomy_family_id": "Minisci",
        "category": "C-C Bond Forming",
        "description": "Minisci radical addition to heteroarenes",
        "difficulty": 0.45,
    },
    "Wittig_reaction": {
        "taxonomy_family_id": "Wittig_reaction",
        "category": "C-C Bond Forming",
        "description": "Wittig olefination (phosphonium ylide + aldehyde)",
        "difficulty": 0.3,
    },
    "Julia-Kocienski_Olefination": {
        "taxonomy_family_id": "Julia_Kocienski_olefination",
        "category": "C-C Bond Forming",
        "description": "Julia-Kocienski olefination (sulfone + aldehyde)",
        "difficulty": 0.35,
    },
    "Olefin_metathesis": {
        "taxonomy_family_id": "Olefin_metathesis",
        "category": "C-C Bond Forming",
        "description": "Olefin metathesis (Grubbs catalyst)",
        "difficulty": 0.35,
    },
    "Mannich__Petasis_boronoMannich": {
        "taxonomy_family_id": "Mannich_Petasis",
        "category": "C-C Bond Forming",
        "description": "Mannich / Petasis borono-Mannich reaction",
        "difficulty": 0.35,
    },
    # ---- C-O Bond Forming ----
    "C_O_Coupling": {
        "taxonomy_family_id": "C_O_Coupling",
        "category": "C-O Bond Forming",
        "description": "Pd-catalyzed C-O coupling / O-arylation",
        "difficulty": 0.35,
    },
    "Williamson_ether_synthesis": {
        "taxonomy_family_id": "Alkyl_Nucleophilic_Substitution",
        "category": "C-O Bond Forming",
        "description": "Williamson ether synthesis",
        "difficulty": 0.15,
    },
    "Fischer_esterification": {
        "taxonomy_family_id": "Esterification",
        "category": "C-O Bond Forming",
        "description": "Fischer esterification (acid + alcohol)",
        "difficulty": 0.1,
    },
    "Steglich_esterification": {
        "taxonomy_family_id": "Esterification",
        "category": "C-O Bond Forming",
        "description": "Steglich esterification (DCC/DMAP mediated)",
        "difficulty": 0.15,
    },
    "Mitsunobu_etherificationinversion": {
        "taxonomy_family_id": "Mitsunobu_reaction",
        "category": "C-O Bond Forming",
        "description": "Mitsunobu etherification with stereoinversion",
        "difficulty": 0.3,
    },
    "Acetal_ketal_formation": {
        "taxonomy_family_id": "Acetal_ketal_formation",
        "category": "C-O Bond Forming",
        "description": "Acetal/ketal formation (carbonyl protection)",
        "difficulty": 0.15,
    },
    "Carbamate_formation": {
        "taxonomy_family_id": "Carbamate_formation",
        "category": "C-O Bond Forming",
        "description": "Carbamate formation (amine + Boc/Cbz/Fmoc reagent)",
        "difficulty": 0.1,
    },
    # ---- C-S Bond Forming ----
    "C_S_Coupling": {
        "taxonomy_family_id": "C_S_Coupling",
        "category": "C-S Bond Forming",
        "description": "C-S coupling (thiol + aryl halide)",
        "difficulty": 0.35,
    },
    "Thioether_formation": {
        "taxonomy_family_id": "Thioether_formation",
        "category": "C-S Bond Forming",
        "description": "Thioether formation (thiol + electrophile)",
        "difficulty": 0.2,
    },
    "Thiolene_and_thiolyne": {
        "taxonomy_family_id": "Thiol_ene",
        "category": "C-S Bond Forming",
        "description": "Thiol-ene / thiol-yne click reaction",
        "difficulty": 0.2,
    },
    "SuFEx": {
        "taxonomy_family_id": "SuFEx",
        "category": "S-F Bond Forming",
        "description": "Sulfur(VI) fluoride exchange (SuFEx click)",
        "difficulty": 0.25,
    },
    # ---- Click / Cycloaddition ----
    "CuAAC_azidealkyne": {
        "taxonomy_family_id": "Click_azide_alkyne_cycloaddition",
        "category": "Click / Cycloaddition",
        "description": "CuAAC azide-alkyne cycloaddition (click)",
        "difficulty": 0.15,
    },
    "SPAAC_strainpromoted_click": {
        "taxonomy_family_id": "Click_azide_alkyne_cycloaddition",
        "category": "Click / Cycloaddition",
        "description": "Strain-promoted azide-alkyne cycloaddition (SPAAC)",
        "difficulty": 0.2,
    },
    "Strainpromoted_click": {
        "taxonomy_family_id": "Click_azide_alkyne_cycloaddition",
        "category": "Click / Cycloaddition",
        "description": "Strain-promoted click (general)",
        "difficulty": 0.2,
    },
    "TetrazineTCO_iEDDA": {
        "taxonomy_family_id": "TetrazineTCO_iEDDA",
        "category": "Click / Cycloaddition",
        "description": "Tetrazine-TCO inverse electron-demand Diels-Alder",
        "difficulty": 0.25,
    },
    # ---- Reduction ----
    "NaBH4_carbonyl_reductions": {
        "taxonomy_family_id": "Reduction_carbonyl_to_alcohol",
        "category": "Reduction",
        "description": "NaBH4 carbonyl reduction to alcohol",
        "difficulty": 0.1,
    },
    "Reduction_LAH_LiAlH4": {
        "taxonomy_family_id": "Reduction_ester_to_alcohol",
        "category": "Reduction",
        "description": "LiAlH4 reduction (ester/acid → alcohol)",
        "difficulty": 0.2,
    },
    "DIBALH_Partial_reductions": {
        "taxonomy_family_id": "DIBALH_reduction",
        "category": "Reduction",
        "description": "DIBAL-H partial reduction (ester → aldehyde)",
        "difficulty": 0.3,
    },
    "Reduction-Borane": {
        "taxonomy_family_id": "Borane_reduction",
        "category": "Reduction",
        "description": "Borane reduction (amide/acid → amine/alcohol)",
        "difficulty": 0.3,
    },
    "Catalytic_hydrogenation": {
        "taxonomy_family_id": "Catalytic_hydrogenation",
        "category": "Reduction",
        "description": "Catalytic hydrogenation (Pd/C, H2)",
        "difficulty": 0.15,
    },
    "Transfer_hydrogenation": {
        "taxonomy_family_id": "Transfer_hydrogenation",
        "category": "Reduction",
        "description": "Transfer hydrogenation (HCO2NH4, cyclohexadiene)",
        "difficulty": 0.2,
    },
    "Clemmensen_reduction": {
        "taxonomy_family_id": "Clemmensen_reduction",
        "category": "Reduction",
        "description": "Clemmensen reduction (carbonyl → methylene, Zn/Hg)",
        "difficulty": 0.35,
    },
    "WolffKishner": {
        "taxonomy_family_id": "WolffKishner_reduction",
        "category": "Reduction",
        "description": "Wolff-Kishner reduction (carbonyl → methylene, N2H4/KOH)",
        "difficulty": 0.35,
    },
    "Staudinger_reduction": {
        "taxonomy_family_id": "Staudinger_reduction",
        "category": "Reduction",
        "description": "Staudinger reduction (azide → amine, PPh3)",
        "difficulty": 0.15,
    },
    "Birch_reduction": {
        "taxonomy_family_id": "Birch_reduction",
        "category": "Reduction",
        "description": "Birch reduction (arene → 1,4-cyclohexadiene, Na/NH3)",
        "difficulty": 0.45,
    },
    "BartonMcCombie_deoxygenation": {
        "taxonomy_family_id": "BartonMcCombie_deoxygenation",
        "category": "Reduction",
        "description": "Barton-McCombie radical deoxygenation",
        "difficulty": 0.4,
    },
    "Rosenmund_reduction": {
        "taxonomy_family_id": "Rosenmund_reduction",
        "category": "Reduction",
        "description": "Rosenmund reduction (acyl chloride → aldehyde)",
        "difficulty": 0.35,
    },
    # ---- Oxidation ----
    "DessMartin_periodinane_DMP_Alcohols__aldehydesketones": {
        "taxonomy_family_id": "Oxidation_primary_alcohol_to_aldehyde",
        "category": "Oxidation",
        "description": "Dess-Martin periodinane oxidation",
        "difficulty": 0.15,
    },
    "Swern_oxidation": {
        "taxonomy_family_id": "Oxidation_primary_alcohol_to_aldehyde",
        "category": "Oxidation",
        "description": "Swern oxidation (DMSO/oxalyl chloride)",
        "difficulty": 0.2,
    },
    "TPAPNMO_Catalytic_Ru_oxidations": {
        "taxonomy_family_id": "Oxidation_secondary_alcohol_to_ketone",
        "category": "Oxidation",
        "description": "TPAP/NMO catalytic Ru oxidation",
        "difficulty": 0.2,
    },
    "TEMPObleach-Primary_alcohols_to_aldehydesacids": {
        "taxonomy_family_id": "TEMPO_oxidation",
        "category": "Oxidation",
        "description": "TEMPO/bleach oxidation (primary alcohol → aldehyde/acid)",
        "difficulty": 0.15,
    },
    "TEMPObleach_Oxidations_Primary_alcohols__aldehydesacids": {
        "taxonomy_family_id": "TEMPO_oxidation",
        "category": "Oxidation",
        "description": "TEMPO/bleach oxidation of primary alcohols",
        "difficulty": 0.15,
    },
    "Jones_oxidation": {
        "taxonomy_family_id": "Jones_oxidation",
        "category": "Oxidation",
        "description": "Jones oxidation (CrO3/H2SO4)",
        "difficulty": 0.2,
    },
    "Riley_oxidation": {
        "taxonomy_family_id": "Riley_oxidation",
        "category": "Oxidation",
        "description": "Riley oxidation (SeO2, allylic/benzylic oxidation)",
        "difficulty": 0.35,
    },
    "BAIB": {
        "taxonomy_family_id": "BAIB_oxidation",
        "category": "Oxidation",
        "description": "BAIB (PhI(OAc)2) oxidation",
        "difficulty": 0.2,
    },
    "Oxidation-BAIB": {
        "taxonomy_family_id": "BAIB_oxidation",
        "category": "Oxidation",
        "description": "BAIB oxidative cleavage",
        "difficulty": 0.2,
    },
    "IBX": {
        "taxonomy_family_id": "IBX_oxidation",
        "category": "Oxidation",
        "description": "IBX oxidation (alcohol → aldehyde/ketone)",
        "difficulty": 0.15,
    },
    "Oxidation_IBX": {
        "taxonomy_family_id": "IBX_oxidation",
        "category": "Oxidation",
        "description": "IBX oxidation",
        "difficulty": 0.15,
    },
    "Benzylic_oxidation": {
        "taxonomy_family_id": "Benzylic_oxidation",
        "category": "Oxidation",
        "description": "Benzylic oxidation (KMnO4 / PDC)",
        "difficulty": 0.25,
    },
    "Oxidation_of_Aromatic_Side_Chains": {
        "taxonomy_family_id": "Benzylic_oxidation",
        "category": "Oxidation",
        "description": "Oxidation of aromatic side chains",
        "difficulty": 0.25,
    },
    "Oxidation_of_Arylmethanes_to_Aldehydes": {
        "taxonomy_family_id": "Benzylic_oxidation",
        "category": "Oxidation",
        "description": "Oxidation of arylmethanes to aldehydes",
        "difficulty": 0.3,
    },
    "Wacker_oxidation_Alkene__ketone_PdCl2_CuCl_O2": {
        "taxonomy_family_id": "Wacker_oxidation",
        "category": "Oxidation",
        "description": "Wacker oxidation (alkene → ketone, PdCl2/CuCl/O2)",
        "difficulty": 0.35,
    },
    "Wacker_oxidation": {
        "taxonomy_family_id": "Wacker_oxidation",
        "category": "Oxidation",
        "description": "Wacker oxidation (alkene → ketone)",
        "difficulty": 0.35,
    },
    "mCPBA_epoxidation_Alkene__epoxide": {
        "taxonomy_family_id": "mCPBA_epoxidation",
        "category": "Oxidation",
        "description": "mCPBA epoxidation (alkene → epoxide)",
        "difficulty": 0.15,
    },
    "Sharpless_dihydroxylation": {
        "taxonomy_family_id": "Sharpless_dihydroxylation",
        "category": "Oxidation",
        "description": "Sharpless asymmetric dihydroxylation (OsO4/AD-mix)",
        "difficulty": 0.3,
    },
    # ---- Halogenation ----
    "Electrophilic_fluorination": {
        "taxonomy_family_id": "Electrophilic_fluorination",
        "category": "Halogenation",
        "description": "Electrophilic fluorination (Selectfluor / NFSI)",
        "difficulty": 0.3,
    },
    "Deoxy_fluorination": {
        "taxonomy_family_id": "Alcohol_to_Alkyl_Halide",
        "category": "Halogenation",
        "description": "Deoxy-fluorination (DAST / Deoxo-Fluor)",
        "difficulty": 0.3,
    },
    "Difluoromethylation": {
        "taxonomy_family_id": "Difluoromethylation",
        "category": "Halogenation",
        "description": "Difluoromethylation (CHF2 installation)",
        "difficulty": 0.4,
    },
    "Trifluoromethylation": {
        "taxonomy_family_id": "Trifluoromethylation",
        "category": "Halogenation",
        "description": "Trifluoromethylation (CF3 installation)",
        "difficulty": 0.4,
    },
    "RuppertPrakash_TMSCF3": {
        "taxonomy_family_id": "Trifluoromethylation",
        "category": "Halogenation",
        "description": "Ruppert-Prakash trifluoromethylation (TMSCF3)",
        "difficulty": 0.35,
    },
    "Sandmeyer_reactions": {
        "taxonomy_family_id": "Sandmeyer",
        "category": "Halogenation",
        "description": "Sandmeyer reaction (diazonium → halide/CN/etc.)",
        "difficulty": 0.3,
    },
    "BalzSchiemann": {
        "taxonomy_family_id": "Sandmeyer",
        "category": "Halogenation",
        "description": "Balz-Schiemann fluorination (diazonium + HBF4)",
        "difficulty": 0.35,
    },
    "Appel_halogenation": {
        "taxonomy_family_id": "Appel_halogenation",
        "category": "Halogenation",
        "description": "Appel halogenation (alcohol → alkyl halide, PPh3/CX4)",
        "difficulty": 0.2,
    },
    "Aliphatic_Halide_Exchange": {
        "taxonomy_family_id": "Aliphatic_halide_exchange",
        "category": "Halogenation",
        "description": "Aliphatic halide exchange (Finkelstein)",
        "difficulty": 0.1,
    },
    "Aromatic_Halogen_Exchange": {
        "taxonomy_family_id": "Aromatic_halogen_exchange",
        "category": "Halogenation",
        "description": "Aromatic halogen exchange (halex)",
        "difficulty": 0.25,
    },
    "NBS_NCS_NIS_halogenation": {
        "taxonomy_family_id": "NBS_NCS_NIS_halogenation",
        "category": "Halogenation",
        "description": "NBS/NCS/NIS halogenation (benzylic, allylic, aromatic)",
        "difficulty": 0.2,
    },
    "Chlorination_SOCl2_oxalyl_chloride": {
        "taxonomy_family_id": "Chlorination_SOCl2",
        "category": "Halogenation",
        "description": "Chlorination with SOCl2 / oxalyl chloride",
        "difficulty": 0.15,
    },
    "Acyl_Halides_formation": {
        "taxonomy_family_id": "Acyl_halide_formation",
        "category": "Halogenation",
        "description": "Acyl halide formation (acid → acid chloride)",
        "difficulty": 0.1,
    },
    "Addition_of_Halogens_to_Double_or_Triple_Bonds": {
        "taxonomy_family_id": "Halogen_addition_to_alkene",
        "category": "Halogenation",
        "description": "Halogen addition to double/triple bonds",
        "difficulty": 0.1,
    },
    "Allylic_Benzylic_and_Vinylic_Halogenations": {
        "taxonomy_family_id": "NBS_NCS_NIS_halogenation",
        "category": "Halogenation",
        "description": "Allylic/benzylic/vinylic halogenation",
        "difficulty": 0.2,
    },
    # ---- Hydroboration ----
    "Brown_Hydroboration": {
        "taxonomy_family_id": "Brown_hydroboration",
        "category": "Hydroboration",
        "description": "Brown hydroboration (anti-Markovnikov)",
        "difficulty": 0.2,
    },
    "Ircatalyzed_CH_borylation": {
        "taxonomy_family_id": "Ircatalyzed_CH_borylation",
        "category": "C-B Bond Forming",
        "description": "Ir-catalyzed C-H borylation",
        "difficulty": 0.35,
    },
    # ---- Heterocycle Formation ----
    "Fischer_indole_synthesis": {
        "taxonomy_family_id": "Fischer_indole_synthesis",
        "category": "Heterocycle Formation",
        "description": "Fischer indole synthesis (arylhydrazine + ketone)",
        "difficulty": 0.35,
    },
    "PaalKnorr_synthesis": {
        "taxonomy_family_id": "PaalKnorr_synthesis",
        "category": "Heterocycle Formation",
        "description": "Paal-Knorr pyrrole/furan/thiophene synthesis",
        "difficulty": 0.25,
    },
    "Knorr_pyrrole_synthesis": {
        "taxonomy_family_id": "Knorr_pyrrole_synthesis",
        "category": "Heterocycle Formation",
        "description": "Knorr pyrrole synthesis (beta-ketoester + amine)",
        "difficulty": 0.3,
    },
    "Hantzsch_synthesis": {
        "taxonomy_family_id": "Hantzsch_synthesis",
        "category": "Heterocycle Formation",
        "description": "Hantzsch dihydropyridine / thiazole synthesis",
        "difficulty": 0.3,
    },
    "Skraup_DoebnerMiller": {
        "taxonomy_family_id": "Skraup_DoebnerMiller",
        "category": "Heterocycle Formation",
        "description": "Skraup / Doebner-Miller quinoline synthesis",
        "difficulty": 0.4,
    },
    "PictetSpengler": {
        "taxonomy_family_id": "PictetSpengler",
        "category": "Heterocycle Formation",
        "description": "Pictet-Spengler tetrahydroisoquinoline synthesis",
        "difficulty": 0.3,
    },
    "Gewald_reaction": {
        "taxonomy_family_id": "Gewald_reaction",
        "category": "Heterocycle Formation",
        "description": "Gewald 2-aminothiophene synthesis",
        "difficulty": 0.35,
    },
    "Beckmann_rearrangement": {
        "taxonomy_family_id": "Beckmann_rearrangement",
        "category": "Rearrangement",
        "description": "Beckmann rearrangement (oxime → lactam/amide)",
        "difficulty": 0.3,
    },
    # ---- Multi-component ----
    "Biginelli_condensation": {
        "taxonomy_family_id": "Biginelli_condensation",
        "category": "Multi-Component",
        "description": "Biginelli DHPM condensation (aldehyde + urea + ketoester)",
        "difficulty": 0.3,
    },
    "Passerini_3component": {
        "taxonomy_family_id": "Passerini_3CR",
        "category": "Multi-Component",
        "description": "Passerini 3-component reaction (acid + aldehyde + isocyanide)",
        "difficulty": 0.35,
    },
    "Ugi_4-component_reaction": {
        "taxonomy_family_id": "Ugi_4CR",
        "category": "Multi-Component",
        "description": "Ugi 4-component reaction (amine + aldehyde + acid + isocyanide)",
        "difficulty": 0.4,
    },
    "GroebkeBlackburnBienaym": {
        "taxonomy_family_id": "GroebkeBlackburnBienayme",
        "category": "Multi-Component",
        "description": "Groebke-Blackburn-Bienayme imidazo[1,2-a]heterocycle",
        "difficulty": 0.4,
    },
    "KabachnikFields": {
        "taxonomy_family_id": "KabachnikFields",
        "category": "Multi-Component",
        "description": "Kabachnik-Fields alpha-aminophosphonate synthesis",
        "difficulty": 0.3,
    },
    # ---- Protecting Group / Functional Group Interconversion ----
    "Formation_of_Oximes": {
        "taxonomy_family_id": "Oxime_formation",
        "category": "Functional Group Interconversion",
        "description": "Oxime formation (carbonyl + NH2OH)",
        "difficulty": 0.1,
    },
    "Oximehy_drazone_formation": {
        "taxonomy_family_id": "Oxime_formation",
        "category": "Functional Group Interconversion",
        "description": "Oxime / hydrazone formation",
        "difficulty": 0.1,
    },
    "Addition_of_Hydrazine_Derivatives_to_Carbonyl_Compounds": {
        "taxonomy_family_id": "Hydrazone_formation",
        "category": "Functional Group Interconversion",
        "description": "Hydrazine derivative addition to carbonyls (hydrazone)",
        "difficulty": 0.15,
    },
    "Hydrazinolysis": {
        "taxonomy_family_id": "Hydrazinolysis",
        "category": "Functional Group Interconversion",
        "description": "Hydrazinolysis (phthalimide deprotection)",
        "difficulty": 0.15,
    },
    "Nitrile_hydration_Nitrile__amide": {
        "taxonomy_family_id": "Nitrile_hydration",
        "category": "Functional Group Interconversion",
        "description": "Nitrile hydration (nitrile → amide)",
        "difficulty": 0.2,
    },
    "BuchererBergs": {
        "taxonomy_family_id": "BuchererBergs",
        "category": "Heterocycle Formation",
        "description": "Bucherer-Bergs hydantoin synthesis",
        "difficulty": 0.3,
    },
    # ---- Miscellaneous / general ----
    "Diazotransfer": {
        "taxonomy_family_id": "Diazotransfer",
        "category": "Functional Group Interconversion",
        "description": "Diazotransfer (amine → azide)",
        "difficulty": 0.25,
    },
    "protocols_all_v2_hte": {
        "taxonomy_family_id": "_general_hte",
        "category": "General HTE",
        "description": "General HTE protocol templates",
        "difficulty": 0.5,
    },
}

# ---------------------------------------------------------------------------
# Quality filters
# ---------------------------------------------------------------------------
MIN_QUALITY = 0.7
MIN_COUNT = 10
REQUIRE_RETRO_VALIDATED = True


def _normalize_source(source: str) -> str:
    """Strip _canonical.csv suffix from source file name."""
    return source.replace("_canonical.csv", "").replace(".csv", "")


def _reverse_smarts(forward_smarts: str) -> str | None:
    """Reverse a forward SMARTS (reactants>>products) to retro (product>>reactants).

    Simply swaps sides around >>. For complex multi-fragment products this
    may need manual curation, but for the vast majority of 2-component
    reactions it works correctly.
    """
    if ">>" not in forward_smarts:
        return None
    parts = forward_smarts.split(">>")
    if len(parts) != 2:
        return None
    reactants, products = parts
    return f"{products}>>{reactants}"


def _make_name(family_id: str, index: int) -> str:
    """Generate a unique template name from family ID and index."""
    # Sanitize: lowercase, collapse non-alpha to underscore
    base = re.sub(r"[^a-zA-Z0-9]+", "_", family_id).strip("_").lower()
    return f"{base}_ext{index}"


def _count_mapped_atoms(smarts: str) -> int:
    """Count unique atom-map numbers in a SMARTS string."""
    return len(set(re.findall(r":([0-9]+)", smarts)))


def _smarts_has_heteroatom(smarts: str) -> bool:
    """Check if SMARTS reactant side contains heteroatom references.

    Templates that only match C/c atoms on the reactant side are much more
    likely to produce false-positive matches on simple molecules.
    """
    # Only check the reactants side (left of >>)
    reactant_side = smarts.split(">>")[0] if ">>" in smarts else smarts
    heteroatoms = re.findall(r'\[?(?:#[0-9]+|[NOSFBPInKZ]|Cl|Br|Se|Si)', reactant_side)
    return len(heteroatoms) > 0


def _determine_n_precursors(retro_smarts: str) -> int:
    """Count the number of precursor fragments (products side of retro SMARTS)."""
    if ">>" not in retro_smarts:
        return 1
    product_side = retro_smarts.split(">>")[1]
    return product_side.count(".") + 1


def _primary_source(template: dict) -> str:
    """Return the primary (most specific) source family for a template."""
    sources = template.get("source_files", [])
    # Prefer non-general sources
    for s in sources:
        norm = _normalize_source(s)
        if norm != "protocols_all_v2_hte":
            return norm
    return _normalize_source(sources[0]) if sources else ""


def _primary_source_from_hte_families(template: dict) -> str:
    """Get the primary source key from a converted template's hte_families."""
    fams = template.get("hte_families", [])
    return fams[0] if fams else ""


def integrate() -> tuple[int, int]:
    """Main integration routine.

    Reads ``results/extracted_templates.json``, filters for quality, maps
    source families to taxonomy IDs, and merges into production JSON files.

    Returns:
        Tuple of (n_new_templates, n_new_retrons) added.
    """
    # Load inputs
    with open(EXTRACTED_FILE, "r", encoding="utf-8") as f:
        extracted = json.load(f)
    with open(HTE_FILE, "r", encoding="utf-8") as f:
        hte_data = json.load(f)
    with open(RETRON_FILE, "r", encoding="utf-8") as f:
        retron_data = json.load(f)

    existing_templates = hte_data.get("templates", [])
    existing_retrons = retron_data.get("retrons", [])

    # Collect existing template names and retro_smarts for dedup
    existing_names = {t["name"] for t in existing_templates}
    existing_retro_smarts = {t.get("retro_smarts", "") for t in existing_templates}

    # Filter high-quality templates
    candidates = [
        t for t in extracted["templates"]
        if t.get("retro_validated", False)
        and t.get("quality_score", 0) >= MIN_QUALITY
        and t.get("count", 0) >= MIN_COUNT
    ]
    print(f"Candidates after quality filter: {len(candidates)}")

    # Skip templates from _general_hte source exclusively
    candidates = [
        t for t in candidates
        if _primary_source(t) != "protocols_all_v2_hte"
    ]
    print(f"Candidates after removing general-only: {len(candidates)}")

    # Group by (cluster_id, primary_source) → pick best quality_score
    cluster_key = {}
    for t in candidates:
        primary = _primary_source(t)
        key = (t.get("cluster_id", id(t)), primary)
        if key not in cluster_key or t.get("quality_score", 0) > cluster_key[key].get("quality_score", 0):
            cluster_key[key] = t
    deduped = list(cluster_key.values())
    print(f"After cluster dedup (best per cluster+family): {len(deduped)}")

    # Sort by quality descending
    deduped.sort(key=lambda t: -t.get("quality_score", 0))

    # Convert to production format & integrate
    new_templates = []
    family_counters: dict[str, int] = defaultdict(int)
    skipped_no_mapping = 0
    skipped_dup_smarts = 0

    for ext_t in deduped:
        primary = _primary_source(ext_t)
        mapping = SOURCE_TO_TAXONOMY.get(primary)
        if not mapping:
            skipped_no_mapping += 1
            continue

        # Reverse the forward SMARTS to get retro SMARTS
        forward_smarts = ext_t["reaction_smarts"]
        retro_smarts = _reverse_smarts(forward_smarts)
        if not retro_smarts:
            continue

        # Skip templates with too few mapped atoms (overly general patterns)
        # Carbon-only templates need more mapped atoms to be selective
        n_mapped = _count_mapped_atoms(forward_smarts)
        has_hetero = _smarts_has_heteroatom(forward_smarts)
        if has_hetero and n_mapped < 6:
            continue
        if not has_hetero and n_mapped < 10:
            continue

        # Skip if we already have this exact retro SMARTS
        if retro_smarts in existing_retro_smarts:
            skipped_dup_smarts += 1
            continue

        fam_id = mapping["taxonomy_family_id"]
        family_counters[fam_id] += 1
        idx = family_counters[fam_id]

        name = _make_name(fam_id, idx)
        while name in existing_names:
            idx += 1
            family_counters[fam_id] = idx
            name = _make_name(fam_id, idx)

        n_prec = _determine_n_precursors(retro_smarts)
        # Cap n_reactants at 2 for forward reactor compatibility
        n_reactants = min(n_prec, 2)

        # Compute difficulty from yield data (higher yield → lower difficulty)
        avg_yield = ext_t.get("avg_yield", 50)
        yield_difficulty = round(max(0.05, 1.0 - avg_yield / 100.0), 2)
        # Blend with family default difficulty
        base_difficulty = mapping["difficulty"]
        difficulty = round((base_difficulty + yield_difficulty) / 2, 2)

        # Build hte_families list from source files
        hte_families = []
        for s in ext_t.get("source_files", []):
            norm = _normalize_source(s)
            if norm not in hte_families and norm != "protocols_all_v2_hte":
                hte_families.append(norm)

        template = {
            "category": mapping["category"],
            "name": name,
            "hte_families": hte_families,
            "taxonomy_id": fam_id,
            "retro_smarts": retro_smarts,
            "forward_smarts": forward_smarts,
            "description": mapping["description"],
            "difficulty": difficulty,
            "n_precursors": n_prec,
            "n_reactants": n_reactants,
            "taxonomy_family_id": fam_id,
            "source": "extracted",
            "quality_score": ext_t.get("quality_score", 0),
            "example_count": ext_t.get("count", 0),
            "avg_yield": avg_yield,
        }
        new_templates.append(template)
        existing_names.add(name)
        existing_retro_smarts.add(retro_smarts)

    print(f"\nSkipped (no source mapping): {skipped_no_mapping}")
    print(f"Skipped (duplicate SMARTS): {skipped_dup_smarts}")
    print(f"New templates to add: {len(new_templates)}")

    # -- Summary by family --
    family_counts = defaultdict(int)
    for t in new_templates:
        family_counts[t["taxonomy_family_id"]] += 1
    print(f"\nNew templates by taxonomy_family_id:")
    for fam, cnt in sorted(family_counts.items()):
        is_new = fam not in {t.get("taxonomy_family_id") for t in existing_templates}
        marker = " [NEW FAMILY]" if is_new else ""
        print(f"  {fam:45s}: {cnt:4d}{marker}")

    # Backup originals
    shutil.copy2(HTE_FILE, HTE_FILE.with_suffix(".json.bak"))
    shutil.copy2(RETRON_FILE, RETRON_FILE.with_suffix(".json.bak"))
    print(f"\nBackups created: *.json.bak")

    # -- Write updated hte_templates.json --
    merged_templates = existing_templates + new_templates
    hte_data["templates"] = merged_templates
    hte_data["description"] = (
        "HTE-backed retrosynthetic template library. "
        "Contains hand-curated templates + high-quality extracted templates (source='extracted'). "
        "Each template applies AllChem.RunReactants() to generate validated precursor SMILES pairs."
    )
    with open(HTE_FILE, "w", encoding="utf-8") as f:
        json.dump(hte_data, f, indent=2, ensure_ascii=False)
    print(f"Updated {HTE_FILE.name}: {len(existing_templates)} → {len(merged_templates)} templates")

    # -- Add retron patterns for NEW taxonomy families --
    existing_retron_fams = {r.get("taxonomy_family_id") for r in existing_retrons}
    new_retron_fams = set()
    for t in new_templates:
        fam = t["taxonomy_family_id"]
        if fam not in existing_retron_fams and fam not in new_retron_fams:
            new_retron_fams.add(fam)

    new_retrons = []
    for fam in sorted(new_retron_fams):
        # Find a representative template for this family
        rep = next(t for t in new_templates if t["taxonomy_family_id"] == fam)
        mapping = SOURCE_TO_TAXONOMY.get(_primary_source_from_hte_families(rep))
        if not mapping:
            # Fallback: find mapping from any source matching the fam
            for src, m in SOURCE_TO_TAXONOMY.items():
                if m["taxonomy_family_id"] == fam:
                    mapping = m
                    break
        if not mapping:
            continue

        # Extract a product_smarts from the retro_smarts (left side of >>)
        retro = rep["retro_smarts"]
        product_smarts = retro.split(">>")[0] if ">>" in retro else retro

        retron = {
            "category": mapping["category"],
            "name": re.sub(r"[^a-zA-Z0-9]+", "_", fam).strip("_").lower(),
            "product_smarts": product_smarts,
            "reaction_name": fam.lower(),
            "taxonomy_id": fam,
            "difficulty": mapping["difficulty"],
            "description": mapping["description"],
            "precursor_hints": [],
            "notes": f"Auto-generated retron from extracted templates (quality >= {MIN_QUALITY}).",
            "taxonomy_family_id": fam,
            "retro_transform_id": fam.lower(),
            "source": "extracted",
        }
        new_retrons.append(retron)

    if new_retrons:
        merged_retrons = existing_retrons + new_retrons
        retron_data["retrons"] = merged_retrons
        with open(RETRON_FILE, "w", encoding="utf-8") as f:
            json.dump(retron_data, f, indent=2, ensure_ascii=False)
        print(f"Updated {RETRON_FILE.name}: {len(existing_retrons)} → {len(merged_retrons)} retrons (+{len(new_retrons)} new families)")
    else:
        print(f"No new retron families to add.")

    print(f"\nIntegration complete!")
    return len(new_templates), len(new_retrons)


if __name__ == "__main__":
    n_templates, n_retrons = integrate()
