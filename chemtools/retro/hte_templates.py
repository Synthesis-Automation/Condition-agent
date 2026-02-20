"""
HTE-backed Retrosynthetic Template Library.

Each template is a retrosynthetic SMARTS in the form:
    [product_pattern] >> [precursor_1] . [precursor_2]

Applied via AllChem.RunReactants((target_mol,)) — the target is the "reactant"
and the output tuples are the precursor pairs.  Atom-mapped atoms (:1, :2, ...)
route the correct molecular context to each precursor fragment; the surrounding
un-mapped atoms follow their mapped neighbours automatically.

Template design rules:
  • Atom mapping must be minimal but unambiguous — only map the atoms at the
    bond being disconnected (1–2 atoms at the reaction centre).
  • Un-mapped substituents in the product SMARTS become NEW atoms (caps like
    Br, F, OH).  They should accurately model the precursor functional group.
  • Templates are stored with their HTE family name so they can be cross-
    referenced against the ~231k-reaction HTE database for conditions/yield.
  • Two-component disconnections only (1 target → 2 precursors).  One-component
    retro-transforms (e.g., reduction) return a single precursor.
  • difficulty follows the retron_patterns.py scale: 0 = trivial, 1 = heroic.

Coverage map (HTE families → template names):
  New families not in the 46 hardcoded retrons are marked with ★.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Template library
# Each entry:
#   name          : unique identifier
#   hte_families  : list of HTE CSV family stem(s) this template covers
#   retro_smarts  : retrosynthetic SMARTS (product >> precursor(s))
#   description   : human-readable
#   difficulty    : 0.0–1.0
#   n_precursors  : 1 or 2 (output fragments from RunReactants)
# ---------------------------------------------------------------------------

HTE_TEMPLATES: List[Dict[str, Any]] = [

    # ── C–N Bond Forming ────────────────────────────────────────────────────

    {
        "name": "buchwald_hartwig",
        "hte_families": ["C_N_Coupling"],
        "retro_smarts": "[c:1]-[NX3;!$(NC=O):2]>>[c:1][Br].[N:2]",
        "description": "Aryl amine ← aryl bromide + amine (Buchwald-Hartwig/Pd)",
        "difficulty": 0.35,
        "n_precursors": 2,
    },
    {   # ★ not in retrons
        "name": "snar_amination",
        "hte_families": ["SNAr_amination"],
        "retro_smarts": "[c:1]-[NX3;!$(NC=O):2]>>[c:1][F].[N:2]",
        "description": "Aryl amine ← electron-poor aryl fluoride + amine (SNAr)",
        "difficulty": 0.20,
        "n_precursors": 2,
        "notes": "Works best on electron-poor arenes (pyridines, nitro-aryls). "
                 "Typically no metal catalyst needed. K2CO3 or Cs2CO3 base, DMSO solvent.",
    },
    {   # ★
        "name": "chan_lam_n_arylation",
        "hte_families": ["ChanLam_Narylation"],
        "retro_smarts": "[c:1]-[NR:2]>>[c:1]B(O)O.[N:2]",
        "description": "N-aryl amine ← arylboronic acid + amine (Chan-Lam, Cu-catalysed)",
        "difficulty": 0.30,
        "n_precursors": 2,
        "notes": "Cu(OAc)2, pyridine or Et3N, open air or O2. Mild conditions, room temp. "
                 "Secondary amines give better yields than primary.",
    },
    {   # ★
        "name": "chan_lam_sulfonamide",
        "hte_families": ["ChanLam_sulfonamide_converted"],
        "retro_smarts": "[c:1]-[NX3:2]-[SX4](=O)(=O)>>[c:1]B(O)O.[N:2][S](=O)(=O)",
        "description": "N-aryl sulfonamide ← arylboronic acid + sulfonamide (Chan-Lam)",
        "difficulty": 0.30,
        "n_precursors": 2,
    },
    {
        "name": "amide_coupling",
        "hte_families": ["Amide_formation", "CDI-mediated_amidation", "Anhydride_coupling",
                         "Weinreb_amide"],
        "retro_smarts": "[C:1](=[O:2])-[NX3:3]>>[C:1](=[O:2])O.[N:3]",
        "description": "Amide ← carboxylic acid + amine (HATU/EDC/DIC coupling)",
        "difficulty": 0.15,
        "n_precursors": 2,
    },
    {
        "name": "reductive_amination",
        "hte_families": ["Reductive_amination", "Reduction_NaBHOAc3_NaBH3CN"],
        "retro_smarts": "[CH:1]([NX3;!$(NC=O):2])>>[C:1]=O.[N:2]",
        "description": "Secondary/tertiary amine ← aldehyde/ketone + amine (reductive amination)",
        "difficulty": 0.20,
        "n_precursors": 2,
        "notes": "NaBH(OAc)3 or NaBH3CN as reductant. pH 5–7 (AcOH). Works for primary and "
                 "secondary amines. Steric hindrance on ketone reduces yield.",
    },
    {   # ★
        "name": "urea_formation",
        "hte_families": ["Urea_and_thiourea_formation"],
        "retro_smarts": "[NX3:1]-[C:2](=[O,S:3])-[NX3:4]>>[N:1].[C:2](=[O,S:3])=[N:4]",
        "description": "Urea/thiourea ← amine + isocyanate/isothiocyanate",
        "difficulty": 0.15,
        "n_precursors": 2,
        "notes": "The isocyanate fragment may be from CDI activation or phosgene-derived.",
    },
    {   # ★
        "name": "sulfonamide_from_sulfonyl_chloride",
        "hte_families": ["Sulfonamide_formation"],
        "retro_smarts": "[NX3:1]-[SX4:2](=O)(=O)>>[N:1].[S:2](=O)(=O)Cl",
        "description": "Sulfonamide ← sulfonyl chloride + amine",
        "difficulty": 0.10,
        "n_precursors": 2,
    },
    {   # ★
        "name": "gabriel_amine",
        "hte_families": ["Gabriel_amine_synthesis"],
        "retro_smarts": "[CX4:1]-[NH2:2]>>[C:1][Br].[N:2]",
        "description": "Primary amine ← alkyl bromide + phthalimide (Gabriel synthesis)",
        "difficulty": 0.25,
        "n_precursors": 2,
    },
    {   # ★
        "name": "curtius_rearrangement",
        "hte_families": ["Curtius_rearrangement", "Hofmann_rearrangement"],
        "retro_smarts": "[NX3:1][C:2](=O)>>[N:1]=[C:2]=O",
        "description": "Carbamate/amine ← isocyanate intermediate (Curtius/Hofmann)",
        "difficulty": 0.35,
        "n_precursors": 1,
    },

    # ── C–C Bond Forming ────────────────────────────────────────────────────

    {
        "name": "suzuki_miyaura",
        "hte_families": ["suzuki_miyaura"],
        "retro_smarts": "[c:1]-[c:2]>>[c:1][Br].[c:2]B(O)O",
        "description": "Biaryl bond ← aryl bromide + arylboronic acid (Suzuki)",
        "difficulty": 0.15,
        "n_precursors": 2,
    },
    {   # ★
        "name": "decarboxylative_arylation",
        "hte_families": ["Decarboxylative_arylation"],
        "retro_smarts": "[c:1]-[CX4:2]>>[c:1][Br].[C:2]C(=O)O",
        "description": "Aryl-alkyl bond ← aryl halide + carboxylic acid (decarboxylative)",
        "difficulty": 0.50,
        "n_precursors": 2,
    },
    {   # ★
        "name": "heck_mizoroki",
        "hte_families": ["HeckMizoroki_coupling"],
        "retro_smarts": "[c:1]-[CH:2]=[CH2:3]>>[c:1][Br].[CH2:2]=[CH2:3]",
        "description": "Aryl vinyl ← aryl bromide + alkene (Heck/Mizoroki)",
        "difficulty": 0.30,
        "n_precursors": 2,
    },
    {   # ★
        "name": "stille_coupling",
        "hte_families": ["Stille_coupling"],
        "retro_smarts": "[c:1]-[c:2]>>[c:1][Br].[c:2][Sn](CC)(CC)CC",
        "description": "Biaryl bond ← aryl bromide + arylstannane (Stille)",
        "difficulty": 0.40,
        "n_precursors": 2,
    },
    {   # ★
        "name": "negishi_coupling",
        "hte_families": ["Negishi_coupling"],
        "retro_smarts": "[c:1]-[CX4:2]>>[c:1][Br].[C:2][Zn]Cl",
        "description": "Aryl–alkyl bond ← aryl bromide + alkylzinc (Negishi)",
        "difficulty": 0.40,
        "n_precursors": 2,
    },
    {   # ★
        "name": "kumada_coupling",
        "hte_families": ["Kumada_coupling"],
        "retro_smarts": "[c:1]-[c:2]>>[c:1][Br].[c:2][Mg]Br",
        "description": "Biaryl bond ← aryl halide + aryl Grignard (Kumada)",
        "difficulty": 0.45,
        "n_precursors": 2,
    },
    {   # ★
        "name": "hiyama_coupling",
        "hte_families": ["Hiyama_coupling"],
        "retro_smarts": "[c:1]-[c:2]>>[c:1][Br].[c:2][Si](OCC)(OCC)OCC",
        "description": "Biaryl bond ← aryl halide + arylsilane (Hiyama)",
        "difficulty": 0.45,
        "n_precursors": 2,
    },
    {   # ★
        "name": "liebeskind_srogl",
        "hte_families": ["LiebeskindSrogl_coupling"],
        "retro_smarts": "[C:1](=O)-[c:2]>>[C:1](=O)[S]CC.[c:2]B(O)O",
        "description": "Ketone ← thioester + arylboronic acid (Liebeskind-Srögl)",
        "difficulty": 0.50,
        "n_precursors": 2,
    },
    {   # ★
        "name": "giese_radical",
        "hte_families": ["Giese_radical_additions"],
        "retro_smarts": "[C:1]-[CH2:2]-[C:3](=[O:4])>>[C:1][Br].[CH2:2]=[C:3](=[O:4])",
        "description": "Radical addition product ← alkyl bromide + Michael acceptor (Giese)",
        "difficulty": 0.45,
        "n_precursors": 2,
    },
    {   # ★ — Michael addition
        "name": "michael_addition",
        "hte_families": ["Michael_addition"],
        "retro_smarts": "[#6:1]-[CH2:2]-[CH2:3]-[C:4](=[O:5])>>[#6:1].[C:2]=[C:3]-[C:4]=[O:5]",
        "description": "Michael adduct ← nucleophile + α,β-unsaturated carbonyl (Michael)",
        "difficulty": 0.30,
        "n_precursors": 2,
        "notes": "Nu must have acidic alpha-H (malonate, nitrile, beta-ketoester). "
                 "Base catalyst (NaH, K2CO3, DBU). Solvent: THF or DMF.",
    },
    {   # ★
        "name": "knoevenagel_condensation",
        "hte_families": ["Knoevenagel_condensation"],
        "retro_smarts": "[CX3:1]=[CX3:2]-[C:3](=[O:4])>>[C:1]=O.[CH:2]-[C:3]=[O:4]",
        "description": "α,β-Unsaturated carbonyl ← aldehyde + active methylene (Knoevenagel)",
        "difficulty": 0.20,
        "n_precursors": 2,
    },
    {   # ★
        "name": "aldol_classic",
        "hte_families": ["Aldol_classic__Mukaiyama"],
        "retro_smarts": "[C:1]([OH:2])-[CH2:3]-[C:4]=O>>[C:1]=O.[CH2:3]-[C:4]=O",
        "description": "β-Hydroxy carbonyl ← two carbonyls (aldol addition)",
        "difficulty": 0.30,
        "n_precursors": 2,
    },
    {   # ★
        "name": "horner_wadsworth_emmons",
        "hte_families": ["HornerWadsworthEmmons"],
        "retro_smarts": "[C:1]=[CH:2]-[C:3](=O)>>[C:1]=O.[CH2:2]-[C:3](=O)",
        "description": "α,β-Unsaturated ester ← aldehyde + phosphonate (HWE)",
        "difficulty": 0.20,
        "n_precursors": 2,
    },
    {   # ★
        "name": "nozaki_hiyama_kishi",
        "hte_families": ["NozakiHiyamaKishi"],
        "retro_smarts": "[C:1]([OH:2])-[CH:3]=[CH2:4]>>[C:1]=O.[CH:3]=[CH2:4]",
        "description": "Allylic alcohol ← aldehyde + vinyl halide (NHK, Cr/Ni)",
        "difficulty": 0.50,
        "n_precursors": 2,
    },

    # ── C–O Bond Forming ────────────────────────────────────────────────────

    {
        "name": "c_o_coupling_pd",
        "hte_families": ["C_O_Coupling"],
        "retro_smarts": "[c:1]-[OX2:2]-[CX4:3]>>[c:1][Br].[O:2][C:3]",
        "description": "Aryl ether ← aryl bromide + alcohol (Pd-catalysed C–O coupling)",
        "difficulty": 0.40,
        "n_precursors": 2,
    },
    {
        "name": "williamson_ether",
        "hte_families": ["Williamson_ether_synthesis"],
        "retro_smarts": "[CX4:1]-[OX2:2]-[CX4:3]>>[C:1][Br].[O:2][C:3]",
        "description": "Ether ← alkyl halide + alkoxide (Williamson)",
        "difficulty": 0.15,
        "n_precursors": 2,
    },
    {
        "name": "fischer_esterification",
        "hte_families": ["Fischer_esterification", "Steglich_esterification"],
        "retro_smarts": "[C:1](=[O:2])-[OX2:3]-[CX4:4]>>[C:1](=[O:2])O.[O:3][C:4]",
        "description": "Ester ← carboxylic acid + alcohol (Fischer/Steglich)",
        "difficulty": 0.10,
        "n_precursors": 2,
    },
    {   # ★
        "name": "mitsunobu_etherification",
        "hte_families": ["Mitsunobu_etherificationinversion"],
        "retro_smarts": "[c:1]-[OX2:2]-[CX4:3]>>[c:1][OH].[O:2][C:3]",
        "description": "Aryl ether ← phenol + alcohol with inversion (Mitsunobu)",
        "difficulty": 0.25,
        "n_precursors": 2,
    },

    # ── C–S Bond Forming ────────────────────────────────────────────────────

    {   # ★
        "name": "c_s_coupling",
        "hte_families": ["C_S_Coupling", "Thioether_formation"],
        "retro_smarts": "[c:1]-[SX2:2]>>[c:1][Br].[S:2]",
        "description": "Aryl thioether ← aryl bromide + thiol (Pd or Cu catalysed)",
        "difficulty": 0.40,
        "n_precursors": 2,
        "notes": "Pd(OAc)2/Xantphos or CuI/diamine. Base: Cs2CO3. Watch for disulfide oxidation.",
    },
    {   # ★
        "name": "thiol_ene",
        "hte_families": ["Thiolene_and_thiolyne"],
        "retro_smarts": "[C:1]-[CH2:2]-[SX2:3]>>[C:1]=[CH2:2].[S:3]",
        "description": "Thioether ← alkene + thiol (radical thiol-ene)",
        "difficulty": 0.20,
        "n_precursors": 2,
    },

    # ── Click Chemistry ─────────────────────────────────────────────────────

    {   # ★
        "name": "cuaac_triazole",
        "hte_families": ["CuAAC_azidealkyne"],
        "retro_smarts": "[n:1]1[n:2][n:3][c:4][c:5]1>>[N-:3]=[N+:2]=[N:1].[C:4]#[C:5]",
        "description": "1,2,3-Triazole ← organic azide + terminal alkyne (CuAAC click)",
        "difficulty": 0.10,
        "n_precursors": 2,
        "notes": "CuI or Cu(OAc)2/sodium ascorbate. Water/t-BuOH or DMSO. Room temperature. "
                 "Highly regioselective → 1,4-substituted triazole. Bioorthogonal.",
    },

    # ── Carbonyl Reductions (1-precursor retro-transforms) ──────────────────

    {   # ★ — retro of NaBH4 reduction
        "name": "ketone_from_nabh4",
        "hte_families": ["NaBH4_carbonyl_reductions"],
        "retro_smarts": "[CH:1]([OH:2])>>[C:1]=O",
        "description": "Secondary alcohol came from NaBH4 reduction of ketone",
        "difficulty": 0.05,
        "n_precursors": 1,
        "notes": "NaBH4 in MeOH/EtOH or THF. Selective for ketones over esters at −78 °C. "
                 "For aldehydes, same conditions at 0 °C.",
    },
    {   # ★ — retro of LAH/DIBAL reduction
        "name": "ester_from_lah_reduction",
        "hte_families": ["Reduction_LAH_LiAlH4", "DIBALH_Partial_reductions"],
        "retro_smarts": "[CX4:1]([OH:2])-[CX4:3]>>[C:1](=O)-[C:3]",
        "description": "Primary/secondary alcohol came from LAH reduction of ester/acid",
        "difficulty": 0.15,
        "n_precursors": 1,
        "notes": "LiAlH4 in Et2O or THF at 0 °C. Not compatible with protic solvents. "
                 "DIBAL-H at −78 °C can stop at aldehyde stage.",
    },

    # ── Halogenation/Fluorination ────────────────────────────────────────────

    {   # ★
        "name": "deoxy_fluorination",
        "hte_families": ["Deoxy_fluorination"],
        "retro_smarts": "[CX4:1][F:2]>>[C:1][OH:2]",
        "description": "Alkyl fluoride came from deoxyfluorination of alcohol (DAST/Deoxofluor)",
        "difficulty": 0.30,
        "n_precursors": 1,
        "notes": "DAST or Deoxofluor in DCM, −78 °C → RT. Inverts configuration at secondary C. "
                 "Not compatible with acid-sensitive groups. Watch for elimination side products.",
    },
    {   # ★
        "name": "electrophilic_fluorination",
        "hte_families": ["Electrophilic_fluorination"],
        "retro_smarts": "[c:1][F:2]>>[c:1][H]",
        "description": "Aryl fluoride came from electrophilic fluorination of arene (Selectfluor)",
        "difficulty": 0.40,
        "n_precursors": 1,
        "notes": "Selectfluor or NFSI. Directed by electron-rich arenes. Regioselectivity "
                 "can be challenging without directing group.",
    },
    {   # ★
        "name": "sandmeyer_bromide",
        "hte_families": ["Sandmeyer_reactions"],
        "retro_smarts": "[c:1][Br:2]>>[c:1][NH2]",
        "description": "Aryl bromide came from Sandmeyer reaction on ArNH2 via diazonium",
        "difficulty": 0.35,
        "n_precursors": 1,
        "notes": "NaNO2/HBr then CuBr. Requires aniline precursor. "
                 "Reaction via diazonium salt at 0–5 °C.",
    },
    {   # ★
        "name": "balz_schiemann",
        "hte_families": ["BalzSchiemann"],
        "retro_smarts": "[c:1][F:2]>>[c:1][NH2]",
        "description": "Aryl fluoride came from Balz-Schiemann (ArNH2 → ArF via diazonium BF4)",
        "difficulty": 0.40,
        "n_precursors": 1,
    },
    {   # ★
        "name": "trifluoromethylation",
        "hte_families": ["Trifluoromethylation", "RuppertPrakash_TMSCF3"],
        "retro_smarts": "[c:1][CF3:2]>>[c:1][Br]",
        "description": "Aryl-CF3 ← aryl bromide (Cu-mediated or Togni reagent)",
        "difficulty": 0.50,
        "n_precursors": 1,
    },

    # ── Oxidations (1-precursor retro) ──────────────────────────────────────

    {   # ★ — Swern/DMP/IBX oxidation
        "name": "aldehyde_from_oxidation",
        "hte_families": ["DessMartin_periodinane_DMP_Alcohols__aldehydesketones",
                         "Swern_oxidation", "IBX", "BAIB"],
        "retro_smarts": "[CH:1]=O>>[CH2:1][OH]",
        "description": "Aldehyde came from oxidation of primary alcohol (DMP/Swern/IBX)",
        "difficulty": 0.10,
        "n_precursors": 1,
    },
    {   # ★
        "name": "ketone_from_oxidation",
        "hte_families": ["Jones_oxidation", "Riley_oxidation",
                         "TPAPNMO_Catalytic_Ru_oxidations",
                         "TEMPObleach_Oxidations_Primary_alcohols__aldehydesacids"],
        "retro_smarts": "[C:1](=O)-[CX4:2]>>[C:1]([OH])-[C:2]",
        "description": "Ketone came from oxidation of secondary alcohol",
        "difficulty": 0.10,
        "n_precursors": 1,
    },

    # ── Miscellaneous ────────────────────────────────────────────────────────

    {   # ★
        "name": "carbamate_formation",
        "hte_families": ["Carbamate_formation"],
        "retro_smarts": "[NX3:1]-[C:2](=O)-[OX2:3]-[CX4:4]>>[N:1].[C:2](=O)=[N:5].[O:3][C:4]",
        "description": "Carbamate ← amine + chloroformate (or isocyanate + alcohol)",
        "difficulty": 0.15,
        "n_precursors": 2,
        "notes": "Use simpler 2-fragment template for amine + chloroformate: "
                 "[NX3:1][C:2](=O)[OX2:3]>>[N:1].[C:2](=O)Cl",
    },
    {   # ★ — simplified carbamate (amine + chloroformate)
        "name": "carbamate_from_amine_chloroformate",
        "hte_families": ["Carbamate_formation"],
        "retro_smarts": "[NX3:1]-[C:2](=O)-[OX2:3]>>[N:1].[C:2](=O)[O:3]Cl",
        "description": "Carbamate ← amine + chloroformate",
        "difficulty": 0.15,
        "n_precursors": 2,
    },
    {   # ★
        "name": "wacker_oxidation",
        "hte_families": ["Wacker_oxidation", "Wacker_oxidation_Alkene__ketone_PdCl2_CuCl_O2"],
        "retro_smarts": "[CH2:1]-[C:2](=O)-[CH3:3]>>[CH2:1]=[CH:2]-[CH3:3]",
        "description": "Methyl ketone ← terminal alkene (Wacker oxidation, Pd/Cu)",
        "difficulty": 0.35,
        "n_precursors": 1,
    },
    {   # ★
        "name": "diazotransfer_azide",
        "hte_families": ["CH_amination_and_azidation"],
        "retro_smarts": "[CX4:1][N3:2]>>[C:1][NH2]",
        "description": "Alkyl azide came from diazotransfer on primary amine",
        "difficulty": 0.30,
        "n_precursors": 1,
    },
    {   # ★
        "name": "borylation_miyaura",
        "hte_families": ["Ircatalyzed_CH_borylation"],
        "retro_smarts": "[c:1][B:2](O)O>>[c:1][H]",
        "description": "Arylboronic acid from C–H borylation (Ir-catalysed)",
        "difficulty": 0.55,
        "n_precursors": 1,
    },
]


# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def get_template_by_name(name: str) -> Optional[Dict[str, Any]]:
    """Return a template dict by name, or None if not found."""
    for t in HTE_TEMPLATES:
        if t["name"] == name:
            return t
    return None


def get_templates_for_family(hte_family: str) -> List[Dict[str, Any]]:
    """Return all templates whose hte_families list contains the given family stem."""
    family_lower = hte_family.lower()
    return [
        t for t in HTE_TEMPLATES
        if any(f.lower() == family_lower for f in t["hte_families"])
    ]


def get_all_template_names() -> List[str]:
    return [t["name"] for t in HTE_TEMPLATES]
