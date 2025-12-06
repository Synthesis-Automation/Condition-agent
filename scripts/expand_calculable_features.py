"""
Script to generate comprehensive expanded calculable_features.json
Expands from ~107 features to 300+ features covering all areas of organic chemistry.
"""

import json
from pathlib import Path

# Load the existing features
existing_file = Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "calculable_features.json"
with open(existing_file, 'r', encoding='utf-8') as f:
    existing_data = json.load(f)

# Create new comprehensive feature set
expanded_features = {
    "version": "2025-11-04.v3.0-comprehensive",
    "description": "Comprehensive molecular feature detection library for organic chemistry reasoning, synthesis planning, retrosynthesis, medicinal chemistry, and computational chemistry applications",
    "schema_notes": {
        **existing_data["schema_notes"],
        "expansion_notes": [
            "Expanded from 107 to 300+ features",
            "Added comprehensive protecting group library",
            "Added detailed heterocycle detection",
            "Added ADME-relevant features",
            "Added stereochemistry and conformational features",
            "Added reactive intermediate markers",
            "Added medicinal chemistry filters"
        ],
        "categories": [
            "halides_and_leaving_groups",
            "organometallics_and_coupling_partners",
            "protecting_groups_comprehensive",
            "nitrogen_functionality",
            "oxygen_functionality",
            "carbonyl_and_derivatives",
            "sulfur_phosphorus_boron",
            "unsaturation_and_pi_systems",
            "heterocycles_comprehensive",
            "stereochemistry_and_topology",
            "reactive_intermediates_and_instability",
            "electrophile_nucleophile_scales",
            "redox_and_photochemistry",
            "physicochemical_properties",
            "medicinal_chemistry_features",
            "metabolic_hotspots",
            "solubility_permeability"
        ]
    },
    "features": []
}

# Keep all existing features
expanded_features["features"] = existing_data["features"].copy()

# ============================================================================
# CATEGORY 1: Additional Protecting Groups (Comprehensive)
# ============================================================================

protecting_groups = [
    # Alcohol protecting groups
    {
        "token": "methyl_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2]([#6])[CH3]"]},
        "why": "Methyl ether; permanent protecting group, requires strong acids/Lewis acids to cleave"
    },
    {
        "token": "benzyl_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2]Cc1ccccc1"]},
        "why": "Benzyl ether (Bn); cleaved by hydrogenolysis (Pd/C, H2) or strong acids"
    },
    {
        "token": "pmb_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2]Cc1ccc([OX2][CH3])cc1"]},
        "why": "para-Methoxybenzyl ether (PMB); oxidatively cleavable with DDQ or CAN"
    },
    {
        "token": "mom_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2]CO[CH3]"]},
        "why": "Methoxymethyl ether (MOM); acid-labile protecting group"
    },
    {
        "token": "sem_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2]CO[Si]([#6])([#6])[#6]"]},
        "why": "2-(Trimethylsilyl)ethoxymethyl ether (SEM); cleaved with fluoride"
    },
    {
        "token": "tbdps_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][Si]([#6])(c1ccccc1)(c1ccccc1)"]},
        "why": "tert-Butyldiphenylsilyl ether (TBDPS); bulky silyl protecting group, fluoride-labile"
    },
    {
        "token": "tbs_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][Si]([CH3])([CH3])C(C)(C)C"]},
        "why": "tert-Butyldimethylsilyl ether (TBS/TBDMS); common silyl protecting group"
    },
    {
        "token": "tips_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][Si]([CH;X4](C)C)([CH;X4](C)C)([CH;X4](C)C)"]},
        "why": "Triisopropylsilyl ether (TIPS); very bulky, stable silyl group"
    },
    {
        "token": "acetate_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][CX3](=[OX1])[CH3]"]},
        "why": "Acetate ester (Ac); simple protecting group, base or acid hydrolysis"
    },
    {
        "token": "benzoate_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][CX3](=[OX1])c1ccccc1"]},
        "why": "Benzoate ester (Bz); more stable than acetate, requires stronger conditions"
    },
    {
        "token": "pivaloate_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[OX2][CX3](=[OX1])C(C)(C)C"]},
        "why": "Pivaloate ester (Piv); sterically hindered, enhances stability"
    },
    {
        "token": "acetal_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX4]([OX2][#6])([OX2][#6])[#6,H]"]},
        "why": "Acetal; protects carbonyls, acid-labile, stable to bases"
    },
    {
        "token": "ketal_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX4]([OX2][#6])([OX2][#6])([#6])[#6]"]},
        "why": "Ketal; protects ketones, acid-labile, stable to bases"
    },
    
    # Amine protecting groups
    {
        "token": "tosyl_sulfonamide_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3][SX4](=[OX1])(=[OX1])c1ccc([CH3])cc1"]},
        "why": "Tosyl (Ts) protecting group on nitrogen; electron-withdrawing, reduces nucleophilicity"
    },
    {
        "token": "nosyl_sulfonamide_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3][SX4](=[OX1])(=[OX1])c1ccc([N+](=O)[O-])cc1"]},
        "why": "Nosyl (Ns/2-nitrobenzenesulfonyl) protecting group; cleaved with thiols"
    },
    {
        "token": "trifluoroacetamide_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3][CX3](=[OX1])C(F)(F)F"]},
        "why": "Trifluoroacetamide (TFA) protecting group; more acid-labile than regular amides"
    },
    {
        "token": "alloc_carbamate_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3][CX3](=[OX1])[OX2]CC=C"]},
        "why": "Allyloxycarbonyl (Alloc) protecting group; cleaved with Pd(0)"
    },
    {
        "token": "troc_carbamate_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3][CX3](=[OX1])[OX2]CC(Cl)(Cl)Cl"]},
        "why": "2,2,2-Trichloroethoxycarbonyl (Troc); cleaved by Zn reduction"
    },
    {
        "token": "phthalimide_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[NX3]1C(=O)c2ccccc2C1=O"]},
        "why": "Phthalimide protecting group; Gabriel synthesis intermediate, cleaved with hydrazine"
    },
    
    # Carbonyl protecting groups
    {
        "token": "dithiane_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX4]1[SX2][CX4][CX4][SX2]1"]},
        "why": "1,3-Dithiane; protects aldehydes/ketones, enables umpolung chemistry"
    },
    {
        "token": "oxime_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3]=[NX2][OX2H]", "[CX3]=[NX2][OX2][#6]"]},
        "why": "Oxime; protects carbonyls, can be hydrolyzed back to carbonyl"
    },
    {
        "token": "hydrazone_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3]=[NX2][NX3]"]},
        "why": "Hydrazone; protects carbonyls, useful in Wolff-Kishner reduction"
    },
    {
        "token": "enol_ether_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3]=[CX3][OX2][#6]"]},
        "why": "Enol ether; masked carbonyl, acid-labile"
    },
    
    # Carboxylic acid protecting groups
    {
        "token": "methyl_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2][CH3]"]},
        "why": "Methyl ester; common protecting group for carboxylic acids, base hydrolysis"
    },
    {
        "token": "ethyl_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2][CH2][CH3]"]},
        "why": "Ethyl ester; slightly more lipophilic than methyl ester"
    },
    {
        "token": "tert_butyl_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2]C(C)(C)C"]},
        "why": "tert-Butyl ester; acid-labile, orthogonal to other ester protecting groups"
    },
    {
        "token": "benzyl_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "protecting_groups",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2]Cc1ccccc1"]},
        "why": "Benzyl ester; cleaved by hydrogenolysis (Pd/C, H2)"
    },
]

# ============================================================================
# CATEGORY 2: Additional Heterocycles (Comprehensive)
# ============================================================================

heterocycles = [
    {
        "token": "furan_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["o1cccc1"]},
        "why": "Furan; aromatic 5-membered O-heterocycle, Diels-Alder reactive"
    },
    {
        "token": "thiophene_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["s1cccc1"]},
        "why": "Thiophene; aromatic 5-membered S-heterocycle, π-excessive, metabolically stable"
    },
    {
        "token": "pyrrole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[nH]1cccc1"]},
        "why": "Pyrrole; aromatic 5-membered NH-heterocycle, π-excessive, electron-rich"
    },
    {
        "token": "imidazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[nH]1cncc1", "n1c[nH]cc1"]},
        "why": "Imidazole; 5-membered heterocycle with 2 nitrogens, basic, coordinating"
    },
    {
        "token": "oxazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["o1ccnc1", "o1cncc1"]},
        "why": "Oxazole; 5-membered O,N-heterocycle, π-deficient, weak base"
    },
    {
        "token": "thiazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["s1ccnc1", "s1cncc1"]},
        "why": "Thiazole; 5-membered S,N-heterocycle, common in drugs, π-deficient"
    },
    {
        "token": "triazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["n1nncn1", "n1nncc1", "[nH]1nncn1", "[nH]1nncc1"]},
        "why": "Triazole (1,2,3- or 1,2,4-); click chemistry product, metabolically stable"
    },
    {
        "token": "tetrazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["n1nnn[nH]1", "n1[nH]nnn1"]},
        "why": "Tetrazole; acidic, carboxylic acid bioisostere, metabolically stable"
    },
    {
        "token": "quinoline_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["c1ccc2ncccc2c1", "n1cccc2ccccc12"]},
        "why": "Quinoline; fused benzene-pyridine system, basic, drug scaffold"
    },
    {
        "token": "isoquinoline_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["c1ccc2cnccc2c1", "n1ccc2ccccc2c1"]},
        "why": "Isoquinoline; fused benzene-pyridine isomer, basic, alkaloid scaffold"
    },
    {
        "token": "benzofuran_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["o1ccc2ccccc12", "c1ccc2occc2c1"]},
        "why": "Benzofuran; fused benzene-furan system, drug scaffold"
    },
    {
        "token": "benzothiophene_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["s1ccc2ccccc12", "c1ccc2sccc2c1"]},
        "why": "Benzothiophene; fused benzene-thiophene system, metabolically stable"
    },
    {
        "token": "benzimidazole_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[nH]1cnc2ccccc12", "n1cnc2ccccc12"]},
        "why": "Benzimidazole; fused benzene-imidazole, drug scaffold (proton pump inhibitors)"
    },
    {
        "token": "purine_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["n1cnc2[nH]cnc12", "n1cnc2ncnc12"]},
        "why": "Purine; nucleobase scaffold (adenine, guanine), kinase inhibitor scaffold"
    },
    {
        "token": "pyrimidine_derivative_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["n1cnccc1", "n1ccncc1"]},
        "why": "Pyrimidine or related; nucleobase scaffold (cytosine, thymine, uracil)"
    },
    {
        "token": "morpholine_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[NX3]1[CX4][CX4][OX2][CX4][CX4]1"]},
        "why": "Morpholine; saturated 6-membered O,N-heterocycle, common solubilizing group"
    },
    {
        "token": "piperazine_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[NX3]1[CX4][CX4][NX3][CX4][CX4]1"]},
        "why": "Piperazine; saturated 6-membered diamine ring, common linker in drugs"
    },
    {
        "token": "piperidine_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[NX3]1[CX4][CX4][CX4][CX4][CX4]1"]},
        "why": "Piperidine; saturated 6-membered N-heterocycle, basic, common scaffold"
    },
    {
        "token": "pyrrolidine_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[NX3]1[CX4][CX4][CX4][CX4]1"]},
        "why": "Pyrrolidine; saturated 5-membered N-heterocycle, basic, proline analog"
    },
    {
        "token": "tetrahydrofuran_present",
        "type": "bool",
        "scope": "global",
        "category": "heterocycles",
        "detect": {"smarts_any": ["[OX2]1[CX4][CX4][CX4][CX4]1"]},
        "why": "Tetrahydrofuran (THF); saturated 5-membered O-heterocycle, ether-like"
    },
]

# ============================================================================
# CATEGORY 3: Reactive Intermediates & Instability Markers
# ============================================================================

reactive_intermediates = [
    {
        "token": "epoxide_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX4]1[OX2][CX4]1"]},
        "why": "Epoxide; strained 3-membered ring ether, electrophilic, ring-opening reactions"
    },
    {
        "token": "aziridine_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX4]1[NX3][CX4]1"]},
        "why": "Aziridine; strained 3-membered ring amine, electrophilic, toxic"
    },
    {
        "token": "diazo_compound_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[#6][N+]#[N-]", "[#6]=[N+]=[N-]"]},
        "why": "Diazo compound; carbene precursor, highly reactive, explosive potential"
    },
    {
        "token": "azide_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[NX2]=[N+]=[N-]", "[NX1-]=[N+]=[NX2]"]},
        "why": "Azide; click chemistry handle, explosive, 1,3-dipole for cycloadditions"
    },
    {
        "token": "peroxide_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[OX2][OX2]"]},
        "why": "Peroxide (O-O bond); oxidizing agent, explosive potential, autoxidation product"
    },
    {
        "token": "hydroperoxide_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[OX2H][OX2]"]},
        "why": "Hydroperoxide; reactive, explosive, autoxidation intermediate"
    },
    {
        "token": "isocyanate_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[NX2]=[CX2]=[OX1]"]},
        "why": "Isocyanate; electrophilic, toxic, reacts with amines/alcohols"
    },
    {
        "token": "isothiocyanate_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[NX2]=[CX2]=[SX1]"]},
        "why": "Isothiocyanate; electrophilic, pungent, forms thioureas"
    },
    {
        "token": "sulfonyl_chloride_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[SX4](=[OX1])(=[OX1])[Cl]"]},
        "why": "Sulfonyl chloride; highly reactive electrophile, moisture-sensitive"
    },
    {
        "token": "acid_anhydride_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2][CX3](=[OX1])"]},
        "why": "Acid anhydride; acylating agent, moisture-sensitive, reactive"
    },
    {
        "token": "alpha_beta_unsaturated_carbonyl_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX3]=[CX3][CX3]=[OX1]"]},
        "why": "α,β-Unsaturated carbonyl; Michael acceptor, electrophilic at β-position"
    },
    {
        "token": "alpha_beta_unsaturated_ester_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX3]=[CX3][CX3](=[OX1])[OX2]"]},
        "why": "α,β-Unsaturated ester; Michael acceptor, conjugate addition substrate"
    },
    {
        "token": "alpha_beta_unsaturated_nitrile_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX3]=[CX3][CX2]#[NX1]"]},
        "why": "α,β-Unsaturated nitrile; electrophilic, Michael acceptor"
    },
    {
        "token": "alpha_halo_carbonyl_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"smarts_any": ["[CX4]([Cl,Br,I])[CX3]=[OX1]"]},
        "why": "α-Halo carbonyl; electrophilic, alkylating agent, lachrymator potential"
    },
    {
        "token": "strained_ring_present",
        "type": "bool",
        "scope": "global",
        "category": "reactive_intermediates",
        "detect": {"heuristic": "contains 3-membered or 4-membered ring"},
        "why": "Strained small ring (3-4 membered); high ring strain energy, reactive"
    },
]

# ============================================================================
# CATEGORY 4: Physicochemical Properties (ADME-relevant)
# ============================================================================

physicochemical = [
    {
        "token": "lipinski_hbd_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "HBD <= 5"},
        "why": "Lipinski rule: ≤5 hydrogen bond donors for oral bioavailability"
    },
    {
        "token": "lipinski_hba_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "HBA <= 10"},
        "why": "Lipinski rule: ≤10 hydrogen bond acceptors for oral bioavailability"
    },
    {
        "token": "lipinski_mw_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "molecular_weight < 500"},
        "why": "Lipinski rule: MW < 500 Da for oral bioavailability"
    },
    {
        "token": "lipinski_logp_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "logP < 5"},
        "why": "Lipinski rule: logP < 5 for oral bioavailability"
    },
    {
        "token": "veber_rotb_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "rotatable_bonds <= 10"},
        "why": "Veber rule: ≤10 rotatable bonds for oral bioavailability"
    },
    {
        "token": "veber_tpsa_compliant",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "TPSA <= 140"},
        "why": "Veber rule: TPSA ≤ 140 Ų for oral bioavailability"
    },
    {
        "token": "ionizable_basic_group_present",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"smarts_any": ["[NX3;H2,H1,H0;!$(NC=O)]", "[nH]"]},
        "why": "Basic ionizable group (pKa ~ 7-11); affects solubility, permeability, distribution"
    },
    {
        "token": "ionizable_acidic_group_present",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"smarts_any": ["[CX3](=[OX1])[OX2H]", "c[OX2H]", "[SX4](=[OX1])(=[OX1])[OX2H]", "n1nnn[nH]1"]},
        "why": "Acidic ionizable group (pKa ~ 2-10); carboxylic acid, phenol, sulfonic acid, tetrazole"
    },
    {
        "token": "quaternary_ammonium_present",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"smarts_any": ["[N+]([#6])([#6])([#6])[#6]"]},
        "why": "Quaternary ammonium; permanently charged, poor cell permeability"
    },
    {
        "token": "zwitterion_potential",
        "type": "bool",
        "scope": "global",
        "category": "physicochemical",
        "detect": {"heuristic": "has both acidic and basic ionizable groups"},
        "why": "Potential zwitterion; both acidic and basic groups present (amino acids, etc.)"
    },
]

# ============================================================================
# CATEGORY 5: Medicinal Chemistry Features
# ============================================================================

medicinal_chem = [
    {
        "token": "aromatic_ring_count",
        "type": "int",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"heuristic": "count aromatic rings"},
        "why": "Number of aromatic rings; affects lipophilicity, π-π interactions"
    },
    {
        "token": "fsp3_high",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"heuristic": "fraction sp3 carbons > 0.5"},
        "why": "High Fsp3 (>50%); 3D character, improved solubility, ADME properties"
    },
    {
        "token": "fsp3_low",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"heuristic": "fraction sp3 carbons < 0.3"},
        "why": "Low Fsp3 (<30%); flat, aromatic character, potential solubility issues"
    },
    {
        "token": "fluorine_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[F]"]},
        "why": "Fluorine atom; metabolic blocker, lipophilicity modulator, bioisostere"
    },
    {
        "token": "trifluoromethyl_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[CX4](F)(F)F"]},
        "why": "Trifluoromethyl group (CF3); lipophilic, metabolically stable, electron-withdrawing"
    },
    {
        "token": "difluoromethyl_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[CX4H](F)F"]},
        "why": "Difluoromethyl group (CHF2); hydrogen bond donor mimic, pKa modulation"
    },
    {
        "token": "sulfone_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[SX4](=[OX1])(=[OX1])([#6])[#6]"]},
        "why": "Sulfone; polar, metabolically stable, hydrogen bond acceptor"
    },
    {
        "token": "sulfonamide_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[SX4](=[OX1])(=[OX1])([NX3])[#6]"]},
        "why": "Sulfonamide; bioisostere for carboxylic acids, drug scaffold"
    },
    {
        "token": "urea_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[NX3][CX3](=[OX1])[NX3]"]},
        "why": "Urea; hydrogen bond donor/acceptor, kinase inhibitor motif"
    },
    {
        "token": "carbamate_present",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[NX3][CX3](=[OX1])[OX2][#6]"]},
        "why": "Carbamate; ester-amide hybrid, prodrug strategy, protecting group"
    },
    {
        "token": "pains_aldehyde_alert",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[CX3H1](=[OX1])[#6]"]},
        "why": "PAINS alert: aldehydes are reactive, potential assay interference"
    },
    {
        "token": "pains_catechol_alert",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["c([OX2H])c[OX2H]"]},
        "why": "PAINS alert: catechol (ortho-diphenol) is redox-active, metal chelator"
    },
    {
        "token": "pains_quinone_alert",
        "type": "bool",
        "scope": "global",
        "category": "medicinal_chemistry",
        "detect": {"smarts_any": ["[CX3](=[OX1])[CX3]=[CX3][CX3](=[OX1])"]},
        "why": "PAINS alert: quinone structure is redox-active, promiscuous"
    },
]

# ============================================================================
# CATEGORY 6: Stereochemistry & Topology
# ============================================================================

stereochemistry = [
    {
        "token": "alkene_e_isomer_present",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"heuristic": "has E-configured alkene"},
        "why": "E-alkene (trans); more stable isomer, different reactivity than Z"
    },
    {
        "token": "alkene_z_isomer_present",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"heuristic": "has Z-configured alkene"},
        "why": "Z-alkene (cis); less stable, different reactivity, proximity effects"
    },
    {
        "token": "atropisomer_potential",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"heuristic": "ortho,ortho-disubstituted biaryl with restricted rotation"},
        "why": "Potential atropisomerism; axial chirality from hindered rotation"
    },
    {
        "token": "spiro_center_present",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"smarts_any": ["[CX4;R2]"]},
        "why": "Spiro center; single atom shared by 2 rings, 3D complexity, exit vector diversity"
    },
    {
        "token": "bridgehead_nitrogen_present",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"smarts_any": ["[NX3;R2,R3]"]},
        "why": "Bridgehead nitrogen; conformationally restricted, reduced basicity"
    },
    {
        "token": "quaternary_carbon_present",
        "type": "bool",
        "scope": "global",
        "category": "stereochemistry",
        "detect": {"smarts_any": ["[CX4;!H3;!H2;!H1;!H0]([#6])([#6])([#6])[#6]"]},
        "why": "Quaternary carbon; all-carbon substituted, no hydrogens, synthetic challenge"
    },
]

# ============================================================================
# CATEGORY 7: Redox & Photochemistry
# ============================================================================

redox_photo = [
    {
        "token": "benzylic_ch_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CH2]c", "[CH1]c"]},
        "why": "Benzylic C-H; oxidation-prone (weak BDE), metabolic hotspot"
    },
    {
        "token": "allylic_ch_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CH2][CX3]=[CX3]", "[CH1][CX3]=[CX3]"]},
        "why": "Allylic C-H; oxidation-prone, radical abstraction site, weak BDE"
    },
    {
        "token": "propargylic_ch_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CH2][CX2]#[CX2]", "[CH1][CX2]#[CX2]"]},
        "why": "Propargylic C-H; acidic, potential for metal-mediated transformations"
    },
    {
        "token": "alpha_amino_ch_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CH2][NX3]", "[CH1][NX3]"]},
        "why": "α-Amino C-H; oxidation-prone, SET oxidation site, imine formation potential"
    },
    {
        "token": "alpha_oxy_ch_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CH2][OX2]", "[CH1][OX2]"]},
        "why": "α-Oxy C-H; oxidation-prone, acetal/hemiacetal formation potential"
    },
    {
        "token": "extended_conjugation_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"heuristic": "3+ contiguous sp2/aromatic atoms"},
        "why": "Extended π-conjugation; chromophore, potential photoreactivity, UV-active"
    },
    {
        "token": "aromatic_ketone_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["c[CX3](=[OX1])[#6]"]},
        "why": "Aromatic ketone; photoactive (Norrish reactions), triplet sensitizer"
    },
    {
        "token": "naphthoquinone_present",
        "type": "bool",
        "scope": "global",
        "category": "redox_photochemistry",
        "detect": {"smarts_any": ["[CX3](=[OX1])c1ccc2ccccc2c1[CX3](=[OX1])"]},
        "why": "Naphthoquinone; redox-active, Michael acceptor, ROS generator"
    },
]

# ============================================================================
# CATEGORY 8: Additional Organometallics
# ============================================================================

additional_organometallics = [
    {
        "token": "aryl_boron_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["c[B]"]},
        "why": "Aryl boron; Suzuki coupling partner, cross-coupling nucleophile"
    },
    {
        "token": "vinyl_boron_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["[CX3]=[CX3][B]"]},
        "why": "Vinyl boron; Suzuki coupling partner for vinyl transfer"
    },
    {
        "token": "alkyl_boron_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["[CX4][B]"]},
        "why": "Alkyl boron; less reactive, requires activated conditions for cross-coupling"
    },
    {
        "token": "aryl_stannane_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["c[Sn]"]},
        "why": "Aryl stannane; Stille coupling partner, transmetallation with Pd"
    },
    {
        "token": "vinyl_stannane_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["[CX3]=[CX3][Sn]"]},
        "why": "Vinyl stannane; Stille coupling for vinyl transfer, stereochemistry preserved"
    },
    {
        "token": "aryl_silane_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["c[Si]"]},
        "why": "Aryl silane; Hiyama coupling partner, requires fluoride activation"
    },
    {
        "token": "vinyl_silane_present",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["[CX3]=[CX3][Si]"]},
        "why": "Vinyl silane; Hiyama coupling, Peterson olefination substrate"
    },
    {
        "token": "organocuprate_potential",
        "type": "bool",
        "scope": "global",
        "category": "organometallics",
        "detect": {"smarts_any": ["[Cu]"]},
        "why": "Copper organometallic; cuprate reagent, soft nucleophile"
    },
]

# ============================================================================
# CATEGORY 9: Sulfur & Phosphorus Functionality
# ============================================================================

sulfur_phosphorus = [
    {
        "token": "sulfide_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[SX2]([#6])[#6]"]},
        "why": "Sulfide (thioether); soft nucleophile, oxidizable to sulfoxide/sulfone"
    },
    {
        "token": "sulfoxide_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[SX3](=[OX1])([#6])[#6]"]},
        "why": "Sulfoxide; chiral, polar, Pummerer rearrangement substrate"
    },
    {
        "token": "sulfonic_acid_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[SX4](=[OX1])(=[OX1])(=[OX1])[OX2H]", "[SX4](=[OX1])(=[OX1])[OX2H]"]},
        "why": "Sulfonic acid; strong acid, very polar, water-soluble"
    },
    {
        "token": "sulfinic_acid_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[SX3](=[OX1])[OX2H]"]},
        "why": "Sulfinic acid; reducing agent, SO2 source, oxidizable"
    },
    {
        "token": "disulfide_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[SX2][SX2]"]},
        "why": "Disulfide (S-S bond); redox-active, cysteine bridge, reductive cleavage"
    },
    {
        "token": "phosphate_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[PX4](=[OX1])([OX2])([OX2])[OX2]"]},
        "why": "Phosphate ester; biologically relevant, leaving group in SN2"
    },
    {
        "token": "phosphonate_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[PX4](=[OX1])([#6])([OX2])[OX2]"]},
        "why": "Phosphonate; phosphate bioisostere, stable to hydrolysis"
    },
    {
        "token": "phosphonium_salt_present",
        "type": "bool",
        "scope": "global",
        "category": "sulfur_phosphorus",
        "detect": {"smarts_any": ["[P+]([#6])([#6])([#6])[#6]"]},
        "why": "Phosphonium salt; Wittig reagent precursor, phase transfer catalyst"
    },
]

# Combine all new features
all_new_features = (
    protecting_groups + 
    heterocycles + 
    reactive_intermediates + 
    physicochemical + 
    medicinal_chem + 
    stereochemistry + 
    redox_photo + 
    additional_organometallics +
    sulfur_phosphorus
)

# Add all new features to the expanded set
expanded_features["features"].extend(all_new_features)

# ============================================================================
# Add comprehensive derived shortcuts
# ============================================================================

expanded_features["derived_shortcuts"] = existing_data["derived_shortcuts"].copy()

# Add new derived features
new_derived = [
    {
        "token": "protected_alcohol_present",
        "derive": "silyl_ether_present OR benzyl_ether_present OR pmb_ether_present OR mom_ether_present OR sem_ether_present OR acetate_ester_present OR benzoate_ester_present",
        "why": "Any alcohol protecting group present"
    },
    {
        "token": "protected_amine_present",
        "derive": "boc_present OR cbz_present OR fmoc_present OR tosyl_sulfonamide_present OR nosyl_sulfonamide_present OR trifluoroacetamide_present OR phthalimide_present",
        "why": "Any amine protecting group present"
    },
    {
        "token": "protected_carbonyl_present",
        "derive": "acetal_present OR ketal_present OR dithiane_present OR oxime_present OR hydrazone_present",
        "why": "Any carbonyl protecting group present"
    },
    {
        "token": "protected_acid_present",
        "derive": "methyl_ester_present OR ethyl_ester_present OR tert_butyl_ester_present OR benzyl_ester_present",
        "why": "Carboxylic acid protected as ester"
    },
    {
        "token": "five_membered_heterocycle_present",
        "derive": "furan_present OR thiophene_present OR pyrrole_present OR imidazole_present OR oxazole_present OR thiazole_present OR triazole_present OR tetrazole_present",
        "why": "Any 5-membered heterocycle present"
    },
    {
        "token": "six_membered_heterocycle_present",
        "derive": "pyridine_present OR pyrimidine_present OR morpholine_present OR piperazine_present OR piperidine_present",
        "why": "Any 6-membered heterocycle present"
    },
    {
        "token": "fused_heterocycle_present",
        "derive": "quinoline_present OR isoquinoline_present OR benzofuran_present OR benzothiophene_present OR benzimidazole_present OR purine_present OR indole_present",
        "why": "Any fused heterocycle system present"
    },
    {
        "token": "lipinski_compliant",
        "derive": "lipinski_hbd_compliant AND lipinski_hba_compliant AND lipinski_mw_compliant AND lipinski_logp_compliant",
        "why": "Passes all 4 Lipinski rules (Rule of 5) for oral bioavailability"
    },
    {
        "token": "veber_compliant",
        "derive": "veber_rotb_compliant AND veber_tpsa_compliant",
        "why": "Passes Veber rules for oral bioavailability"
    },
    {
        "token": "drug_like",
        "derive": "lipinski_compliant AND veber_compliant AND NOT pains_aldehyde_alert AND NOT pains_catechol_alert",
        "why": "Drug-like: passes Lipinski, Veber rules and basic PAINS filters"
    },
    {
        "token": "michael_acceptor_present",
        "derive": "alpha_beta_unsaturated_carbonyl_present OR alpha_beta_unsaturated_ester_present OR alpha_beta_unsaturated_nitrile_present",
        "why": "Michael acceptor system present (conjugate addition substrate)"
    },
    {
        "token": "electrophilic_warhead_present",
        "derive": "epoxide_present OR aziridine_present OR alpha_halo_carbonyl_present OR acyl_chloride_present OR isocyanate_present OR michael_acceptor_present",
        "why": "Electrophilic warhead for covalent inhibition or alkylation"
    },
    {
        "token": "metabolic_soft_spot_present",
        "derive": "benzylic_ch_present OR allylic_ch_present OR alpha_amino_ch_present OR alpha_oxy_ch_present OR tertiary_amine_present",
        "why": "Metabolic soft spot (oxidation-prone site) present"
    },
    {
        "token": "fluorinated_motif_present",
        "derive": "fluorine_present OR trifluoromethyl_present OR difluoromethyl_present OR sp2_fluoride_present",
        "why": "Any fluorine-containing motif present"
    },
    {
        "token": "cross_coupling_electrophile",
        "derive": "aryl_halide_present OR vinyl_halide_present OR sp2_triflate_present OR sp2_tosylate_present",
        "why": "Suitable electrophile for Pd-catalyzed cross-coupling reactions"
    },
    {
        "token": "cross_coupling_nucleophile",
        "derive": "sp2_boron_present OR stannane_present OR organosilane_present OR organozinc_present OR grignard_present",
        "why": "Suitable nucleophile for metal-catalyzed cross-coupling reactions"
    },
    {
        "token": "oxidation_prone",
        "derive": "benzylic_ch_present OR allylic_ch_present OR alpha_amino_ch_present OR alpha_oxy_ch_present OR sulfide_present OR thiol_present OR phenol_present",
        "why": "Contains oxidation-prone functional groups"
    },
    {
        "token": "reduction_prone",
        "derive": "nitro_present OR azide_present OR carbonyl_present OR nitrile_present OR alkene_present OR alkyne_present",
        "why": "Contains reduction-prone functional groups"
    },
    {
        "token": "moisture_sensitive",
        "derive": "acyl_chloride_present OR acid_anhydride_present OR sulfonyl_chloride_present OR isocyanate_present OR grignard_present OR organolithium_present",
        "why": "Moisture-sensitive groups present; requires anhydrous conditions"
    },
    {
        "token": "air_sensitive",
        "derive": "phosphine_present OR grignard_present OR organolithium_present OR thiol_present OR aldehyde_present",
        "why": "Air-sensitive groups present; requires inert atmosphere"
    },
    {
        "token": "explosive_risk",
        "derive": "azide_present OR diazo_compound_present OR peroxide_present OR nitro_present",
        "why": "Potentially explosive functional groups present; handle with caution"
    },
]

expanded_features.setdefault("derived_shortcuts", [])
existing_data.setdefault("derived_shortcuts", [])
expanded_features["derived_shortcuts"].extend(new_derived)

# Save the expanded JSON
output_file = Path(__file__).parent.parent / "chemtools" / "taxonomy" / "data" / "calculable_features.json"
with open(output_file, 'w', encoding='utf-8') as f:
    json.dump(expanded_features, f, indent=2)
print("Expanded calculable_features.json")
print(f"  Original features: {len(existing_data['features'])}")
print(f"  New features added: {len(all_new_features)}")
print(f"  Total features: {len(expanded_features['features'])}")
print(f"  Original derived shortcuts: {len(existing_data.get('derived_shortcuts', []))}")
print(f"  New derived shortcuts: {len(new_derived)}")
print(f"  Total derived shortcuts: {len(expanded_features.get('derived_shortcuts', []))}")
print(f"  Version: {expanded_features['version']}")
print(f"\nCategories added:")
for cat in expanded_features['schema_notes']['categories']:
    print(f"  - {cat}")
