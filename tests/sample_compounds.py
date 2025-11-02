"""
Sample Compounds Data
====================

This file contains comprehensive sample compound SMILES with known reaction types
for testing the calculable features detection system and reaction classification.

Each compound is annotated with:
- SMILES string
- Common name or description
- Typical reaction role (electrophile, nucleophile, coupling partner)
- Expected reaction types
- Key structural features
"""

from typing import Dict, List, Any


# ═══════════════════════════════════════════════════════════
# ELECTROPHILES - Aryl Halides
# ═══════════════════════════════════════════════════════════

ARYL_HALIDES = [
    {
        "smiles": "Brc1ccccc1",
        "name": "Bromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Ullmann", "Sonogashira", "Heck", "Stille", "Negishi", "Kumada"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Standard aryl bromide for all cross-coupling reactions"
    },
    {
        "smiles": "Clc1ccccc1",
        "name": "Chlorobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Ullmann", "Sonogashira", "Heck"],
        "features": ["sp2_chloride_present", "aryl_halide_present", "ArCl_present"],
        "notes": "Less reactive; requires activated catalyst systems"
    },
    {
        "smiles": "Ic1ccccc1",
        "name": "Iodobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira", "Heck", "Stille", "Negishi", "Kumada"],
        "features": ["sp2_iodide_present", "aryl_halide_present", "ArI_present"],
        "notes": "Most reactive aryl halide"
    },
    {
        "smiles": "Fc1ccccc1",
        "name": "Fluorobenzene",
        "role": "electrophile",
        "reaction_types": ["SNAr"],
        "features": ["sp2_fluoride_present", "aryl_halide_present", "ArF_present"],
        "notes": "Requires strong nucleophiles for SNAr; not for standard cross-coupling"
    },
    {
        "smiles": "Brc1ccc(C#N)cc1",
        "name": "4-Bromobenzonitrile",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Electron-poor aryl bromide; more reactive"
    },
    {
        "smiles": "Brc1ccc(OC)cc1",
        "name": "4-Bromoanisole",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira", "Heck"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Electron-rich aryl bromide; standard substrate"
    },
    {
        "smiles": "Brc1ccc([N+](=O)[O-])cc1",
        "name": "4-Bromonitrobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "SNAr"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Highly electron-poor; very reactive"
    },
    {
        "smiles": "Brc1ccc(C(F)(F)F)cc1",
        "name": "4-Bromobenzotrifluoride",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Electron-poor CF3 substituent"
    },
    {
        "smiles": "Brc1cc(Br)cc(Br)c1",
        "name": "1,3,5-Tribromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Multiple coupling sites; sp2_halide_site_count = 3"
    },
    {
        "smiles": "Clc1ccc(Cl)cc1",
        "name": "1,4-Dichlorobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_chloride_present", "aryl_halide_present", "ArCl_present"],
        "notes": "Di-functionalization substrate"
    },
    {
        "smiles": "Brc1ccccc1C",
        "name": "2-Bromotoluene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Ortho-substituted; sterically hindered"
    },
    {
        "smiles": "Brc1cc(C)cc(C)c1",
        "name": "3,5-Dimethylbromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Meta-substituted; moderate sterics"
    },
    {
        "smiles": "Brc1c(C)c(C)c(C)c(C)c1C",
        "name": "Pentamethylbromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Highly sterically hindered"
    },
    {
        "smiles": "Brc1ccc2ccccc2c1",
        "name": "2-Bromonaphthalene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Fused aromatic system"
    },
    {
        "smiles": "Clc1ccc2ccccc2c1",
        "name": "2-Chloronaphthalene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_chloride_present", "aryl_halide_present", "ArCl_present"],
        "notes": "Activated chloride on naphthalene"
    },
    {
        "smiles": "Brc1cccc2c1cccc2",
        "name": "1-Bromonaphthalene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Peri-position naphthyl bromide"
    },


    # ═══════════════════════════════════════════════════════════
    # ELECTROPHILES - Heteroaryl Halides
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "Brc1ccncc1",
        "name": "4-Bromopyridine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira", "Stille", "Negishi"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present", "pyridine_present", "pyridine_poison_risk"],
        "notes": "Heteroaryl halide; pyridine can poison some catalysts"
    },
    {
        "smiles": "Brc1cccnc1",
        "name": "3-Bromopyridine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present", "pyridine_present", "pyridine_poison_risk"],
        "notes": "Meta-bromo pyridine"
    },
    {
        "smiles": "Brc1ccccn1",
        "name": "2-Bromopyridine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present", "pyridine_present", "pyridine_poison_risk"],
        "notes": "Ortho-bromo pyridine; strong chelation"
    },
    {
        "smiles": "Clc1ccncc1",
        "name": "4-Chloropyridine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_chloride_present", "heteroaryl_halide_present", "pyridine_present"],
        "notes": "Activated heterocyclic chloride"
    },
    {
        "smiles": "Brc1cnccn1",
        "name": "5-Bromopyrimidine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present", "pyrimidine_present"],
        "notes": "Diazine substrate"
    },
    {
        "smiles": "Brc1cccs1",
        "name": "2-Bromothiophene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira", "Stille"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "Five-membered heterocycle"
    },
    {
        "smiles": "Brc1ccco1",
        "name": "2-Bromofuran",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "Oxygen heterocycle"
    },
    {
        "smiles": "Brc1cc[nH]c1",
        "name": "3-Bromopyrrole",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "NH heterocycle; acidic proton"
    },
    {
        "smiles": "Brc1ccc2[nH]ccc2c1",
        "name": "5-Bromoindole",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present", "indole_present"],
        "notes": "Indole scaffold; NH can coordinate"
    },
    {
        "smiles": "Brc1nc2ccccc2s1",
        "name": "2-Bromobenzothiazole",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "Fused bicyclic heterocycle"
    },
    {
        "smiles": "Ic1nc2ccccc2o1",
        "name": "2-Iodobenzoxazole",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_iodide_present", "heteroaryl_halide_present"],
        "notes": "Oxygen-containing fused heterocycle"
    },
    {
        "smiles": "Brc1ccnc2ccccc12",
        "name": "3-Bromoquinoline",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "Quinoline scaffold"
    },
    {
        "smiles": "Brc1cnc2ccccc2n1",
        "name": "3-Bromoquinoxaline",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "heteroaryl_halide_present"],
        "notes": "Diaza-naphthalene"
    },
    {
        "smiles": "Clc1ccc(Cl)nc1",
        "name": "2,5-Dichloropyridine",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_chloride_present", "heteroaryl_halide_present", "pyridine_present"],
        "notes": "Selective coupling possible at different positions"
    },


    # ═══════════════════════════════════════════════════════════
    # ELECTROPHILES - Aryl Sulfonates (Pseudohalides)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "c1ccc(OS(=O)(=O)C(F)(F)F)cc1",
        "name": "Phenyl triflate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira", "Heck"],
        "features": ["sp2_triflate_present", "aryl_halide_present", "ArOTf_present", "sp2_pseudohalide_present"],
        "notes": "Excellent leaving group; very reactive"
    },
    {
        "smiles": "COc1ccc(OS(=O)(=O)C(F)(F)F)cc1",
        "name": "4-Methoxyphenyl triflate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_triflate_present", "aryl_halide_present", "ArOTf_present"],
        "notes": "Electron-rich triflate"
    },
    {
        "smiles": "c1ccc(OS(=O)(=O)c2ccc(C)cc2)cc1",
        "name": "Phenyl tosylate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_tosylate_present", "aryl_halide_present", "ArOTs_present", "sp2_pseudohalide_present"],
        "notes": "Good leaving group; cheaper than triflate"
    },
    {
        "smiles": "c1ccc(OS(=O)(=O)C)cc1",
        "name": "Phenyl mesylate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_mesylate_present", "aryl_halide_present", "ArOMs_present", "sp2_pseudohalide_present"],
        "notes": "Mesylate leaving group"
    },
    {
        "smiles": "N#Cc1ccc(OS(=O)(=O)C(F)(F)F)cc1",
        "name": "4-Cyanophenyl triflate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_triflate_present", "aryl_halide_present", "ArOTf_present"],
        "notes": "Electron-poor triflate; highly reactive"
    },


    # ═══════════════════════════════════════════════════════════
    # ELECTROPHILES - Vinyl Halides
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "C=CBr",
        "name": "Vinyl bromide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira", "Heck", "Stille", "Negishi", "Kumada"],
        "features": ["sp2_bromide_present", "vinyl_halide_present", "VinylBr_present", "alkene_present", "terminal_alkene_present"],
        "notes": "Simple vinyl halide"
    },
    {
        "smiles": "C=CCl",
        "name": "Vinyl chloride",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Heck"],
        "features": ["sp2_chloride_present", "vinyl_halide_present", "VinylCl_present", "alkene_present"],
        "notes": "Less reactive vinyl halide"
    },
    {
        "smiles": "C=CI",
        "name": "Vinyl iodide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira", "Heck", "Stille"],
        "features": ["sp2_iodide_present", "vinyl_halide_present", "VinylI_present", "alkene_present"],
        "notes": "Most reactive vinyl halide"
    },
    {
        "smiles": "C/C=C/Br",
        "name": "(E)-1-Bromopropene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira", "Heck"],
        "features": ["sp2_bromide_present", "vinyl_halide_present", "VinylBr_present", "alkene_present"],
        "notes": "Stereodefined vinyl bromide"
    },
    {
        "smiles": "C=C(Br)C",
        "name": "2-Bromopropene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "vinyl_halide_present", "alkene_present"],
        "notes": "Geminal disubstituted vinyl bromide"
    },
    {
        "smiles": "BrC(=C)c1ccccc1",
        "name": "α-Bromostyrene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "vinyl_halide_present", "alkene_present"],
        "notes": "Aryl-substituted vinyl bromide"
    },


    # ═══════════════════════════════════════════════════════════
    # ELECTROPHILES - Alkyl Halides (sp3)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "CCBr",
        "name": "Ethyl bromide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Negishi", "Kumada", "SN2"],
        "features": ["sp3_bromide_present", "beta_hydride_possible"],
        "notes": "Simple primary alkyl bromide; β-hydride elimination risk"
    },
    {
        "smiles": "CCCl",
        "name": "Ethyl chloride",
        "role": "electrophile",
        "reaction_types": ["Negishi", "Kumada", "SN2"],
        "features": ["sp3_chloride_present", "beta_hydride_possible"],
        "notes": "Primary alkyl chloride"
    },
    {
        "smiles": "CCI",
        "name": "Ethyl iodide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Negishi", "Kumada", "SN2"],
        "features": ["sp3_iodide_present", "beta_hydride_possible"],
        "notes": "Most reactive primary alkyl halide"
    },
    {
        "smiles": "CC(C)Br",
        "name": "Isopropyl bromide",
        "role": "electrophile",
        "reaction_types": ["Negishi", "Kumada"],
        "features": ["sp3_bromide_present", "beta_hydride_possible"],
        "notes": "Secondary alkyl bromide; high β-hydride risk"
    },
    {
        "smiles": "CC(C)(C)Br",
        "name": "tert-Butyl bromide",
        "role": "electrophile",
        "reaction_types": ["SN1"],
        "features": ["sp3_bromide_present"],
        "notes": "Tertiary alkyl bromide; no β-hydride; SN1 mechanism"
    },
    {
        "smiles": "BrCc1ccccc1",
        "name": "Benzyl bromide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "SN2", "Negishi"],
        "features": ["sp3_bromide_present", "benzylic_halide_present"],
        "notes": "Benzylic position; activated sp3 halide"
    },
    {
        "smiles": "C(Br)Cc1ccccc1",
        "name": "2-Phenylethyl bromide",
        "role": "electrophile",
        "reaction_types": ["SN2", "Negishi", "Kumada"],
        "features": ["sp3_bromide_present", "beta_hydride_possible"],
        "notes": "Primary alkyl bromide with aryl substituent"
    },
    {
        "smiles": "C=CCCBr",
        "name": "4-Bromo-1-butene",
        "role": "electrophile",
        "reaction_types": ["SN2", "Negishi"],
        "features": ["sp3_bromide_present", "alkene_present", "beta_hydride_possible"],
        "notes": "Homoallylic bromide"
    },
    {
        "smiles": "C=C(C)CBr",
        "name": "3-Bromo-2-methylpropene (allylic)",
        "role": "electrophile",
        "reaction_types": ["SN2", "Allylic-substitution"],
        "features": ["sp3_bromide_present", "allylic_halide_present", "alkene_present"],
        "notes": "Allylic bromide; can undergo allylic substitution"
    },
    {
        "smiles": "CC(C)(C)CBr",
        "name": "Neopentyl bromide",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Negishi"],
        "features": ["sp3_bromide_present"],
        "notes": "No β-hydrogens; minimal elimination risk"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Boronic Acids and Esters (Suzuki)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "c1ccc(B(O)O)cc1",
        "name": "Phenylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present"],
        "notes": "Standard boronic acid coupling partner"
    },
    {
        "smiles": "COc1ccc(B(O)O)cc1",
        "name": "4-Methoxyphenylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present"],
        "notes": "Electron-rich boronic acid"
    },
    {
        "smiles": "N#Cc1ccc(B(O)O)cc1",
        "name": "4-Cyanophenylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present"],
        "notes": "Electron-poor boronic acid"
    },
    {
        "smiles": "c1ccc(B2OC(C)(C)C(C)(C)O2)cc1",
        "name": "Phenylboronic acid pinacol ester",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present", "boron_bpin_present"],
        "notes": "Protected boronic acid; more stable"
    },
    {
        "smiles": "c1ccc([B-](F)(F)F)cc1.[K+]",
        "name": "Potassium phenyltrifluoroborate",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present", "boron_bf3k_present"],
        "notes": "BF3K salt; air-stable boronate"
    },
    {
        "smiles": "C=CB(O)O",
        "name": "Vinylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present", "alkene_present", "terminal_alkene_present"],
        "notes": "Vinyl boronic acid"
    },
    {
        "smiles": "C/C=C/B(O)O",
        "name": "(E)-Propenylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present", "alkene_present"],
        "notes": "Stereodefined vinyl boronic acid"
    },
    {
        "smiles": "CCB(O)O",
        "name": "Ethylboronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp3_boron_present"],
        "notes": "Alkyl boronic acid; prone to β-hydride issues"
    },
    {
        "smiles": "c1ccncc1B(O)O",
        "name": "Pyridine-4-boronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present", "pyridine_present", "pyridine_poison_risk"],
        "notes": "Heteroaryl boronic acid; pyridine coordination risk"
    },
    {
        "smiles": "c1coc(B(O)O)c1",
        "name": "Furan-2-boronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present"],
        "notes": "Five-membered heterocyclic boronic acid"
    },
    {
        "smiles": "c1csc(B(O)O)c1",
        "name": "Thiophene-2-boronic acid",
        "role": "nucleophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_boron_present"],
        "notes": "Thiophene boronic acid; sulfur can coordinate"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Amines (Buchwald-Hartwig, Ullmann)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "c1ccc(N)cc1",
        "name": "Aniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann", "SNAr"],
        "features": ["aniline_present", "aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Aromatic primary amine; standard C-N coupling partner"
    },
    {
        "smiles": "CCN",
        "name": "Ethylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann", "SNAr"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Primary aliphatic amine"
    },
    {
        "smiles": "CC(C)N",
        "name": "Isopropylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Branched primary amine; α-branching"
    },
    {
        "smiles": "CC(C)(C)N",
        "name": "tert-Butylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Highly hindered primary amine"
    },
    {
        "smiles": "CCNC",
        "name": "N-Methylethylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Secondary aliphatic amine"
    },
    {
        "smiles": "c1ccc(NC)cc1",
        "name": "N-Methylaniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Secondary aromatic amine"
    },
    {
        "smiles": "c1ccc(N(C)C)cc1",
        "name": "N,N-Dimethylaniline",
        "role": "nucleophile",
        "reaction_types": [],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Tertiary amine; not reactive in C-N coupling"
    },
    {
        "smiles": "COc1ccc(N)cc1",
        "name": "4-Methoxyaniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Electron-rich aniline"
    },
    {
        "smiles": "N#Cc1ccc(N)cc1",
        "name": "4-Aminobenzonitrile",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Electron-poor aniline; less nucleophilic"
    },
    {
        "smiles": "Nc1ccncc1",
        "name": "4-Aminopyridine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig"],
        "features": ["aliphatic_amine_present", "pyridine_present", "pyridine_poison_risk"],
        "notes": "Heteroaryl amine; pyridine coordination"
    },
    {
        "smiles": "c1ccc2[nH]ccc2c1",
        "name": "Indole",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "C-H activation"],
        "features": ["indole_present", "acidic_proton_present"],
        "notes": "NH heterocycle; can act as nucleophile"
    },
    {
        "smiles": "c1c[nH]c(N)c1",
        "name": "3-Aminopyrrole",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Heteroaryl amine with NH"
    },
    {
        "smiles": "Nc1ccccc1C",
        "name": "2-Methylaniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Ortho-substituted aniline; sterically hindered"
    },
    {
        "smiles": "Nc1cc(C)cc(C)c1",
        "name": "3,5-Dimethylaniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Meta-substituted aniline"
    },
    {
        "smiles": "NCc1ccccc1",
        "name": "Benzylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aliphatic_amine_present", "acidic_proton_present"],
        "notes": "Benzylic amine; aliphatic but stabilized"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Alcohols and Phenols (C-O Coupling)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "c1ccc(O)cc1",
        "name": "Phenol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "SNAr", "Mitsunobu"],
        "features": ["phenol_present", "acidic_proton_present"],
        "notes": "Aromatic alcohol; C-O coupling partner"
    },
    {
        "smiles": "CCO",
        "name": "Ethanol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "Mitsunobu"],
        "features": ["alcohol_present", "acidic_proton_present"],
        "notes": "Simple aliphatic alcohol"
    },
    {
        "smiles": "CC(C)O",
        "name": "Isopropanol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "Mitsunobu"],
        "features": ["alcohol_present", "acidic_proton_present"],
        "notes": "Secondary alcohol"
    },
    {
        "smiles": "CC(C)(C)O",
        "name": "tert-Butanol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann"],
        "features": ["alcohol_present", "acidic_proton_present"],
        "notes": "Tertiary alcohol; sterically hindered"
    },
    {
        "smiles": "COc1ccc(O)cc1",
        "name": "4-Methoxyphenol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "SNAr"],
        "features": ["phenol_present", "acidic_proton_present"],
        "notes": "Electron-rich phenol"
    },
    {
        "smiles": "N#Cc1ccc(O)cc1",
        "name": "4-Hydroxybenzonitrile",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "SNAr"],
        "features": ["phenol_present", "acidic_proton_present"],
        "notes": "Electron-poor phenol; more acidic"
    },
    {
        "smiles": "OCc1ccccc1",
        "name": "Benzyl alcohol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "Mitsunobu"],
        "features": ["alcohol_present", "acidic_proton_present"],
        "notes": "Benzylic alcohol"
    },
    {
        "smiles": "Oc1ccncc1",
        "name": "4-Hydroxypyridine",
        "role": "nucleophile",
        "reaction_types": ["Ullmann"],
        "features": ["phenol_present", "pyridine_present", "acidic_proton_present"],
        "notes": "Heteroaryl phenol"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Thiols (C-S Coupling)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "c1ccc(S)cc1",
        "name": "Thiophenol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "SNAr"],
        "features": ["thiol_present", "thiol_poison_risk", "acidic_proton_present"],
        "notes": "Aromatic thiol; strong nucleophile; catalyst poison risk"
    },
    {
        "smiles": "CCS",
        "name": "Ethanethiol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann", "SNAr"],
        "features": ["thiol_present", "thiol_poison_risk", "acidic_proton_present"],
        "notes": "Aliphatic thiol; catalyst poison"
    },
    {
        "smiles": "CC(C)S",
        "name": "2-Propanethiol",
        "role": "nucleophile",
        "reaction_types": ["Ullmann"],
        "features": ["thiol_present", "thiol_poison_risk", "acidic_proton_present"],
        "notes": "Secondary thiol"
    },
    {
        "smiles": "CSc1ccccc1",
        "name": "Thioanisole",
        "role": "none",
        "reaction_types": [],
        "features": ["thiol_poison_risk"],
        "notes": "Sulfide; not nucleophilic but can poison catalysts"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Terminal Alkynes (Sonogashira)
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "C#C",
        "name": "Acetylene",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present"],
        "notes": "Simplest terminal alkyne"
    },
    {
        "smiles": "C#Cc1ccccc1",
        "name": "Phenylacetylene",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present"],
        "notes": "Aryl-substituted terminal alkyne"
    },
    {
        "smiles": "C#CC(C)(C)C",
        "name": "3,3-Dimethylbut-1-yne",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present"],
        "notes": "Sterically hindered terminal alkyne"
    },
    {
        "smiles": "C#C[Si](C)(C)C",
        "name": "Trimethylsilylacetylene",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present", "organosilane_present", "base_sensitive"],
        "notes": "TMS-protected alkyne; deprotectable"
    },
    {
        "smiles": "C#CCOC",
        "name": "Methyl propargyl ether",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present"],
        "notes": "Ether-substituted alkyne"
    },
    {
        "smiles": "C#Cc1ccncc1",
        "name": "4-Ethynylpyridine",
        "role": "nucleophile",
        "reaction_types": ["Sonogashira"],
        "features": ["alkyne_present", "terminal_alkyne_present", "pyridine_present"],
        "notes": "Heteroaryl terminal alkyne"
    },


    # ═══════════════════════════════════════════════════════════
    # NUCLEOPHILES - Organometallic Reagents
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "[Mg]Brc1ccccc1",
        "name": "Phenylmagnesium bromide",
        "role": "nucleophile",
        "reaction_types": ["Kumada"],
        "features": ["grignard_present"],
        "notes": "Grignard reagent for Kumada coupling"
    },
    {
        "smiles": "CC[Mg]Br",
        "name": "Ethylmagnesium bromide",
        "role": "nucleophile",
        "reaction_types": ["Kumada"],
        "features": ["grignard_present"],
        "notes": "Alkyl Grignard"
    },
    {
        "smiles": "[Zn]Clc1ccccc1",
        "name": "Phenylzinc chloride",
        "role": "nucleophile",
        "reaction_types": ["Negishi"],
        "features": ["organozinc_present"],
        "notes": "Organozinc for Negishi coupling"
    },
    {
        "smiles": "CC[Zn]Br",
        "name": "Ethylzinc bromide",
        "role": "nucleophile",
        "reaction_types": ["Negishi"],
        "features": ["organozinc_present"],
        "notes": "Alkylzinc reagent"
    },
    {
        "smiles": "c1ccc([Sn](C)(C)C)cc1",
        "name": "Phenyltrimethylstannane",
        "role": "nucleophile",
        "reaction_types": ["Stille"],
        "features": ["stannane_present"],
        "notes": "Organotin for Stille coupling"
    },
    {
        "smiles": "C=C[Sn](C)(C)C",
        "name": "Vinyltrimethylstannane",
        "role": "nucleophile",
        "reaction_types": ["Stille"],
        "features": ["stannane_present", "alkene_present"],
        "notes": "Vinyl stannane"
    },


    # ═══════════════════════════════════════════════════════════
    # SPECIAL CASES - Multifunctional Compounds
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "Brc1ccc(B(O)O)cc1",
        "name": "4-Bromophenylboronic acid",
        "role": "bifunctional",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "sp2_boron_present", "aryl_halide_present"],
        "notes": "Can act as both electrophile and nucleophile; selective coupling"
    },
    {
        "smiles": "Brc1ccc(N)cc1",
        "name": "4-Bromoaniline",
        "role": "bifunctional",
        "reaction_types": ["Buchwald-Hartwig", "Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "aniline_present", "aryl_halide_present", "aliphatic_amine_present"],
        "notes": "Both electrophile (Br) and nucleophile (NH2); orthogonal protection needed"
    },
    {
        "smiles": "Brc1ccc(O)cc1",
        "name": "4-Bromophenol",
        "role": "bifunctional",
        "reaction_types": ["Ullmann", "Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "phenol_present", "aryl_halide_present"],
        "notes": "Both electrophile (Br) and nucleophile (OH)"
    },
    {
        "smiles": "Brc1ccc(C#C)cc1",
        "name": "4-Bromophenylacetylene",
        "role": "bifunctional",
        "reaction_types": ["Sonogashira", "Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "alkyne_present", "terminal_alkyne_present", "aryl_halide_present"],
        "notes": "Both electrophile (Br) and nucleophile (alkyne)"
    },
    {
        "smiles": "Brc1ccc(Br)cc1",
        "name": "1,4-Dibromobenzene",
        "role": "di-electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "ArBr_present"],
        "notes": "Two identical electrophilic sites; sp2_halide_site_count = 2"
    },
    {
        "smiles": "Nc1ccc(N)cc1",
        "name": "1,4-Phenylenediamine",
        "role": "di-nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Two identical nucleophilic sites"
    },


    # ═══════════════════════════════════════════════════════════
    # PROTECTED COMPOUNDS
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "Brc1ccc(NC(=O)OC(C)(C)C)cc1",
        "name": "4-Bromo-N-Boc-aniline",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Boc-protected amine; prevents C-N coupling"
    },
    {
        "smiles": "Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1",
        "name": "4-Bromo-TBS-phenol",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "organosilane_present", "base_sensitive"],
        "notes": "TBS-protected phenol"
    },
    {
        "smiles": "Brc1ccc(C(=O)OC)cc1",
        "name": "Methyl 4-bromobenzoate",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Ester-protected carboxylic acid"
    },
    {
        "smiles": "Brc1ccc(C=O)cc1",
        "name": "4-Bromobenzaldehyde",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Aldehyde functional group; can be reduced/oxidized"
    },


    # ═══════════════════════════════════════════════════════════
    # STERICALLY DEMANDING SUBSTRATES
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "Brc1c(C)c(C)c(C)c(C)c1C",
        "name": "Pentamethylbromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Extremely sterically hindered; requires bulky ligands"
    },
    {
        "smiles": "Brc1ccc(C(C)(C)C)cc1",
        "name": "4-Bromo-tert-butylbenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Para-tert-butyl substituent"
    },
    {
        "smiles": "CC(C)(C)N",
        "name": "tert-Butylamine",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig"],
        "features": ["aliphatic_amine_present"],
        "notes": "Highly hindered primary amine; challenging substrate"
    },
    {
        "smiles": "CC(C)(C)c1ccc(N)cc1",
        "name": "4-tert-Butylaniline",
        "role": "nucleophile",
        "reaction_types": ["Buchwald-Hartwig", "Ullmann"],
        "features": ["aniline_present", "aliphatic_amine_present"],
        "notes": "Sterically hindered aniline"
    },


    # ═══════════════════════════════════════════════════════════
    # ELECTRON-RICH/POOR SUBSTRATES
    # ═══════════════════════════════════════════════════════════

    {
        "smiles": "Brc1ccc([N+](=O)[O-])cc1",
        "name": "4-Bromonitrobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Buchwald-Hartwig", "SNAr"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Highly electron-poor; very reactive in SNAr"
    },
    {
        "smiles": "Brc1ccc(N(C)C)cc1",
        "name": "4-Bromo-N,N-dimethylaniline",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "Sonogashira"],
        "features": ["sp2_bromide_present", "aryl_halide_present", "aniline_present"],
        "notes": "Highly electron-rich; tertiary amine substituent"
    },
    {
        "smiles": "COc1cc(OC)c(Br)c(OC)c1OC",
        "name": "1-Bromo-2,3,4,5-tetramethoxybenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura"],
        "features": ["sp2_bromide_present", "aryl_halide_present"],
        "notes": "Extremely electron-rich; multiple methoxy groups"
    },
    {
        "smiles": "Fc1c(F)c(F)c(Br)c(F)c1F",
        "name": "Pentafluorobromobenzene",
        "role": "electrophile",
        "reaction_types": ["Suzuki-Miyaura", "SNAr"],
        "features": ["sp2_bromide_present", "sp2_fluoride_present", "aryl_halide_present"],
        "notes": "Extremely electron-poor; multiple fluorines"
    },
]


# Convenience function to get compounds by role
def get_compounds_by_role(role: str) -> List[Dict[str, Any]]:
    """
    Get all compounds matching a specific role.
    
    Args:
        role: 'electrophile', 'nucleophile', 'bifunctional', or 'di-electrophile'
    
    Returns:
        List of compound dictionaries
    """
    return [c for c in ARYL_HALIDES if c.get("role") == role]


# Convenience function to get compounds by reaction type
def get_compounds_by_reaction(reaction_type: str) -> List[Dict[str, Any]]:
    """
    Get all compounds suitable for a specific reaction type.
    
    Args:
        reaction_type: e.g., 'Suzuki-Miyaura', 'Buchwald-Hartwig', etc.
    
    Returns:
        List of compound dictionaries
    """
    return [c for c in ARYL_HALIDES if reaction_type in c.get("reaction_types", [])]


# Convenience function to get compounds by feature
def get_compounds_by_feature(feature: str) -> List[Dict[str, Any]]:
    """
    Get all compounds with a specific calculable feature.
    
    Args:
        feature: e.g., 'sp2_bromide_present', 'ArBr_present', etc.
    
    Returns:
        List of compound dictionaries
    """
    return [c for c in ARYL_HALIDES if feature in c.get("features", [])]


# Export all compounds as flat list
ALL_SAMPLE_COMPOUNDS = ARYL_HALIDES


if __name__ == "__main__":
    print(f"Total sample compounds: {len(ALL_SAMPLE_COMPOUNDS)}")
    print(f"\nBreakdown by role:")
    print(f"  Electrophiles: {len(get_compounds_by_role('electrophile'))}")
    print(f"  Nucleophiles: {len(get_compounds_by_role('nucleophile'))}")
    print(f"  Bifunctional: {len(get_compounds_by_role('bifunctional'))}")
    print(f"\nSample compounds for Suzuki-Miyaura: {len(get_compounds_by_reaction('Suzuki-Miyaura'))}")
    print(f"Sample compounds for Buchwald-Hartwig: {len(get_compounds_by_reaction('Buchwald-Hartwig'))}")
