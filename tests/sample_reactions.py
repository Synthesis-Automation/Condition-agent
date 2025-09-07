"""
Sample Reactions Data
====================

This file contains comprehensive sample reaction SMILES for testing
and demonstration purposes - over 80 diverse examples covering
various reaction types, with extensive C-N coupling coverage.
"""

# Comprehensive sample reactions for testing - over 80 diverse examples
SAMPLE_REACTIONS = [
    "Select a sample reaction...",
    
    # ═══════════════════════════════════════════════════════════
    # C-C COUPLING REACTIONS (Cross-coupling)
    # ═══════════════════════════════════════════════════════════
    
    # Suzuki-Miyaura Coupling
    "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Simple Ph-Ph)",
    "Clc1ccc(C#N)cc1.c1ccc(B(O)O)cc1>>N#Cc1ccc(-c2ccccc2)cc1 (Suzuki - Electron-poor ArCl)",
    "Brc1ccc(OC)cc1.c1ccc(B(O)O)cc1>>COc1ccc(-c2ccccc2)cc1 (Suzuki - Electron-rich ArBr)",
    "Ic1ccncc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccncc2)cc1 (Suzuki - Heteroaryl pyridine)",
    "Brc1ccc(C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>FC(F)(F)c1ccc(-c2ccccc2)cc1 (Suzuki - CF3 substrate)",
    "Clc1ccc2ccccc2c1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc3ccccc3c2)cc1 (Suzuki - Naphthyl chloride)",
    "Brc1cc(C)ccc1C.c1ccc(B(O)O)cc1>>Cc1ccc(C)c(-c2ccccc2)c1 (Suzuki - Sterically hindered)",
    
    # Stille Coupling
    "Brc1ccccc1.c1ccc([Sn](C)(C)C)cc1>>c1ccc(-c2ccccc2)cc1 (Stille - Ph-Ph)",
    "Ic1ccncc1.C=C[Sn](C)(C)C>>C=Cc1ccncc1 (Stille - Vinyl pyridine)",
    
    # Sonogashira Coupling
    "Brc1ccccc1.C#Cc1ccccc1>>c1ccc(C#Cc2ccccc2)cc1 (Sonogashira - Diphenylacetylene)",
    "Ic1ccncc1.C#CC>>C#Cc1ccncc1 (Sonogashira - Pyridine acetylene)",
    "Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1 (Sonogashira - tert-butyl)",
    
    # Heck Reaction
    "Brc1ccccc1.C=C>>c1ccc(/C=C/)cc1 (Heck - Simple styrene)",
    "Ic1ccc(C=O)cc1.C=CC(=O)OCC>>CCOC(=O)/C=C/c1ccc(C=O)cc1 (Heck - Acrylate)",
    "Brc1ccncc1.C=CC#N>>N#C/C=C/c1ccncc1 (Heck - Acrylonitrile)",
    
    # Negishi Coupling
    "Brc1ccccc1.c1ccc([Zn]Cl)cc1>>c1ccc(-c2ccccc2)cc1 (Negishi - Ph-Ph)",
    
    # ═══════════════════════════════════════════════════════════
    # C-N COUPLING REACTIONS (Comprehensive Test Set)
    # Both Ullmann (Cu-catalyzed) and Buchwald-Hartwig (Pd-catalyzed)
    # ═══════════════════════════════════════════════════════════
    
    # Classic Buchwald-Hartwig Examples (Pd-catalyzed)
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (Buchwald-Hartwig - Diphenylamine)",
    "Clc1ccncc1.NCC>>CCNc1ccncc1 (Buchwald-Hartwig - Pyridine ethylamine)",
    "Brc1ccc(C(F)(F)F)cc1.NC1CCCCC1>>FC(F)(F)c1ccc(NC2CCCCC2)cc1 (B-H - Cyclohexylamine)",
    "Ic1ccc(C=O)cc1.Nc1ccccc1>>O=Cc1ccc(Nc2ccccc2)cc1 (B-H - Aldehyde substrate)",
    "Brc1ccc2ccccc2c1.NCC>>CCNc1ccc2ccccc2c1 (B-H - Naphthylamine)",
    "Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3[nH]2)cc1 (B-H - Benzimidazole)",

    # ═══════════════════════════════════════════════════════════
    # COMPREHENSIVE C-N COUPLING TEST SET (Ullmann & Buchwald)
    # 40+ diverse examples covering substrate scope and selectivity
    # ═══════════════════════════════════════════════════════════
    
    # Primary anilines with various aryl halides
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + aniline → diphenylamine)",
    "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (C-N - Ph-Cl + aniline → diphenylamine)",
    "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (C-N - Ph-I + aniline → diphenylamine)",
    "Brc1ccc(OC)cc1.Nc1ccccc1>>COc1ccc(Nc2ccccc2)cc1 (C-N - 4-MeO-Ph-Br + aniline)",
    "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1 (C-N - 4-Me-Ph-Br + aniline)",
    "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - 4-CF3-Ph-Br + aniline)",
    "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1 (C-N - 4-CN-Ph-Br + aniline)",
    "Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>[O-][N+](=O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-NO2-Ph-Br + aniline)",
    "Brc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC(=O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-acetyl-Ph-Br + aniline)",
    "Brc1ccc(F)cc1.Nc1ccccc1>>Fc1ccc(Nc2ccccc2)cc1 (C-N - 4-F-Ph-Br + aniline)",
    
    # Heteroaryl halides with aniline
    "Clc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1 (C-N - 4-Cl-pyridine + aniline)",
    "Brc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1 (C-N - 4-Br-pyridine + aniline)",
    "Clc1cccnc1.Nc1ccccc1>>c1ccc(Nc2cccnc2)cc1 (C-N - 3-Cl-pyridine + aniline)",
    "Brc1cnccn1.Nc1ccccc1>>c1ccc(Nc2cnccn2)cc1 (C-N - 2-Br-pyrimidine + aniline)",
    "Clc1cccc2ncccc12.Nc1ccccc1>>c1ccc(Nc2cccc3ncccc23)cc1 (C-N - 4-Cl-quinoline + aniline)",
    "Brc1ccc2[nH]ccc2c1.Nc1ccccc1>>c1ccc(Nc2ccc3[nH]ccc3c2)cc1 (C-N - 5-Br-indole + aniline)",
    
    # Substituted anilines with aryl halides
    "Brc1ccccc1.Nc1ccc(C)cc1>>Cc1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methylaniline)",
    "Brc1ccccc1.Nc1ccc(OC)cc1>>COc1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methoxyaniline)",
    "Brc1ccccc1.Nc1ccc(F)cc1>>Fc1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-fluoroaniline)",
    "Brc1ccccc1.Nc1ccc(C(F)(F)F)cc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - Ph-Br + 4-CF3-aniline)",
    "Brc1ccccc1.Nc1cccc(C)c1>>Cc1cccc(Nc2ccccc2)c1 (C-N - Ph-Br + 3-methylaniline)",
    "Brc1ccccc1.Nc1cc(C)cc(C)c1>>Cc1cc(C)cc(Nc2ccccc2)c1 (C-N - Ph-Br + 3,5-dimethylaniline)",
    
    # Primary aliphatic amines
    "Brc1ccccc1.CN>>CNc1ccccc1 (C-N - Ph-Br + methylamine)",
    "Brc1ccccc1.CCN>>CCNc1ccccc1 (C-N - Ph-Br + ethylamine)",
    "Brc1ccccc1.CCCN>>CCCNc1ccccc1 (C-N - Ph-Br + propylamine)",
    "Brc1ccccc1.CC(C)N>>CC(C)Nc1ccccc1 (C-N - Ph-Br + isopropylamine)",
    "Brc1ccccc1.CC(C)(C)N>>CC(C)(C)Nc1ccccc1 (C-N - Ph-Br + tert-butylamine)",
    "Brc1ccccc1.NCc1ccccc1>>c1ccc(CNc2ccccc2)cc1 (C-N - Ph-Br + benzylamine)",
    "Brc1ccccc1.NCCOC>>COCCNc1ccccc1 (C-N - Ph-Br + 2-methoxyethylamine)",
    
    # Secondary aliphatic amines
    "Brc1ccccc1.CCN(CC)CC>>CCN(CC)c1ccccc1 (C-N - Ph-Br + diethylamine)",
    "Brc1ccccc1.CN(C)Cc1ccccc1>>CN(Cc1ccccc1)c2ccccc2 (C-N - Ph-Br + N,N-dimethylbenzylamine)",
    
    # Cyclic amines (heterocycles)
    "Brc1ccccc1.N1CCCC1>>c1ccc(N2CCCC2)cc1 (C-N - Ph-Br + pyrrolidine)",
    "Brc1ccccc1.N1CCCCC1>>c1ccc(N2CCCCC2)cc1 (C-N - Ph-Br + piperidine)",
    "Brc1ccccc1.N1CCOCC1>>c1ccc(N2CCOCC2)cc1 (C-N - Ph-Br + morpholine)",
    "Brc1ccccc1.N1CCN(C)CC1>>CN1CCN(c2ccccc2)CC1 (C-N - Ph-Br + N-methylpiperazine)",
    "Brc1ccccc1.N1CCC(O)CC1>>OC1CCN(c2ccccc2)CC1 (C-N - Ph-Br + 4-hydroxypiperidine)",
    "Brc1ccccc1.N1Cc2ccccc2C1>>c1ccc(N2Cc3ccccc3C2)cc1 (C-N - Ph-Br + tetrahydroisoquinoline)",
    
    # Dihalide substrates (potential bis-amination)
    "Brc1ccc(Br)cc1.Nc1ccccc1>>c1ccc(Nc2ccc(Nc3ccccc3)cc2)cc1 (C-N - 1,4-dibromobenzene + aniline)",
    
    # Ortho-substituted challenging substrates
    "Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1 (C-N - 2-methylbromobenzene + aniline)",
    "Brc1ccccc1OC.Nc1ccccc1>>COc1ccccc1Nc1ccccc1 (C-N - 2-bromoanisole + aniline)",
    "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1 (C-N - 2-bromoacetophenone + aniline)",
    
    # Pharmaceutical-relevant scaffolds
    "Clc1ccc2nc(Cl)ccc2c1.Nc1ccccc1>>c1ccc(Nc2ccc3cc(Nc4ccccc4)ccc3n2)cc1 (C-N - 4,7-dichloroquinoline + aniline)",
    "Brc1ccc2c(c1)OCO2.Nc1ccccc1>>c1ccc(Nc2ccc3c(c2)OCO3)cc1 (C-N - 5-bromobenzo[d][1,3]dioxole + aniline)",
    "Clc1nc2ccccc2s1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3s2)cc1 (C-N - 2-chlorobenzothiazole + aniline)",
    
    # Chan-Lam Coupling (C-N)
    "c1ccccc1B(O)O.Nc1ccccc1>>[O]>>c1ccc(Nc2ccccc2)cc1 (Chan-Lam - Oxidative)",
    
    # ═══════════════════════════════════════════════════════════
    # C-O COUPLING REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Ullmann Ether Synthesis
    "Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1 (Ullmann Ether - Diphenyl ether)",
    "Ic1ccncc1.OCC>>CCOc1ccncc1 (Ullmann Ether - Ethyl pyridyl ether)",
    
    # ═══════════════════════════════════════════════════════════
    # ESTERIFICATION & AMIDATION
    # ═══════════════════════════════════════════════════════════
    
    "CC(=O)O.CCO>>CC(=O)OCC (Esterification - Ethyl acetate)",
    "c1ccc(C(=O)O)cc1.OCC>>c1ccc(C(=O)OCC)cc1 (Esterification - Ethyl benzoate)",
    "CC(=O)Cl.CCO>>CC(=O)OCC (Acyl chloride esterification)",
    "CC(=O)O.Nc1ccccc1>>CC(=O)Nc1ccccc1 (Amidation - Acetanilide)",
    "c1ccc(C(=O)O)cc1.NCC>>c1ccc(C(=O)NCC)cc1 (Amidation - N-ethylbenzamide)",
    
    # ═══════════════════════════════════════════════════════════
    # REDUCTION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Hydrogenation
    "C=Cc1ccccc1>>[H][H]>>CCc1ccccc1 (Hydrogenation - Ethylbenzene)",
    "c1ccc(C=O)cc1>>[H][H]>>c1ccc(CO)cc1 (Hydrogenation - Benzyl alcohol)",
    "C#Cc1ccccc1>>[H][H]>>CCc1ccccc1 (Hydrogenation - Complete alkyne)",
    "c1ccc([N+](=O)[O-])cc1>>[H][H]>>c1ccc(N)cc1 (Hydrogenation - Nitro to amine)",
    
    # Metal Hydride Reductions
    "c1ccc(C=O)cc1.[BH4-].[Na+]>>c1ccc(CO)cc1 (NaBH4 - Benzaldehyde reduction)",
    "CC(=O)c1ccccc1.[BH4-].[Na+]>>CC(O)c1ccccc1 (NaBH4 - Acetophenone)",
    "c1ccc(C(=O)C)cc1.[AlH4-].[Li+]>>c1ccc(C(O)C)cc1 (LiAlH4 - Ketone reduction)",
    "CC(=O)OCC.[AlH4-].[Li+]>>CCO (LiAlH4 - Ester to alcohol)",
    "c1ccc(C#N)cc1.[AlH4-].[Li+]>>c1ccc(CN)cc1 (LiAlH4 - Nitrile reduction)",
    
    # ═══════════════════════════════════════════════════════════
    # OXIDATION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    "c1ccc(CO)cc1.[O]>>c1ccc(C=O)cc1 (Oxidation - Benzyl alcohol to aldehyde)",
    "CC(O)c1ccccc1.[O]>>CC(=O)c1ccccc1 (Oxidation - Secondary alcohol to ketone)",
    "CCO.[O]>>CC=O (Oxidation - Ethanol to acetaldehyde)",
    "c1ccc(CO)cc1.Cl[Cr](=O)(=O)OC(C)(C)C>>c1ccc(C=O)cc1 (PCC Oxidation)",
    "CCO.OS(=O)(=O)O>>CC=O (Swern-type oxidation)",
    
    # ═══════════════════════════════════════════════════════════
    # CYCLOADDITION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Diels-Alder
    "C=CC=C.C=C>>C1C=CCC1 (Diels-Alder - Simple)",
    "C=CC(=C)C.C=CC(=O)OCC>>CC1=CC(C(=O)OCC)CC(C)=C1 (Diels-Alder - Substituted)",
    "c1ccc2ccccc2c1C=CC=C.C=CC(=O)C>>O=C1CCc2c(ccc3ccccc23)C1 (Diels-Alder - Naphthyl)",
    
    # 1,3-Dipolar Cycloaddition (Click Chemistry)
    "C#CCO.[N-]=[N+]=Nc1ccccc1>>c1ccc(n2cc(CO)nn2)cc1 (CuAAC Click - Triazole)",
    "C#Cc1ccccc1.N=[N+]=[N-]C>>Cc1nn(-c2ccccc2)cc1-c1ccccc1 (Click - Phenyl triazole)",
    
    # ═══════════════════════════════════════════════════════════
    # SUBSTITUTION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # SN2 Reactions
    "CCCCBr.[OH-]>>CCCCO (SN2 - Primary alcohol)",
    "CC(C)CH2Br.[OH-]>>CC(C)CO (SN2 - Primary branched)",
    "CCBr.N#C>>CCC#N (SN2 - Nitrile formation)",
    "CCCCBr.Oc1ccccc1>>CCCCOc1ccccc1 (SN2 - Phenoxide)",
    "CH3I.SC>>CCSC (SN2 - Thioether formation)",
    
    # ═══════════════════════════════════════════════════════════
    # C-S COUPLING REACTIONS (Thioether formation via coupling)
    # ═══════════════════════════════════════════════════════════
    
    # Pd-catalyzed arylation of thiols (generic examples)
    "Brc1ccccc1.SCc1ccccc1>>c1ccc(S-c2ccccc2)cc1 (C-S Coupling - Thioether Formation)",
    "Clc1ccc(C)cc1.SCc1ccccc1>>Cc1ccc(S-c2ccccc2)cc1 (C-S Coupling - Thioether Formation)",
    
    # SN1 Reactions  
    "CC(C)(C)Br.O>>CC(C)(C)O (SN1 - tert-Butyl)",
    "c1ccc(C(C)(C)Cl)cc1.O>>c1ccc(C(C)(C)O)cc1 (SN1 - Benzyl tertiary)",
    
    # ═══════════════════════════════════════════════════════════
    # ELIMINATION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # E2 Eliminations
    "CC(Br)C>>C=CC (E2 - Simple alkene)",
    "CC(Br)C(C)C>>CC=C(C)C (E2 - Zaitsev product)",
    "c1ccc(C(C)Br)cc1>>c1ccc(/C=C/)cc1 (E2 - Styrene formation)",
    "CCCBr>>CC=C (E2 - Propene)",
    
    # E1 Eliminations
    "CC(C)(C)Br>>CC(C)=C (E1 - tert-Butyl)",
    
    # ═══════════════════════════════════════════════════════════
    # CARBONYL CHEMISTRY
    # ═══════════════════════════════════════════════════════════
    
    # Aldol Reactions
    "CC=O.CC=O>>CC(O)CC=O (Aldol - Aldol product)",
    "CC(=O)C.CC=O>>CC(O)C(C)C=O (Aldol - Mixed aldol)",
    "c1ccc(C=O)cc1.CC(=O)C>>c1ccc(C(O)C(C)C=O)cc1 (Aldol - Benzaldehyde)",
    
    # Wittig Reactions
    "CC=O.C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[base]>>C=CC (Wittig - Propene)",
    "c1ccc(C=O)cc1.C=C[P+](c1ccccc1)(c1ccccc1)c1ccccc1>>c1ccc(/C=C/C=C)cc1 (Wittig - Diene)",
    
    # Mannich Reaction
    "CC=O.NC.C=O>>CC(=O)CN (Mannich - Simple)",
    
    # ═══════════════════════════════════════════════════════════
    # ORGANOMETALLIC REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Grignard Reactions
    "CC=O.[Mg]Br.CC>>CC(O)(C)C (Grignard - Tertiary alcohol)",
    "c1ccc(C=O)cc1.[Mg]Br.C>>c1ccc(C(O)C)cc1 (Grignard - Secondary alcohol)",
    "CC(=O)C.[Mg]Br.c1ccccc1>>CC(O)(c1ccccc1)C (Grignard - Phenyl addition)",
    
    # ═══════════════════════════════════════════════════════════
    # METATHESIS REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    "C=CCCCC=C>>C=C.C1CCC=CC1 (Ring-Closing Metathesis)",
    "C=CC=C.C=C>>C=CCC=C (Cross-Metathesis)",
    
    # ═══════════════════════════════════════════════════════════
    # REARRANGEMENT REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    "CC(C)C(C([O-])=O)C(C)C>>CC(C)CC(C)(C)C(=O)O (Rearrangement)",
    
    # ═══════════════════════════════════════════════════════════
    # MULTICOMPONENT REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    "CC=O.NC.CC(=O)O>>CC1=C(C)OC(=O)C(C)(N)C1 (Ugi-type reaction)",
    
    # ═══════════════════════════════════════════════════════════
    # HETEROCYCLE SYNTHESIS
    # ═══════════════════════════════════════════════════════════
    
    "c1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>c1ccc2[nH]c3ccccc3c2c1 (Carbazole synthesis)",
    "Nc1ccccc1.c1ccc(C=O)cc1>>c1ccc2nc(-c3ccccc3)cc(-c3ccccc3)c2c1 (Quinoline synthesis)",
    
    # ═══════════════════════════════════════════════════════════
    # AMIDE FORMATION REACTIONS (Acid + Amine → Amide)
    # ═══════════════════════════════════════════════════════════
    
    # Simple carboxylic acid + amine amidations
    "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1 (Amide: Benzoic acid + aniline)",
    "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1 (Amide: Acetic acid + benzylamine)",
    "O=C(O)CCc1ccccc1.NC1CCCCC1>>O=C(N1CCCCC1)CCc1ccccc1 (Amide: Phenylacetic acid + cyclohexylamine)",
    "CC(C)(C)C(=O)O.Nc1ccc(C)cc1>>CC(C)(C)C(=O)Nc1ccc(C)cc1 (Amide: Pivalic acid + p-toluidine)",
    "O=C(O)c1ccc(F)cc1.NCCc1ccccc1>>O=C(NCCc1ccccc1)c1ccc(F)cc1 (Amide: 4-Fluorobenzoic acid + phenethylamine)",
    
    # Aromatic acids with different amines
    "O=C(O)c1ccc(OC)cc1.NC1CC1>>COc1ccc(C(=O)NC2CC2)cc1 (Amide: p-Methoxybenzoic acid + cyclopropylamine)",
    "O=C(O)c1ccc(C#N)cc1.NCCCC>>N#Cc1ccc(C(=O)NCCCC)cc1 (Amide: 4-Cyanobenzoic acid + butylamine)",
    "O=C(O)c1ccc([N+](=O)[O-])cc1.NC(C)(C)C>>CC(C)(C)NC(=O)c1ccc([N+](=O)[O-])cc1 (Amide: 4-Nitrobenzoic acid + tert-butylamine)",
    "O=C(O)c1cccc(Cl)c1.NCc1ccncc1>>O=C(NCc1ccncc1)c1cccc(Cl)c1 (Amide: 3-Chlorobenzoic acid + 4-picolylamine)",
    "O=C(O)c1ccc2ccccc2c1.NC1CCNCC1>>O=C(N1CCC(N)CC1)c1ccc2ccccc2c1 (Amide: 2-Naphthoic acid + 4-aminopiperidine)",
    
    # Heteroaromatic carboxylic acids
    "O=C(O)c1ccncc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccncc1 (Amide: Isonicotinic acid + aniline)",
    "O=C(O)c1cnccn1.NCc1ccc(F)cc1>>O=C(NCc1ccc(F)cc1)c1cnccn1 (Amide: Pyrimidine-2-carboxylic acid + 4-fluorobenzylamine)",
    "O=C(O)c1cccs1.NC1CCOCC1>>O=C(N1CCOCC1)c1cccs1 (Amide: Thiophene-3-carboxylic acid + morpholine)",
    "O=C(O)c1cc[nH]c1.NCCOc1ccccc1>>O=C(NCCOc1ccccc1)c1cc[nH]c1 (Amide: Pyrrole-3-carboxylic acid + 2-phenoxyethylamine)",
    "O=C(O)c1ccc[nH]1.NC1CCC(C)(C)CC1>>O=C(NC1CCC(C)(C)CC1)c1ccc[nH]1 (Amide: Pyrrole-2-carboxylic acid + 4,4-dimethylcyclohexylamine)",
    
    # Aliphatic carboxylic acids  
    "CCCCCCCC(=O)O.Nc1ccc(OC)cc1>>CCCCCCCC(=O)Nc1ccc(OC)cc1 (Amide: Octanoic acid + p-anisidine)",
    "CC(C)CC(=O)O.NCc1ccc(Cl)cc1>>CC(C)CC(=O)NCc1ccc(Cl)cc1 (Amide: Isovaleric acid + 4-chlorobenzylamine)",
    "CCC(C)C(=O)O.NC1CCN(C)CC1>>CCC(C)C(=O)N1CCN(C)CC1 (Amide: 2-Methylbutyric acid + N-methylpiperazine)",
    "O=C(O)C1CCCCC1.NCCc1ccc[nH]1>>O=C(NCCc1ccc[nH]1)C1CCCCC1 (Amide: Cyclohexanecarboxylic acid + tryptamine)",
    
    # Amino acids as starting materials
    "CC(C)CC(N)C(=O)O.NCc1ccccc1>>CC(C)CC(N)C(=O)NCc1ccccc1 (Amide: Leucine + benzylamine)",
    "N[C@@H](Cc1ccccc1)C(=O)O.NC1CCCCC1>>N[C@@H](Cc1ccccc1)C(=O)N1CCCCC1 (Amide: Phenylalanine + cyclohexylamine)",
    
    # Secondary amines (N-substituted amides)
    "O=C(O)c1ccccc1.CN(C)c1ccccc1>>CN(C)c1ccccc1.O=C(N(C)c1ccccc1)c1ccccc1 (Amide: Benzoic acid + N,N-dimethylaniline → tertiary amide)",
    "CC(=O)O.N1CCCC1>>CC(=O)N1CCCC1 (Amide: Acetic acid + pyrrolidine)",
    "O=C(O)c1ccc(Br)cc1.N1CCOCC1>>O=C(N1CCOCC1)c1ccc(Br)cc1 (Amide: 4-Bromobenzoic acid + morpholine)",
    
    # Coupling reagent-specific examples (indicated in labels)
    "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1 (Amide: EDC/HOBt coupling)",
    "CC(C)(C)C(=O)O.NCc1ccccc1>>CC(C)(C)C(=O)NCc1ccccc1 (Amide: DCC/NHS coupling)",
    "O=C(O)c1ccc(OC)cc1.NC1CCCCC1>>COc1ccc(C(=O)N1CCCCC1)cc1 (Amide: HATU/DIPEA coupling)"
]

# ═══════════════════════════════════════════════════════════
# COMPREHENSIVE BUCHWALD-HARTWIG AMINATION COLLECTION
# ═══════════════════════════════════════════════════════════

BUCHWALD_HARTWIG_REACTIONS = [
    "Select a Buchwald-Hartwig reaction...",
    
    # ═══════ PRIMARY AMINES ═══════
    
    # Aniline derivatives with aryl bromides
    "Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (B-H: PhBr + aniline → diphenylamine)",
    "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1 (B-H: 4-MeBr + aniline → 4-methyldiphenylamine)",
    "Brc1ccc(OC)cc1.Nc1ccccc1>>COc1ccc(Nc2ccccc2)cc1 (B-H: 4-MeOBr + aniline → 4-methoxydiphenylamine)",
    "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1 (B-H: 4-CF3Br + aniline → electron-poor)",
    "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1 (B-H: 4-CNBr + aniline → cyano diphenylamine)",
    "Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>O=[N+]([O-])c1ccc(Nc2ccccc2)cc1 (B-H: 4-NO2Br + aniline → nitro)",
    
    # Substituted anilines with bromobenzene
    "Brc1ccccc1.Nc1ccc(C)cc1>>Cc1ccc(Nc2ccccc2)cc1 (B-H: PhBr + 4-toluidine → N-phenyl-4-toluidine)",
    "Brc1ccccc1.Nc1ccc(OC)cc1>>COc1ccc(Nc2ccccc2)cc1 (B-H: PhBr + 4-anisidine → N-phenyl-4-anisidine)",
    "Brc1ccccc1.Nc1ccc(F)cc1>>Fc1ccc(Nc2ccccc2)cc1 (B-H: PhBr + 4-fluoroaniline → fluoro diphenylamine)",
    "Brc1ccccc1.Nc1ccc(Cl)cc1>>Clc1ccc(Nc2ccccc2)cc1 (B-H: PhBr + 4-chloroaniline → chloro diphenylamine)",
    "Brc1ccccc1.Nc1ccc(C(F)(F)F)cc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1 (B-H: PhBr + 4-CF3-aniline → CF3)",
    
    # Ortho-substituted aromatics (challenging substrates)
    "Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1 (B-H: 2-bromotoluene + aniline → sterically hindered)",
    "Brc1c(C)cccc1C.Nc1ccccc1>>Cc1cccc(C)c1Nc1ccccc1 (B-H: 2,6-dimethylbromobenzene → highly hindered)",
    "Brc1ccccc1OC.Nc1ccccc1>>COc1ccccc1Nc1ccccc1 (B-H: 2-bromoanisole + aniline → ortho methoxy)",
    
    # Heteroaryl bromides
    "Brc1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1 (B-H: 4-bromopyridine + aniline → 4-phenylaminopyridine)",
    "Brc1cccnc1.Nc1ccccc1>>c1ccc(Nc2cccnc2)cc1 (B-H: 3-bromopyridine + aniline → 3-phenylaminopyridine)",
    "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1 (B-H: 2-bromopyridine + aniline → 2-phenylaminopyridine)",
    "Brc1cncc(C)c1.Nc1ccccc1>>Cc1cncc(Nc2ccccc2)c1 (B-H: 5-bromo-2-methylpyrimidine)",
    "Brc1nc2ccccc2s1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3s2)cc1 (B-H: 2-bromobenzothiazole)",
    
    # ═══════ SECONDARY AMINES ═══════
    
    # Dialkyl amines
    "Brc1ccccc1.CN(C)C>>CN(C)c1ccccc1 (B-H: PhBr + dimethylamine → N,N-dimethylaniline)",
    "Brc1ccc(C)cc1.CN(C)C>>CN(C)c1ccc(C)cc1 (B-H: 4-MeBr + dimethylamine → N,N-dimethyl-4-toluidine)",
    "Brc1ccc(OC)cc1.CN(C)C>>CN(C)c1ccc(OC)cc1 (B-H: 4-MeOBr + dimethylamine → N,N-dimethyl-4-anisidine)",
    "Brc1ccc(C(F)(F)F)cc1.CN(C)C>>CN(C)c1ccc(C(F)(F)F)cc1 (B-H: 4-CF3Br + dimethylamine → electron-poor)",
    "Brc1ccccc1.C1CCNCC1>>C1CCN(c2ccccc2)CC1 (B-H: PhBr + piperidine → N-phenylpiperidine)",
    "Brc1ccccc1.C1COCCN1>>C1COCN(c2ccccc2)C1 (B-H: PhBr + morpholine → N-phenylmorpholine)",
    "Brc1ccccc1.c1ccncc1>>c1ccc(N2CCNCC2)cc1 (B-H: PhBr + piperazine → N-phenylpiperazine)",
    
    # Ethyl substituted amines
    "Brc1ccccc1.CCN(CC)C>>CCN(CC)c1ccccc1 (B-H: PhBr + diethylamine → N,N-diethylaniline)",
    "Brc1ccccc1.CCNC>>CCNc1ccccc1 (B-H: PhBr + ethylamine → N-ethylaniline)",
    "Brc1ccccc1.CCCNC>>CCCNc1ccccc1 (B-H: PhBr + propylamine → N-propylaniline)",
    
    # Cyclical secondary amines
    "Brc1ccccc1.C1CCCC1N>>C1CCCC1Nc1ccccc1 (B-H: PhBr + cyclopentylamine → N-cyclopentylaniline)",
    "Brc1ccccc1.C1CCC(N)CC1>>C1CCC(Nc2ccccc2)CC1 (B-H: PhBr + cyclohexylamine → N-cyclohexylaniline)",
    
    # ═══════ CHLORO SUBSTRATES ═══════
    
    # Aryl chlorides (more challenging than bromides)
    "Clc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (B-H: PhCl + aniline → diphenylamine, difficult)",
    "Clc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1 (B-H: 4-ClCN + aniline → activated ArCl)",
    "Clc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>O=[N+]([O-])c1ccc(Nc2ccccc2)cc1 (B-H: 4-ClNO2 + aniline)",
    "Clc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC(=O)c1ccc(Nc2ccccc2)cc1 (B-H: 4-ClCOMe + aniline → activated)",
    "Clc1ccncc1.CN(C)C>>CN(C)c1ccncc1 (B-H: 4-chloropyridine + dimethylamine → electron-poor)",
    "Clc1nc2ccccc2n1.Nc1ccccc1>>c1ccc(Nc2nc3ccccc3n2)cc1 (B-H: 2-chloroquinoxaline + aniline)",
    
    # ═══════ IODO SUBSTRATES ═══════
    
    # Aryl iodides (most reactive)
    "Ic1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1 (B-H: PhI + aniline → diphenylamine, fast)",
    "Ic1ccc(C)cc1.CN(C)C>>CN(C)c1ccc(C)cc1 (B-H: 4-iodotoluene + dimethylamine → fast)",
    "Ic1ccncc1.Nc1ccccc1>>c1ccc(Nc2ccncc2)cc1 (B-H: 4-iodopyridine + aniline → heteroaryl)",
    
    # ═══════ BENZYL AMINES ═══════
    
    # Primary benzyl amines
    "Brc1ccccc1.NCc1ccccc1>>c1ccc(CNc2ccccc2)cc1 (B-H: PhBr + benzylamine → N-benzylaniline)",
    "Brc1ccc(C)cc1.NCc1ccccc1>>Cc1ccc(NCc2ccccc2)cc1 (B-H: 4-MeBr + benzylamine)",
    "Brc1ccc(OC)cc1.NCc1ccccc1>>COc1ccc(NCc2ccccc2)cc1 (B-H: 4-MeOBr + benzylamine)",
    
    # Substituted benzyl amines
    "Brc1ccccc1.NCc1ccc(C)cc1>>Cc1ccc(CNc2ccccc2)cc1 (B-H: PhBr + 4-methylbenzylamine)",
    "Brc1ccccc1.NCc1ccc(OC)cc1>>COc1ccc(CNc2ccccc2)cc1 (B-H: PhBr + 4-methoxybenzylamine)",
    "Brc1ccccc1.NCc1ccc(F)cc1>>Fc1ccc(CNc2ccccc2)cc1 (B-H: PhBr + 4-fluorobenzylamine)",
    
    # ═══════ INTRAMOLECULAR CYCLIZATIONS ═══════
    
    # Formation of carbazoles and related heterocycles
    "Brc1ccccc1Nc1ccccc1Br>>c1ccc2[nH]c3ccccc3c2c1 (B-H: intramolecular → carbazole formation)",
    "Brc1ccc2ccccc2c1Nc1ccccc1>>c1ccc2c(c1)c1ccccc1n2c1ccccc1 (B-H: naphthyl cyclization)",
    
    # ═══════ AMINO ACID AND PEPTIDE SUBSTRATES ═══════
    
    # Amino acid derivatives
    "Brc1ccccc1.NCCC(=O)O>>O=C(O)CCNc1ccccc1 (B-H: PhBr + β-alanine → N-phenyl-β-alanine)",
    "Brc1ccccc1.NCC(=O)OC>>COC(=O)CNc1ccccc1 (B-H: PhBr + glycine ester → N-phenylglycine ester)",
    
    # ═══════ CHALLENGING STERICALLY HINDERED CASES ═══════
    
    # Highly substituted aromatics
    "Brc1c(C)c(C)c(C)c(C)c1C.Nc1ccccc1>>Cc1c(C)c(C)c(C)c(Nc2ccccc2)c1C (B-H: pentamethylbromobenzene)",
    "Brc1ccc2c(c1)cccc2C(C)(C)C.Nc1ccccc1>>CC(C)(C)c1cccc2cc(Nc3ccccc3)ccc12 (B-H: bulky naphthyl)",
    
    # ═══════ PHARMACEUTICAL INTERMEDIATES ═══════
    
    # Drug-like substrates
    "Brc1ccc(S(=O)(=O)N)cc1.Nc1ccccc1>>NS(=O)(=O)c1ccc(Nc2ccccc2)cc1 (B-H: sulfonamide derivative)",
    "Brc1ccc(C(=O)N)cc1.CN(C)C>>CN(C)c1ccc(C(=O)N)cc1 (B-H: benzamide derivative)",
    "Clc1ccc2nc(Cl)nc(N)c2c1.Nc1ccccc1>>Nc1nc(Nc2ccccc2)nc2ccc(Cl)cc12 (B-H: quinazoline scaffold)",
    
    # ═══════ COMPLEX FUNCTIONAL GROUP TOLERANCE ═══════
    
    # Multiple functional groups
    "Brc1cc(C(=O)OC)cc(OC)c1OC.Nc1ccccc1>>COC(=O)c1cc(Nc2ccccc2)c(OC)c(OC)c1 (B-H: complex ester)",
    "Brc1ccc(CCC(=O)O)cc1.CN(C)C>>O=C(O)CCc1ccc(N(C)C)cc1 (B-H: carboxylic acid tolerance)",
    "Brc1ccc(c2ccccc2)cc1.Nc1ccccc1>>c1ccc(-c2ccc(Nc3ccccc3)cc2)cc1 (B-H: biphenyl system)",
    
    # ═══════ HETEROCYCLE-CONTAINING AMINES ═══════
    
    # Pyridyl amines
    "Brc1ccccc1.Nc1ccccn1>>c1ccc(Nc2ccccn2)cc1 (B-H: PhBr + 2-aminopyridine → 2-phenylaminopyridine)",
    "Brc1ccccc1.Nc1cccnc1>>c1ccc(Nc2cccnc2)cc1 (B-H: PhBr + 3-aminopyridine → 3-phenylaminopyridine)",
    "Brc1ccccc1.Nc1ccncc1>>c1ccc(Nc2ccncc2)cc1 (B-H: PhBr + 4-aminopyridine → 4-phenylaminopyridine)",
    
    # Thiazole and other heteroaryl amines
    "Brc1ccccc1.Nc1nccs1>>c1ccc(Nc2nccs2)cc1 (B-H: PhBr + 2-aminothiazole)",
    "Brc1ccccc1.Nc1ccco1>>c1ccc(Nc2ccco2)cc1 (B-H: PhBr + 3-aminofuran)",
    "Brc1ccccc1.Nc1ccc[nH]1>>c1ccc(Nc2ccc[nH]2)cc1 (B-H: PhBr + 3-aminopyrrole)",
    
    # ═══════ MACROCYCLE PRECURSORS ═══════
    
    # Long-chain substrates for macrocyclization
    "Brc1ccccc1CCCCCCCCC(=O)Nc1ccccc1Br>>O=C1CCCCCCCc2ccccc2N1c1ccccc1 (B-H: macrocycle formation)",
    
    # ═══════════════════════════════════════════════════════════
    # AMIDE FORMATION REACTIONS (Acid + Amine → Amide)
    # Various coupling agents and conditions
    # ═══════════════════════════════════════════════════════════
    
    # Simple carboxylic acid + amine amidations
    "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1 (Amide: Benzoic acid + aniline)",
    "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1 (Amide: Acetic acid + benzylamine)",
    "O=C(O)CCc1ccccc1.NC1CCCCC1>>O=C(N1CCCCC1)CCc1ccccc1 (Amide: Phenylacetic acid + cyclohexylamine)",
    "CC(C)(C)C(=O)O.Nc1ccc(C)cc1>>CC(C)(C)C(=O)Nc1ccc(C)cc1 (Amide: Pivalic acid + p-toluidine)",
    "O=C(O)c1ccc(F)cc1.NCCc1ccccc1>>O=C(NCCc1ccccc1)c1ccc(F)cc1 (Amide: 4-Fluorobenzoic acid + phenethylamine)",
    
    # Aromatic acids with different amines
    "O=C(O)c1ccc(OC)cc1.NC1CC1>>COc1ccc(C(=O)NC2CC2)cc1 (Amide: p-Methoxybenzoic acid + cyclopropylamine)",
    "O=C(O)c1ccc(C#N)cc1.NCCCC>>N#Cc1ccc(C(=O)NCCCC)cc1 (Amide: 4-Cyanobenzoic acid + butylamine)",
    "O=C(O)c1ccc([N+](=O)[O-])cc1.NC(C)(C)C>>CC(C)(C)NC(=O)c1ccc([N+](=O)[O-])cc1 (Amide: 4-Nitrobenzoic acid + tert-butylamine)",
    "O=C(O)c1cccc(Cl)c1.NCc1ccncc1>>O=C(NCc1ccncc1)c1cccc(Cl)c1 (Amide: 3-Chlorobenzoic acid + 4-picolylamine)",
    "O=C(O)c1ccc2ccccc2c1.NC1CCNCC1>>O=C(N1CCC(N)CC1)c1ccc2ccccc2c1 (Amide: 2-Naphthoic acid + 4-aminopiperidine)",
    
    # Heteroaromatic carboxylic acids
    "O=C(O)c1ccncc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccncc1 (Amide: Isonicotinic acid + aniline)",
    "O=C(O)c1cnccn1.NCc1ccc(F)cc1>>O=C(NCc1ccc(F)cc1)c1cnccn1 (Amide: Pyrimidine-2-carboxylic acid + 4-fluorobenzylamine)",
    "O=C(O)c1cccs1.NC1CCOCC1>>O=C(N1CCOCC1)c1cccs1 (Amide: Thiophene-3-carboxylic acid + morpholine)",
    "O=C(O)c1cc[nH]c1.NCCOc1ccccc1>>O=C(NCCOc1ccccc1)c1cc[nH]c1 (Amide: Pyrrole-3-carboxylic acid + 2-phenoxyethylamine)",
    "O=C(O)c1ccc[nH]1.NC1CCC(C)(C)CC1>>O=C(NC1CCC(C)(C)CC1)c1ccc[nH]1 (Amide: Pyrrole-2-carboxylic acid + 4,4-dimethylcyclohexylamine)",
    
    # Aliphatic carboxylic acids
    "CCCCCCCC(=O)O.Nc1ccc(OC)cc1>>CCCCCCCC(=O)Nc1ccc(OC)cc1 (Amide: Octanoic acid + p-anisidine)",
    "CC(C)CC(=O)O.NCc1ccc(Cl)cc1>>CC(C)CC(=O)NCc1ccc(Cl)cc1 (Amide: Isovaleric acid + 4-chlorobenzylamine)",
    "CCC(C)C(=O)O.NC1CCN(C)CC1>>CCC(C)C(=O)N1CCN(C)CC1 (Amide: 2-Methylbutyric acid + N-methylpiperazine)",
    "O=C(O)C1CCCCC1.NCCc1ccc[nH]1>>O=C(NCCc1ccc[nH]1)C1CCCCC1 (Amide: Cyclohexanecarboxylic acid + tryptamine)",
    "CC1CC(C(=O)O)CC(C)(C)C1.Nc1cccc(C)c1>>CC1CC(C(=O)Nc2cccc(C)c2)CC(C)(C)C1 (Amide: Camphoric acid derivative + m-toluidine)",
    
    # Amino acids as starting materials
    "CC(C)CC(N)C(=O)O.NCc1ccccc1>>CC(C)CC(N)C(=O)NCc1ccccc1 (Amide: Leucine + benzylamine)",
    "N[C@@H](Cc1ccccc1)C(=O)O.NC1CCCCC1>>N[C@@H](Cc1ccccc1)C(=O)N1CCCCC1 (Amide: Phenylalanine + cyclohexylamine)",
    "N[C@@H](CC(=O)O)C(=O)O.NCc1ccc(OC)cc1>>N[C@@H](CC(=O)O)C(=O)NCc1ccc(OC)cc1 (Amide: Aspartic acid + 4-methoxybenzylamine)",
    "N[C@@H](CCCCN)C(=O)O.Nc1ccc(C(F)(F)F)cc1>>N[C@@H](CCCCN)C(=O)Nc1ccc(C(F)(F)F)cc1 (Amide: Lysine + 4-trifluoromethylaniline)",
    
    # Dicarboxylic acids (mono-amide formation)
    "O=C(O)CCCCCC(=O)O.Nc1ccccc1>>O=C(O)CCCCCC(=O)Nc1ccccc1 (Amide: Adipic acid mono-amide + aniline)",
    "O=C(O)c1ccc(C(=O)O)cc1.NCc1ccncc1>>O=C(O)c1ccc(C(=O)NCc2ccncc2)cc1 (Amide: Terephthalic acid mono-amide + 4-picolylamine)",
    "O=C(O)CCC(=O)O.NC1CCOCC1>>O=C(O)CCC(=O)N1CCOCC1 (Amide: Succinic acid mono-amide + morpholine)",
    
    # Secondary amines (N-substituted amides)
    "O=C(O)c1ccccc1.CN(C)c1ccccc1>>CN(C)c1ccccc1.O=C(N(C)c1ccccc1)c1ccccc1 (Amide: Benzoic acid + N,N-dimethylaniline → tertiary amide)",
    "CC(=O)O.N1CCCC1>>CC(=O)N1CCCC1 (Amide: Acetic acid + pyrrolidine)",
    "O=C(O)c1ccc(Br)cc1.N1CCOCC1>>O=C(N1CCOCC1)c1ccc(Br)cc1 (Amide: 4-Bromobenzoic acid + morpholine)",
    "O=C(O)CCc1ccccc1.N1CCN(c2ccccc2)CC1>>O=C(N1CCN(c2ccccc2)CC1)CCc1ccccc1 (Amide: Phenylacetic acid + N-phenylpiperazine)",
    
    # Substituted anilines
    "O=C(O)c1ccccc1.Nc1ccc(Cl)c(Cl)c1>>O=C(Nc1ccc(Cl)c(Cl)c1)c1ccccc1 (Amide: Benzoic acid + 3,4-dichloroaniline)",
    "CC(=O)O.Nc1cc(C)cc(C)c1>>CC(=O)Nc1cc(C)cc(C)c1 (Amide: Acetic acid + 3,5-dimethylaniline)",
    "O=C(O)c1ccc(F)cc1.Nc1ccc(OCC(F)(F)F)cc1>>O=C(Nc1ccc(OCC(F)(F)F)cc1)c1ccc(F)cc1 (Amide: 4-Fluorobenzoic acid + 4-trifluoroethoxyaniline)",
    
    # Coupling reagent-specific examples (indicated in labels)
    "O=C(O)c1ccccc1.Nc1ccccc1>>O=C(Nc1ccccc1)c1ccccc1 (Amide: EDC/HOBt coupling)",
    "CC(C)(C)C(=O)O.NCc1ccccc1>>CC(C)(C)C(=O)NCc1ccccc1 (Amide: DCC/NHS coupling)",
    "O=C(O)c1ccc(OC)cc1.NC1CCCCC1>>COc1ccc(C(=O)N1CCCCC1)cc1 (Amide: HATU/DIPEA coupling)",
    "O=C(O)CCc1ccccc1.Nc1ccncc1>>O=C(Nc1ccncc1)CCc1ccccc1 (Amide: PyBOP/NMM coupling)",
    "O=C(O)c1cccs1.NCCc1ccccc1>>O=C(NCCc1ccccc1)c1cccs1 (Amide: T3P/Et3N coupling)",
    
    # Mixed aliphatic/aromatic systems
    "O=C(O)C(c1ccccc1)c1ccccc1.NCCCCCC>>O=C(NCCCCCC)C(c1ccccc1)c1ccccc1 (Amide: Diphenylacetic acid + hexylamine)",
    "CCCCCCC(=O)O.Nc1ccc2[nH]c3ccccc3c2c1>>CCCCCCC(=O)Nc1ccc2[nH]c3ccccc3c2c1 (Amide: Heptanoic acid + 2-aminocarbazole)",
    "O=C(O)c1ccc(Oc2ccccc2)cc1.NC1CC1C(F)(F)F>>O=C(NC1CC1C(F)(F)F)c1ccc(Oc2ccccc2)cc1 (Amide: 4-Phenoxybenzoic acid + 2-trifluoromethylcyclopropylamine)",
]

def get_sample_reactions():
    """Get the list of all sample reactions"""
    return SAMPLE_REACTIONS

def get_buchwald_hartwig_reactions():
    """Get comprehensive Buchwald-Hartwig amination reaction collection"""
    return BUCHWALD_HARTWIG_REACTIONS

def get_coupling_reactions():
    """Get only coupling reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(coupling in r for coupling in 
            ["Suzuki", "Stille", "Sonogashira", "Heck", "Negishi", "Buchwald-Hartwig", "Chan-Lam", "Ullmann"])]

def get_amide_formation_reactions():
    """Get amide formation reaction examples (acid + amine → amide)"""
    return [r for r in SAMPLE_REACTIONS if "(Amide:" in r]

def get_cc_coupling_reactions():
    """Get C-C coupling examples (Suzuki, Stille, Sonogashira, Heck, Negishi, Kumada)"""
    tokens = ["Suzuki", "Stille", "Sonogashira", "Heck", "Negishi", "Kumada"]
    return [r for r in SAMPLE_REACTIONS if any(t in r for t in tokens)]

def get_cn_coupling_reactions():
    """Get C-N coupling examples (Buchwald-Hartwig, Ullmann C-N, Chan-Lam)"""
    tokens = ["Buchwald-Hartwig", "Ullmann C-N", "Chan-Lam"]
    return [r for r in SAMPLE_REACTIONS if any(t in r for t in tokens)]

def get_co_coupling_reactions():
    """Get C-O coupling examples (Ullmann Ether, Mitsunobu)"""
    tokens = ["Ullmann Ether", "Mitsunobu"]
    return [r for r in SAMPLE_REACTIONS if any(t in r for t in tokens)]

def get_cs_coupling_reactions():
    """Get C-S coupling examples (Thioether Formation)"""
    tokens = ["C-S Coupling", "Thioether Formation"]
    return [r for r in SAMPLE_REACTIONS if any(t in r for t in tokens)]

def get_reduction_reactions():
    """Get only reduction reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(reduction in r for reduction in 
            ["Hydrogenation", "NaBH4", "LiAlH4", "reduction"])]

def get_oxidation_reactions():
    """Get only oxidation reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(oxidation in r for oxidation in 
            ["Oxidation", "PCC", "Swern"])]

def get_substitution_reactions():
    """Get only substitution reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(sub in r for sub in ["SN1", "SN2"])]

def get_elimination_reactions():
    """Get only elimination reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(elim in r for elim in ["E1", "E2"])]

def get_cycloaddition_reactions():
    """Get only cycloaddition reaction examples"""
    return [r for r in SAMPLE_REACTIONS if any(cyclo in r for cyclo in 
            ["Diels-Alder", "Click", "Cycloaddition"])]

def search_reactions(query):
    """Search for reactions containing the query string"""
    query_lower = query.lower()
    return [r for r in SAMPLE_REACTIONS if query_lower in r.lower()]
