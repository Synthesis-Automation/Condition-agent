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
    "Fc1ccc(OS(=O)(=O)C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>Fc1ccc(-c2ccccc2)cc1 (Suzuki - Aryl triflate electrophile)",
    "C/C=C/Br.C=CB(O)O>>C/C=C/C=C (Suzuki - Vinyl bromide + vinyl boronic acid)",
    "Brc1ccc(OC)cc1.[K+].[B-](F)(F)c2ccccc2>>COc1ccc(-c2ccccc2)cc1 (Suzuki - Potassium aryltrifluoroborate)",
    "Clc1ccncc1.c1ccncc1B1OC(=O)CN(CC(=O)O)C(=O)CN(CC(=O)O)C(=O)O1>>c1ccc(-c2ccccn2)nc1 (Suzuki - 2-pyridyl MIDA slow-release)",
    "Brc1ccc(Cl)nc1.c2ccc(B(O)O)cc2>>c1ccc(-c2ccc(Cl)nc2)cc1 (Suzuki - Dichloropyridine)",
    "Brc1cc(Cl)nc(Cl)c1.c2ccc(B(O)O)nc2>>c1ccc(-c2cc(Cl)nc(Cl)c2)nc1 (Suzuki - Mixed halide pyridine)",
    "Brc1cccc2[nH]ccc12.c3ccc(B(O)O)cc3>>c1ccc(-c2ccc3[nH]ccc3c2)cc1 (Suzuki - Indole heterobiaryl)",
    "FC(F)(F)c1ccc(Cl)cc1.Nc1ccnc(B(O)O)n1>>FC(F)(F)c1ccc(-c2cc(N)ncn2)cc1 (Suzuki - Trifluoromethyl aryl chloride + pyrimidine boronate)",
    "Clc1cncc(B(O)O)c1.Brc2ccc(F)cc2>>Fc1ccc(-c2cc(Cl)cnc2)cc1 (Suzuki - Chloropyridyl boronic acid with aryl bromide)",
    "Brc1ccc(Br)cc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccc(-c3ccccc3)cc2)cc1 (Suzuki - Bis-coupling to biphenyl)",
    "Brc1cccc(Br)c1.c1ccc(B(O)O)cc1>>c1ccc(-c2cccc(-c3ccccc3)c2)cc1 (Suzuki - Meta-dibromobenzene)",
    "Brc1ccccc1.c1ncc(B(O)O)cn1>>c1ccc(-c2cnccn2)cc1 (Suzuki - Pyrimidine-5-boronic acid)",
    "Clc1ccc(C#N)cc1.c1coc(B(O)O)c1>>N#Cc1ccc(-c2ccoc2)cc1 (Suzuki - Furan-2-boronic acid)",
    "Brc1ccncc1.c1csc(B(O)O)c1>>c1cc(-c2cccs2)ccn1 (Suzuki - Thiophene-2-boronic acid)",
    "Ic1ccc(C(=O)OC)cc1.c1c[nH]c(B(O)O)c1>>COC(=O)c1ccc(-c2cc[nH]c2)cc1 (Suzuki - Pyrrole-3-boronic acid)",
    "Brc1cnc2ccccc2n1.c1ccc(B(O)O)cc1>>c1ccc(-c2cnc3ccccc3n2)cc1 (Suzuki - 3-Bromoquinoxaline)",
    "Clc1cccc2c1cc[nH]2.c1ccc(B(O)O)nc1>>c1ccc(-c2cccc3[nH]ccc23)nc1 (Suzuki - 4-Chloroindole + pyridine boronic acid)",
    "Brc1nc2ccccc2s1.Cc1ccc(B(O)O)cc1>>Cc1ccc(-c2nc3ccccc3s2)cc1 (Suzuki - 2-Bromobenzothiazole)",
    "Ic1nc2ccccc2o1.COc1ccc(B(O)O)cc1>>COc1ccc(-c2nc3ccccc3o2)cc1 (Suzuki - 2-Iodobenzoxazole)",
    "Brc1cccc2nsnc12.c1ccc(B(O)O)cc1>>c1ccc(-c2cccc3nsnc23)cc1 (Suzuki - Benzothiadiazole)",
    "Brc1ccccc1C(C)C.c1ccc(B(O)O)cc1>>CC(C)c1ccccc1-c1ccccc1 (Suzuki - 2-Isopropyl-bromobenzene)",
    "Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1 (Suzuki - 2,6-Dimethyl-chlorobenzene)",
    "Brc1ccccc1OCC.c1ccc(B(O)O)cc1C>>CCOc1ccccc1-c1ccccc1C (Suzuki - Ortho-ethoxy + ortho-methyl)",
    "Fc1c(F)c(F)c(Br)c(F)c1F.c1ccc(B(O)O)cc1>>Fc1c(F)c(F)c(-c2ccccc2)c(F)c1F (Suzuki - Pentafluorobromobenzene)",
    "Clc1cc(Cl)cc(Br)c1.c1ccc(B(O)O)cc1>>Clc1cc(Cl)cc(-c2ccccc2)c1 (Suzuki - 3,5-Dichloro-bromobenzene)",
    "Brc1ccc([N+](=O)[O-])cc1[N+](=O)[O-].c1ccc(B(O)O)cc1>>c1ccc(-c2ccc([N+](=O)[O-])cc2[N+](=O)[O-])cc1 (Suzuki - 2,5-Dinitro-bromobenzene)",
    "Brc1ccc(C(=O)OCC)cc1.c1ccc(B(O)O)cc1>>CCOC(=O)c1ccc(-c2ccccc2)cc1 (Suzuki - Ethyl ester protected)",
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.c1ccc(B(O)O)cc1>>CC(C)(C)OC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - Boc-protected amine)",
    "Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1.c1ccc(B(O)O)cc1>>CC(C)(C)[Si](C)(C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TBS-protected phenol)",
    "Brc1ccccc1.C=CB(O)O>>C=Cc1ccccc1 (Suzuki - Vinylboronic acid to styrene)",
    "Ic1ccc(C=O)cc1.C/C=C/B(O)O>>C/C=C/c1ccc(C=O)cc1 (Suzuki - (E)-Propenylboronic acid)",
    "Brc1ccncc1.C=C(C)B(O)O>>C=C(C)c1ccncc1 (Suzuki - Isopropenylboronic acid)",
    "Brc1ccc(OC)cc1.C#CB(O)O>>C#Cc1ccc(OC)cc1 (Suzuki - Ethynylboronic acid)",
    "Brc1ccc(Br)cc1CCCCCCCC(=O)O.c1ccc(B(O)O)cc1>>O=C(O)CCCCCCCc1ccc(-c2ccccc2)cc1-c1ccccc1 (Suzuki - Macrocyclization precursor)",
    "Brc1ccccc1.[K+].[B-](F)(F)(F)c1ccccc1>>c1ccc(-c2ccccc2)cc1 (Suzuki - Potassium phenyltrifluoroborate)",
    "Clc1ccc(C(F)(F)F)cc1.[K+].[B-](F)(F)(F)c1ccc(OC)cc1>>COc1ccc(-c2ccc(C(F)(F)F)cc2)cc1 (Suzuki - BF3K + ArCl)",
    "[O-][n+]1ccccc1Br.c1ccc(B(O)O)cc1>>[O-][n+]1ccccc1-c1ccccc1 (Suzuki - Pyridine N-oxide)",
    "BrC1(c2ccccc2)CC1.c1ccc(B(O)O)cc1>>c1ccc(C2(c3ccccc3)CC2)cc1 (Suzuki - Cyclopropyl bromide)",
    
    # Stille Coupling
    "Brc1ccccc1.c1ccc([Sn](C)(C)C)cc1>>c1ccc(-c2ccccc2)cc1 (Stille - Ph-Ph)",
    "Ic1ccncc1.C=C[Sn](C)(C)C>>C=Cc1ccncc1 (Stille - Vinyl pyridine)",
    
    # Sonogashira Coupling
    "Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1 (Sonogashira - Diphenylacetylene)",
    "Ic1ccncc1.C#CC>>C#Cc1ccncc1 (Sonogashira - Pyridine acetylene)",
    "Clc1ccc(C#N)cc1.C#CC(C)(C)C>>CC(C)(C)C#Cc1ccc(C#N)cc1 (Sonogashira - tert-butyl)",
    "Brc1ccc(OC)cc1.C#Cc1ccccc1>>COc1ccc(C#Cc2ccccc2)cc1 (Sonogashira - Electron-rich ArBr)",
    "Ic1ccc(C(F)(F)F)cc1.C#C[Si](C)(C)C>>CC(C)(C)[Si]C#Cc1ccc(C(F)(F)F)cc1 (Sonogashira - TMS-acetylene)",
    "Brc1cccs1.C#Cc1ccccc1>>c1cc(C#Cc2ccccc2)sc1 (Sonogashira - Thiophene coupling)",
    "Clc1ccncc1.C#CCOC>>COC(C)C#Cc1ccncc1 (Sonogashira - Ether-substituted alkyne)",
    "Brc1ccc2[nH]ccc2c1.C#CC>>C#Cc1ccc2[nH]ccc2c1 (Sonogashira - Indole substrate)",
    
    # Heck Reaction
    "Brc1ccccc1.C=C>>C(=Cc1ccccc1) (Heck - Simple styrene)",
    "Ic1ccc(C=O)cc1.C=CC(=O)OCC>>CCOC(=O)/C=C/c1ccc(C=O)cc1 (Heck - Acrylate)",
    "Brc1ccncc1.C=CC#N>>N#C/C=C/c1ccncc1 (Heck - Acrylonitrile)",
    "Brc1ccc(OC)cc1.C=CC(=O)OC>>COc1ccc(/C=C/C(=O)OC)cc1 (Heck - Methyl acrylate)",
    "Ic1ccc2ccccc2c1.C=C>>C(=Cc1ccc2ccccc2c1) (Heck - Naphthyl iodide)",
    "Brc1cccs1.C=CC(=O)N>>NC(=O)/C=C/c1cccs1 (Heck - Thiophene + acrylamide)",
    "Clc1ccc(C#N)cc1.C=C>>C=Cc1ccc(C#N)cc1 (Heck - Electron-poor chloride)",
    
    # Negishi Coupling
    "Brc1ccccc1.c1ccc([Zn]Cl)cc1>>c1ccc(-c2ccccc2)cc1 (Negishi - Ph-Ph)",
    "Ic1ccncc1.C=C[Zn]Br>>C=Cc1ccncc1 (Negishi - Vinyl pyridine)",
    "Brc1ccc(C(F)(F)F)cc1.c1ccc([Zn]Cl)cc1>>FC(F)(F)c1ccc(-c2ccccc2)cc1 (Negishi - CF3 substrate)",
    "Clc1ccc(C#N)cc1.CC[Zn]Br>>CCc1ccc(C#N)cc1 (Negishi - Alkyl zinc + ArCl)",
    
    # Kumada Coupling
    "Brc1ccccc1.c1ccc([Mg]Br)cc1>>c1ccc(-c2ccccc2)cc1 (Kumada - Ph-Ph)",
    "Ic1ccncc1.C=C[Mg]Cl>>C=Cc1ccncc1 (Kumada - Vinyl pyridine)",
    "Brc1ccc(OC)cc1.CC[Mg]Br>>CCc1ccc(OC)cc1 (Kumada - Alkyl Grignard)",
    "Clc1ccc(C#N)cc1.c1ccc([Mg]Br)cc1>>N#Cc1ccc(-c2ccccc2)cc1 (Kumada - Activated ArCl)",
    
    # ═══════════════════════════════════════════════════════════
    # C-N COUPLING REACTIONS (Comprehensive Test Set)
    # Both Ullmann (Cu-catalyzed) and Buchwald-Hartwig (Pd-catalyzed)
    # ═══════════════════════════════════════════════════════════
    
    # Classic Buchwald-Hartwig Examples (Pd-catalyzed)
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (Buchwald-Hartwig - Diphenylamine)",
    "Clc1ccncc1.NCC>>CCNc1ccncc1 (Buchwald-Hartwig - Pyridine ethylamine)",
    "Brc1ccc(C(F)(F)F)cc1.NC1CCCCC1>>FC(F)(F)c1ccc(NC2CCCCC2)cc1 (B-H - Cyclohexylamine)",
    "Ic1ccc(C=O)cc1.Nc1ccccc1>>O=Cc1ccc(Nc2ccccc2)cc1 (B-H - Aldehyde substrate)",
    "Brc1ccc2ccccc2c1.NCC>>CCNc1ccc2ccccc2c1 (B-H - Naphthylamine)",
    "Clc1nc2ccccc2[nH]1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2[nH]1 (B-H - Benzimidazole)",

    # ═══════════════════════════════════════════════════════════
    # COMPREHENSIVE C-N COUPLING TEST SET (Ullmann & Buchwald)
    # 40+ diverse examples covering substrate scope and selectivity
    # ═══════════════════════════════════════════════════════════
    
    # Primary anilines with various aryl halides
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (C-N - Ph-Br + aniline → diphenylamine)",
    "Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (C-N - Ph-Cl + aniline → diphenylamine)",
    "Ic1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (C-N - Ph-I + aniline → diphenylamine)",
    "Brc1ccc(OC)cc1.Nc1ccccc1>>COc1ccc(Nc2ccccc2)cc1 (C-N - 4-MeO-Ph-Br + aniline)",
    "Brc1ccc(C)cc1.Nc1ccccc1>>Cc1ccc(Nc2ccccc2)cc1 (C-N - 4-Me-Ph-Br + aniline)",
    "Brc1ccc(C(F)(F)F)cc1.Nc1ccccc1>>FC(F)(F)c1ccc(Nc2ccccc2)cc1 (C-N - 4-CF3-Ph-Br + aniline)",
    "Brc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1 (C-N - 4-CN-Ph-Br + aniline)",
    "Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>[O-][N+](=O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-NO2-Ph-Br + aniline)",
    "Brc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC(=O)c1ccc(Nc2ccccc2)cc1 (C-N - 4-acetyl-Ph-Br + aniline)",
    "Brc1ccc(F)cc1.Nc1ccccc1>>Fc1ccc(Nc2ccccc2)cc1 (C-N - 4-F-Ph-Br + aniline)",
    
    # Heteroaryl halides with aniline
    "Clc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1 (C-N - 4-Cl-pyridine + aniline)",
    "Brc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1 (C-N - 4-Br-pyridine + aniline)",
    "Clc1cccnc1.Nc1ccccc1>>c1ccccc1Nc1cccnc1 (C-N - 3-Cl-pyridine + aniline)",
    "Brc1cnccn1.Nc1ccccc1>>c1ccccc1Nc1cnccn1 (C-N - 2-Br-pyrimidine + aniline)",
    "Clc1cccc2ncccc12.Nc1ccccc1>>c1ccccc1Nc1cccc2ncccc12 (C-N - 4-Cl-quinoline + aniline)",
    "Brc1ccc2[nH]ccc2c1.Nc1ccccc1>>c1ccccc1Nc1ccc2[nH]ccc2c1 (C-N - 5-Br-indole + aniline)",
    "Brc1ccco1.Nc1ccccc1>>c1ccccc1Nc1ccco1 (C-N - 3-Br-furan + aniline)",
    "Clc1cccs1.Nc1ccccc1>>c1ccccc1Nc1cccs1 (C-N - 2-Cl-thiophene + aniline)",
    "Brc1cnc(C)cn1.Nc1ccccc1>>Cc1cnc(Nc2ccccc2)cn1 (C-N - 2-Br-4-methylpyrimidine + aniline)",
    
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
    "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1 (C-N - Ph-Br + benzylamine)",
    "Brc1ccccc1.NCCOC>>COCCNc1ccccc1 (C-N - Ph-Br + 2-methoxyethylamine)",
    "Brc1ccncc1.NCCC>>CCCNc1ccncc1 (C-N - 4-Br-pyridine + propylamine)",
    "Clc1ccc(C#N)cc1.CC(C)N>>CC(C)Nc1ccc(C#N)cc1 (C-N - 4-Cl-benzonitrile + isopropylamine)",
    
    # Secondary aliphatic amines
    "Brc1ccccc1.CCN(CC)CC>>CCN(CC)c1ccccc1 (C-N - Ph-Br + diethylamine)",
    
    # Cyclic amines (heterocycles)
    "Brc1ccccc1.N1CCCC1>>c1ccccc1N1CCCC1 (C-N - Ph-Br + pyrrolidine)",
    "Brc1ccccc1.N1CCCCC1>>c1ccccc1N1CCCCC1 (C-N - Ph-Br + piperidine)",
    "Brc1ccccc1.N1CCOCC1>>c1ccccc1N1CCOCC1 (C-N - Ph-Br + morpholine)",
    "Brc1ccccc1.N1CCN(C)CC1>>CN1CCN(c2ccccc2)CC1 (C-N - Ph-Br + N-methylpiperazine)",
    "Brc1ccccc1.N1CCC(O)CC1>>OC1CCN(c2ccccc2)CC1 (C-N - Ph-Br + 4-hydroxypiperidine)",
    "Brc1ccccc1.N1Cc2ccccc2C1>>c1ccccc1N1Cc2ccccc2C1 (C-N - Ph-Br + tetrahydroisoquinoline)",
    "Brc1ccncc1.N1CCCC1>>c1ccncc1N1CCCC1 (C-N - 4-Br-pyridine + pyrrolidine)",
    "Clc1ccc(C(F)(F)F)cc1.N1CCOCC1>>FC(F)(F)c1ccc(N2CCOCC2)cc1 (C-N - 4-Cl-benzotrifluoride + morpholine)",
    "Brc1ccc(OC)cc1.N1CCCCC1>>COc1ccc(N2CCCCC2)cc1 (C-N - 4-Br-anisole + piperidine)",
    
    # Dihalide substrates (potential bis-amination)
    "Brc1ccc(Br)cc1.Nc1ccccc1>>c1ccc(Nc2ccc(Nc3ccccc3)cc2)cc1 (C-N - 1,4-dibromobenzene + aniline)",
    
    # Ortho-substituted challenging substrates
    "Brc1ccccc1C.Nc1ccccc1>>Cc1ccccc1Nc1ccccc1 (C-N - 2-methylbromobenzene + aniline)",
    "Brc1ccccc1OC.Nc1ccccc1>>COc1ccccc1Nc1ccccc1 (C-N - 2-bromoanisole + aniline)",
    "Brc1ccccc1C(=O)C.Nc1ccccc1>>CC(=O)c1ccccc1Nc1ccccc1 (C-N - 2-bromoacetophenone + aniline)",
    
    # Pharmaceutical-relevant scaffolds
    "Clc1ccc2nc(Cl)ccc2c1.Nc1ccccc1>>c1ccc(Nc2ccc3cc(Nc4ccccc4)ccc3n2)cc1 (C-N - 4,7-dichloroquinoline + aniline)",
    "Brc1ccc2c(c1)OCO2.Nc1ccccc1>>c1ccc(Nc2ccc3c(c2)OCO3)cc1 (C-N - 5-bromobenzo[d][1,3]dioxole + aniline)",
    "Clc1nc2ccccc2s1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2s1 (C-N - 2-chlorobenzothiazole + aniline)",
    "Brc1nc2ccccc2o1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2o1 (C-N - 2-bromobenzoxazole + aniline)",
    "Clc1cnc2ccccc2n1.Nc1ccccc1>>c1ccccc1Nc1cnc2ccccc2n1 (C-N - 2-chloroquinoxaline + aniline)",
    
    # Aromatic N-H nucleophiles (Indoles, Pyrazoles, Carbazoles, etc.)
    "Brc1ccccc1.c1cc[nH]c1>>c1ccccc1n1cccc1 (C-N - Ph-Br + pyrrole → N-phenylpyrrole)",
    "Brc1ccc(C)cc1.c1cc[nH]c1>>Cc1ccc(n2cccc2)cc1 (C-N - 4-bromotoluene + pyrrole)",
    "Clc1ccc(C#N)cc1.c1cc[nH]c1>>N#Cc1ccc(n2cccc2)cc1 (C-N - 4-chlorobenzonitrile + pyrrole)",
    "Brc1ccccc1.c1ccc2[nH]ccc2c1>>c1ccc2c(c1)ccn2c1ccccc1 (C-N - Ph-Br + indole → N-phenylindole)",
    "Brc1ccc(OC)cc1.c1ccc2[nH]ccc2c1>>COc1ccc(n2ccc3ccccc32)cc1 (C-N - 4-bromoanisole + indole)",
    "Clc1ccncc1.c1ccc2[nH]ccc2c1>>c1ccncc1n1ccc2ccccc21 (C-N - 4-chloropyridine + indole)",
    "Brc1ccc(C(F)(F)F)cc1.c1ccc2[nH]ccc2c1>>FC(F)(F)c1ccc(n2ccc3ccccc32)cc1 (C-N - 4-bromobenzotrifluoride + indole)",
    "Ic1ccncc1.c1ccc2[nH]ccc2c1>>c1ccncc1n1ccc2ccccc21 (C-N - 4-iodopyridine + indole)",
    "Brc1ccccc1.c1c[nH]nc1>>c1ccccc1n1ccnc1 (C-N - Ph-Br + pyrazole → N-phenylpyrazole)",
    "Brc1ccc(Cl)cc1.c1c[nH]nc1>>Clc1ccc(n2ccnc2)cc1 (C-N - 4-bromochlorobenzene + pyrazole)",
    "Brc1ccncc1.c1c[nH]nc1>>c1ccncc1n1ccnc1 (C-N - 4-bromopyridine + pyrazole)",
    "Brc1ccc([N+](=O)[O-])cc1.c1c[nH]nc1>>[O-][N+](=O)c1ccc(n2ccnc2)cc1 (C-N - 4-bromonitrobenzene + pyrazole)",
    "Brc1ccccc1.c1c[nH]nc1C>>Cc1nn(c2ccccc2)cc1 (C-N - Ph-Br + 3-methylpyrazole)",
    "Clc1ccc(OC)cc1.c1c[nH]nc1C>>Cc1nn(c2ccc(OC)cc2)cc1 (C-N - 4-chloroanisole + 3-methylpyrazole)",
    "Brc1ccccc1.c1cnc[nH]1>>c1ccccc1n1cncc1 (C-N - Ph-Br + imidazole → N-phenylimidazole)",
    "Brc1ccc(F)cc1.c1cnc[nH]1>>Fc1ccc(n2cncc2)cc1 (C-N - 4-bromofluorobenzene + imidazole)",
    "Clc1ccncc1.c1cnc[nH]1>>c1ccncc1n1cncc1 (C-N - 4-chloropyridine + imidazole)",
    "Brc1ccccc1.c1ccc2c(c1)c1ccccc1[nH]2>>c1ccc2c(c1)c1ccccc1n2c1ccccc1 (C-N - Ph-Br + carbazole → N-phenylcarbazole)",
    "Brc1ccc(C)cc1.c1ccc2c(c1)c1ccccc1[nH]2>>Cc1ccc(n2c3ccccc3c3ccccc32)cc1 (C-N - 4-bromotoluene + carbazole)",
    "Ic1ccc(OC)cc1.c1ccc2c(c1)c1ccccc1[nH]2>>COc1ccc(n2c3ccccc3c3ccccc32)cc1 (C-N - 4-iodoanisole + carbazole)",
    "Brc1ccccc1.c1nc[nH]n1>>c1ccccc1n1ncnc1 (C-N - Ph-Br + 1,2,4-triazole → N-phenyltriazole)",
    "Brc1ccc(C#N)cc1.c1nc[nH]n1>>N#Cc1ccc(n2ncnc2)cc1 (C-N - 4-bromobenzonitrile + 1,2,4-triazole)",
    "Brc1ccccc1.c1ncncc1N>>Nc1cnccn1c1ccccc1 (C-N - Ph-Br + 2-aminopyrimidine)",
    "Clc1ccncc1.c1ncncc1N>>Nc1cnccn1c1ccncc1 (C-N - 4-chloropyridine + 2-aminopyrimidine)",
    "Brc1ccc(OC)cc1.c1cn[nH]c1>>COc1ccc(n2nccn2)cc1 (C-N - 4-bromoanisole + 1,2,4-triazole)",
    "Brc1ccccc1C.c1ccc2[nH]ccc2c1>>Cc1ccccc1n1ccc2ccccc21 (C-N - 2-bromotoluene + indole, ortho-substituted)",
    "Brc1c(C)cccc1C.c1ccc2[nH]ccc2c1>>Cc1cccc(C)c1n1ccc2ccccc21 (C-N - 2,6-dimethylbromobenzene + indole, sterically hindered)",
    "Brc1ccccc1OC.c1cc[nH]c1>>COc1ccccc1n1cccc1 (C-N - 2-bromoanisole + pyrrole, ortho-substituted)",
    "Ic1ccc2ccccc2c1.c1ccc2[nH]ccc2c1>>c1ccc(n2ccc3ccccc32)c2ccccc12 (C-N - 2-bromonaphthalene + indole, fused ring)",
    "Brc1ccc(Br)cc1.c1ccc2[nH]ccc2c1>>c1ccc(n2ccc3ccccc32)cc1n1ccc2ccccc21 (C-N - 1,4-dibromobenzene + indole, bis-coupling)",
    "Brc1cccc(Br)c1.c1cc[nH]c1>>c1cccc(n2cccc2)c1n1cccc1 (C-N - 1,3-dibromobenzene + pyrrole, bis-coupling)",
    "Brc1cnc2ccccc2n1.c1ccc2[nH]ccc2c1>>c1ccc2c(c1)ccn2c1cnc2ccccc2n1 (C-N - 3-bromoquinoxaline + indole)",
    "Clc1cccc2c1cc[nH]2.c1ccccc1Br>>c1ccccc1n1ccc2c(Cl)cccc21 (C-N - 4-chloroindole + bromobenzene)",
    
    # Amides as nucleophiles (Primary and Secondary)
    "Brc1ccccc1.NC(=O)C>>CC(=O)Nc1ccccc1 (C-N - Ph-Br + acetamide → N-phenylacetamide)",
    "Brc1ccc(C)cc1.NC(=O)C>>CC(=O)Nc1ccc(C)cc1 (C-N - 4-bromotoluene + acetamide)",
    "Clc1ccc(OC)cc1.NC(=O)C>>CC(=O)Nc1ccc(OC)cc1 (C-N - 4-chloroanisole + acetamide)",
    "Brc1ccncc1.NC(=O)C>>CC(=O)Nc1ccncc1 (C-N - 4-bromopyridine + acetamide)",
    "Ic1ccc(C(F)(F)F)cc1.NC(=O)C>>CC(=O)Nc1ccc(C(F)(F)F)cc1 (C-N - 4-iodobenzotrifluoride + acetamide)",
    "Brc1ccccc1.NC(=O)CC>>CCC(=O)Nc1ccccc1 (C-N - Ph-Br + propionamide)",
    "Brc1ccc(Cl)cc1.NC(=O)CC>>CCC(=O)Nc1ccc(Cl)cc1 (C-N - 4-bromochlorobenzene + propionamide)",
    "Brc1ccccc1.NC(=O)c1ccccc1>>O=C(Nc1ccccc1)c1ccccc1 (C-N - Ph-Br + benzamide → N-phenylbenzamide)",
    "Brc1ccc(OC)cc1.NC(=O)c1ccccc1>>COc1ccc(NC(=O)c2ccccc2)cc1 (C-N - 4-bromoanisole + benzamide)",
    "Clc1ccncc1.NC(=O)c1ccccc1>>O=C(Nc1ccncc1)c1ccccc1 (C-N - 4-chloropyridine + benzamide)",
    "Brc1ccc(C#N)cc1.NC(=O)c1ccccc1>>N#Cc1ccc(NC(=O)c2ccccc2)cc1 (C-N - 4-bromobenzonitrile + benzamide)",
    "Brc1ccccc1.NC(=O)c1ccc(OC)cc1>>COc1ccc(C(=O)Nc2ccccc2)cc1 (C-N - Ph-Br + 4-methoxybenzamide)",
    "Brc1ccc(F)cc1.NC(=O)c1ccc(F)cc1>>Fc1ccc(C(=O)Nc2ccc(F)cc2)cc1 (C-N - 4-bromofluorobenzene + 4-fluorobenzamide)",
    "Brc1ccccc1.NC(=O)CC(C)C>>CC(C)CC(=O)Nc1ccccc1 (C-N - Ph-Br + isobutyramide)",
    "Brc1ccccc1.NC(=O)C(C)(C)C>>CC(C)(C)C(=O)Nc1ccccc1 (C-N - Ph-Br + pivalamide)",
    "Clc1ccc(C(F)(F)F)cc1.NC(=O)CC(C)C>>CC(C)CC(=O)Nc1ccc(C(F)(F)F)cc1 (C-N - 4-chlorobenzotrifluoride + isobutyramide)",
    "Brc1ccccc1.NC(=O)CCc1ccccc1>>O=C(CCc1ccccc1)Nc1ccccc1 (C-N - Ph-Br + hydrocinnamamide)",
    "Brc1ccc(OC)cc1.NC(=O)CCc1ccccc1>>COc1ccc(NC(=O)CCc2ccccc2)cc1 (C-N - 4-bromoanisole + hydrocinnamamide)",
    "Brc1ccccc1.C1CCNC1=O>>O=C1CCCN1c1ccccc1 (C-N - Ph-Br + 2-pyrrolidinone → N-phenyl-2-pyrrolidinone)",
    "Brc1ccc(C)cc1.C1CCNC1=O>>Cc1ccc(N2CCCC2=O)cc1 (C-N - 4-bromotoluene + 2-pyrrolidinone)",
    "Clc1ccncc1.C1CCNC1=O>>O=C1CCCN1c1ccncc1 (C-N - 4-chloropyridine + 2-pyrrolidinone)",
    "Brc1ccccc1.C1CCCCNC1=O>>O=C1CCCCN1c1ccccc1 (C-N - Ph-Br + caprolactam → N-phenylcaprolactam)",
    "Brc1ccc(OC)cc1.C1CCCCNC1=O>>COc1ccc(N2CCCCC2=O)cc1 (C-N - 4-bromoanisole + caprolactam)",
    "Ic1ccc(C(F)(F)F)cc1.C1CCCCNC1=O>>FC(F)(F)c1ccc(N2CCCCC2=O)cc1 (C-N - 4-iodobenzotrifluoride + caprolactam)",
    "Brc1ccccc1.NC(C)=O>>CC(=O)Nc1ccccc1 (C-N - Ph-Br + acetamide, alternative)",
    "Brc1ccccc1.NC(=O)c1ccncc1>>O=C(Nc1ccccc1)c1ccncc1 (C-N - Ph-Br + isonicotinamide)",
    "Clc1ccc(C#N)cc1.NC(=O)c1ccncc1>>N#Cc1ccc(NC(=O)c2ccncc2)cc1 (C-N - 4-chlorobenzonitrile + isonicotinamide)",
    "Brc1ccccc1.NC(=O)c1cccs1>>O=C(Nc1ccccc1)c1cccs1 (C-N - Ph-Br + thiophene-3-carboxamide)",
    "Brc1ccc(OC)cc1.NC(=O)c1cccs1>>COc1ccc(NC(=O)c2cccs2)cc1 (C-N - 4-bromoanisole + thiophene-3-carboxamide)",
    "Brc1ccccc1.NC(=O)c1ccco1>>O=C(Nc1ccccc1)c1ccco1 (C-N - Ph-Br + furan-3-carboxamide)",
    "Brc1ccccc1C.NC(=O)C>>CC(=O)Nc1ccccc1C (C-N - 2-bromotoluene + acetamide, ortho-substituted)",
    "Brc1c(C)cccc1C.NC(=O)c1ccccc1>>Cc1cccc(C)c1NC(=O)c1ccccc1 (C-N - 2,6-dimethylbromobenzene + benzamide, sterically hindered)",
    "Ic1ccc2ccccc2c1.NC(=O)C>>CC(=O)Nc1ccc2ccccc2c1 (C-N - 2-iodonaphthalene + acetamide, fused ring)",
    "Brc1ccc(Br)cc1.NC(=O)C>>CC(=O)Nc1ccc(NC(C)=O)cc1 (C-N - 1,4-dibromobenzene + acetamide, bis-coupling)",
    
    # Chan-Lam Coupling (C-N)
    "c1ccccc1B(O)O.Nc1ccccc1.[O]>>c1ccccc1Nc1ccccc1 (Chan-Lam - Oxidative)",
    "c1ccc(OC)cc1B(O)O.Nc1ccccc1.[O]>>COc1ccc(Nc2ccccc2)cc1 (Chan-Lam - Methoxy boronic acid)",
    "c1ccncc1B(O)O.Nc1ccccc1.[O]>>c1ccccc1Nc1ccncc1 (Chan-Lam - Pyridine boronic acid)",
    
    # ═══════════════════════════════════════════════════════════
    # C-O COUPLING REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Ullmann Ether Synthesis
    "Brc1ccccc1.Oc1ccccc1>>c1ccccc1Oc1ccccc1 (Ullmann Ether - Diphenyl ether)",
    "Ic1ccncc1.OCC>>CCOc1ccncc1 (Ullmann Ether - Ethyl pyridyl ether)",
    "Brc1ccc(C(F)(F)F)cc1.Oc1ccccc1>>FC(F)(F)c1ccc(Oc2ccccc2)cc1 (Ullmann Ether - CF3 substrate)",
    "Clc1ccc(C#N)cc1.OC>>COc1ccc(C#N)cc1 (Ullmann Ether - Methoxy benzonitrile)",
    "Brc1cccs1.Oc1ccccc1>>c1cc(Oc2ccccc2)sc1 (Ullmann Ether - Thiophene ether)",
    
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
    "C=Cc1ccccc1.[H][H]>>CCc1ccccc1 (Hydrogenation - Ethylbenzene)",
    "c1ccc(C=O)cc1.[H][H]>>c1ccc(CO)cc1 (Hydrogenation - Benzyl alcohol)",
    "C#Cc1ccccc1.[H][H]>>CCc1ccccc1 (Hydrogenation - Complete alkyne)",
    "c1ccc([N+](=O)[O-])cc1.[H][H]>>c1ccc(N)cc1 (Hydrogenation - Nitro to amine)",
    "C=CC=C.[H][H]>>CCCC (Hydrogenation - 1,3-butadiene)",
    "c1ccc(C#N)cc1.[H][H]>>c1ccc(CN)cc1 (Hydrogenation - Nitrile to benzylamine)",
    "O=C1CCCCC1.[H][H]>>OC1CCCCC1 (Hydrogenation - Cyclohexanone to cyclohexanol)",
    "C1=CC=C(C=C1)C=CC(=O)O.[H][H]>>c1ccc(CCC(=O)O)cc1 (Hydrogenation - Cinnamic acid)",
    
    # Metal Hydride Reductions
    "c1ccc(C=O)cc1.[BH4-].[Na+]>>c1ccc(CO)cc1 (NaBH4 - Benzaldehyde reduction)",
    "CC(=O)c1ccccc1.[BH4-].[Na+]>>CC(O)c1ccccc1 (NaBH4 - Acetophenone)",
    "c1ccc(C(=O)C)cc1.[AlH4-].[Li+]>>c1ccc(C(O)C)cc1 (LiAlH4 - Ketone reduction)",
    "CC(=O)OCC.[AlH4-].[Li+]>>CCO (LiAlH4 - Ester to alcohol)",
    "c1ccc(C#N)cc1.[AlH4-].[Li+]>>c1ccc(CN)cc1 (LiAlH4 - Nitrile reduction)",
    "c1ccc(C(=O)Nc2ccccc2)cc1.[AlH4-].[Li+]>>c1ccccc1CNc1ccccc1 (LiAlH4 - Amide to amine)",
    "O=C1CCC(=O)O1.[AlH4-].[Li+]>>OCCCO (LiAlH4 - Lactone to diol)",
    "CC(=O)OC(=O)C.[BH4-].[Na+]>>CC(=O)CO (NaBH4 - Selective reduction)",
    
    # Transfer Hydrogenation
    "CC(=O)c1ccccc1.C1CCNC1>>[Ir]>>CC(O)c1ccccc1 (Transfer Hydrogenation - Acetophenone)",
    "c1ccc(C=O)cc1.CCOC(=O)C=O>>[Ru]>>c1ccc(CO)cc1 (Transfer Hydrogenation - Benzaldehyde)",
    
    # Dissolving Metal Reductions
    "c1ccc([N+](=O)[O-])cc1.[Zn].[HCl]>>c1ccc(N)cc1 (Zn/HCl - Nitro reduction)",
    "C1=CC=CC=C1>>[Na].[NH3]>>C1CCC=CC1 (Birch Reduction - Benzene to cyclohexadiene)",
    "c1ccc(OC)cc1>>[Na].[NH3]>>COC1C=CCC=C1 (Birch Reduction - Anisole)",
    
    # Selective Reductions
    "c1ccc(C(=O)C=C)cc1.[NaBH(OAc)3]>>c1ccc(C(O)C=C)cc1 (NaBH(OAc)3 - Selective ketone)",
    "CC(=O)CC=C.[DIBAL]>>CC(O)CC=C (DIBAL - Ketone selective)",
    "CC(=O)OCC.[DIBAL]>>CC=O (DIBAL - Ester to aldehyde)",
    
    # OXIDATION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    "c1ccc(CO)cc1.[O]>>c1ccc(C=O)cc1 (Oxidation - Benzyl alcohol to aldehyde)",
    "CC(O)c1ccccc1.[O]>>CC(=O)c1ccccc1 (Oxidation - Secondary alcohol to ketone)",
    "CCO.[O]>>CC=O (Oxidation - Ethanol to acetaldehyde)",
    "c1ccc(CO)cc1.Cl[Cr](=O)(=O)OC(C)(C)C>>c1ccc(C=O)cc1 (PCC Oxidation)",
    "CCO.OS(=O)(=O)O>>CC=O (Swern-type oxidation)",
    "c1ccc(CO)cc1.[MnO2]>>c1ccc(C=O)cc1 (MnO2 - Benzyl alcohol oxidation)",
    "CC(O)C.[CrO3]>>CC(=O)C (Jones Oxidation - Acetone)",
    "c1ccc(C(O)C)cc1.[Dess-Martin]>>c1ccc(C(=O)C)cc1 (Dess-Martin - Acetophenone)",
    "CC(O)CCO.[IBX]>>CC(=O)CC=O (IBX - Diol to diketone)",
    "c1ccc(CO)cc1.[TEMPO].[NaOCl]>>c1ccc(C(=O)O)cc1 (TEMPO - Alcohol to acid)",
    "C1CCC(CO)CC1.[PCC]>>C1CCC(C=O)CC1 (PCC - Cyclohexylmethanol)",
    
    # Epoxidations
    "C=Cc1ccccc1.[O]>>C1OC1c1ccccc1 (mCPBA - Styrene epoxidation)",
    "C=C(C)C.[O]>>CC1(C)CO1 (mCPBA - Isobutylene oxide)",
    "C=CCCCCC=C.[O]>>C1OC1CCCCC1CO1 (mCPBA - Diepoxide)",
    
    # Baeyer-Villiger Oxidation
    "CC(=O)c1ccccc1.[O]>>CC(=O)Oc1ccccc1 (Baeyer-Villiger - Phenyl acetate)",
    "O=C1CCCCC1.[O]>>O=C1OCCCC1 (Baeyer-Villiger - Caprolactone)",
    
    # Sulfide Oxidations
    "CSc1ccccc1.[O]>>CS(=O)c1ccccc1 (Sulfide to sulfoxide)",
    "CSc1ccccc1.[mCPBA].[O]>>CS(=O)(=O)c1ccccc1 (Sulfide to sulfone)",
    
    # ═══════════════════════════════════════════════════════════
    # CYCLOADDITION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # Diels-Alder
    "C=CC=C.C=C>>C1C=CCC1 (Diels-Alder - Simple)",
    "C=CC(=C)C.C=CC(=O)OCC>>CC1=CC(C(=O)OCC)CC(C)=C1 (Diels-Alder - Substituted)",
    "c1ccc2ccccc2c1C=CC=C.C=CC(=O)C>>O=C1CCc2c(ccc3ccccc23)C1 (Diels-Alder - Naphthyl)",
    
    # 1,3-Dipolar Cycloaddition (Click Chemistry)
    "C#CCO.c1ccc([N-][N+]#N)cc1>>c1ccccc1n1cc(CO)nn1 (CuAAC Click - Triazole)",
    "C#Cc1ccccc1.C[N-][N+]#N>>Cc1nnn(-c2ccccc2)c1 (Click - Phenyl triazole)",
    "C#CCCOC.CC[N-][N+]#N>>CCn1ncc(CCOC)n1 (Click - Ether-substituted triazole)",
    "C#CC(C)(C)O.c1ccc(C[N-][N+]#N)cc1>>CC(C)(O)Cn1cc(Cc2ccccc2)nn1 (Click - Benzyl azide + propargyl alcohol)",
    
    # ═══════════════════════════════════════════════════════════
    # SUBSTITUTION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # SN2 Reactions
    "CCCCBr.[OH-]>>CCCCO (SN2 - Primary alcohol)",
    "CC(Br)C.[OH-]>>CC(O)C (SN2 - isopropyl alcohol)",
    "CCBr.N#C>>CCC#N (SN2 - Nitrile formation)",
    "CCCCBr.Oc1ccccc1>>CCCCOc1ccccc1 (SN2 - Phenoxide)",
    "CI.SC>>CSC (SN2 - Thioether formation)",
    "CCCBr.[I-]>>CCCI (SN2 - Finkelstein halide exchange)",
    "BrCCCCCl.[C-]#N>>N#CCCCCCN (SN2 - Dinitrile)",
    "ICc1ccccc1.[OH-]>>OCc1ccccc1 (SN2 - Benzyl alcohol)",
    "BrCC(=O)OC.[N-]=[N+]=[N-]>>[N-]=[N+]=[N-]CC(=O)OC (SN2 - Azide substitution)",
    
    # C-S COUPLING REACTIONS (Thioether formation via coupling)
    # ═══════════════════════════════════════════════════════════
    
    # Pd-catalyzed arylation of thiols (generic examples)
    "Brc1ccccc1.Sc1ccccc1>>c1ccccc1Sc1ccccc1 (C-S Coupling - Thioether Formation)",
    "Clc1ccc(C)cc1.Sc1ccccc1>>Cc1ccc(Sc2ccccc2)cc1 (C-S Coupling - Thioether Formation)",
    "Brc1ccncc1.SCC>>CCSc1ccncc1 (C-S Coupling - Pyridine thioether)",
    "Ic1ccc(C(F)(F)F)cc1.SCc1ccccc1>>FC(F)(F)c1ccc(SCc2ccccc2)cc1 (C-S Coupling - Benzyl thioether)",
    "Brc1ccc(OC)cc1.SC>>CSc1ccc(OC)cc1 (C-S Coupling - Methyl thioether)",
    
    # SN1 Reactions  
    "CC(C)(C)Br.O>>CC(C)(C)O (SN1 - tert-Butyl)",
    "c1ccc(C(C)(C)Cl)cc1.O>>c1ccc(C(C)(C)O)cc1 (SN1 - Benzyl tertiary)",
    
    # SNAr Reactions (Nucleophilic Aromatic Substitution)
    "Fc1ccc([N+](=O)[O-])cc1.[OH-]>>[O-][N+](=O)c1ccc(O)cc1 (SNAr - Nitrofluorobenzene)",
    "Clc1nc(Cl)nc(Cl)n1.OC>>COc1nc(Cl)nc(Cl)n1 (SNAr - Cyanuric chloride)",
    "Fc1cc([N+](=O)[O-])cc([N+](=O)[O-])c1.NC>>[O-][N+](=O)c1cc([N+](=O)[O-])cc(NC)c1 (SNAr - Dinitrofluorobenzene)",
    "Fc1ccc(C(F)(F)F)cc1.Oc1ccccc1>>FC(F)(F)c1ccc(Oc2ccccc2)cc1 (SNAr - Trifluoromethyl activation)",
    "Clc1ccncc1.OC>>COc1ccncc1 (SNAr - Pyridine)",
    
    # Allylic Substitution
    "C=CCC(Br)C.[OH-]>>C=CCC(C)O (Allylic SN2)",
    "C=CCCC(Cl)C.SC>>C=CCCC(C)SC (Allylic substitution - Thiol)",
    
    # Finkelstein Reaction
    "CCCCBr.[I-]>>CCCCI (Finkelstein - Iodide exchange)",
    "c1ccc(CCl)cc1.[I-]>>c1ccc(CI)cc1 (Finkelstein - Benzyl)",
    
    # Appel Reaction
    "c1ccc(CO)cc1.ClP(c1ccccc1)(c1ccccc1)c1ccccc1.C(Cl)(Cl)(Cl)Cl>>c1ccc(CCl)cc1 (Appel - Alcohol to chloride)",
    "CCCCO.BrP(c1ccccc1)(c1ccccc1)c1ccccc1>>CCCCBr (Appel - Alcohol to bromide)",
    
    # Mitsunobu Reaction
    "c1ccc(CO)cc1.Oc1ccccc1.N=NC(C(=O)OCC)C(=O)OCC.P(c1ccccc1)(c1ccccc1)c1ccccc1>>c1ccc(COc2ccccc2)cc1 (Mitsunobu - Ether)",
    "CCCO.NC(=O)c1ccccc1.N=NC(C(=O)OCC)C(=O)OCC.P(c1ccccc1)(c1ccccc1)c1ccccc1>>CCCNC(=O)c1ccccc1 (Mitsunobu - Amide)",
    
    # ═══════════════════════════════════════════════════════════
    # ELIMINATION REACTIONS
    # ═══════════════════════════════════════════════════════════
    
    # E2 Eliminations
    "CC(Br)C>>CC=C (E2 - Propene from 2-bromopropane)",
    "CC(Br)C(C)C>>CC=C(C)C (E2 - Zaitsev product)",
    "c1ccc(C(C)Br)cc1>>C=Cc1ccccc1 (E2 - Styrene formation)",
    "CCCBr>>CC=C (E2 - Propene)",
    
    # E1 Eliminations
    "CC(C)(C)Br>>CC(C)=C (E1 - tert-Butyl)",
    
    # CARBONYL CHEMISTRY
    # ═══════════════════════════════════════════════════════════
    
    # Aldol Reactions
    "CC=O.CC=O>>CC(O)CC=O (Aldol - Aldol product)",
    "CC(=O)C.CC=O>>CC(O)C(C)C=O (Aldol - Mixed aldol)",
    "c1ccc(C=O)cc1.CC(=O)C>>c1ccc(C(O)C(C)C=O)cc1 (Aldol - Benzaldehyde)",
    "CCC=O.CC(=O)C>>CCC(O)C(C)C=O (Aldol - Propionaldehyde)",
    "c1ccc(C=O)cc1.CC=O>>c1ccc(/C=C/C=O)cc1 (Aldol - Dehydration to cinnamaldehyde)",
    
    # Wittig Reactions
    "CC=O.C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[base]>>C=CC (Wittig - Propene)",
    "c1ccc(C=O)cc1.C=C[P+](c1ccccc1)(c1ccccc1)c1ccccc1>>c1ccc(/C=C/C=C)cc1 (Wittig - Diene)",
    "CC(=O)C.CCOC(=O)C=P(c1ccccc1)(c1ccccc1)c1ccccc1>>CCOC(=O)C=C(C)C (Wittig - α,β-Unsaturated ester)",
    "c1ccc(C=O)cc1.[Ph3P]=CHC>>c1ccc(C=CC)cc1 (Wittig - Propiophenone)",
    
    # Mannich Reaction
    "CC=O.NC.C=O>>CC(=O)CN (Mannich - Simple)",
    "c1ccc(C=O)cc1.C1CCNC1.C=O>>c1ccc(C(=O)CN1CCCC1)cc1 (Mannich - Pyrrolidine)",
    
    # Michael Addition
    "C=CC(=O)C.CC(=O)CC(=O)C>>CC(=O)C(CC(=O)C)CC(=O)C (Michael - 1,5-diketone)",
    "C=CC(=O)OCC.CC(=O)CC(=O)OCC>>CCOC(=O)C(CC(=O)OCC)CC(=O)OCC (Michael - Malonate)",
    "c1ccc(C=CC(=O)C)cc1.CC(=O)CC>>c1ccc(C(CC(C)=O)CC(=O)CC)cc1 (Michael - Chalcone)",
    
    # Claisen Condensation
    "CC(=O)OCC.CC(=O)OCC>>CC(=O)CC(=O)OCC (Claisen - Acetoacetate)",
    "c1ccc(C(=O)OC)cc1.CC(=O)OC>>c1ccc(C(=O)CC(=O)OC)cc1 (Claisen - Benzoylacetate)",
    
    # Knoevenagel Condensation
    "c1ccc(C=O)cc1.CC(C#N)C#N>>c1ccc(/C=C(C#N)C#N)cc1 (Knoevenagel - Malononitrile)",
    "CC=O.CC(C(=O)OCC)C(=O)OCC>>C/C=C(C(=O)OCC)C(=O)OCC (Knoevenagel - Malonate)",
    
    # Henry Reaction
    "c1ccc(C=O)cc1.C[N+](=O)[O-]>>c1ccc(C(O)C[N+](=O)[O-])cc1 (Henry - Nitroaldol)",
    "CC=O.CC[N+](=O)[O-]>>CC(O)CC[N+](=O)[O-] (Henry - Nitropropanol)",
    
    # Reformatsky Reaction
    "c1ccc(C=O)cc1.BrCC(=O)OCC.[Zn]>>c1ccc(C(O)CC(=O)OCC)cc1 (Reformatsky)",
    "CC=O.BrCC(=O)OCC.[Zn]>>CC(O)CC(=O)OCC (Reformatsky - Acetaldehyde)",
    
    # Horner-Wadsworth-Emmons
    "c1ccc(C=O)cc1.CCOP(=O)(OCC)CC(=O)OCC>>c1ccc(/C=C/C(=O)OCC)cc1 (HWE - Cinnamate)",
    "CC=O.CCOP(=O)(OCC)CC(=O)OCC>>C/C=C/C(=O)OCC (HWE - Crotonate)",
    
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
    
    # HETEROCYCLE SYNTHESIS
    # ═══════════════════════════════════════════════════════════
    
    "c1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>c1ccc2[nH]c3ccccc3c2c1 (Carbazole synthesis)",
    "Nc1ccccc1.c1ccc(C=O)cc1>>c1ccc2nc(-c3ccccc3)cc(-c3ccccc3)c2c1 (Quinoline synthesis)",
    
    # Paal-Knorr Synthesis
    "CC(=O)CC(=O)C.NC>>Cc1ccc(C)n1C (Paal-Knorr - Pyrrole)",
    "CC(=O)CCC(=O)C.N>>Cc1ccc(C)[nH]1 (Paal-Knorr - Simple pyrrole)",
    "c1ccc(C(=O)CCC(=O)c2ccccc2)cc1.NC>>c1ccc(-c2ccc(-c3ccccc3)n2C)cc1 (Paal-Knorr - Diphenylpyrrole)",
    
    # Hantzsch Synthesis
    "CC(=O)CC(=O)OCC.CC(=O)CC(=O)OCC.CC=O.N>>CCOC(=O)C1=C(C)NC(C)=C(C(=O)OCC)C1 (Hantzsch - Dihydropyridine)",
    
    # Biginelli Reaction
    "O=C(N)N.CC(=O)CC(=O)OCC.c1ccc(C=O)cc1>>CCOC(=O)C1=C(C)NC(=O)NC1c1ccccc1 (Biginelli - Dihydropyrimidinone)",
    "O=C(N)N.CC(=O)CC(=O)OCC.CC=O>>CCOC(=O)C1=C(C)NC(=O)NC1C (Biginelli - Acetaldehyde)",
    
    # Fischer Indole Synthesis
    "c1ccc(C(C)=O)cc1.c1ccc(NN)cc1>>Cc1c2ccccc2[nH]c1-c1ccccc1 (Fischer Indole)",
    "CC(=O)CC.c1ccccc1NN>>CCc1c2ccccc2[nH]c1C (Fischer Indole - 3-Ethyl)",
    
    # Pictet-Spengler Reaction
    "NCCc1ccc(O)cc1.CC=O>>Cc1nccc2ccc(O)cc12 (Pictet-Spengler - Tetrahydroisoquinoline)",
    "NCCc1c[nH]c2ccccc12.C=O>>c1ccc2[nH]cc3c2c1CCN3 (Pictet-Spengler - β-Carboline)",
    
    # Benzimidazole Synthesis
    "c1ccc(c(N)c1)N.c1ccc(C=O)cc1>>c1ccc(-c2nc3ccccc3[nH]2)cc1 (Benzimidazole from o-phenylenediamine)",
    "c1ccc(c(N)c1)N.CC(=O)O>>Cc1nc2ccccc2[nH]1 (Benzimidazole - Acetic acid)",
    
    # Benzoxazole Synthesis
    "c1ccc(c(N)c1)O.c1ccc(C=O)cc1>>c1ccc(-c2nc3ccccc3o2)cc1 (Benzoxazole)",
    "c1ccc(c(N)c1)O.CC(=O)O>>Cc1nc2ccccc2o1 (Benzoxazole - Acetic acid)",
    
    # Benzothiazole Synthesis
    "c1ccc(c(N)c1)S.c1ccc(C=O)cc1>>c1ccc(-c2nc3ccccc3s2)cc1 (Benzothiazole)",
    
    # Imidazole Synthesis
    "CC(=O)C=O.NC=O.N>>Cc1nccn1C (Imidazole)",
    
    # ═══════════════════════════════════════════════════════════
    # PROTECTING GROUP CHEMISTRY
    # ═══════════════════════════════════════════════════════════
    
    # Boc Protection/Deprotection
    "c1ccc(N)cc1.CC(C)(C)OC(=O)OC(=O)OC(C)(C)C>>CC(C)(C)OC(=O)Nc1ccccc1 (Boc Protection - Aniline)",
    "NCCc1ccccc1.CC(C)(C)OC(=O)OC(=O)OC(C)(C)C>>CC(C)(C)OC(=O)NCCc1ccccc1 (Boc Protection - Phenethylamine)",
    "CC(C)(C)OC(=O)Nc1ccccc1.[HCl]>>c1ccc(N)cc1 (Boc Deprotection - TFA or HCl)",
    
    # Cbz Protection/Deprotection
    "c1ccc(N)cc1.O=C(OCc1ccccc1)OCc1ccccc1>>O=C(Nc1ccccc1)OCc1ccccc1 (Cbz Protection - Aniline)",
    "NCC(=O)O.O=C(OCc1ccccc1)OCc1ccccc1>>O=C(NCC(=O)O)OCc1ccccc1 (Cbz Protection - Glycine)",
    "O=C(Nc1ccccc1)OCc1ccccc1.[H][H].[Pd]>>c1ccc(N)cc1 (Cbz Deprotection - Hydrogenolysis)",
    
    # FMOC Protection/Deprotection
    "c1ccc(N)cc1.O=C1OC(=O)c2ccccc21>>O=C(Nc1ccccc1)OCC1c2ccccc2-c2ccccc21 (FMOC Protection)",
    "O=C(Nc1ccccc1)OCC1c2ccccc2-c2ccccc21.[Base]>>c1ccc(N)cc1 (FMOC Deprotection - Base)",
    
    # TBS/TBDMS Protection/Deprotection (Silyl Ethers)
    "c1ccc(CO)cc1.CC(C)(C)[Si](C)(C)Cl>>CC(C)(C)[Si](C)(C)OCc1ccccc1 (TBS Protection - Benzyl alcohol)",
    "OCCO.CC(C)(C)[Si](C)(C)Cl>>CC(C)(C)[Si](C)(C)OCCO[Si](C)(C)C(C)(C)C (TBS Protection - Diol)",
    "CC(C)(C)[Si](C)(C)OCc1ccccc1.[F-]>>c1ccc(CO)cc1 (TBS Deprotection - TBAF)",
    
    # PMB Protection/Deprotection (p-Methoxybenzyl)
    "c1ccc(CO)cc1.BrCc1ccc(OC)cc1>>COc1ccc(COCc2ccccc2)cc1 (PMB Protection)",
    "COc1ccc(COCc2ccccc2)cc1.[DDQ]>>c1ccc(CO)cc1 (PMB Deprotection - DDQ)",
    
    # Acetonide Protection/Deprotection
    "OCC(O)CO.CC(=O)C>>CC1(C)OCC(CO)O1 (Acetonide Protection - Glycerol)",
    "CC1(C)OCC(CO)O1.[H+]>>OCC(O)CO (Acetonide Deprotection - Acid)",
    
    # Benzyl Ether Protection/Deprotection
    "c1ccc(CO)cc1.BrCc1ccccc1>>c1ccc(COCc2ccccc2)cc1 (Bn Protection)",
    "c1ccc(COCc2ccccc2)cc1.[H][H].[Pd]>>c1ccc(CO)cc1 (Bn Deprotection - Hydrogenolysis)",
    
    # Trityl Protection/Deprotection
    "c1ccc(CO)cc1.ClC(c1ccccc1)(c1ccccc1)c1ccccc1>>c1ccc(COC(c2ccccc2)(c2ccccc2)c2ccccc2)cc1 (Trityl Protection)",
    "c1ccc(COC(c2ccccc2)(c2ccccc2)c2ccccc2)cc1.[H+]>>c1ccc(CO)cc1 (Trityl Deprotection)",
    
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
    "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (B-H: PhBr + aniline → diphenylamine)",
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
    "Brc1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1 (B-H: 4-bromopyridine + aniline → 4-phenylaminopyridine)",
    "Brc1cccnc1.Nc1ccccc1>>c1ccccc1Nc1cccnc1 (B-H: 3-bromopyridine + aniline → 3-phenylaminopyridine)",
    "Brc1ccccn1.Nc1ccccc1>>c1ccc(Nc2ccccn2)cc1 (B-H: 2-bromopyridine + aniline → 2-phenylaminopyridine)",
    "Brc1cncc(C)c1.Nc1ccccc1>>Cc1cncc(Nc2ccccc2)c1 (B-H: 5-bromo-2-methylpyrimidine)",
    "Brc1nc2ccccc2s1.Nc1ccccc1>>c1ccccc1Nc1nc2ccccc2s1 (B-H: 2-bromobenzothiazole)",
    
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
    "Clc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (B-H: PhCl + aniline → diphenylamine, difficult)",
    "Clc1ccc(C#N)cc1.Nc1ccccc1>>N#Cc1ccc(Nc2ccccc2)cc1 (B-H: 4-ClCN + aniline → activated ArCl)",
    "Clc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>O=[N+]([O-])c1ccc(Nc2ccccc2)cc1 (B-H: 4-ClNO2 + aniline)",
    "Clc1ccc(C(=O)C)cc1.Nc1ccccc1>>CC(=O)c1ccc(Nc2ccccc2)cc1 (B-H: 4-ClCOMe + aniline → activated)",
    "Clc1ccncc1.CN(C)C>>CN(C)c1ccncc1 (B-H: 4-chloropyridine + dimethylamine → electron-poor)",
    "Clc1cnc2ccccc2n1.Nc1ccccc1>>c1ccc(Nc2cnc3ccccc3n2)cc1 (B-H: 2-chloroquinoxaline + aniline)",
    
    # ═══════ IODO SUBSTRATES ═══════
    
    # Aryl iodides (most reactive)
    "Ic1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1 (B-H: PhI + aniline → diphenylamine, fast)",
    "Ic1ccncc1.Nc1ccccc1>>c1ccccc1Nc1ccncc1 (B-H: 4-iodopyridine + aniline → heteroaryl)",
    
    # ═══════ BENZYL AMINES ═══════
    
    # Primary benzyl amines
    "Brc1ccccc1.NCc1ccccc1>>c1ccccc1CNc1ccccc1 (B-H: PhBr + benzylamine → N-benzylaniline)",
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
    "Brc1ccccc1.Nc1cccnc1>>c1ccccc1Nc1cccnc1 (B-H: PhBr + 3-aminopyridine → 3-phenylaminopyridine)",
    "Brc1ccccc1.Nc1ccncc1>>c1ccccc1Nc1ccncc1 (B-H: PhBr + 4-aminopyridine → 4-phenylaminopyridine)",
    
    # Thiazole and other heteroaryl amines
    "Brc1ccccc1.Nc1nccs1>>c1ccc(Nc2nccs2)cc1 (B-H: PhBr + 2-aminothiazole)",
    "Brc1ccccc1.Nc1ccco1>>c1ccccc1Nc1ccco1 (B-H: PhBr + 3-aminofuran)",
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
    
    # ═══════════════════════════════════════════════════════════
    # PHASE 3/4 FEATURE TEST REACTIONS - Advanced Substrates
    # ═══════════════════════════════════════════════════════════
    
    # Polyhalogenated substrates
    "Brc1cc(Br)cc(Br)c1.OB(O)c1ccccc1>>Brc1cc(Br)cc(-c2ccccc2)c1 (Suzuki - 1,3,5-Tribromobenzene selective coupling)",
    "Fc1c(F)c(F)c(Br)c(F)c1F.OB(O)c1ccccc1>>Fc1c(F)c(F)c(-c2ccccc2)c(F)c1F (Suzuki - Pentafluorobromobenzene)",
    "Clc1ccc(Cl)nc1.OB(O)c1ccccc1>>Clc1ccc(-c2ccccc2)nc1 (Suzuki - 2,5-Dichloropyridine)",
    "Brc1c(Cl)cccc1Cl.OB(O)c1ccccc1>>c1ccccc1-c1c(Cl)cccc1Cl (Suzuki - Mixed halide with ortho-chlorines)",
    
    # Sterically hindered substrates  
    "Brc1ccccc1C(C)(C)C.OB(O)c1ccccc1>>CC(C)(C)c1ccccc1-c1ccccc1 (Suzuki - Ortho-tert-butyl sterics)",
    "Clc1cc(C(C)C)ccc1C(C)C.OB(O)c1ccccc1>>CC(C)c1cccc(C(C)C)c1-c1ccccc1 (Suzuki - Diisopropyl sterics)",
    "Brc1c(C)c(C)c(C)c(C)c1C.OB(O)c1ccccc1>>Cc1c(C)c(C)c(C)c(C)c1-c1ccccc1 (Suzuki - Pentamethyl extreme sterics)",
    "Ic1c(OC)cccc1OC.OB(O)c1ccccc1>>COc1cccc(OC)c1-c1ccccc1 (Suzuki - Ortho-dimethoxy with EDG)",
    
    # Protected intermediates (BOC, CBZ, FMOC, silyl)
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.OB(O)c1ccccc1>>CC(C)(C)OC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - BOC-protected aniline)",
    "Brc1ccc(NC(=O)OCc2ccccc2)cc1.OB(O)c1ccccc1>>c1ccccc1COC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - CBZ-protected aniline)",
    "Clc1ccc(NC(=O)OCC2c3ccccc3-c3ccccc32)cc1.OB(O)c1ccccc1>>c1ccccc1-c1ccc(NC(=O)OCC2c3ccccc3-c3ccccc32)cc1 (Suzuki - FMOC-protected aniline)",
    "Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1.OB(O)c1ccccc1>>CC(C)(C)[Si](C)(C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TBS-protected phenol)",
    "Ic1ccc(O[Si](C(C)C)(C(C)C)C(C)C)cc1.OB(O)c1ccccc1>>CC(C)[Si](C(C)C)(C(C)C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TIPS-protected phenol)",
    
    # Strong EWG/EDG substrates
    "Brc1ccc([N+](=O)[O-])cc1.OB(O)c1ccccc1>>[O-][N+](=O)c1ccc(-c2ccccc2)cc1 (Suzuki - p-Nitrobromobenzene strong EWG)",
    "Clc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1 (Suzuki - p-Chlorobenzonitrile strong EWG)",
    "Brc1ccc(N(C)C)cc1.OB(O)c1ccccc1>>CN(C)c1ccc(-c2ccccc2)cc1 (Suzuki - p-Dimethylaminobromobenzene strong EDG)",
    "Ic1ccc(OC)cc1.OB(O)c1ccccc1>>COc1ccc(-c2ccccc2)cc1 (Suzuki - p-Iodoanisole strong EDG)",
    
    # Chiral substrates
    "NC(Cc1ccccc1)C(=O)O.Brc1ccccc1>>Brc1ccccc1.NC(Cc1ccccc1)C(=O)O (C-N - Phenylalanine amination, chiral)",
    "C[C@H](N)C(=O)O.Brc1ccccc1>>Brc1ccccc1.C[C@H](N)C(=O)O (C-N - (S)-Alanine amination, chiral center)",
    
    # Fused ring systems
    "Brc1ccc2ccccc2c1.OB(O)c1ccccc1>>c1ccccc1-c1ccc2ccccc2c1 (Suzuki - 2-Bromonaphthalene fused bicyclic)",
    "Clc1ccc2c3ccccc3ccc2c1.OB(O)c1ccccc1>>c1ccccc1-c1ccc2c3ccccc3ccc2c1 (Suzuki - 9-Chloroanthracene fused tricyclic)",
    
    # Buchwald-Hartwig with Phase 3/4 features
    "Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>[O-][N+](=O)c1ccc(Nc2ccccc2)cc1 (C-N - p-Nitrobromobenzene + aniline, strong EWG)",
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.Nc1ccccc1>>CC(C)(C)OC(=O)Nc1ccc(Nc2ccccc2)cc1 (C-N - BOC-protected with amination)",
    "Brc1c(OC)cccc1OC.Nc1ccccc1>>COc1cccc(OC)c1Nc1ccccc1 (C-N - Ortho-dimethoxy sterically hindered)",
    "Brc1ccc2ccccc2c1.NCCC>>CCCNc1ccc2ccccc2c1 (C-N - 2-Bromonaphthalene fused ring + propylamine)",
    
    # Sonogashira with Phase 3/4 features
    "Brc1cc(Br)cc(Br)c1.C#Cc1ccccc1>>Brc1cc(Br)cc(C#Cc2ccccc2)c1 (Sonogashira - Tribromobenzene polyhalogenated)",
    "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1 (Sonogashira - Ortho-tert-butyl sterics)",
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


# ═══════════════════════════════════════════════════════════
# SUZUKI DATABASE VALIDATION TEST SET
# One test reaction for each rule in suzuki_db.json
# ═══════════════════════════════════════════════════════════

SUZUKI_DB_TEST_REACTIONS = {
    # Scheme entries (18 total)
    "SCDB-SUZ-ARBRI-GENERAL-PPh3": {
        "smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "description": "Simple PhBr + PhB(OH)2 → biphenyl (ArBr, boronic acid)",
        "expected_features": {"lg_class": "Br", "boron_type": "boronic_acid"}
    },
    
    "SCDB-SUZ-ARBRI-GENERAL-SPhos": {
        "smiles": "Ic1ccc(C)cc1.OB(O)c1ccccc1>>Cc1ccc(-c2ccccc2)cc1",
        "description": "4-Iodotoluene + PhB(OH)2 (ArI, SPhos set)",
        "expected_features": {"lg_class": "I", "boron_type": "boronic_acid"}
    },
    
    "SCDB-SUZ-BULKY-NUC-XPHOS": {
        "smiles": "Brc1ccccc1.OB(O)c1c(C)cccc1C>>Cc1cccc(C)c1-c1ccccc1",
        "description": "PhBr + 2,6-dimethylphenylboronic acid (hindered nucleophile)",
        "expected_features": {"lg_class": "Br", "nucleophile_ortho_sub": 2}
    },
    
    "SCDB-SUZ-2HETARYL-SPHOS": {
        "smiles": "Brc1ccncc1.OB(O)c1ccccc1>>c1ccc(-c2ccncc2)cc1",
        "description": "4-Bromopyridine + PhB(OH)2 (2-hetaryl halide)",
        "expected_features": {"lg_class": "Br", "ring_hetero_count": 1}
    },
    
    "SCDB-SUZ-ARI-META-TBUXPHOS": {
        "smiles": "Ic1cc(C)cc(C)c1.OB(O)c1ccccc1>>Cc1cc(C)cc(-c2ccccc2)c1",
        "description": "3,5-Dimethyl iodobenzene + PhB(OH)2 (meta-hindered ArI)",
        "expected_features": {"lg_class": "I", "meta_sub_count": 2}
    },
    
    "SCDB-SUZ-ARBR-ORTHO-XPhos": {
        "smiles": "Brc1ccccc1C.OB(O)c1ccccc1>>Cc1ccccc1-c1ccccc1",
        "description": "2-Bromotoluene + PhB(OH)2 (ortho-hindered ArBr)",
        "expected_features": {"lg_class": "Br", "ortho_sub_count": 1}
    },
    
    "SCDB-SUZ-ARCL-EPoor-XPhos": {
        "smiles": "Clc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1",
        "description": "4-Chlorobenzonitrile + PhB(OH)2 (electron-poor ArCl)",
        "expected_features": {"lg_class": "Cl", "ewg_count": 1}
    },
    
    "SCDB-SUZ-ARCL-ERich-L95": {
        "smiles": "Clc1ccc(OC)cc1.OB(O)c1ccccc1>>COc1ccc(-c2ccccc2)cc1",
        "description": "4-Chloroanisole + PhB(OH)2 (electron-rich ArCl)",
        "expected_features": {"lg_class": "Cl", "edg_count": 1}
    },
    
    "SCDB-SUZ-OTf-DPPF": {
        "smiles": "FC(F)(F)S(=O)(=O)Oc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "description": "Phenyl triflate + PhB(OH)2 (aryl triflate)",
        "expected_features": {"lg_class": "OTf"}
    },
    
    "SCDB-SUZ-HET-2PYRIDYL-SLOWRELEASE": {
        "smiles": "Brc1ccccc1.Fc1ccncc1.FB(F)F.[K+]>>c1ccc(-c2ccncc2)cc1",
        "description": "PhBr + 4-pyridyl-BF3K (slow-release 2-pyridyl partner)",
        "expected_features": {"boron_type": "BF3K", "protodeboronation_prone": True}
    },
    
    "SCDB-SUZ-VINYL-DPPF-RT": {
        "smiles": "BrC=C.OB(O)C=C>>C=CC=C",
        "description": "Vinyl bromide + vinylboronic acid (vinyl-vinyl coupling at RT)",
        "expected_features": {"is_vinyl": True, "lg_class": "Br"}
    },
    
    "SCDB-SUZ-ALKYL-PRIMARYI-9BBN": {
        "smiles": "CCCCI.B1CCCCCCCC1>>CCCCC1CCCCCCCC1",
        "description": "Propyl iodide + 9-BBN (primary alkyl sp3)",
        "expected_features": {"is_primary_alkyl": True, "lg_class": "I"}
    },
    
    "SCDB-SUZ-HET-THIOPHENE-FURAN-TRIBORONATE": {
        "smiles": "Brc1ccccc1.OB(O)c1cccs1>>c1ccc(-c2cccs2)cc1",
        "description": "PhBr + 2-thiophenylboronic acid (protodeboronation-prone)",
        "expected_features": {"boron_protodeboronation_prone": True}
    },
    
    "SCDB-SUZ-MIYAURA-BORYLATION": {
        "smiles": "Brc1ccccc1.B2OC(C)(C)C(C)(C)O2.B2OC(C)(C)C(C)(C)O2>>c1ccccc1B1OC(C)(C)C(C)(C)O1",
        "description": "PhBr + B2pin2 → PhBpin (Miyaura borylation precursor - using correct B2pin2 structure)",
        "expected_features": {"boron_type": "B2pin2", "lg_class": "Br"}
    },
    
    "M-SUZ-OTf-DPPF": {
        "smiles": "FC(F)(F)S(=O)(=O)Oc1ccc(C)cc1.OB(O)c1ccccc1>>Cc1ccc(-c2ccccc2)cc1",
        "description": "4-Methyl phenyl triflate + PhB(OH)2 (scheme-based OTf)",
        "expected_features": {"lg_class": "OTf"}
    },
    
    "M-SUZ-VINYL-RT": {
        "smiles": "IC=C.OB(O)C=Cc1ccccc1>>C=C/C=C/c1ccccc1",
        "description": "Vinyl iodide + styrylboronic acid (scheme vinyl)",
        "expected_features": {"is_vinyl": True}
    },
    
    "M-SUZ-BF3K-GENERAL": {
        "smiles": "Brc1ccc(OC)cc1.OB(O)c1ccccc1>>COc1ccc(-c2ccccc2)cc1",
        "description": "4-Bromoanisole + PhBF3K (simplified - use boronic acid as proxy)",
        "expected_features": {"lg_class": "Br"}
    },
    
    "SCDB-SUZ-HET-AZINE-BORON": {
        "smiles": "Brc1ccccc1.OB(O)c1cnccn1>>c1ccc(-c2cnccn2)cc1",
        "description": "PhBr + pyrimidine-5-boronic acid (azine boron partner)",
        "expected_features": {"boron_hetero_count": 2}
    },
    
    # Default condition entries (5 total) - these should catch general cases
    "SCDB-SUZ-DEFAULT-ArI_ArBr": {
        "smiles": "Brc1ccc(F)cc1.OB(O)c1ccccc1>>Fc1ccc(-c2ccccc2)cc1",
        "description": "4-Fluorobromobenzene + PhB(OH)2 (general ArBr default)",
        "expected_features": {"lg_class": "Br"}
    },
    
    "SCDB-SUZ-DEFAULT-ArCl": {
        "smiles": "Clc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "description": "Chlorobenzene + PhB(OH)2 (general ArCl default, no special activation)",
        "expected_features": {"lg_class": "Cl"}
    },
    
    "SCDB-SUZ-DEFAULT-2HET_prone": {
        "smiles": "Brc1ccccc1.OB(O)c1ccco1>>c1ccc(-c2ccco2)cc1",
        "description": "PhBr + 2-furylboronic acid (protodeboronation-prone default)",
        "expected_features": {"boron_protodeboronation_prone": True}
    },
    
    "SCDB-SUZ-DEFAULT-vinyl": {
        "smiles": "BrC=Cc1ccccc1.OB(O)C=C>>C=C/C=C/c1ccccc1",
        "description": "β-Bromostyrene + vinylboronic acid (vinyl default)",
        "expected_features": {"is_vinyl": True}
    },
    
    "SCDB-SUZ-DEFAULT-OTf": {
        "smiles": "FC(F)(F)S(=O)(=O)Oc1ccc(C(F)(F)F)cc1.OB(O)c1ccccc1>>FC(F)(F)c1ccc(-c2ccccc2)cc1",
        "description": "4-CF3 phenyl triflate + PhB(OH)2 (general OTf default)",
        "expected_features": {"lg_class": "OTf"}
    }
}


# ═══════════════════════════════════════════════════════════
# PHASE 3/4 FEATURE TEST REACTIONS
# ═══════════════════════════════════════════════════════════

PHASE_3_4_REACTIONS = [
    # Polyhalogenated substrates
    "Brc1cc(Br)cc(Br)c1.OB(O)c1ccccc1>>Brc1cc(Br)cc(-c2ccccc2)c1 (Suzuki - 1,3,5-Tribromobenzene selective coupling)",
    "Fc1c(F)c(F)c(Br)c(F)c1F.OB(O)c1ccccc1>>Fc1c(F)c(F)c(-c2ccccc2)c(F)c1F (Suzuki - Pentafluorobromobenzene)",
    "Clc1ccc(Cl)nc1.OB(O)c1ccccc1>>Clc1ccc(-c2ccccc2)nc1 (Suzuki - 2,5-Dichloropyridine)",
    "Brc1c(Cl)cccc1Cl.OB(O)c1ccccc1>>c1ccccc1-c1c(Cl)cccc1Cl (Suzuki - Mixed halide with ortho-chlorines)",
    
    # Sterically hindered substrates
    "Brc1ccccc1C(C)(C)C.OB(O)c1ccccc1>>CC(C)(C)c1ccccc1-c1ccccc1 (Suzuki - Ortho-tert-butyl sterics)",
    "Clc1cc(C(C)C)ccc1C(C)C.OB(O)c1ccccc1>>CC(C)c1cccc(C(C)C)c1-c1ccccc1 (Suzuki - Diisopropyl sterics)",
    "Brc1c(C)c(C)c(C)c(C)c1C.OB(O)c1ccccc1>>Cc1c(C)c(C)c(C)c(C)c1-c1ccccc1 (Suzuki - Pentamethyl extreme sterics)",
    "Ic1c(OC)cccc1OC.OB(O)c1ccccc1>>COc1cccc(OC)c1-c1ccccc1 (Suzuki - Ortho-dimethoxy with EDG)",
    
    # Protected intermediates (BOC, CBZ, FMOC, silyl)
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.OB(O)c1ccccc1>>CC(C)(C)OC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - BOC-protected aniline)",
    "Brc1ccc(NC(=O)OCc2ccccc2)cc1.OB(O)c1ccccc1>>c1ccccc1COC(=O)Nc1ccc(-c2ccccc2)cc1 (Suzuki - CBZ-protected aniline)",
    "Clc1ccc(NC(=O)OCC2c3ccccc3-c3ccccc32)cc1.OB(O)c1ccccc1>>c1ccccc1-c1ccc(NC(=O)OCC2c3ccccc3-c3ccccc32)cc1 (Suzuki - FMOC-protected aniline)",
    "Brc1ccc(O[Si](C)(C)C(C)(C)C)cc1.OB(O)c1ccccc1>>CC(C)(C)[Si](C)(C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TBS-protected phenol)",
    "Ic1ccc(O[Si](C(C)C)(C(C)C)C(C)C)cc1.OB(O)c1ccccc1>>CC(C)[Si](C(C)C)(C(C)C)Oc1ccc(-c2ccccc2)cc1 (Suzuki - TIPS-protected phenol)",
    
    # Strong EWG/EDG substrates
    "Brc1ccc([N+](=O)[O-])cc1.OB(O)c1ccccc1>>[O-][N+](=O)c1ccc(-c2ccccc2)cc1 (Suzuki - p-Nitrobromobenzene strong EWG)",
    "Clc1ccc(C#N)cc1.OB(O)c1ccccc1>>N#Cc1ccc(-c2ccccc2)cc1 (Suzuki - p-Chlorobenzonitrile strong EWG)",
    "Brc1ccc(N(C)C)cc1.OB(O)c1ccccc1>>CN(C)c1ccc(-c2ccccc2)cc1 (Suzuki - p-Dimethylaminobromobenzene strong EDG)",
    "Ic1ccc(OC)cc1.OB(O)c1ccccc1>>COc1ccc(-c2ccccc2)cc1 (Suzuki - p-Iodoanisole strong EDG)",
    
    # Chiral substrates
    "NC(Cc1ccccc1)C(=O)O.Brc1ccccc1>>Brc1ccccc1.NC(Cc1ccccc1)C(=O)O (C-N - Phenylalanine amination, chiral)",
    "C[C@H](N)C(=O)O.Brc1ccccc1>>Brc1ccccc1.C[C@H](N)C(=O)O (C-N - (S)-Alanine amination, chiral center)",
    
    # Fused ring systems
    "Brc1ccc2ccccc2c1.OB(O)c1ccccc1>>c1ccccc1-c1ccc2ccccc2c1 (Suzuki - 2-Bromonaphthalene fused bicyclic)",
    "Clc1ccc2c3ccccc3ccc2c1.OB(O)c1ccccc1>>c1ccccc1-c1ccc2c3ccccc3ccc2c1 (Suzuki - 9-Chloroanthracene fused tricyclic)",
    
    # Buchwald-Hartwig with Phase 3/4 features
    "Brc1ccc([N+](=O)[O-])cc1.Nc1ccccc1>>[O-][N+](=O)c1ccc(Nc2ccccc2)cc1 (C-N - p-Nitrobromobenzene + aniline, strong EWG)",
    "Ic1ccc(NC(=O)OC(C)(C)C)cc1.Nc1ccccc1>>CC(C)(C)OC(=O)Nc1ccc(Nc2ccccc2)cc1 (C-N - BOC-protected with amination)",
    "Brc1c(OC)cccc1OC.Nc1ccccc1>>COc1cccc(OC)c1Nc1ccccc1 (C-N - Ortho-dimethoxy sterically hindered)",
    "Brc1ccc2ccccc2c1.NCCC>>CCCNc1ccc2ccccc2c1 (C-N - 2-Bromonaphthalene fused ring + propylamine)",
    
    # Sonogashira with Phase 3/4 features
    "Brc1cc(Br)cc(Br)c1.C#Cc1ccccc1>>Brc1cc(Br)cc(C#Cc2ccccc2)c1 (Sonogashira - Tribromobenzene polyhalogenated)",
    "Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1 (Sonogashira - Ortho-tert-butyl sterics)",
]


def get_suzuki_db_tests():
    """Get test reactions for suzuki_db.json validation"""
    return SUZUKI_DB_TEST_REACTIONS
