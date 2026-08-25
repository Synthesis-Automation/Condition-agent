# My notes

python -m pip install --upgrade rdkit
python -m app.web_api

cmd /d /c "cd /d C:\Git-softwares\Condition-agent && C:\Users\xubar\AppData\Local\Programs\Python\Python312\python.exe -m chem_coworker"

## 2-step tandem reaction

CC(C)(C)OC(=O)NCCN.O=c1oc2cc(Br)ccc2cc1-c1ccccc1>>NCCNc1ccc2cc(-c3ccccc3)c(=O)oc2c1

## double sites reaction

Brc1ccc(Br)cc1.CC(N)=O>>CC(=O)Nc1ccc(NC(C)=O)cc1





 Step 1                                                                                                                                           │
│ BrCc1ccccc1.OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCOCc3ccccc3)(CC1)CC2                                │
│ verified_signature; L2; score 0.401                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Step 2                                                                                                                                           │
│ O=CC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2 >> OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2                                                    │
│ verified_signature; L2; score 0.375                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Unresolved leaves                                                                                                                                │
│ • O=CC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2                                                                                                  │
│                                                                                                                                                  │
│ LLM review (advisory)                                                                                                                            │
│ downrank (70%): Route-8 avoids the impossible quaternary ammonium alkylation present in many other routes. Step-2 is a simple aldehyde           │
│ reduction, which is plausible. Step-1 is an etherification of a primary alcohol with benzyl bromide; a structural competition warning indicates  │
│ the tertiary alcohol on the diphenylmethanol group could also react, but this may be manageable with appropriate conditions. However, the route  │
│ is incomplete (aldehyde leaf not terminal), and the competition warning introduces uncertainty. No strategic guidance was provided. Downranked   │
│ due to the selectivity risk and incomplete status.                                                                                               │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────────────────────────────── Shown row 2 details · original rank 4 ──────────────────────────────────────────────────────╮
│ Step 1                                                                                                                                           │
│ BrCc1ccccc1.OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCOCc3ccccc3)(CC1)CC2                                │
│ verified_signature; L2; score 0.401                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Step 2                                                                                                                                           │
│ Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCO)(CC1)CC2 >> OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2                                                   │
│ verified_signature; L2; score 0.386                                                                                                              │
│ Conditions: recommended fallback                                                                                                                 │
│                                                                                                                                                  │
│ Unresolved leaves                                                                                                                                │
│ • O=C(c1ccccc1)C12CC[N+](CCO)(CC1)CC2                                                                                                            │
│                                                                                                                                                  │
│ LLM review (advisory)                                                                                                                            │
│ downrank (80%): Route-4 uses a Grignard addition (step-2) to install the diphenylmethanol. The ketone precursor contains a free primary hydroxyl │
│ group (CCOH), which would quench the Grignard reagent unless protected. The route does not include a protection step, and the condition evidence │
│ does not address this incompatibility. Step-1 etherification has a competition warning with the tertiary alcohol. The route is incomplete        │
│ (ketone leaf not terminal). Downranked for chemoselectivity and missing protection.                                                              │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────────────────────────────── Shown row 3 details · original rank 5 ──────────────────────────────────────────────────────╮
│ Step 1                                                                                                                                           │
│ OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCBr)(CC1)CC2.OCc1ccccc1 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCOCc3ccccc3)(CC1)CC2                                │
│ verified_signature; L2; score 0.387                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Step 2                                                                                                                                           │
│ Brc1ccccc1.O=C(c1ccccc1)C12CC[N+](CCBr)(CC1)CC2 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCBr)(CC1)CC2                                                 │
│ verified_signature; L2; score 0.381                                                                                                              │
│ Conditions: recommended fallback                                                                                                                 │
│                                                                                                                                                  │
│ Unresolved leaves                                                                                                                                │
│ • O=C(c1ccccc1)C12CC[N+](CCBr)(CC1)CC2                                                                                                           │
│                                                                                                                                                  │
│ LLM review (advisory)                                                                                                                            │
│ downrank (85%): Route-5 performs a Grignard addition (step-2) on a ketone that also bears an alkyl bromide group. Phenylmagnesium bromide would  │
│ react competitively with both the ketone and the alkyl halide, leading to severe selectivity issues. The route is incomplete (ketone leaf not    │
│ terminal). Downranked for chemoselectivity and lack of a protecting strategy.                                                                    │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────────────────────────────── Shown row 4 details · original rank 7 ──────────────────────────────────────────────────────╮
│ Step 1                                                                                                                                           │
│ BrCc1ccccc1.OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCOCc3ccccc3)(CC1)CC2                                │
│ verified_signature; L2; score 0.401                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Step 2                                                                                                                                           │
│ C1CCOC1.CCOC(=O)C12CC[N+](CCO)(CC1)CC2.[Li][c]1ccccc1 >> OCC[N+]12CCC(C(O)(c3ccccc3)c3ccccc3)(CC1)CC2                                            │
│ verified_signature; L2; score 0.333                                                                                                              │
│ Conditions: recommended fallback                                                                                                                 │
│                                                                                                                                                  │
│ Unresolved leaves                                                                                                                                │
│ • CCOC(=O)C12CC[N+](CCO)(CC1)CC2                                                                                                                 │
│                                                                                                                                                  │
│ LLM review (advisory)                                                                                                                            │
│ downrank (80%): Route-7 uses phenyllithium to add to an ester (step-2). The ester precursor contains a free primary hydroxyl group, which would  │
│ instantly deprotonate and consume the organolithium reagent. Additionally, the quaternary ammonium may be susceptible to Hofmann elimination     │
│ under strong base. The route is incomplete (ester leaf not terminal). Downranked for chemoselectivity and incompatibility.                       │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭───────────────────────────────────────────────────── Shown row 5 details · original rank 1 ──────────────────────────────────────────────────────╮
│ Step 1                                                                                                                                           │
│ OC(c1ccccc1)(c1ccccc1)C1CC[N+](CCBr)(CCOCc2ccccc2)CC1 >> OC(c1ccccc1)(c1ccccc1)C12CC[N+](CCOCc3ccccc3)(CC1)CC2                                   │
│ verified_signature; L2; score 0.392                                                                                                              │
│ Conditions: recommended fallback                                                                                                                 │
│                                                                                                                                                  │
│ Step 2                                                                                                                                           │
│ BrCc1ccccc1.OCC[N+]1(CCBr)CCC(C(O)(c2ccccc2)c2ccccc2)CC1 >> OC(c1ccccc1)(c1ccccc1)C1CC[N+](CCBr)(CCOCc2ccccc2)CC1                                │
│ verified_signature; L2; score 0.417                                                                                                              │
│ Conditions: recommended direct                                                                                                                   │
│                                                                                                                                                  │
│ Unresolved leaves                                                                                                                                │
│ • OCC[N+]1(CCBr)CCC(C(O)(c2ccccc2)c2ccccc2)CC1                                                                                                   │
│                                                                                                                                                  │
│ LLM review (advisory)                                                                                                                            │
│ flag (95%): Step-1 proposes an intramolecular N-alkylation of a quaternary ammonium. The nitrogen is already fully substituted and positively    │
│ charged; it cannot act as a nucleophile. This step is chemically impossible. The forward validation may have passed the signature, but the       │
│ reaction is not feasible. The route is therefore unsound and flagged.                                                                            │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭────────────────────────────────────────────────────────────── LLM multistep review ──────────────────────────────────────────────────────────────╮
│ All routes are partial and unsolved. A critical recurring error is an intramolecular N-alkylation of a quaternary ammonium, which is chemically  │
│ impossible because the nitrogen is already saturated. This affects routes 1, 2, 3, 6, 9, and 10. Among the remaining routes, routes 4, 5, and 7  │
│ propose Grignard or organolithium additions to carbonyl compounds in the presence of unprotected hydroxyl groups or alkyl halides, which would   │
│ cause chemoselectivity or compatibility problems. Route 8 avoids these impossible steps and uses a plausible aldehyde reduction and              │
│ etherification, but it still has a functional‑group competition warning and an unresolved leaf. Therefore, route 8 is ranked highest, followed   │
│ by routes 4, 5, and 7 (all downranked due to chemoselectivity concerns), and the routes containing the impossible quaternary ammonium alkylation │
│ are flagged and ranked lowest.                                                                                                                   │
│                                                                                                                                                  │
│ Questions                                                                                                                                        │
│ • For routes 4, 5, and 7, would the Grignard/organolithium step be feasible if the hydroxyl group were protected? The routes do not show         │
│ protection, so the current proposal is not viable.                                                                                               │
│                                                                                                                                                  │
│ deepseek-v4-pro · advisory · 46977 tokens · 2 provider attempt(s)                                                                                │
╰──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
╭──────────────────────────────────────────────────────────────────── Warnings ────────────────────────────────────────────────────────────────────╮
│ • No fully solved route was found within the configured search bounds; partial routes remain explicitly unresolved.                 