# When you “understand” an organic reaction (especially from a reaction SMILES / scheme), it helps to decompose it into a small set of major aspects that are largely orthogonal. Your “electrophilic carbon + leaving group” lens is one of the most reliable starting points. Here’s a structured checklist I use

## Net transformation (what changed)

Bond changes: which σ/π bonds formed/broken (C–C, C–N, C–O, C–X, etc.)

Functional group interconversion: e.g., alcohol → halide, aryl halide → aryl amine, alkene → epoxide

Redox change (if any): oxidation state changes, hydrogen balance

Practical: compute/visualize the reaction center (changed atoms) and treat everything else as context.

## Electrophile analysis (your point)

When a leaving group departs from carbon (or carbon becomes more electrophilic):

Hybridization / topology: sp, sp² (aryl/vinyl), sp³; benzylic/allylic/propargylic

Substitution level: methyl/1°/2°/3°; neopentyl penalty

Stabilization of cationic character: resonance (benzylic/allylic), adjacent heteroatoms

Electronic activation: EWG (carbonyl, CN, NO₂, sulfone), heteroaryl effects

Leaving group quality: OTf/OMs/OTs > I > Br > Cl >> F (except SNAr); sulfonates; diazonium

Sterics and conformation: backside access for SN2; planarity for SNAr; β-H availability for E2

This largely predicts SN1/SN2/E1/E2/SNAr/oxidative addition compatibility.

## Nucleophile / basicity profile

You often need to decide: is the reagent behaving as a nucleophile, a base, a reductant, or a ligand?

Hard/soft character: thiolates/phosphines/organocuprates vs alkoxides/amides

Charge and aggregation: anions vs neutral nucleophiles; ion pairing

Steric bulk: tBuO⁻ pushes elimination; bulky amines slow SN2

Ambident nucleophiles: enolates, cyanide, nitrite, thiocyanate

Basicity vs nucleophilicity split: e.g., DBU strong base but not great nucleophile at carbon

## Driving forces & thermodynamics

A fast way to sanity-check plausibility:

Leaving group departure (formation of stable anion or neutral molecule)

Formation of strong bonds: C=O, aromaticity gained, σ-bond formation from π

Relief of strain / ring closure: 5-exo/6-endo preferences; Baldwin rules

Precipitation / gas evolution: salts out, CO₂, N₂, etc.

Acid–base neutralization: proton transfers that “lock in” products

## Mechanistic class (the “minimum viable mechanism”)

Pick the simplest mechanistic family consistent with 2–4:

Polar 2e pathways: SN1/SN2/E1/E2, additions to π systems, acyl substitution

Aromatic substitution: SNAr (Meisenheimer), benzyne, SEAr

Radical 1e pathways: SRN1, Giese-type additions, HAT, radical–polar crossover

Organometallic cycles: oxidative addition / transmetalation / reductive elimination; migratory insertion; β-H elimination

Key diagnostic questions:

Is the electrophile aryl/vinyl? If yes, SN2 is out; think SNAr (if activated) or metal-catalyzed OA.

Is there a strong base + β-H? E2 becomes likely.

Are there redox reagents, light, peroxides, AIBN? Radical pathways become plausible.

## Selectivity map

Even for “one-step” reactions, selectivity is often the core:

Chemoselectivity: which functional group reacts first

Regioselectivity: Markovnikov/anti, SNAr site, conjugate vs 1,2-addition

Stereoselectivity: inversion/retention, E/Z control, diastereocontrol

Competing pathways: substitution vs elimination; over-alkylation; rearrangements

## Conditions as constraints (solvent, temp, catalysts)

Conditions are not decoration—they often determine the mechanism.

Solvent: polar aprotic boosts SN2; protic boosts SN1; HFIP effects; coordinating solvents for metals

Temperature: pushes elimination, rearrangements, equilibrium

Catalysts/ligands/additives: dictate metal oxidation state and rate-limiting step

Counterions: Li vs Na vs K; phase-transfer; crown ethers

Concentration / dilution: intramolecular vs intermolecular outcomes

## Side reactions & failure modes

A good “understanding” includes what could go wrong:

Hydrolysis, elimination, polymerization

Overreaction (multiple alkylations), O- vs N-alkylation

Catalyst poisoning (thiols, amines, heteroarenes), β-H elimination in cross-coupling

Rearrangements (carbocation shifts), racemization

A compact workflow you can apply quickly

Identify bond changes (reaction center).

For each broken C–LG bond, classify carbon: sp²/sp³, 1°/2°/3°, benzylic/allylic, activated by EWG?

Classify the “attacking” species: nucleophile vs base vs radical vs organometallic.

Choose the minimal mechanism class consistent with those facts.

Predict selectivity + main competitors.

Practical add-on (if you’re building an agent)

Represent the above as a structured record:

Electrophile: {center_atom, hybridization, degree, leaving_group, activation_tags, steric_tags}

Nucleophile/base: {role_candidates, pKa_proxy, steric, softness}

Mechanism_candidates: ranked list with rules/heuristics

Selectivity_risks: list

Missing info: what condition data would disambiguate (solvent/base/metal/light)

## Practical “highest ROI” defaults (if you must keep it simple)

Atom-map (external) → bond changes

Electrophile center descriptors (hybridization, degree, benzylic/allylic, EWG radius-2)

Nucleophile identity (which atom forms new bond)

Mechanism shortlist: SN2/E2/SNAr/acyl-substitution/OA-based coupling flags

Selectivity risks: substitution vs elimination; overreaction; ambident sites
