# Split Mode Feature - UnifiedRecommender Interactive CLI

## Overview

The interactive CLI now supports **Split Mode**, which displays the top-ranked protocol and top-ranked rule separately with their full detailed recommended conditions. This makes it easy to compare both experimental protocols and general reaction guidelines side-by-side.

## Usage

### Command Line

```bash
# Show top protocol + top rule with detailed conditions
python app/unified_interactive_cli.py --rxn "REACTION_SMILES" --split
```

### Interactive Mode

```bash
# Start interactive mode
python app/unified_interactive_cli.py

# Enable split mode
reaction> /split on

# Query a reaction
reaction> Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1
```

## Example Output

### Suzuki-Miyaura Coupling

```bash
python app/unified_interactive_cli.py --rxn "CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1" --split
```

**Output:**

```
================================================================================
Query: CC(=O)c1ccc(Br)cc1.c1ccc(B(O)O)cc1>>CC(=O)c1ccc(-c2ccccc2)cc1
================================================================================

━━━ TOP PROTOCOL ━━━
📋 Nickel-Catalyzed Suzuki-Miyaura Coupling in t-Amyl Alcohol
   Similarity: 0.220
   Family: Ni_Suzuki_ArylHalide+BoronicAcid_tAmOH
   ID: nickel_catalyzed_suzuki_miyaura_coupling_in_t_amyl

  Protocol Conditions:
    • Catalyst: Bis(tricyclohexylphosphine)nickel(II) chloride
    • Solvent: tert-Amyl alcohol (2-methyl-2-butanol)
    • Base: Potassium phosphate tribasic (anhydrous)
    • Temperature: 120 °C
    • Time: 1.0 h
    • Atmosphere: N₂ (inert atmosphere)

━━━ TOP RULE ━━━
📖 Suzuki-Miyaura Coupling (Improved, HTE-ready)
   Similarity: 0.794
   Family: Suzuki_Miyaura
   ID: suzuki_v2

  Default Conditions:
    (General starter for sp2–sp2 Suzuki; robust across many substrates.)
    • Pd Source: PdCl2(dtbpf)
    • Catalyst Loading: 1.5 mol%
    • Ligand: dtbpf
    • Solvent: 1,4-dioxane/H2O 4:1
    • Base: K3PO4 (aq, 3.25 M)
    • Base Equivalents: 2.0
    • Temperature: 80–100 °C
    • Time: 1–6 h
    • Atmosphere: N2 or Ar

  Alternative Conditions:
    1. dtbpf / K3PO4 (most frequent hit)
       General, scalable; good for broad aryl bromide/iodide and many aryl chloride cases.
       Pd: PdCl2(dtbpf) | Ligand: dtbpf | Base: K3PO4 (aq, 3.25 M) | Temp: 80–100 °C

    2. P(tBu)3–Pd G3 / K3PO4 (hard substrates)
       Excellent for hindered aryl chlorides and electron-poor boronates; fast at elevated T.
       Catalyst: Pd-PEPPSI-type G3 with P(tBu)3 (representative) | Ligand: P(tBu)3 | Base: K3PO4 (aq, 3.25 M) | Temp: 90–110 °C

    3. XPhos–Pd G3 / K3PO4
       Dialkylbiaryl phosphine platform; reliable on aryl bromides and many chlorides.
       Catalyst: XPhos–Pd G3 | Ligand: XPhos | Base: K3PO4 (aq, 3.25 M) | Temp: 90–110 °C
```

### Buchwald-Hartwig C-N Coupling

```bash
python app/unified_interactive_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --split
```

**Output:**

```
================================================================================
Query: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1
================================================================================

No protocol recommendations found.

━━━ TOP RULE ━━━
📖 Buchwald–Hartwig C–N Coupling
   Similarity: 1.000
   Family: Buchwald–Hartwig_C–N
   ID: c_n_coupling_pd_v2

  Default Conditions:
    (General starter for Buchwald–Hartwig C–N (weak aqueous base; ligand precatalyst).)
    • Ligand: match selected base_rule (default tBuBrettPhos‑class)
    • Solvent: MeCN / t‑AmOH / dioxane / toluene (match polarity)
    • Base: K3PO4 (aq, ~3.25 M) 3.0 equiv (or Cs2CO3)
    • Temperature: 80 start (80–110 typical; Arom.NH often 100 start) °C
    • Time: 2–8 (if stalled after ~4 h, add 2–4 mol% Pd) h
    • Atmosphere: N2 or Ar; degas before base addition

  Alternative Conditions:
    1. Primary amines (aliphatic or aniline)
       Primary NH2 nucleophiles (aliphatic or aryl). First shot: t‑AmOH/Cs2CO3.
       Ligand: tBuBrettPhos or BippyPhos | Base: K3PO4 (aq) or Cs2CO3; alt: NaOtBu or KHMDS (strong) | Temp: 80–100 °C

    2. Secondary amines
       Secondary NHR nucleophiles. Weak base often Cs2CO3; strong base if needed.
       Ligand: QPhos, XantPhos, or tBuBrettPhos | Base: Cs2CO3 (weak); alt: KHMDS or NaOtBu (strong) | Temp: 80–105 °C

    3. Heteroaromatic NH (indoles, pyrazoles, etc.)
       Heteroaryl N–H coupling; Arom.NH often needs hotter start (~100 °C).
       Ligand: Josiphos SL‑J009‑1, tBuBrettPhos, or XantPhos | Base: K3PO4 (aq) or Cs2CO3; alt: NaOtBu | Temp: 100–110 (start ~100 for Arom.NH) °C
```

### Sonogashira Coupling

```bash
python app/unified_interactive_cli.py --rxn "Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1" --split
```

**Output:**

```
================================================================================
Query: Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1
================================================================================

No protocol recommendations found.

━━━ TOP RULE ━━━
📖 Sonogashira Coupling Guidelines (Pd-Cu / Cu-free)
   Similarity: 1.000
   Family: Sonogashira
   ID: sonogashira_v2

  Default Conditions:
    (General Sonogashira coupling with CuI cocatalyst and amine base. Thoroughly degassed conditions.)
    • Catalyst Loading: 0.5-2.0 mol%
    • Ligand: Built-in (PPh3 or dppf)
    • Solvent: THF, toluene, or DMF
    • Base: Et3N or DIPEA
    • Base Equivalents: 2.0-3.0
    • Temperature: 40-80 °C
    • Time: 1-8 h
    • Atmosphere: N2 or Ar; thoroughly degassed

  Alternative Conditions:
    1. Aryl iodides/bromides (standard reactivity)
       Fast track with mild bases. CuI typically helpful but optional if Glaser coupling is problematic.
       Base: Et3N or DIPEA | Temp: 40-70 °C

    2. Aryl chlorides or deactivated/ortho-blocked aryls
       Hard electrophiles requiring stronger ligation and higher temperature. Copper-free preferred to avoid side reactions.
       Base: DIPEA or Cs2CO3 | Temp: 80-110 °C

    3. Vinyl halides/triflates
       Clean at lower temperatures. Monitor for isomerization under strongly basic or high-temperature conditions.
       Base: Et3N or DIPEA | Temp: 25-60 °C
```

## Features

### For Protocols
- **Catalyst**: Metal catalyst and ligands used
- **Solvent**: Reaction solvent(s)
- **Base**: Base reagent
- **Temperature**: Reaction temperature in °C
- **Time**: Reaction duration in hours
- **Atmosphere**: Inert atmosphere (N₂/Ar)
- **Reported Yield**: If available from the protocol

### For Rules
- **Default Conditions**: The general starting conditions for the reaction type
  - Catalyst/Ligand system
  - Solvent system
  - Base and equivalents
  - Temperature range
  - Typical reaction time
  - Atmosphere requirements

- **Alternative Conditions** (Top 3): Different condition sets optimized for specific substrates or situations
  - Substrate-specific recommendations
  - Troubleshooting alternatives
  - Special case handling

## Interactive Mode Commands

```bash
/split on       # Enable split mode
/split off      # Disable split mode
/split          # Toggle split mode
/show           # Display current settings including split_mode
/help           # Show all commands
```

## Notes

- **Split mode is mutually exclusive** with compact mode
- If no protocol or rule is found, only the available type will be displayed
- The top result from each category (protocol/rule) is selected based on DRFP similarity
- Alternative conditions for rules provide substrate-specific optimizations based on HTE data and literature precedent

## Benefits

1. **Quick Comparison**: See both experimental protocols and general guidelines at once
2. **Detailed Conditions**: All key reaction parameters displayed clearly
3. **Alternative Strategies**: For rules, see multiple condition sets optimized for different scenarios
4. **Research-Ready**: Formatted output suitable for experimental planning
