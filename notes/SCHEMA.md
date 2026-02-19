# Notes Schema — Contributor Reference

This document describes the format for all note files in the `notes/` system.
Every note file has a YAML front matter block followed by source content sections.

---

## Folder Structure

```
notes/
├── reactions/    ← per-reaction-type chemistry notes
├── mechanisms/   ← mechanistic principle notes
├── substrates/   ← substrate class notes
├── protocols/    ← practical technique / procedure notes
├── _index.json   ← auto-generated; run `python notes/build_index.py` to rebuild
└── SCHEMA.md     ← this file
```

Add new notes to the appropriate subfolder. Run `build_index.py` after adding files.

---

## Universal Front Matter

Every note file **must** start with a YAML front matter block:

```yaml
---
id: snake_case_identifier        # unique, used for lookup
type: reaction|mechanism|substrate|protocol
title: Human-Readable Title
aliases: [alternative_name, another_alias]   # optional
tags: [keyword1, keyword2, keyword3]         # 5–10 snake_case keywords
---
```

Rebuild `_index.json` after editing front matter: `python notes/build_index.py`

---

## Reaction Notes (`notes/reactions/`)

```yaml
---
id: suzuki_miyaura
type: reaction
title: Suzuki-Miyaura Cross-Coupling
aliases: [suzuki, miyaura coupling]
bond_formed: [C-C]
metal: [palladium, copper, nickel]
mechanism: [oxidative_addition, transmetalation, reductive_elimination]
substrates: [aryl_halide, boronic_acid, boronate_ester]
tags: [pd, cross_coupling, boronate, sp2_coupling, c_c_bond]
related_reactions: [negishi_coupling, heck_reaction]
related_mechanisms: [oxidative_addition, reductive_elimination]
related_substrates: [aryl_halides, boronic_acids]
sources_count: 3
last_updated: YYYY-MM-DD
---
```

**Sections** (add after front matter, one `## Source:` block per literature source):

```markdown
## Source: Journal/Procedure Name  ·  YYYY-MM-DD
url: https://...
doi: 10.xxxx/...
journal: Org. Synth.
year: 2025
pages: 102, 86–99
tags: tag1, tag2, tag3

### Reaction Type
Named reaction and key variant.

### Mechanism Overview          ← NEW: brief mechanistic rationale
Key steps, rate-limiting step, why conditions matter electronically/sterically.

### Solvent Notes
✓ Good: "THF/H₂O (3:1) — reason"
✗ Avoid: "DMF — reason"

### Reagent and Catalyst Notes
Catalyst/ligand preferences and incompatibilities.

### Substrate Scope and Limitations
What works, what is problematic.

### Functional Group Tolerance   ← NEW: explicit compatibility
✓ Tolerates: ester, nitrile, ketone
✗ Incompatible: free amine (coordinates catalyst), aldehyde (reduction risk)

### Critical Conditions
Temperature, atmosphere, addition order, concentration effects.

### Side Reactions
Known side reactions + suppression strategies.

### Procedure Hints              ← NEW: practical steps
Addition order, inert atmosphere requirements, key observations.

### Scale-up Notes               ← NEW
Concentration, heat transfer, mixing effects at larger scale.

### Analytical Notes             ← NEW
TLC Rf, NMR monitoring tips, LCMS signatures.

### Yield / Troubleshooting Tips
Practical observations for optimization.

---

### HTE Data                     ← for entries with HTE experimental backing
Source: AbbVie HTE Perspective — doi: ...
| Rank | Catalyst | Base | Solvent | Top-1 Hits | Top-3 Hits |
|------|----------|------|---------|-----------|-----------|
| 1    | ...      | ...  | ...     | ...       | ...       |
```

---

## Mechanism Notes (`notes/mechanisms/`)

```yaml
---
id: oxidative_addition
type: mechanism
title: Oxidative Addition
applies_to_reactions: [suzuki_miyaura, buchwald_hartwig, negishi, heck]
metal: [palladium, nickel, rhodium]
oxidation_change: "M(0) → M(II)"
bond_broken: [C-X, C-H]
tags: [organometallic, d8_d10, ligand_effect, rate_determining, transition_metal]
related: [reductive_elimination, transmetalation]
---
```

**Sections:**

```markdown
### Overview
What the step is, in which reactions it appears, and why it matters.

### Elementary Steps
Stepwise description of the process with key intermediates.

### Electronic Factors
How ligand/substrate electronics affect rate or selectivity.

### Steric Factors
How ligand/substrate sterics affect rate or selectivity.

### Competing Pathways
What can go wrong or compete (e.g. β-H elimination vs. reductive elimination).

### Examples
Named reactions where this step is rate-limiting or selectivity-determining.
```

---

## Substrate Notes (`notes/substrates/`)

```yaml
---
id: aryl_halides
type: substrate
title: Aryl Halides
variants: [aryl_bromide, aryl_chloride, aryl_iodide, aryl_triflate]
used_in: [suzuki_miyaura, buchwald_hartwig, negishi, heck]
tags: [electrophile, aryl, halide, oxidative_addition, leaving_group]
related: [heteroaryl_halides, vinyl_halides, alkyl_halides]
---
```

**Sections:**

```markdown
### Overview
Reactivity, bond character, key properties.

### Reactivity Trends
Relative order (ArI > ArOTf > ArBr >> ArCl), electronic effects.

### Functional Group Compatibility
What FGs on the substrate are tolerated vs. problematic.

### Preparation / Availability
How to make or source; common commercial forms.

### Handling and Stability
Storage, moisture/air sensitivity, shelf life.

### Used In Reactions
Which reactions use this substrate class, with notes on conditions.
```

---

## Protocol Notes (`notes/protocols/`)

```yaml
---
id: inert_atmosphere
type: protocol
title: Inert Atmosphere Setup (Schlenk / Glovebox)
applies_to: [air_sensitive, organometallic, pd_coupling, cu_coupling]
tags: [technique, schlenk, n2, argon, air_sensitive, glove_box, vacuum]
---
```

**Sections:**

```markdown
### Purpose
What problem this protocol solves.

### When to Use
Which reactions or conditions require this protocol.

### Steps
Numbered procedure. Keep principles, not specific quantities.

### Common Mistakes
What goes wrong and how to avoid it.

### Variations
Alternative approaches (Schlenk vs. glovebox, balloon vs. Schlenk line).
```

---

## Naming Conventions

| Pattern | Example |
|---------|---------|
| Reaction IDs | `suzuki_miyaura`, `buchwald_hartwig_amination` |
| Mechanism IDs | `oxidative_addition`, `reductive_elimination` |
| Substrate IDs | `aryl_halides`, `boronic_acids` |
| Protocol IDs | `inert_atmosphere`, `dry_solvents` |
| Tags | `snake_case` only — `beta_hydride_elimination`, `sp3_coupling` |
| Metal values | lowercase full name — `palladium`, `copper` |
| Bond values | hyphenated — `C-C`, `C-N`, `C-O`, `C-B` |

---

## Intake via CLI

```bash
# Add a reaction note from a URL
python chem_coworker/cli.py intake https://orgsyn.org/... --note-type reactions

# Add a mechanism note from a PDF
python chem_coworker/cli.py intake paper.pdf --note-type mechanisms --reaction-type oxidative_addition

# Rebuild index after manual edits
python notes/build_index.py
```
