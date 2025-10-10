# Protocol CLI Tester - User Guide

## 🚀 Quick Start

### Single Query Mode

Search for a specific reaction:

```bash
python test_protocol_cli.py "BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1"
```

### Interactive Mode

Run without arguments for interactive mode:

```bash
python test_protocol_cli.py
```

Then enter reactions one at a time:

```
🔬 Reaction SMILES (or command): CCBr.c1ccccc1B(O)O>>CCc1ccccc1
```

## 📖 Command-Line Options

### Basic Options

```bash
# Get top 3 matches
python test_protocol_cli.py "REACTION_SMILES" -k 3

# Filter by reaction family
python test_protocol_cli.py "REACTION_SMILES" --family Suzuki_Cu_alkyl_halide+aryl_boron

# Filter by tags
python test_protocol_cli.py "REACTION_SMILES" --tags suzuki,palladium

# Hide detailed conditions
python test_protocol_cli.py "REACTION_SMILES" --no-conditions

# Enable debug mode
python test_protocol_cli.py "REACTION_SMILES" --debug
```

### Combined Options

```bash
python test_protocol_cli.py "CCBr.c1ccccc1B(O)O>>CCc1ccccc1" -k 3 --tags suzuki --no-conditions
```

## 🎮 Interactive Mode Commands

### Searching for Reactions

Just type or paste a reaction SMILES:

```
🔬 Reaction SMILES (or command): BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1
```

### Setting Options

```
set k 3                                    # Show top 3 results
set family Suzuki_Cu_alkyl_halide+aryl_boron   # Filter by family
set tags suzuki,palladium                  # Filter by tags (comma-separated)
set conditions on                          # Show conditions (default)
set conditions off                         # Hide conditions
```

### Viewing Settings

```
settings                                   # Show current settings
```

### Clearing Filters

```
clear                                      # Clear all filters (family and tags)
```

### Getting Help

```
help                                       # Show help message
```

### Exiting

```
quit
exit
q
```

Or press `Ctrl+C`

## 📊 Output Format

For each match, you'll see:

```
1. Match Score: 0.8018 (80.2%)
----------------------------------------------------------------------
📄 File: Suzuki_Cu_C(sp3)-C(sp2).json
🔬 Family: Suzuki_Cu_alkyl_halide+aryl_boron
📖 Title: Copper-Catalyzed Suzuki-Miyaura Coupling...
📰 Journal: Organic Syntheses (2025)
🔗 DOI: 10.15227/orgsyn.102.0086
🌐 URL: https://www.orgsyn.org/demo.aspx?prep=v102p0086
⚗️  Reaction: BrC1CCCCC1.c1ccc(Cl)cc1B(...)>>...
🏷️  Tags: Suzuki, Cu, alkyl_halide, Coupling
📝 Notes: CuBr·SMe2 (5 mol%), bathophenanthroline...

🧪 Extracted Conditions:
   Catalyst: CuBr·SMe2
   Ligand: Bathophenanthroline
   Base: NaOtBu
   Solvent: toluene
   Temperature: 80 °C
   Time: 24 h
   Atmosphere: air
```

## 💡 Example Reactions to Try

### Suzuki Coupling (Alkyl-Aryl)

```
BrC1CCCCC1.c1ccc(Cl)cc1B(O)O>>Clc1ccc(C2CCCCC2)cc1
```

### Borylation

```
CCCCCI.B(B1OC(C)(C)C(C)(C)O1)B2OC(C)(C)C(C)(C)O2>>CCCCB1OC(C)(C)C(C)(C)O1
```

### C-N Coupling (Buchwald-Hartwig)

```
CC(C)(C)c1ccc(Br)cc1.CNc2ccccc2>>CC(C)(C)c1ccc(N(C)c2ccccc2)cc1
```

### Suzuki Coupling (Aryl-Aryl)

```
Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccccc1c1ccccc1
```

### Sonogashira Coupling

```
Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1
```

## 🎯 Interactive Session Example

```bash
$ python test_protocol_cli.py

======================================================================
🧪 Interactive Protocol Recommendation CLI
======================================================================

Enter a reaction SMILES to find matching protocols.
Type 'help' for options, 'quit' to exit.

----------------------------------------------------------------------

🔬 Reaction SMILES (or command): set k 3

✅ Set k = 3

----------------------------------------------------------------------

🔬 Reaction SMILES (or command): set tags suzuki

✅ Set tag filter: suzuki

----------------------------------------------------------------------

🔬 Reaction SMILES (or command): settings

----------------------------------------------------------------------
⚙️  Current Settings
----------------------------------------------------------------------
  Number of results (k): 3
  Show conditions: ON
  Family filter: None
  Tag filter: suzuki

----------------------------------------------------------------------

🔬 Reaction SMILES (or command): BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1

🔍 Searching for similar protocols...
   Query: BrC1CCCCC1.c1ccccc1B(O)O>>c1ccccc1C1CCCCC1
   Tag filter: suzuki

✅ Found 2 matches
   (Searched 2 candidates)

[... match details ...]

----------------------------------------------------------------------

🔬 Reaction SMILES (or command): quit

Goodbye! 👋
```

## 🔧 Troubleshooting

### Index Not Found

```
❌ Protocol index not found!

Please build the index first:
  python -m chemtools.protocol.cli build
```

**Solution:** Build the index first
```bash
python -m chemtools.protocol.cli build
```

### Invalid Reaction SMILES

```
❌ Invalid reaction SMILES. Must contain '>>' separator.
   Example: CCBr.c1ccccc1B(O)O>>CCc1ccccc1
```

**Solution:** Make sure your reaction SMILES has the `>>` separator between reactants and products.

### No Matches Found

```
❌ No matching protocols found.
   Searched 16 candidates

💡 Tip: Try removing filters with 'clear' command
```

**Solutions:**
1. Remove filters: type `clear`
2. Try a different reaction
3. Check that your reaction SMILES is correct

## 🎨 Features

- ✅ **Single query mode** for quick searches
- ✅ **Interactive mode** for testing multiple reactions
- ✅ **DRFP similarity scoring** (0.0 to 1.0)
- ✅ **Detailed protocol information** (title, DOI, source, etc.)
- ✅ **Extracted conditions** (catalyst, solvent, temperature, etc.)
- ✅ **Flexible filtering** by family and tags
- ✅ **Adjustable result count** (top-k)
- ✅ **Clean, colorful output** with emojis
- ✅ **User-friendly commands** and help system

## 📚 Related Documentation

- **Full Module Docs:** `docs/PROTOCOL_MODULE.md`
- **Quick Start:** `docs/PROTOCOL_QUICKSTART.md`
- **Module Summary:** `docs/PROTOCOL_MODULE_SUMMARY.md`

## 🎉 Happy Testing!

The CLI makes it easy to test the protocol recommendation system with your own reactions. Try different reactions, adjust filters, and explore the protocol database!
