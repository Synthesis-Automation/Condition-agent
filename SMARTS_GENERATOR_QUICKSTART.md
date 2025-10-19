# SMARTS Generator Quick Start

This guide shows you how to use the Reaction SMARTS Applicability Pattern Generator.

## What is it?

A tool that helps you define the **scope and limitations** of chemical reaction protocols by converting reaction SMILES into structured SMARTS patterns.

## Why use it?

When adding protocols to your database, you need to specify:
- ✅ What substrates the protocol works with
- ❌ What substrates it doesn't work with
- 📋 Any structural requirements

This tool automates the generation of these patterns.

## Quick Start

### Option 1: Interactive Mode (Easiest)

```powershell
# Windows
.\smarts_generator.bat --interactive

# Or use Python directly
python -m chemtools.protocol.smarts_generator_cli --interactive
```

### Option 2: Single Reaction

```powershell
python -m chemtools.protocol.smarts_generator_cli `
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" `
  --output my_protocol.json
```

### Option 3: Batch Processing

```powershell
# Process multiple reactions from a file
python -m chemtools.protocol.smarts_generator_cli `
  --batch examples\example_reactions.txt `
  --output results.json
```

## Example Workflow

Let's say you have the alkyl iodide borylation protocol:

```
Reaction SMILES: CCCCCCCCI.B2pin2>>CCCCCCCCBpin
```

**Step 1:** Run the generator
```powershell
python -m chemtools.protocol.smarts_generator_cli --interactive
```

**Step 2:** Enter your reaction when prompted
```
Enter reaction SMILES: CCCCCCCCI.B2pin2>>CCCCCCCCBpin
```

**Step 3:** Refine the core pattern
```
Enter core SMARTS pattern:
[C:1;H2,H3;X4;!$(C-[a]);!$(C-[C]=[C])]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])
```

This pattern means:
- `[C:1;H2,H3;X4]` = Primary or secondary carbon (2-3 hydrogens)
- `!$(C-[a])` = NOT bonded to aromatic
- `!$(C-[C]=[C])` = NOT allylic
- `[I:2]` = Iodide leaving group
- `[B:3;X3]` = Boron product

**Step 4:** Add forbidden patterns
```
Forbid: [CH]-I         # No secondary iodides
Forbid: [C;H0]-I       # No tertiary iodides  
Forbid: [CH2]-[a]-I    # No benzylic iodides
Forbid: [CH2]-[C]=[C]-I # No allylic iodides
```

**Step 5:** Add notes
```
Enter description: Primary alkyl iodides only; no aromatic/allylic substrates
```

**Step 6:** Save the output
```json
{
  "reaction_smarts_applicability": {
    "core": "[C:1;H2,H3;X4;!$(C-[a]);!$(C-[C]=[C])]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])",
    "guards_forbid": [
      "[CH]-I",
      "[C;H0]-I", 
      "[CH2]-[a]-I",
      "[CH2]-[C]=[C]-I"
    ],
    "notes": "Primary alkyl iodides only; no aromatic/allylic substrates"
  }
}
```

## Using the Output

Copy the generated JSON into your protocol file:

```json
{
  "source": { ... },
  "reaction": {
    "reaction_smiles": "CCCCCCCCI.B2pin2>>CCCCCCCCBpin",
    "family": "Alkyl_Iodide_Borylation",
    "reaction_smarts_applicability": {
      "core": "[C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3]",
      "guards_forbid": ["[CH]-I", "[C;H0]-I"]
    }
  },
  "reaction_setup": [ ... ]
}
```

## Common SMARTS Patterns

### Atom Types
- `[C]` = Any carbon
- `[C;H2]` = Carbon with 2 hydrogens (CH2)
- `[C;H3]` = Carbon with 3 hydrogens (CH3)
- `[C;X4]` = Carbon with 4 connections
- `[a]` = Any aromatic atom
- `[A]` = Any aliphatic atom

### Common Exclusions
- `[CH]-I` = Secondary iodide
- `[C;H0]-I` = Tertiary iodide  
- `[CH2]-[a]-I` = Benzylic iodide
- `[CH2]-[C]=[C]-I` = Allylic iodide
- `!$(C-[a])` = NOT bonded to aromatic
- `!$(C-C=C)` = NOT bonded to double bond

### Atom Mapping
Use `:1`, `:2`, `:3` to track atoms from reactant to product:
```
[C:1]-[I:2]>>[C:1]-[B:3]
  :1 = carbon (stays)
  :2 = iodine (leaves)
  :3 = boron (new)
```

## Tips

1. **Start with auto-generated pattern**, then refine manually
2. **Test your patterns** against known examples
3. **Document edge cases** in the notes field
4. **Use guard patterns** to exclude problematic substrates
5. **Balance specificity** - not too broad, not too narrow

## Need Help?

- Check if RDKit is installed: `python -m chemtools.protocol.smarts_generator_cli --check-rdkit`
- See full documentation: `docs\SMARTS_GENERATOR_GUIDE.md`
- View examples: `examples\example_reactions.txt`

## More Information

For detailed documentation, see:
- [Full SMARTS Generator Guide](docs/SMARTS_GENERATOR_GUIDE.md)
- [SMARTS Syntax Reference](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
- [RDKit SMARTS Documentation](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support)
