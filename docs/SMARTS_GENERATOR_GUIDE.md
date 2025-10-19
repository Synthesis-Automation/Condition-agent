# Reaction SMARTS Applicability Pattern Generator

A CLI tool for generating structured SMARTS patterns that define the scope and limitations of chemical reaction protocols.

## Purpose

When defining a chemical protocol (like the alkyl iodide borylation example), it's critical to specify:
1. **What reactions it applies to** - Core transformation pattern
2. **What structures it excludes** - Guard patterns that forbid certain features
3. **What structures are required** - Guard patterns that must be present

This tool helps convert reaction SMILES into structured SMARTS patterns with these constraints.

## Installation

```bash
# Ensure RDKit is installed (required for automatic pattern generation)
pip install rdkit

# The tool is included in the chemtools package
# No additional installation needed
```

## Quick Start

### Interactive Mode (Recommended)

The easiest way to use the tool is in interactive mode with guided prompts:

```bash
# Windows
smarts_generator.bat --interactive

# Linux/Mac or direct Python
python -m chemtools.protocol.smarts_generator_cli --interactive
```

The tool will ask if you want to generate visualizations of your patterns to verify they match correctly.

Example session:
```
Enter reaction SMILES: CCCCCCCCI.B2pin2>>CCCCCCCCBpin

Step 1: Core Transformation Pattern
------------------------------------
Auto-generated starting point:
  [#6:1]-[#53:2]>>[#6:3]-[#5:4]

Enter core SMARTS pattern: [C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])

Step 2: Guard Patterns - Forbidden Structures
---------------------------------------------
Suggestions:
  1. [C;H0]-I  # Exclude tertiary iodides
  2. [CH]-I    # Exclude secondary iodides

Enter forbidden patterns:
  Forbid: [CH]-I
  Forbid: [C;H0]-I
  Forbid: [CH2]-[a]-I
  Forbid: [CH2]-[C]=[C]-I
  Forbid: 

Step 3: Guard Patterns - Required Structures
--------------------------------------------
(press Enter to skip)

Step 4: Notes
-------------
Enter description: Primary alkyl iodides only; no aromatic/allylic substrates
```

### Single Reaction Mode

Process a single reaction and save to JSON:

```bash
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \
  --output alkyl_borylation_smarts.json
```

### Batch Processing Mode

Process multiple reactions from a text file (one per line):

```bash
# Create input file
echo "CCCCI>>CCCBpin" > reactions.txt
echo "CCCBr>>CCCOH" >> reactions.txt

# Process all reactions
python -m chemtools.protocol.smarts_generator_cli \
  --batch reactions.txt \
  --output results.json
```

## Visualization & Verification

The tool can generate visual depictions of your SMARTS patterns using RDKit to help verify correctness.

### Generate Visualizations

Add the `--visualize` flag to create PNG images of your patterns:

```bash
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \
  --visualize \
  --viz-dir my_patterns
```

This creates:
- `core_transformation.png` - Visual representation of the reaction transformation
- `guard_forbid_N.png` - Images for each forbidden pattern
- `guard_require_N.png` - Images for each required pattern (if any)

### Test Against Example Molecules

Verify your pattern matches the expected substrates and excludes the wrong ones:

```bash
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \
  --visualize \
  --test-smiles "CCCCI,CC(C)I,CC(C)(C)I,c1ccccc1CI"
```

The tool will report:
- ✅ **MATCH**: Molecule matches core pattern and passes all guards
- ❌ **NO MATCH**: Molecule fails (with explanation of why)

Example output:
```
Testing Patterns Against Examples
----------------------------------------------------------------------
  ✅ MATCH: CCCCI                    # Primary iodide - OK
  ❌ NO MATCH: CC(C)I                # Secondary - forbidden
    → Matches forbidden pattern
  ❌ NO MATCH: CC(C)(C)I             # Tertiary - forbidden
    → Matches forbidden pattern
  ❌ NO MATCH: c1ccccc1CI            # Benzylic - forbidden
    → Matches forbidden pattern
```

### In Interactive Mode

The tool will automatically ask if you want visualizations:

```
Generate visualization images? (Y/n): y
✅ Saved reaction SMARTS visualization to smarts_visualizations\core_transformation.png
✅ Saved SMARTS visualization to smarts_visualizations\guard_forbid_1.png

Enter test SMILES to validate (comma-separated, or press Enter to skip): CCCCI,CC(C)I
  ✅ MATCH: CCCCI
  ❌ NO MATCH: CC(C)I
    → Matches forbidden pattern
```

## Output Format

The tool generates JSON in the format expected by protocol database entries:

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
    "guards_require": [],
    "notes": "Primary alkyl iodides only; no aromatic or allylic substrates"
  }
}
```

### Field Descriptions

- **`core`**: SMARTS pattern describing the core transformation with atom mapping
  - Use `:1`, `:2`, `:3` for atom mapping
  - Include specific atom properties: `H2,H3` (hydrogen count), `X4` (total connections)
  - Use negation: `!$(...)` to exclude substructures
  
- **`guards_forbid`**: List of SMARTS patterns where the reaction does NOT work
  - `[CH]-I` = secondary iodide (1 hydrogen on carbon)
  - `[C;H0]-I` = tertiary iodide (0 hydrogens on carbon)
  - `[CH2]-[a]-I` = benzylic iodide (aromatic adjacent)
  - `[CH2]-[C]=[C]-I` = allylic iodide (double bond adjacent)

- **`guards_require`**: List of SMARTS patterns that MUST be present (optional)

- **`notes`**: Human-readable description of scope and limitations

## SMARTS Pattern Syntax

### Basic Patterns

```
[C]          Carbon atom
[C;H2]       Carbon with exactly 2 hydrogens
[C;H2,H3]    Carbon with 2 OR 3 hydrogens
[C;X4]       Carbon with 4 total connections
[a]          Any aromatic atom
[A]          Any aliphatic atom
```

### Atom Mapping (for transformations)

```
[C:1]-[I:2]>>[C:1]-[B:3]
  :1 = first atom (carbon from reactant → product)
  :2 = second atom (iodine, leaving group)
  :3 = third atom (boron, new in product)
```

### Negation Patterns

```
[C;!H0]          Carbon that is NOT primary (H0 = 0 hydrogens)
[C;!$(C-[a])]    Carbon NOT bonded to aromatic
!$(C-C=C)        NOT matching the substructure C-C=C
```

### Common Exclusions

```
[CH]-I              Secondary iodide
[C;H0]-I            Tertiary iodide
[CH2]-[a]-I         Benzylic iodide
[CH2]-[C]=[C]-I     Allylic iodide
[CH2]-[C;H0]        Neopentyl system
[C;r]               Carbon in ring
[C;!r]              Carbon NOT in ring
```

## Integration with Protocol Database

To add the generated SMARTS pattern to your protocol JSON file:

1. Generate the pattern using this tool
2. Copy the output into your protocol JSON under `reaction.reaction_smarts_applicability`
3. Manually refine the patterns based on experimental evidence
4. Test against known substrates to validate

Example protocol structure:
```json
{
  "source": { ... },
  "reaction": {
    "reaction_smiles": "CCCCCCCCI.B2pin2>>CCCCCCCCBpin",
    "family": "Alkyl_Iodide_Borylation",
    "notes": "...",
    "reaction_smarts_applicability": {
      "core": "[C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3]",
      "guards_forbid": ["[CH]-I", "[C;H0]-I"]
    }
  },
  "reaction_setup": [ ... ]
}
```

## Command-Line Options

```
usage: smarts-generator [-h] [--interactive] [--reaction REACTION]
                       [--batch BATCH] [--output OUTPUT] [--visualize]
                       [--viz-dir VIZ_DIR] [--test-smiles TEST_SMILES]
                       [--check-rdkit]

Generate reaction SMARTS applicability patterns from reaction SMILES

optional arguments:
  -h, --help            show this help message and exit
  --interactive, -i     Run in interactive mode with guided input
  --reaction REACTION, -r REACTION
                        Single reaction SMILES to process
  --batch BATCH, -b BATCH
                        Batch process reactions from text file (one per line)
  --output OUTPUT, -o OUTPUT
                        Output JSON file path
  --visualize, -v       Generate visualization images of SMARTS patterns (requires RDKit)
  --viz-dir VIZ_DIR     Directory for visualization images (default: smarts_visualizations)
  --test-smiles TEST_SMILES
                        Comma-separated SMILES to test against the pattern
  --check-rdkit         Check if RDKit is available and exit
```

## Tips for Writing Good SMARTS Patterns

1. **Start Simple, Add Specificity**
   - Begin with basic transformation: `[C]-[I]>>[C]-[B]`
   - Add atom mapping: `[C:1]-[I:2]>>[C:1]-[B:3]`
   - Add constraints: `[C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3]`

2. **Test Against Known Examples**
   - Use RDKit to test patterns: `Chem.MolFromSmarts(pattern)`
   - Verify matches on positive examples
   - Verify non-matches on negative examples

3. **Document Edge Cases**
   - Note any borderline substrates in `notes` field
   - Include literature references if applicable

4. **Balance Specificity vs Generality**
   - Too specific: Won't match valid substrates
   - Too general: Will match incompatible substrates
   - Use guard patterns to exclude problematic cases

## Troubleshooting

### RDKit Not Available

If you see "RDKit required" error:
```bash
pip install rdkit
# or
conda install -c conda-forge rdkit
```

### Invalid SMILES

Ensure your reaction SMILES:
- Uses correct arrow format (`>>` or `>`)
- Contains valid SMILES on both sides
- Separates multiple reactants with `.`

Example: `CCCCI.B2pin2>>CCCBpin` ✅
Not: `CCCCI + B2pin2 → CCCBpin` ❌

### Pattern Too Generic

If the auto-generated pattern is too generic:
- Use interactive mode to manually refine
- Add specific atom properties (H-count, connectivity)
- Add negation patterns for exclusions

## Examples

### Example 1: Alkyl Iodide Borylation

```bash
python -m chemtools.protocol.smarts_generator_cli --interactive
# Enter: CCCCCCCCI.B2pin2>>CCCCCCCCBpin
```

Result:
```json
{
  "core": "[C:1;H2,H3;X4;!$(C-[a]);!$(C-[C]=[C])]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])",
  "guards_forbid": ["[CH]-I", "[C;H0]-I", "[CH2]-[a]-I", "[CH2]-[C]=[C]-I"],
  "notes": "Primary alkyl iodides only"
}
```

### Example 2: Suzuki Coupling

```bash
python -m chemtools.protocol.smarts_generator_cli \
  -r "c1ccccc1Br.c1ccccc1B(O)O>>c1ccccc1-c1ccccc1"
```

This generates a starting pattern for aryl-aryl Suzuki coupling that you can refine.

### Example 3: Batch Processing Literature Examples

Create `literature_reactions.txt`:
```
CCCCI>>CCCBpin
c1ccccc1Br>>c1ccccc1B(O)O
CCCCl>>CCCOH
```

Process:
```bash
python -m chemtools.protocol.smarts_generator_cli --batch literature_reactions.txt
```

## See Also

- [Protocol Database Documentation](../docs/PROTOCOL_DATABASE.md)
- [SMARTS Tutorial](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
- [RDKit SMARTS Documentation](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support)

## Support

For issues or questions:
1. Check if RDKit is properly installed: `smarts_generator.bat --check-rdkit`
2. Review error messages for SMILES parsing issues
3. Consult SMARTS syntax documentation for pattern writing

## License

Part of the Condition-agent project. See main LICENSE file.
