# SMARTS Pattern Visualization

## Overview

The SMARTS Generator includes powerful visualization capabilities to help verify that your reaction patterns match the intended substrates.

## Why Visualize?

SMARTS patterns can be complex and error-prone. Visualization helps you:

- **Catch mistakes early** - See if your pattern matches unexpected structures
- **Verify exclusions** - Confirm forbidden patterns actually block the right molecules  
- **Test scope** - Validate against real substrate examples
- **Communication** - Share visual patterns with colleagues

## Generated Images

The tool creates PNG images for:

### 1. Core Transformation
`core_transformation.png` - Comprehensive visualization showing:

**Top Section:** Complete reaction transformation with arrow
- Full reaction view showing reactant → product

**Bottom Section:** Individual pattern components side-by-side
- **Reactant Pattern (Blue):** Shows the substrate pattern with SMARTS text
- **Product Pattern (Green):** Shows the product pattern with SMARTS text

Example for `[C:1;H2,H3]-[I:2]>>[C:1]-[B:3]`:
- Complete view shows the full transformation
- Reactant section shows the C-I fragment structure
- Product section shows the C-B fragment structure
- SMARTS text displayed below each component for reference

### 2. Forbidden Patterns
`guard_forbid_N.png` - Visual of each forbidden substructure

Examples:
- `[CH]-I` → Shows secondary iodide structure
- `[C;H0]-I` → Shows tertiary iodide structure
- `[CH2]-[a]-I` → Shows benzylic iodide structure

### 3. Required Patterns
`guard_require_N.png` - Visual of each required substructure (if any)

## Usage Examples

### Example 1: Basic Visualization

```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI>>CCCBpin" \
  --visualize
```

Output:
```
✅ Saved reaction SMARTS visualization to smarts_visualizations\core_transformation.png
✅ Saved SMARTS visualization to smarts_visualizations\guard_forbid_1.png
✅ Saved SMARTS visualization to smarts_visualizations\guard_forbid_2.png
```

### Example 2: Test Against Molecules

```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI>>CCCBpin" \
  --visualize \
  --test-smiles "CCCCI,CC(C)I,CC(C)(C)I"
```

Output:
```
Testing Patterns Against Examples
----------------------------------------------------------------------
  ✅ MATCH: CCCCI                    # Primary iodide works
  ❌ NO MATCH: CC(C)I                # Secondary blocked
    → Matches forbidden pattern
  ❌ NO MATCH: CC(C)(C)I             # Tertiary blocked
    → Matches forbidden pattern
```

### Example 3: Custom Output Directory

```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCI>>CCCBpin" \
  --visualize \
  --viz-dir "protocol_images/borylation"
```

## Pattern Validation Workflow

### Step 1: Generate Initial Pattern

```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \
  --visualize \
  --output initial_pattern.json
```

### Step 2: Review Visual Output

Open the generated images and check:
- ✅ Does the core pattern show the right transformation?
- ✅ Do forbidden patterns match what you want to exclude?
- ✅ Are there any unexpected matches?

### Step 3: Test Against Known Substrates

Create a test file `test_substrates.txt`:
```
CCCCI           # Should work - primary
CC(C)I          # Should fail - secondary  
CC(C)(C)I       # Should fail - tertiary
c1ccccc1CI      # Should fail - benzylic
CC=CCI          # Should fail - allylic
```

Run validation:
```powershell
python -m chemtools.protocol.smarts_generator_cli \
  --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \
  --visualize \
  --test-smiles "CCCCI,CC(C)I,CC(C)(C)I,c1ccccc1CI,CC=CCI"
```

### Step 4: Refine Pattern

If tests reveal issues:
1. Run in interactive mode
2. Adjust core pattern or guards
3. Re-visualize and re-test
4. Repeat until all tests pass correctly

## Understanding Match Results

### ✅ MATCH
The molecule:
- Matches the core transformation pattern
- Does NOT match any forbidden patterns
- Matches ALL required patterns (if any)

### ❌ NO MATCH (Core pattern doesn't match)
The core transformation pattern didn't match the molecule structure.

**Action:** Check if core pattern is too specific or molecule is actually incompatible

### ❌ NO MATCH (Matches forbidden pattern)
The molecule matched one of the guard_forbid patterns.

**Action:** This is expected behavior if molecule should be excluded

### ❌ NO MATCH (Missing required pattern)
The molecule didn't match one of the guard_require patterns.

**Action:** Check if required pattern is too restrictive

## Real-World Example

### Alkyl Iodide Borylation Protocol

**Reaction:** Primary alkyl iodides → Alkyl boronic esters

**Generated Pattern:**
```json
{
  "core": "[C:1;H2,H3;X4;!$(C-[a]);!$(C-[C]=[C])]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])",
  "guards_forbid": [
    "[CH]-I",
    "[C;H0]-I",
    "[CH2]-[a]-I",
    "[CH2]-[C]=[C]-I"
  ]
}
```

**Test Results:**
```
✅ MATCH: CCCCI                     # n-Propyl iodide → OK
✅ MATCH: CCCCCCI                   # n-Hexyl iodide → OK
❌ NO MATCH: CC(C)I                # Isopropyl iodide → Blocked (secondary)
❌ NO MATCH: CC(C)(C)I             # tert-Butyl iodide → Blocked (tertiary)
❌ NO MATCH: c1ccccc1CI            # Benzyl iodide → Blocked (benzylic)
❌ NO MATCH: CC=CCI                # Allyl iodide → Blocked (allylic)
```

**Visual Verification:**
- `core_transformation.png` → Shows C-I → C-B transformation
- `guard_forbid_1.png` → Shows secondary iodide (blocked)
- `guard_forbid_2.png` → Shows tertiary iodide (blocked)
- `guard_forbid_3.png` → Shows benzylic iodide (blocked)
- `guard_forbid_4.png` → Shows allylic iodide (blocked)

## Tips for Effective Visualization

### 1. Start Broad, Then Narrow
- Generate initial pattern
- Visualize to see what it matches
- Add guards to exclude unwanted matches
- Re-visualize to confirm

### 2. Test Edge Cases
Always test borderline substrates:
- Minimum/maximum chain lengths
- Functional group tolerance
- Steric hindrance cases
- Electronic effects

### 3. Document Findings
Save visualization images with your protocol for future reference:
```
protocol_docs/
  alkyl_iodide_borylation/
    protocol.json
    core_transformation.png
    guard_patterns/
      forbid_secondary.png
      forbid_tertiary.png
      forbid_benzylic.png
```

### 4. Share With Team
- Include images in protocol documentation
- Use in training materials
- Reference in literature writeups

## Troubleshooting

### Issue: "RDKit not available"

**Solution:**
```powershell
pip install rdkit
# or
conda install -c conda-forge rdkit
```

### Issue: Blank or garbled images

**Cause:** Complex SMARTS pattern that RDKit can't depict

**Solution:**
- Simplify pattern
- Test parts separately
- Check SMARTS syntax

### Issue: Pattern matches too many/few molecules

**Solution:**
1. Visualize current pattern
2. Identify what's wrong in the visual
3. Adjust pattern specificity
4. Re-test with more examples

### Issue: Can't see small details in image

**Solution:** Increase image size in code or zoom in on PNG file

## Integration with Protocol Database

Generated visualizations can be stored alongside protocol JSON files:

```
data/protocol_db/
  Alkyl_Iodide_Borylation.json
  Alkyl_Iodide_Borylation_viz/
    core_transformation.png
    guard_forbid_1.png
    guard_forbid_2.png
```

Reference in documentation:
```markdown
## Scope

This protocol works for primary alkyl iodides only.

![Core Transformation](Alkyl_Iodide_Borylation_viz/core_transformation.png)

### Excluded Substrates

![Secondary Iodides](Alkyl_Iodide_Borylation_viz/guard_forbid_1.png)
![Tertiary Iodides](Alkyl_Iodide_Borylation_viz/guard_forbid_2.png)
```

## Advanced: Custom Visualization

For more control, use the visualization functions directly:

```python
from chemtools.protocol.smarts_generator_cli import (
    visualize_reaction_smarts,
    visualize_pattern_with_examples
)
from pathlib import Path

# Custom visualization
visualize_reaction_smarts(
    "[C:1]-[I:2]>>[C:1]-[B:3]",
    Path("my_custom_image.png"),
    img_size=(1200, 600)
)

# Test with examples
from chemtools.protocol.smarts_generator_cli import ReactionSmartsApplicability

pattern = ReactionSmartsApplicability(
    core="[C:1]-[I:2]>>[C:1]-[B:3]",
    guards_forbid=["[CH]-I"]
)

test_molecules = ["CCCCI", "CC(C)I"]
results = visualize_pattern_with_examples(
    pattern,
    test_molecules,
    Path("output_dir")
)

for smiles, matches in results.items():
    print(f"{smiles}: {'✅' if matches else '❌'}")
```

## Summary

Visualization is a powerful feature that:
- ✅ Catches pattern errors before they cause problems
- ✅ Validates scope against real substrates
- ✅ Improves pattern quality and accuracy
- ✅ Facilitates team communication
- ✅ Documents protocol limitations clearly

**Best Practice:** Always visualize and test your patterns before adding them to production protocol databases!
