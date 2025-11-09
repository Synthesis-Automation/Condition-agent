# Image Generator Quick Start Guide

## Overview
The `chemtools.visualization.rendering` module provides high-quality image generation for molecules and reactions.

## Installation Requirements

```bash
# RDKit is required
conda install -c conda-forge rdkit
```

## Basic Usage

### 1. Render a Molecule

```python
from chemtools.visualization.rendering import render_molecule_image

# Simple molecule
render_molecule_image(
    smiles="c1ccccc1",
    output_path="benzene.png",
    size=(400, 300),
    legend="Benzene"
)
```

### 2. Render a Reaction

```python
from chemtools.visualization.rendering import render_reaction_image

# Suzuki coupling reaction
render_reaction_image(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    output_path="suzuki_reaction.png",
    size=(960, 320)
)
```

### 3. Generate SVG Format

```python
# For molecules
render_molecule_image(
    smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    output_path="caffeine.svg",
    image_format="svg",
    size=(400, 300),
    legend="Caffeine"
)

# For reactions
render_reaction_image(
    reaction_smiles="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    output_path="buchwald_hartwig.svg",
    image_format="svg",
    size=(960, 320)
)
```

## Function Reference

### `render_molecule_image()`

**Parameters:**
- `smiles` (str): Molecular SMILES string
- `output_path` (str | Path): Output file path
- `size` (tuple[int, int]): Image size in pixels (width, height)
  - Default: `(480, 300)`
- `image_format` (str): "png" or "svg"
  - Default: `"png"`
- `kekulize` (bool): Whether to kekulize aromatic rings
  - Default: `True`
- `legend` (str | None): Optional caption below molecule
  - Default: `None`

**Returns:**
- `Path`: Path to the generated image file

**Example:**
```python
render_molecule_image(
    smiles="c1ccc(C(F)(F)F)cc1",
    output_path="trifluorotoluene.png",
    size=(600, 400),
    image_format="png",
    legend="4-Trifluorotoluene"
)
```

### `render_reaction_image()`

**Parameters:**
- `reaction_smiles` (str): Reaction SMILES (reactants>agents>products)
- `output_path` (str | Path): Output file path
- `size` (tuple[int, int]): Image size in pixels (width, height)
  - Default: `(960, 320)`
- `image_format` (str): "png" or "svg"
  - Default: `"png"`
- `kekulize` (bool): Whether to kekulize aromatic rings
  - Default: `True`

**Returns:**
- `Path`: Path to the generated image file

**Example:**
```python
render_reaction_image(
    reaction_smiles="Clc1ccncc1.NCC>>CCNc1ccncc1",
    output_path="cn_coupling.png",
    size=(1200, 400),
    image_format="png"
)
```

## Reaction SMILES Format

### Standard Format
```
reactants>agents>products
```

### Alternative Format (no agents)
```
reactants>>products
```

### Examples

**With agents (catalysts, bases, etc.):**
```python
"Brc1ccccc1.c1ccc(B(O)O)cc1>[Pd].[K2CO3]>c1ccc(-c2ccccc2)cc1"
```

**Without agents (agents moved to reactants):**
```python
"Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"
```

**Multiple reactants (separated by `.`):**
```python
"Brc1ccccc1.c1ccc(B(O)O)cc1.c1ccc(OC)cc1>>product"
```

## Size Recommendations

### Molecules
- **Small icons**: `(200, 150)`
- **Standard**: `(400, 300)` - Default
- **Large**: `(600, 450)`
- **Poster**: `(800, 600)`

### Reactions
- **Compact**: `(640, 240)`
- **Standard**: `(960, 320)` - Default
- **Wide**: `(1200, 400)`
- **Publication**: `(1600, 480)`

## Styling Features

### Automatic Features
- ✅ Kekulization of aromatic systems
- ✅ Explicit hydrogens on heteroatoms
- ✅ Clean bond rendering (2.5 pt width)
- ✅ White background
- ✅ Proper padding (4% of canvas)
- ✅ Reaction arrow between reactants and products
- ✅ Multi-component layout (multiple reactants/products)

### Style Constants
```python
DEFAULT_MOLECULE_SIZE = (480, 300)
DEFAULT_REACTION_SIZE = (960, 320)
SUPPORTED_FORMATS = {"png", "svg"}
```

## Error Handling

```python
from chemtools.visualization.rendering import rdkit_available

# Check RDKit availability
if not rdkit_available():
    print("RDKit is not installed")
    exit(1)

# Handle invalid SMILES
try:
    render_molecule_image(
        smiles="invalid_smiles",
        output_path="test.png"
    )
except ValueError as e:
    print(f"Invalid SMILES: {e}")

# Handle invalid reaction SMILES
try:
    render_reaction_image(
        reaction_smiles=">>invalid",
        output_path="test.png"
    )
except ValueError as e:
    print(f"Invalid reaction: {e}")
```

## Batch Processing Example

```python
from pathlib import Path
from chemtools.visualization.rendering import render_molecule_image

# List of molecules to render
molecules = {
    "benzene": "c1ccccc1",
    "pyridine": "c1ccncc1",
    "indole": "c1ccc2[nH]ccc2c1",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
}

output_dir = Path("molecules")
output_dir.mkdir(exist_ok=True)

for name, smiles in molecules.items():
    render_molecule_image(
        smiles=smiles,
        output_path=output_dir / f"{name}.png",
        legend=name.title()
    )
    print(f"✓ Generated {name}.png")
```

## Integration with Sample Reactions

```python
from tests.sample_reactions import SAMPLE_REACTIONS
from chemtools.visualization.rendering import render_reaction_image

# Get Suzuki reactions
suzuki_reactions = [
    rxn for rxn in SAMPLE_REACTIONS[1:]  # Skip label
    if "Suzuki" in rxn
]

# Render first 5 Suzuki reactions
for i, rxn_str in enumerate(suzuki_reactions[:5], 1):
    smiles = rxn_str.split("(")[0].strip()
    description = rxn_str.split("(")[1].rstrip(")")
    
    output_path = f"suzuki_{i:02d}.png"
    render_reaction_image(
        reaction_smiles=smiles,
        output_path=output_path
    )
    print(f"✓ {description}")
```

## Advanced Usage

### Custom Directory Creation
```python
from pathlib import Path

def render_categorized_reactions(reactions_dict):
    """Render reactions organized by category."""
    for category, reactions in reactions_dict.items():
        output_dir = Path("reactions") / category
        output_dir.mkdir(parents=True, exist_ok=True)
        
        for i, (name, smiles) in enumerate(reactions, 1):
            render_reaction_image(
                reaction_smiles=smiles,
                output_path=output_dir / f"{name}.png"
            )
```

### Quality Settings
```python
# High-quality publication figure
render_reaction_image(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    output_path="publication_figure.png",
    size=(1600, 480),  # Large size
    image_format="png",
    kekulize=True      # Clean aromatic representation
)

# Vector graphic for scaling
render_reaction_image(
    reaction_smiles="Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    output_path="vector_figure.svg",
    size=(1200, 400),
    image_format="svg"  # Scalable
)
```

## Troubleshooting

### Issue: "RDKit is required for rendering"
**Solution:** Install RDKit
```bash
conda install -c conda-forge rdkit
```

### Issue: "Invalid SMILES string"
**Solution:** Verify SMILES syntax
- Check for typos
- Ensure valid atom symbols
- Check bracket notation for atoms with charges or explicit hydrogens
- Use a SMILES validator

### Issue: "Invalid reaction SMILES string"
**Solution:** Check reaction format
- Must contain `>>` separator
- Format: `reactants>>products` or `reactants>agents>products`
- Multiple components separated by `.`

### Issue: "Images look small/blurry"
**Solution:** Increase image size
```python
# Instead of default (960, 320)
render_reaction_image(
    reaction_smiles=smiles,
    output_path=path,
    size=(1600, 480)  # Larger size
)
```

### Issue: "Directory not found"
**Solution:** Directories are created automatically
```python
# This is handled automatically
render_molecule_image(
    smiles="c1ccccc1",
    output_path="deep/nested/path/benzene.png"
)
# Creates: deep/nested/path/ automatically
```

## Performance Tips

1. **Batch processing**: Process multiple structures in a loop
2. **Format selection**: Use PNG for raster, SVG for vector graphics
3. **Size optimization**: Choose appropriate size for use case
4. **Caching**: Pre-generate commonly used structures

## Use Cases

### Research Documentation
```python
# Document a synthesis route
steps = [
    ("Step 1", "Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1"),
    ("Step 2", "c1ccc(-c2ccccc2)cc1.[H][H]>>c1ccc(CCc2ccccc2)cc1"),
]

for name, smiles in steps:
    render_reaction_image(smiles, f"{name}.png")
```

### Web Application
```python
from flask import Flask, send_file
from chemtools.visualization.rendering import render_molecule_image

app = Flask(__name__)

@app.route('/molecule/<smiles>')
def get_molecule(smiles):
    output_path = f"temp/{smiles}.png"
    render_molecule_image(smiles, output_path)
    return send_file(output_path)
```

### Automated Reports
```python
def generate_report(reactions):
    """Generate a report with reaction images."""
    for i, rxn in enumerate(reactions, 1):
        img_path = f"report/reaction_{i}.png"
        render_reaction_image(rxn['smiles'], img_path)
        print(f"Added reaction {i} to report")
```

## See Also

- `tests/sample_reactions.py` - Comprehensive reaction database
- `chemtools/util/rdkit_helpers.py` - RDKit utility functions
- RDKit documentation: https://www.rdkit.org/docs/
