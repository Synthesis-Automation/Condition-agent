"""
Demonstration of image generation capabilities using sample reactions.

This script showcases the rendering.py module's ability to generate
high-quality images for molecules and reactions from the sample database.
"""

import sys
from pathlib import Path

# Add the project root to the path
sys.path.insert(0, str(Path(__file__).parent))

from chemtools.visualization.rendering import (
    render_molecule_image,
    render_reaction_image,
    rdkit_available,
)


def demo_molecules():
    """Demonstrate molecule rendering."""
    print("\n" + "="*70)
    print("MOLECULE RENDERING DEMO")
    print("="*70 + "\n")
    
    molecules = [
        ("c1ccccc1", "benzene", "Simple aromatic ring"),
        ("CCO", "ethanol", "Simple alcohol"),
        ("c1ccc(Br)cc1", "bromobenzene", "Aryl halide"),
        ("Nc1ccccc1", "aniline", "Primary aromatic amine"),
        ("c1ccc(C(F)(F)F)cc1", "trifluorotoluene", "Electron-withdrawing group"),
        ("c1ccc(OC)cc1", "anisole", "Electron-donating group"),
        ("c1ccncc1", "pyridine", "Heterocycle"),
        ("c1ccc2[nH]ccc2c1", "indole", "Fused heterocycle"),
        ("c1ccc(B(O)O)cc1", "phenylboronic_acid", "Boronic acid"),
        ("C#Cc1ccccc1", "phenylacetylene", "Terminal alkyne"),
    ]
    
    output_dir = Path("demo_images/molecules")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating {len(molecules)} molecule images in: {output_dir}\n")
    
    for smiles, name, description in molecules:
        try:
            output_path = output_dir / f"{name}.png"
            render_molecule_image(
                smiles=smiles,
                output_path=output_path,
                size=(400, 300),
                legend=name.replace('_', ' ').title()
            )
            print(f"✓ {name:25s} - {description}")
        except Exception as e:
            print(f"✗ {name:25s} - Error: {e}")


def demo_reactions():
    """Demonstrate reaction rendering with key reaction types."""
    print("\n" + "="*70)
    print("REACTION RENDERING DEMO")
    print("="*70 + "\n")
    
    reactions = [
        # C-C Coupling
        ("Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1", 
         "suzuki_coupling", "Suzuki-Miyaura Coupling"),
        
        ("Brc1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1", 
         "sonogashira_coupling", "Sonogashira Coupling"),
        
        ("Brc1ccccc1.C=C>>C(=Cc1ccccc1)", 
         "heck_reaction", "Heck Reaction"),
        
        # C-N Coupling
        ("Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1", 
         "buchwald_hartwig", "Buchwald-Hartwig Amination"),
        
        ("Clc1ccncc1.NCC>>CCNc1ccncc1", 
         "buchwald_hartwig_pyridine", "B-H with Pyridine"),
        
        # Reduction
        ("C=Cc1ccccc1.[H][H]>>CCc1ccccc1", 
         "hydrogenation", "Catalytic Hydrogenation"),
        
        ("c1ccc(C=O)cc1.[BH4-].[Na+]>>c1ccc(CO)cc1", 
         "nabh4_reduction", "NaBH4 Reduction"),
        
        # Oxidation
        ("c1ccc(CO)cc1.[O]>>c1ccc(C=O)cc1", 
         "alcohol_oxidation", "Alcohol to Aldehyde"),
        
        # Cycloaddition
        ("C=CC=C.C=C>>C1C=CCC1", 
         "diels_alder", "Diels-Alder Cycloaddition"),
        
        ("C#CCO.c1ccc([N-][N+]#N)cc1>>c1ccccc1n1cc(CO)nn1", 
         "click_reaction", "Click Chemistry (CuAAC)"),
        
        # Carbonyl chemistry
        ("c1ccc(C=O)cc1.C[P+](c1ccccc1)(c1ccccc1)c1ccccc1.[base]>>C=Cc1ccccc1", 
         "wittig_reaction", "Wittig Olefination"),
        
        ("c1ccc(C=O)cc1.Nc1ccccc1>>c1ccc(CNc2ccccc2)cc1", 
         "reductive_amination", "Reductive Amination"),
    ]
    
    output_dir = Path("demo_images/reactions")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating {len(reactions)} reaction images in: {output_dir}\n")
    
    for smiles, name, description in reactions:
        try:
            output_path = output_dir / f"{name}.png"
            render_reaction_image(
                reaction_smiles=smiles,
                output_path=output_path,
                size=(960, 320)
            )
            print(f"✓ {name:30s} - {description}")
        except Exception as e:
            print(f"✗ {name:30s} - Error: {str(e)[:50]}")


def demo_svg_format():
    """Demonstrate SVG format generation."""
    print("\n" + "="*70)
    print("SVG FORMAT DEMO")
    print("="*70 + "\n")
    
    output_dir = Path("demo_images/svg")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating SVG images in: {output_dir}\n")
    
    # Molecule SVG
    try:
        mol_path = output_dir / "caffeine.svg"
        render_molecule_image(
            smiles="CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
            output_path=mol_path,
            size=(400, 300),
            image_format="svg",
            legend="Caffeine"
        )
        print(f"✓ Molecule SVG: caffeine.svg")
    except Exception as e:
        print(f"✗ Molecule SVG: Error: {e}")
    
    # Reaction SVG
    try:
        rxn_path = output_dir / "suzuki_reaction.svg"
        render_reaction_image(
            reaction_smiles="Brc1ccc(C(F)(F)F)cc1.c1ccc(B(O)O)cc1>>FC(F)(F)c1ccc(-c2ccccc2)cc1",
            output_path=rxn_path,
            size=(960, 320),
            image_format="svg"
        )
        print(f"✓ Reaction SVG: suzuki_reaction.svg")
    except Exception as e:
        print(f"✗ Reaction SVG: Error: {e}")


def demo_custom_sizes():
    """Demonstrate different image sizes."""
    print("\n" + "="*70)
    print("CUSTOM SIZE DEMO")
    print("="*70 + "\n")
    
    output_dir = Path("demo_images/sizes")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating images at different sizes in: {output_dir}\n")
    
    smiles = "c1ccc2c(c1)c1ccccc1[nH]2"  # Carbazole
    sizes = [
        (300, 200, "small"),
        (600, 400, "medium"),
        (900, 600, "large"),
    ]
    
    for width, height, label in sizes:
        try:
            output_path = output_dir / f"carbazole_{label}.png"
            render_molecule_image(
                smiles=smiles,
                output_path=output_path,
                size=(width, height),
                legend=f"Carbazole ({label})"
            )
            print(f"✓ {label:10s} ({width}×{height}) - carbazole_{label}.png")
        except Exception as e:
            print(f"✗ {label:10s} ({width}×{height}) - Error: {e}")


def main():
    """Run all demonstrations."""
    print("\n" + "="*70)
    print("IMAGE GENERATION DEMO")
    print("Testing chemtools.visualization.rendering module")
    print("="*70)
    
    # Check RDKit availability
    if not rdkit_available():
        print("\n❌ RDKit is not available. Cannot generate images.")
        print("Install RDKit to use image generation features:")
        print("  conda install -c conda-forge rdkit")
        return 1
    
    print("\n✓ RDKit is available - Ready to generate images")
    
    # Run all demonstrations
    demo_molecules()
    demo_reactions()
    demo_svg_format()
    demo_custom_sizes()
    
    # Summary
    print("\n" + "="*70)
    print("DEMO COMPLETE")
    print("="*70)
    print("\nGenerated images can be found in:")
    print("  - demo_images/molecules/  (10 molecular structures)")
    print("  - demo_images/reactions/  (12 reaction schemes)")
    print("  - demo_images/svg/        (2 SVG format examples)")
    print("  - demo_images/sizes/      (3 size variants)")
    print("\n✅ Image generation demo completed successfully!")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
