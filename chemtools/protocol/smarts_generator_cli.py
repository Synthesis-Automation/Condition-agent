"""CLI tool for generating reaction SMARTS applicability patterns from reaction SMILES.

This tool helps define the scope and limitations of chemical protocols by creating
structured SMARTS patterns that specify:
1. Core transformation pattern (reactant→product with atom mapping)
2. Guard patterns that forbid specific structural features

Usage:
    python -m chemtools.protocol.smarts_generator_cli --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin"
    python -m chemtools.protocol.smarts_generator_cli --interactive
    python -m chemtools.protocol.smarts_generator_cli --batch input.txt --output output.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field

# Add parent directory to path if running as script (not as module)
if __name__ == "__main__" and __package__ is None:
    _script_dir = Path(__file__).resolve().parent
    _project_root = _script_dir.parent.parent
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))

try:
    from rdkit import Chem  # type: ignore[import-not-found]
    from rdkit.Chem import AllChem, Draw  # type: ignore[import-not-found]
    from rdkit.Chem.Draw import rdMolDraw2D  # type: ignore[import-not-found]
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    Chem = None  # type: ignore[assignment]

# Import our reusable chemistry utilities
from chemtools.util.substrate_classifier import SubstrateClassifier, classify_substrate
from chemtools.util.smarts_builders import SmartsBuilder, build_smarts_with_guards


# ============================================================================
# Data Models
# ============================================================================

@dataclass
class ReactionSmartsApplicability:
    """Structured representation of reaction SMARTS applicability pattern."""
    
    core: str  # Core SMARTS transformation with atom mapping
    guards_forbid: List[str] = field(default_factory=list)  # SMARTS patterns that invalidate the reaction
    guards_require: List[str] = field(default_factory=list)  # SMARTS patterns that must be present
    notes: str = ""  # Human-readable description
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary format matching protocol JSON schema."""
        result = {"core": self.core}
        if self.guards_forbid:
            result["guards_forbid"] = self.guards_forbid
        if self.guards_require:
            result["guards_require"] = self.guards_require
        if self.notes:
            result["notes"] = self.notes
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> ReactionSmartsApplicability:
        """Create from dictionary."""
        return cls(
            core=data["core"],
            guards_forbid=data.get("guards_forbid", []),
            guards_require=data.get("guards_require", []),
            notes=data.get("notes", "")
        )


# ============================================================================
# SMARTS Generation Logic
# ============================================================================

class SmartsGenerator:
    """Generate reaction SMARTS patterns from reaction SMILES."""
    
    def __init__(self, reaction_smiles: str):
        """Initialize with reaction SMILES.
        
        Args:
            reaction_smiles: Reaction in SMILES format (reactants>>products or reactants>reagents>products)
        """
        self.reaction_smiles = reaction_smiles.strip()
        self.reactants_smiles = ""
        self.products_smiles = ""
        self.reagents_smiles = ""
        self._parse_reaction()
        
    def _parse_reaction(self):
        """Parse reaction SMILES into components."""
        # Check for invalid triple arrow
        if ">>>" in self.reaction_smiles:
            raise ValueError(f"Invalid reaction SMILES format (too many arrows): {self.reaction_smiles}")
        
        # Handle both >> and > separators
        if ">>" in self.reaction_smiles:
            parts = self.reaction_smiles.split(">>")
            if len(parts) != 2:
                raise ValueError(f"Invalid reaction SMILES format: {self.reaction_smiles}")
            self.reactants_smiles = parts[0].strip()
            self.products_smiles = parts[1].strip()
        elif ">" in self.reaction_smiles:
            parts = self.reaction_smiles.split(">")
            if len(parts) == 2:
                self.reactants_smiles = parts[0].strip()
                self.products_smiles = parts[1].strip()
            elif len(parts) == 3:
                self.reactants_smiles = parts[0].strip()
                self.reagents_smiles = parts[1].strip()
                self.products_smiles = parts[2].strip()
            else:
                raise ValueError(f"Invalid reaction SMILES format: {self.reaction_smiles}")
        else:
            raise ValueError(f"No reaction arrow found in: {self.reaction_smiles}")
    
    def generate_core_smarts(self, 
                            reactant_pattern: Optional[str] = None,
                            product_pattern: Optional[str] = None,
                            atom_mapping: Optional[Dict[int, int]] = None) -> str:
        """Generate core SMARTS transformation pattern.
        
        Args:
            reactant_pattern: Optional custom SMARTS pattern for reactant side
            product_pattern: Optional custom SMARTS pattern for product side
            atom_mapping: Optional atom mapping dict (reactant_idx -> product_idx)
            
        Returns:
            SMARTS string with format "[reactant_pattern]>>[product_pattern]"
        """
        if not RDKIT_AVAILABLE:
            # Fallback: Use provided patterns or simple conversion
            if reactant_pattern and product_pattern:
                return f"{reactant_pattern}>>{product_pattern}"
            raise RuntimeError("RDKit required for automatic SMARTS generation. Install with: pip install rdkit")
        
        # If patterns provided, use them directly
        if reactant_pattern and product_pattern:
            return f"{reactant_pattern}>>{product_pattern}"
        
        # Parse molecules
        reactant_mols = [Chem.MolFromSmiles(smi.strip()) for smi in self.reactants_smiles.split('.') if smi.strip()]
        product_mols = [Chem.MolFromSmiles(smi.strip()) for smi in self.products_smiles.split('.') if smi.strip()]
        
        if None in reactant_mols:
            raise ValueError(f"Invalid SMILES in reactants: {self.reactants_smiles}")
        if None in product_mols:
            raise ValueError(f"Invalid SMILES in products: {self.products_smiles}")
        
        # Generate generic patterns focusing on functional groups
        # For substrates (first reactant), focus on the reactive site
        reactant_smarts = self._mol_to_generic_smarts(reactant_mols[0])
        
        # For products, also focus on the functional group
        product_smarts = self._mol_to_generic_smarts(product_mols[0])
        
        return f"{reactant_smarts}>>{product_smarts}"
    
    def _mol_to_generic_smarts(self, mol: Any, add_mapping: bool = False) -> str:  # mol is Chem.Mol when rdkit available
        """Convert molecule to chemistry-aware SMARTS pattern.
        
        Uses SubstrateClassifier and SmartsBuilder for context-aware pattern generation.
        This creates patterns that reflect actual chemical context (primary vs secondary,
        aryl vs alkyl, aniline vs aliphatic amine, etc.).
        
        Users should refine the generated pattern to add:
        - Atom mapping (:1, :2, etc.)
        - Additional constraints if needed
        
        Returns:
            Chemistry-aware SMARTS pattern
        """
        # Convert RDKit mol to SMILES
        smiles = Chem.MolToSmiles(mol)
        
        # Use SmartsBuilder for chemistry-aware pattern generation
        builder = SmartsBuilder()
        try:
            smarts = builder.build_from_smiles(smiles)
            return smarts
        except Exception as e:
            # Fallback: use full SMARTS if builder fails
            print(f"  ⚠ Warning: SmartsBuilder failed ({e}), using fallback")
            return Chem.MolToSmarts(mol)
    
    def suggest_guard_patterns(self, verbose: bool = False) -> List[str]:
        """Suggest context-aware guard patterns based on substrate classification.
        
        Uses SubstrateClassifier and SmartsBuilder to generate chemistry-aware
        exclusion patterns based on the actual substrate type detected.
        
        Returns:
            List of SMARTS patterns that should be forbidden
        """
        if not RDKIT_AVAILABLE:
            return []
        
        all_guards = []
        
        # Analyze each reactant
        reactant_smiles_list = [smi.strip() for smi in self.reactants_smiles.split('.') if smi.strip()]
        
        for smiles in reactant_smiles_list:
            try:
                # Use build_smarts_with_guards for context-aware guard generation
                result = build_smarts_with_guards(smiles)
                
                if result.get('guards_forbid'):
                    all_guards.extend(result['guards_forbid'])
                    
                    if verbose:
                        print(f"  ℹ Substrate detected as: {result['substrate_class']} ({result['substrate_family']})")
                        print(f"    Generated {len(result['guards_forbid'])} guard patterns")
                        
            except Exception as e:
                if verbose:
                    print(f"  ⚠ Warning: Could not analyze {smiles}: {e}")
                continue
        
        return all_guards
    
    def generate_interactive(self) -> ReactionSmartsApplicability:
        """Interactive mode to build SMARTS pattern with user guidance."""
        print("\n" + "="*70)
        print("  Reaction SMARTS Applicability Pattern Generator")
        print("="*70)
        print(f"\nReaction SMILES: {self.reaction_smiles}")
        print(f"  Reactants: {self.reactants_smiles}")
        print(f"  Products:  {self.products_smiles}")
        if self.reagents_smiles:
            print(f"  Reagents:  {self.reagents_smiles}")
        print()
        
        # Guide user through pattern creation
        print("Step 1: Core Transformation Pattern")
        print("-" * 70)
        print("Define the core SMARTS pattern with atom mapping.")
        print("Example: [C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])")
        print()
        print("Tips:")
        print("  - Use :1, :2, :3 for atom mapping")
        print("  - H2,H3 = 2 or 3 hydrogens")
        print("  - X4 = 4 total connections")
        print("  - !$(...) = NOT matching subpattern")
        print()
        
        default_core = self.generate_core_smarts()
        print(f"Auto-generated starting point:\n  {default_core}\n")
        
        core_pattern = input("Enter core SMARTS pattern (or press Enter to use above): ").strip()
        if not core_pattern:
            core_pattern = default_core
        
        # Guard patterns (forbid)
        print("\nStep 2: Guard Patterns - Forbidden Structures")
        print("-" * 70)
        print("Define SMARTS patterns that INVALIDATE the reaction.")
        print("These are structural features where the reaction does NOT work.\n")
        
        suggested_guards = self.suggest_guard_patterns(verbose=True)
        if suggested_guards:
            print("Suggestions based on your reaction:")
            for i, guard in enumerate(suggested_guards, 1):
                print(f"  {i}. {guard}")
            print()
        
        guards_forbid = []
        print("Enter forbidden patterns one per line (empty line to finish):")
        while True:
            pattern = input("  Forbid: ").strip()
            if not pattern:
                break
            guards_forbid.append(pattern)
        
        # Guard patterns (require)
        print("\nStep 3: Guard Patterns - Required Structures (Optional)")
        print("-" * 70)
        print("Define SMARTS patterns that MUST be present (optional).\n")
        
        guards_require = []
        print("Enter required patterns one per line (empty line to finish):")
        while True:
            pattern = input("  Require: ").strip()
            if not pattern:
                break
            guards_require.append(pattern)
        
        # Notes
        print("\nStep 4: Notes (Optional)")
        print("-" * 70)
        notes = input("Enter description of applicability scope: ").strip()
        
        result = ReactionSmartsApplicability(
            core=core_pattern,
            guards_forbid=guards_forbid,
            guards_require=guards_require,
            notes=notes
        )
        
        print("\n" + "="*70)
        print("  Generated Pattern")
        print("="*70)
        print(json.dumps(result.to_dict(), indent=2))
        print()
        
        return result


# ============================================================================
# Visualization Functions
# ============================================================================

def visualize_smarts_pattern(pattern_str: str, output_path: Optional[Path] = None, 
                            label: str = "", img_size: tuple = (800, 400)) -> bool:
    """Visualize a SMARTS pattern by drawing example molecules that match it.
    
    Args:
        pattern_str: SMARTS pattern string
        output_path: Path to save image (PNG format)
        label: Label to add to the image
        img_size: Image size as (width, height)
        
    Returns:
        True if visualization succeeded, False otherwise
    """
    if not RDKIT_AVAILABLE:
        print("⚠️  RDKit not available - cannot generate visualization")
        return False
    
    try:
        # Try to parse the SMARTS pattern
        pattern = Chem.MolFromSmarts(pattern_str)
        if pattern is None:
            print(f"⚠️  Invalid SMARTS pattern: {pattern_str}")
            return False
        
        # Create a drawing of the SMARTS pattern
        drawer = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])
        drawer.DrawMolecule(pattern)
        drawer.FinishDrawing()
        
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'wb') as f:
                f.write(drawer.GetDrawingText())
            print(f"✅ Saved SMARTS visualization to {output_path}")
        
        return True
        
    except Exception as e:
        print(f"⚠️  Error visualizing SMARTS: {e}")
        return False


def visualize_reaction_smarts(core_smarts: str, output_path: Optional[Path] = None,
                              img_size: tuple = (1200, 500)) -> bool:
    """Visualize a reaction SMARTS pattern showing both full reaction and individual patterns.
    
    Args:
        core_smarts: Reaction SMARTS string (format: reactant>>product)
        output_path: Path to save image (PNG format)
        img_size: Image size as (width, height)
        
    Returns:
        True if visualization succeeded, False otherwise
    """
    if not RDKIT_AVAILABLE:
        print("⚠️  RDKit not available - cannot generate visualization")
        return False
    
    try:
        from PIL import Image, ImageDraw, ImageFont
        import io
        
        # Parse reaction SMARTS and split into reactant/product
        if ">>" not in core_smarts:
            print(f"⚠️  Invalid reaction SMARTS format (missing >>): {core_smarts}")
            return False
        
        reactant_smarts, product_smarts = core_smarts.split(">>")
        
        # Parse reaction
        rxn = AllChem.ReactionFromSmarts(core_smarts)
        if rxn is None:
            print(f"⚠️  Invalid reaction SMARTS: {core_smarts}")
            return False
        
        # Create images for each component
        images = []
        labels = []
        
        # 1. Full reaction view
        drawer = rdMolDraw2D.MolDraw2DCairo(1000, 300)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        img_bytes = drawer.GetDrawingText()
        img_reaction = Image.open(io.BytesIO(img_bytes))
        images.append(img_reaction)
        labels.append("Complete Reaction Transformation")
        
        # 2. Reactant pattern
        reactant_mol = Chem.MolFromSmarts(reactant_smarts)
        if reactant_mol:
            drawer = rdMolDraw2D.MolDraw2DCairo(500, 300)
            drawer.DrawMolecule(reactant_mol)
            drawer.FinishDrawing()
            img_bytes = drawer.GetDrawingText()
            img_reactant = Image.open(io.BytesIO(img_bytes))
            images.append(img_reactant)
            labels.append(f"Reactant Pattern\n{reactant_smarts}")
        
        # 3. Product pattern
        product_mol = Chem.MolFromSmarts(product_smarts)
        if product_mol:
            drawer = rdMolDraw2D.MolDraw2DCairo(500, 300)
            drawer.DrawMolecule(product_mol)
            drawer.FinishDrawing()
            img_bytes = drawer.GetDrawingText()
            img_product = Image.open(io.BytesIO(img_bytes))
            images.append(img_product)
            labels.append(f"Product Pattern\n{product_smarts}")
        
        # Combine images into a single figure
        total_width = max(img_reaction.width, img_reactant.width + img_product.width + 20)
        total_height = img_reaction.height + max(img_reactant.height, img_product.height) + 120
        
        combined = Image.new('RGB', (total_width, total_height), 'white')
        draw = ImageDraw.Draw(combined)
        
        # Try to use a nice font, fallback to default
        try:
            font_large = ImageFont.truetype("arial.ttf", 16)
            font_small = ImageFont.truetype("arial.ttf", 12)
        except:
            font_large = ImageFont.load_default()
            font_small = ImageFont.load_default()
        
        # Add title
        draw.text((10, 10), "Reaction SMARTS Pattern Visualization", fill='black', font=font_large)
        
        # Add full reaction
        y_offset = 40
        draw.text((10, y_offset), labels[0], fill='black', font=font_large)
        y_offset += 25
        combined.paste(img_reaction, ((total_width - img_reaction.width) // 2, y_offset))
        y_offset += img_reaction.height + 20
        
        # Add individual patterns side by side
        draw.text((10, y_offset), "Individual Pattern Components:", fill='black', font=font_large)
        y_offset += 25
        
        x_offset = 20
        if len(images) > 1:
            # Reactant
            label_lines = labels[1].split('\n')
            draw.text((x_offset, y_offset), label_lines[0], fill='blue', font=font_large)
            if len(label_lines) > 1:
                draw.text((x_offset, y_offset + 20), label_lines[1], fill='gray', font=font_small)
            combined.paste(images[1], (x_offset, y_offset + 40))
            x_offset += images[1].width + 40
        
        if len(images) > 2:
            # Product
            label_lines = labels[2].split('\n')
            draw.text((x_offset, y_offset), label_lines[0], fill='green', font=font_large)
            if len(label_lines) > 1:
                draw.text((x_offset, y_offset + 20), label_lines[1], fill='gray', font=font_small)
            combined.paste(images[2], (x_offset, y_offset + 40))
        
        # Save combined image
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            combined.save(output_path, 'PNG')
            print(f"✅ Saved reaction SMARTS visualization to {output_path}")
        
        return True
        
    except ImportError:
        print("⚠️  PIL/Pillow not available - using simple reaction drawing")
        # Fallback to simple reaction drawing
        try:
            rxn = AllChem.ReactionFromSmarts(core_smarts)
            if rxn is None:
                return False
            
            drawer = rdMolDraw2D.MolDraw2DCairo(1000, 400)
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            
            if output_path:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                with open(output_path, 'wb') as f:
                    f.write(drawer.GetDrawingText())
                print(f"✅ Saved reaction SMARTS visualization to {output_path}")
            
            return True
        except Exception as e:
            print(f"⚠️  Error visualizing reaction SMARTS: {e}")
            return False
        
    except Exception as e:
        print(f"⚠️  Error visualizing reaction SMARTS: {e}")
        return False


def visualize_pattern_with_examples(applicability: ReactionSmartsApplicability,
                                    test_smiles: List[str],
                                    output_dir: Path) -> Dict[str, bool]:
    """Visualize SMARTS patterns and test against example molecules.
    
    Args:
        applicability: ReactionSmartsApplicability object with patterns
        test_smiles: List of SMILES strings to test matching
        output_dir: Directory to save visualization images
        
    Returns:
        Dictionary mapping test SMILES to match results
    """
    if not RDKIT_AVAILABLE:
        print("⚠️  RDKit not available - cannot generate visualizations")
        return {}
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Visualize the core reaction transformation
    print("\n" + "="*70)
    print("  Generating Visualizations")
    print("="*70)
    
    core_path = output_dir / "core_transformation.png"
    visualize_reaction_smarts(applicability.core, core_path)
    
    # Visualize guard patterns
    for i, guard in enumerate(applicability.guards_forbid, 1):
        # Clean up the pattern (remove comments)
        clean_pattern = guard.split('#')[0].strip()
        guard_path = output_dir / f"guard_forbid_{i}.png"
        visualize_smarts_pattern(clean_pattern, guard_path, f"Forbidden: {clean_pattern}")
    
    for i, guard in enumerate(applicability.guards_require, 1):
        clean_pattern = guard.split('#')[0].strip()
        guard_path = output_dir / f"guard_require_{i}.png"
        visualize_smarts_pattern(clean_pattern, guard_path, f"Required: {clean_pattern}")
    
    # Test patterns against example molecules
    results = {}
    if test_smiles:
        print("\n" + "-"*70)
        print("  Testing Patterns Against Examples")
        print("-"*70)
        
        # Parse the reaction SMARTS to get reactant pattern
        if ">>" in applicability.core:
            reactant_pattern = applicability.core.split(">>")[0]
            pattern_mol = Chem.MolFromSmarts(reactant_pattern)
            
            for smi in test_smiles:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    print(f"  ⚠️  Invalid SMILES: {smi}")
                    results[smi] = False
                    continue
                
                # Check if molecule matches core pattern
                core_match = mol.HasSubstructMatch(pattern_mol) if pattern_mol else False
                
                # Check forbidden patterns
                forbidden_match = False
                for guard in applicability.guards_forbid:
                    clean_guard = guard.split('#')[0].strip()
                    guard_mol = Chem.MolFromSmarts(clean_guard)
                    if guard_mol and mol.HasSubstructMatch(guard_mol):
                        forbidden_match = True
                        break
                
                # Check required patterns
                required_match = True
                for guard in applicability.guards_require:
                    clean_guard = guard.split('#')[0].strip()
                    guard_mol = Chem.MolFromSmarts(clean_guard)
                    if guard_mol and not mol.HasSubstructMatch(guard_mol):
                        required_match = False
                        break
                
                # Overall match: core matches AND not forbidden AND all required
                matches = core_match and not forbidden_match and required_match
                results[smi] = matches
                
                status = "✅ MATCH" if matches else "❌ NO MATCH"
                print(f"  {status}: {smi}")
                if not matches:
                    if not core_match:
                        print(f"    → Core pattern doesn't match")
                    if forbidden_match:
                        print(f"    → Matches forbidden pattern")
                    if not required_match:
                        print(f"    → Missing required pattern")
    
    print()
    return results


# ============================================================================
# CLI Interface
# ============================================================================

def run_interactive_mode():
    """Run interactive SMARTS pattern generator."""
    print("\n🧪 Reaction SMARTS Generator - Interactive Mode\n")
    
    while True:
        reaction_smiles = input("Enter reaction SMILES (or 'quit' to exit): ").strip()
        if reaction_smiles.lower() in ('quit', 'exit', 'q'):
            break
        
        if not reaction_smiles:
            print("❌ Please enter a reaction SMILES\n")
            continue
        
        try:
            generator = SmartsGenerator(reaction_smiles)
            result = generator.generate_interactive()
            
            # Ask about visualization
            if RDKIT_AVAILABLE:
                visualize = input("\nGenerate visualization images? (Y/n): ").strip().lower()
                if visualize != 'n':
                    viz_dir = Path("smarts_visualizations")
                    visualize_pattern_with_examples(result, [], viz_dir)
                    
                    # Ask for test molecules
                    test_molecules = input("\nEnter test SMILES to validate (comma-separated, or press Enter to skip): ").strip()
                    if test_molecules:
                        test_list = [s.strip() for s in test_molecules.split(',') if s.strip()]
                        visualize_pattern_with_examples(result, test_list, viz_dir)
            
            # Ask if user wants to save
            save = input("\nSave to file? (y/N): ").strip().lower()
            if save == 'y':
                filename = input("Enter filename (e.g., my_protocol.json): ").strip()
                if filename:
                    output_path = Path(filename)
                    output_path.write_text(json.dumps(result.to_dict(), indent=2), encoding='utf-8')
                    print(f"✅ Saved to {output_path}\n")
        
        except Exception as e:
            print(f"❌ Error: {e}\n")
        
        another = input("Generate another pattern? (Y/n): ").strip().lower()
        if another == 'n':
            break
    
    print("👋 Goodbye!\n")


def run_batch_mode(input_file: Path, output_file: Optional[Path] = None):
    """Process multiple reactions from file."""
    if not input_file.exists():
        print(f"❌ Input file not found: {input_file}")
        return 1
    
    lines = input_file.read_text(encoding='utf-8').strip().split('\n')
    results = []
    
    print(f"\n🧪 Processing {len(lines)} reactions from {input_file}\n")
    
    for i, line in enumerate(lines, 1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        print(f"  [{i}/{len(lines)}] {line[:50]}...")
        try:
            generator = SmartsGenerator(line)
            core = generator.generate_core_smarts()
            guards = generator.suggest_guard_patterns()
            
            result = ReactionSmartsApplicability(
                core=core,
                guards_forbid=guards,
                notes=f"Auto-generated from: {line}"
            )
            results.append({
                "reaction_smiles": line,
                "smarts_applicability": result.to_dict()
            })
        except Exception as e:
            print(f"    ⚠️  Error: {e}")
            results.append({
                "reaction_smiles": line,
                "error": str(e)
            })
    
    # Save results
    if output_file:
        output_path = output_file
    else:
        output_path = input_file.with_suffix('.json')
    
    output_path.write_text(json.dumps(results, indent=2), encoding='utf-8')
    print(f"\n✅ Saved {len(results)} results to {output_path}\n")
    
    return 0


def run_single_reaction(reaction_smiles: str, output_file: Optional[Path] = None, 
                       visualize: bool = False, viz_dir: Optional[Path] = None):
    """Process single reaction and print result."""
    try:
        generator = SmartsGenerator(reaction_smiles)
        
        print(f"\n🧪 Analyzing reaction: {reaction_smiles}")
        
        # Show substrate classification for each reactant
        print("\n📊 Substrate Analysis:")
        print("-" * 70)
        reactant_smiles_list = [smi.strip() for smi in generator.reactants_smiles.split('.') if smi.strip()]
        
        for i, smiles in enumerate(reactant_smiles_list, 1):
            try:
                info = classify_substrate(smiles)
                print(f"  Reactant {i}: {smiles}")
                print(f"    └─ Class:  {info.substrate_class}")
                print(f"    └─ Family: {info.substrate_family}")
                if info.special_positions.benzylic:
                    print(f"    └─ Special: Benzylic position(s) detected")
                if info.special_positions.allylic:
                    print(f"    └─ Special: Allylic position(s) detected")
                if info.special_positions.propargylic:
                    print(f"    └─ Special: Propargylic position(s) detected")
            except Exception as e:
                print(f"  Reactant {i}: {smiles} (classification failed: {e})")
        
        print()
        
        core = generator.generate_core_smarts()
        guards = generator.suggest_guard_patterns(verbose=False)
        
        result = ReactionSmartsApplicability(
            core=core,
            guards_forbid=guards,
            notes="Auto-generated pattern - please review and refine"
        )
        
        print("\n" + "="*70)
        print("  Generated SMARTS Applicability Pattern")
        print("="*70)
        output = result.to_dict()
        print(json.dumps(output, indent=2))
        print("\n⚠️  Note: This is a starting point. Please refine the pattern to match")
        print("    the exact scope of your protocol, especially:")
        print("    - Add atom mapping (:1, :2, :3)")
        print("    - Add specific atom properties (H2, H3, X4)")
        print("    - Add negation patterns (!$(...)) for excluded structures")
        print()
        
        # Generate visualization if requested
        if visualize and RDKIT_AVAILABLE:
            if viz_dir is None:
                viz_dir = Path("smarts_visualizations")
            visualize_pattern_with_examples(result, [], viz_dir)
        
        if output_file:
            output_file.write_text(json.dumps(output, indent=2), encoding='utf-8')
            print(f"✅ Saved to {output_file}\n")
        
        return 0
    
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1
        print("    the exact scope of your protocol, especially:")
        print("    - Add atom mapping (:1, :2, :3)")
        print("    - Add specific atom properties (H2, H3, X4)")
        print("    - Add negation patterns (!$(...)) for excluded structures")
        print()
        
        if output_file:
            output_file.write_text(json.dumps(output, indent=2), encoding='utf-8')
            print(f"✅ Saved to {output_file}\n")
        
        return 0
    
    except Exception as e:
        print(f"❌ Error: {e}")
        return 1


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="smarts-generator",
        description="Generate reaction SMARTS applicability patterns from reaction SMILES",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode with guided input
  python -m chemtools.protocol.smarts_generator_cli --interactive
  
  # Single reaction
  python -m chemtools.protocol.smarts_generator_cli \\
    --reaction "CCCCCCCCI.B2pin2>>CCCCCCCCBpin" \\
    --output my_protocol_smarts.json
  
  # Batch processing
  python -m chemtools.protocol.smarts_generator_cli \\
    --batch reactions.txt \\
    --output results.json

Format for protocol JSON:
  "reaction_smarts_applicability": {
    "core": "[C:1;H2,H3;X4]-[I:2]>>[C:1]-[B:3;X3](-[O;H0])(-[O;H0])",
    "guards_forbid": ["[CH]-I", "[C;H0]-I", "[CH2]-[a]-I"],
    "guards_require": [],
    "notes": "Primary alkyl iodides only, no aromatic/allylic"
  }
        """
    )
    
    parser.add_argument(
        '--interactive', '-i',
        action='store_true',
        help="Run in interactive mode with guided input"
    )
    
    parser.add_argument(
        '--reaction', '-r',
        type=str,
        help="Single reaction SMILES to process"
    )
    
    parser.add_argument(
        '--batch', '-b',
        type=Path,
        help="Batch process reactions from text file (one per line)"
    )
    
    parser.add_argument(
        '--output', '-o',
        type=Path,
        help="Output JSON file path"
    )
    
    parser.add_argument(
        '--visualize', '-v',
        action='store_true',
        help="Generate visualization images of SMARTS patterns (requires RDKit)"
    )
    
    parser.add_argument(
        '--viz-dir',
        type=Path,
        default=Path("smarts_visualizations"),
        help="Directory for visualization images (default: smarts_visualizations)"
    )
    
    parser.add_argument(
        '--test-smiles',
        type=str,
        help="Comma-separated SMILES to test against the pattern"
    )
    
    parser.add_argument(
        '--check-rdkit',
        action='store_true',
        help="Check if RDKit is available and exit"
    )
    
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """Main CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args(argv)
    
    # Check RDKit availability
    if args.check_rdkit:
        if RDKIT_AVAILABLE:
            print("✅ RDKit is available")
            return 0
        else:
            print("❌ RDKit is not available")
            print("   Install with: pip install rdkit")
            return 1
    
    # Determine mode
    if args.interactive:
        run_interactive_mode()
        return 0
    elif args.batch:
        return run_batch_mode(args.batch, args.output)
    elif args.reaction:
        # Parse test SMILES if provided
        test_smiles = []
        if args.test_smiles:
            test_smiles = [s.strip() for s in args.test_smiles.split(',') if s.strip()]
        
        result = run_single_reaction(args.reaction, args.output, args.visualize, args.viz_dir)
        
        # If visualization requested and we have test molecules, visualize them
        if args.visualize and test_smiles and result == 0:
            try:
                generator = SmartsGenerator(args.reaction)
                core = generator.generate_core_smarts()
                guards = generator.suggest_guard_patterns()
                applicability = ReactionSmartsApplicability(
                    core=core,
                    guards_forbid=guards
                )
                visualize_pattern_with_examples(applicability, test_smiles, args.viz_dir)
            except Exception as e:
                print(f"⚠️  Error during visualization: {e}")
        
        return result
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
