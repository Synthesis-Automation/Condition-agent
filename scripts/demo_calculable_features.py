"""
Demo: Calculable Feature Detection

Demonstrates the usage of the calculable features module for detecting
structural motifs and chemical properties in molecules.

Usage:
    python scripts/demo_calculable_features.py
    python scripts/demo_calculable_features.py "c1ccc(Br)cc1"
    python scripts/demo_calculable_features.py --batch molecules.txt
"""

import sys
import json
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemtools.featurizers import calculable
from chemtools.util.rdkit_helpers import rdkit_available


def demo_single_molecule(smiles: str):
    """Demonstrate feature detection for a single molecule."""
    print(f"\n{'='*70}")
    print(f"Molecule: {smiles}")
    print(f"{'='*70}\n")
    
    # Detect all features
    features = calculable.detect_all_features(smiles)
    
    # Get only present features
    present = calculable.get_present_features(smiles)
    
    print(f"Total features detected: {len(present)}/{len(features)}\n")
    
    # Organize by category
    categories = {
        "Halides & Leaving Groups": [],
        "Boron & Organometallics": [],
        "Nucleophiles": [],
        "Unsaturation": [],
        "Heterocycles": [],
        "Polarity & Reactivity": [],
        "Counts": [],
        "Derived Features": [],
    }
    
    for token in sorted(present):
        value = features[token]
        
        # Categorize
        if any(x in token for x in ["halide", "bromide", "chloride", "iodide", "fluoride", 
                                     "triflate", "tosylate", "mesylate", "pseudohalide"]):
            if "count" in token:
                categories["Counts"].append((token, value))
            else:
                categories["Halides & Leaving Groups"].append((token, value))
        elif any(x in token for x in ["boron", "grignard", "zinc", "lithium", "stannane", "silane"]):
            categories["Boron & Organometallics"].append((token, value))
        elif any(x in token for x in ["amine", "aniline", "phenol", "alcohol", "thiol"]):
            categories["Nucleophiles"].append((token, value))
        elif any(x in token for x in ["alkene", "alkyne", "allylic", "benzylic"]):
            categories["Unsaturation"].append((token, value))
        elif any(x in token for x in ["pyridine", "pyrimidine", "azole", "indole", "heteroaryl"]):
            categories["Heterocycles"].append((token, value))
        elif any(x in token for x in ["polarity", "sensitive", "acidic", "beta_hydride", "poison"]):
            categories["Polarity & Reactivity"].append((token, value))
        elif "count" in token:
            categories["Counts"].append((token, value))
        elif any(x in token for x in ["Ar", "Vinyl", "internal"]):
            categories["Derived Features"].append((token, value))
    
    # Print organized results
    for category, items in categories.items():
        if items:
            print(f"{category}:")
            for token, value in items:
                if isinstance(value, bool):
                    print(f"  ✓ {token}")
                else:
                    print(f"  • {token}: {value}")
            print()


def demo_comparison():
    """Compare features across different molecule types."""
    print("\n" + "="*70)
    print("Comparison: Different Electrophile Types")
    print("="*70 + "\n")
    
    molecules = {
        "Aryl Bromide": "c1ccc(Br)cc1",
        "Aryl Chloride": "c1ccc(Cl)cc1",
        "Aryl Triflate": "c1ccc(OS(=O)(=O)C(F)(F)F)cc1",
        "Vinyl Bromide": "C=CBr",
        "Alkyl Bromide": "CCCBr",
    }
    
    comparison = {}
    for name, smiles in molecules.items():
        features = calculable.detect_all_features(smiles)
        comparison[name] = {
            "SMILES": smiles,
            "sp2_halide": features.get("sp2_bromide_present") or features.get("sp2_chloride_present"),
            "sp3_halide": features.get("sp3_bromide_present") or features.get("sp3_chloride_present"),
            "aryl": features.get("aryl_halide_present"),
            "vinyl": features.get("vinyl_halide_present"),
            "triflate": features.get("sp2_triflate_present"),
            "beta_hydride_risk": features.get("beta_hydride_possible"),
        }
    
    # Print as table
    print(f"{'Molecule':<20} {'SMILES':<30} {'sp2':<5} {'sp3':<5} {'Aryl':<5} {'Vinyl':<6} {'OTf':<5} {'β-H':<5}")
    print("-" * 90)
    for name, data in comparison.items():
        print(f"{name:<20} {data['SMILES']:<30} "
              f"{'✓' if data['sp2_halide'] else '-':<5} "
              f"{'✓' if data['sp3_halide'] else '-':<5} "
              f"{'✓' if data['aryl'] else '-':<5} "
              f"{'✓' if data['vinyl'] else '-':<6} "
              f"{'✓' if data['triflate'] else '-':<5} "
              f"{'✓' if data['beta_hydride_risk'] else '-':<5}")


def demo_nucleophiles():
    """Demonstrate nucleophile classification."""
    print("\n" + "="*70)
    print("Nucleophile Classification")
    print("="*70 + "\n")
    
    nucleophiles = {
        "Primary Amine": "CCN",
        "Secondary Amine": "CCNC",
        "Aniline": "c1ccc(N)cc1",
        "Indole": "c1ccc2[nH]ccc2c1",
        "Phenol": "c1ccc(O)cc1",
        "Alcohol": "CCO",
        "Thiol": "CCS",
    }
    
    print(f"{'Nucleophile':<20} {'SMILES':<25} {'Features'}")
    print("-" * 80)
    
    for name, smiles in nucleophiles.items():
        present = calculable.get_present_features(smiles)
        relevant = [f for f in present if any(x in f for x in 
                   ["amine", "aniline", "phenol", "alcohol", "thiol", "indole"])]
        print(f"{name:<20} {smiles:<25} {', '.join(relevant[:3])}")


def demo_suzuki_coupling():
    """Demonstrate Suzuki coupling substrate analysis."""
    print("\n" + "="*70)
    print("Suzuki-Miyaura Coupling: Substrate Analysis")
    print("="*70 + "\n")
    
    print("Electrophile (Aryl Bromide):")
    ar_br = "c1ccc(Br)c(C(=O)OC)c1"
    print(f"  SMILES: {ar_br}")
    features = calculable.detect_all_features(ar_br)
    print(f"  ArBr: {features['ArBr_present']}")
    print(f"  sp2 halide count: {features['sp2_halide_site_count']}")
    print(f"  Ester present: {features.get('ester_present', 'N/A')}")
    
    print("\nNucleophile (Aryl Boronic Acid):")
    ar_b = "c1ccc(B(O)O)cc1"
    print(f"  SMILES: {ar_b}")
    features = calculable.detect_all_features(ar_b)
    print(f"  sp2 boron: {features['sp2_boron_present']}")
    print(f"  Aromatic: Present")


def demo_batch_analysis(file_path: str):
    """Analyze multiple molecules from a file."""
    print("\n" + "="*70)
    print(f"Batch Analysis: {file_path}")
    print("="*70 + "\n")
    
    # Read SMILES from file
    smiles_list = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                smiles_list.append(line)
    
    # Batch detection
    results = calculable.detect_features_batch(smiles_list)
    
    # Summary statistics
    print(f"Analyzed {len(results)} molecules\n")
    
    # Count feature occurrences
    feature_counts = {}
    for result in results:
        for token, value in result.items():
            if (isinstance(value, bool) and value) or (isinstance(value, int) and value > 0):
                feature_counts[token] = feature_counts.get(token, 0) + 1
    
    # Print top 20 most common features
    print("Most common features:")
    for token, count in sorted(feature_counts.items(), key=lambda x: x[1], reverse=True)[:20]:
        pct = (count / len(results)) * 100
        print(f"  {token:<30} {count:>4} ({pct:>5.1f}%)")


def main():
    """Main entry point."""
    if not rdkit_available():
        print("ERROR: RDKit is not available. Please install RDKit to use this module.")
        return 1
    
    print("\n" + "="*70)
    print("Calculable Features Demo")
    print("="*70)
    
    # Check command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == "--batch" and len(sys.argv) > 2:
            demo_batch_analysis(sys.argv[2])
        else:
            # Single molecule from command line
            demo_single_molecule(sys.argv[1])
    else:
        # Run all demos
        demo_single_molecule("c1ccc(Br)cc1")
        demo_comparison()
        demo_nucleophiles()
        demo_suzuki_coupling()
    
    print("\n" + "="*70)
    print("Demo complete!")
    print("="*70 + "\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
