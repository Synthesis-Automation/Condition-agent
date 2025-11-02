"""
Quick Start: Calculable Features

This example shows how to use the calculable features system.
"""

from chemtools.featurizers import calculable
from chemtools.featurizers.molecular import featurize

# Example 1: Detect all features for a single molecule
print("="*70)
print("Example 1: Single Molecule Feature Detection")
print("="*70)

aryl_bromide = "c1ccc(Br)cc1"
features = calculable.detect_all_features(aryl_bromide)

print(f"\nMolecule: {aryl_bromide}")
print(f"sp2_bromide_present: {features['sp2_bromide_present']}")
print(f"aryl_halide_present: {features['aryl_halide_present']}")
print(f"ArBr_present: {features['ArBr_present']}")
print(f"sp2_halide_site_count: {features['sp2_halide_site_count']}")
print(f"polarity_low: {features['polarity_low']}")

# Example 2: Get only present features
print("\n" + "="*70)
print("Example 2: Get Present Features Only")
print("="*70)

present = calculable.get_present_features(aryl_bromide)
print(f"\nMolecule: {aryl_bromide}")
print(f"Present features ({len(present)}):")
for feature in present:
    print(f"  - {feature}")

# Example 3: Compare multiple molecules
print("\n" + "="*70)
print("Example 3: Batch Comparison")
print("="*70)

molecules = {
    "Aryl Bromide": "c1ccc(Br)cc1",
    "Aryl Chloride": "c1ccc(Cl)cc1", 
    "Alkyl Bromide": "CCBr",
    "Vinyl Bromide": "C=CBr",
}

print(f"\n{'Molecule':<20} {'ArBr':<8} {'ArCl':<8} {'sp3-Br':<8} {'Vinyl':<8}")
print("-" * 60)

for name, smiles in molecules.items():
    f = calculable.detect_all_features(smiles)
    print(f"{name:<20} "
          f"{'✓' if f['ArBr_present'] else '-':<8} "
          f"{'✓' if f['ArCl_present'] else '-':<8} "
          f"{'✓' if f['sp3_bromide_present'] else '-':<8} "
          f"{'✓' if f['vinyl_halide_present'] else '-':<8}")

# Example 4: Integration with molecular featurizer
print("\n" + "="*70)
print("Example 4: Integration with Molecular Featurizer")
print("="*70)

electrophile = "c1ccc(Br)cc1"
nucleophile = "CCN"

result = featurize(
    electrophile=electrophile,
    nucleophile=nucleophile,
    include_calculable=True
)

print(f"\nElectrophile: {electrophile}")
print(f"Nucleophile: {nucleophile}")
print(f"\nBase features:")
print(f"  LG: {result['LG']}")
print(f"  elec_class: {result['elec_class']}")
print(f"  nuc_class: {result['nuc_class']}")

print(f"\nCalculable features (electrophile):")
elec_calc = result['calculable']['electrophile']
print(f"  ArBr_present: {elec_calc['ArBr_present']}")
print(f"  sp2_halide_site_count: {elec_calc['sp2_halide_site_count']}")

print(f"\nCalculable features (nucleophile):")
nuc_calc = result['calculable']['nucleophile']
print(f"  aliphatic_amine_present: {nuc_calc['aliphatic_amine_present']}")
print(f"  aniline_present: {nuc_calc['aniline_present']}")

# Example 5: Feature summary (human-readable)
print("\n" + "="*70)
print("Example 5: Human-Readable Summary")
print("="*70)

print(calculable.feature_summary("c1ccc(Br)c(C(=O)O)c1"))

print("\n" + "="*70)
print("Examples Complete!")
print("="*70)
