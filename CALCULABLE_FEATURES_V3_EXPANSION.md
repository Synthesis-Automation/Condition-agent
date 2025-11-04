# Calculable Features v3.0 - Comprehensive Expansion

## Overview

The `calculable_features.json` has been expanded from **107 features to 244+ tokens** (206 direct features + 38 derived shortcuts), creating a comprehensive molecular feature detection library for organic chemistry reasoning, synthesis planning, retrosynthesis, and medicinal chemistry applications.

## Expansion Summary

- **Version**: 2025-11-04.v3.0-comprehensive
- **Original features**: 91 base features
- **New features added**: 115 features
- **Total features**: 206 direct detection features
- **Derived shortcuts**: 38 logical combinations
- **Total tokens available**: 244 tokens
- **Categories**: 17 organized domains

## New Feature Categories

### 1. Protecting Groups (Comprehensive) - 27 features

Added comprehensive protecting group library covering all major functional groups:

#### Alcohol Protecting Groups
- `methyl_ether_present` - Permanent protection, strong acid cleavage
- `benzyl_ether_present` - Hydrogenolysis cleavable (Pd/C, H2)
- `pmb_ether_present` - DDQ/CAN oxidatively cleavable
- `mom_ether_present` - Acid-labile methoxymethyl ether
- `sem_ether_present` - Fluoride-cleavable silyl ether
- `tbdps_ether_present` - Bulky silyl group
- `tbs_ether_present` - Common TBDMS protection
- `tips_ether_present` - Very bulky, stable silyl group
- `acetate_ester_present` - Simple, base/acid hydrolyzable
- `benzoate_ester_present` - More stable than acetate
- `pivaloate_ester_present` - Sterically hindered
- `acetal_present` - Carbonyl protection, acid-labile
- `ketal_present` - Ketone protection, acid-labile

#### Amine Protecting Groups
- `tosyl_sulfonamide_present` - Electron-withdrawing, reduces nucleophilicity
- `nosyl_sulfonamide_present` - Thiol-cleavable
- `trifluoroacetamide_present` - More acid-labile than regular amides
- `alloc_carbamate_present` - Pd(0)-cleavable
- `troc_carbamate_present` - Zn-reduction cleavable
- `phthalimide_present` - Gabriel synthesis, hydrazine-cleavable

#### Carbonyl Protecting Groups
- `dithiane_present` - Enables umpolung chemistry
- `oxime_present` - Reversible carbonyl protection
- `hydrazone_present` - Wolff-Kishner precursor
- `enol_ether_present` - Masked carbonyl, acid-labile

#### Carboxylic Acid Protecting Groups
- `methyl_ester_present` - Base hydrolysis
- `ethyl_ester_present` - Slightly more lipophilic
- `tert_butyl_ester_present` - Acid-labile, orthogonal
- `benzyl_ester_present` - Hydrogenolysis cleavable

### 2. Heterocycles (Comprehensive) - 20 features

Expanded heterocycle detection covering common drug scaffolds:

#### 5-Membered Heterocycles
- `furan_present` - Aromatic O-heterocycle, Diels-Alder reactive
- `thiophene_present` - Aromatic S-heterocycle, metabolically stable
- `pyrrole_present` - Aromatic NH-heterocycle, π-excessive
- `imidazole_present` - 2 nitrogens, basic, coordinating
- `oxazole_present` - O,N-heterocycle, π-deficient
- `thiazole_present` - S,N-heterocycle, drug scaffold
- `triazole_present` - Click chemistry product
- `tetrazole_present` - Carboxylic acid bioisostere

#### 6-Membered Heterocycles
- `morpholine_present` - Saturated O,N-heterocycle, solubilizing
- `piperazine_present` - Saturated diamine, drug linker
- `piperidine_present` - Saturated N-heterocycle, basic
- `pyrrolidine_present` - Saturated 5-membered, proline analog
- `tetrahydrofuran_present` - Saturated O-heterocycle

#### Fused Heterocycles
- `quinoline_present` - Benzene-pyridine fusion, drug scaffold
- `isoquinoline_present` - Isomeric fusion, alkaloid scaffold
- `benzofuran_present` - Benzene-furan fusion
- `benzothiophene_present` - Benzene-thiophene fusion
- `benzimidazole_present` - PPI scaffold
- `purine_present` - Nucleobase, kinase inhibitor scaffold
- `pyrimidine_derivative_present` - Nucleobase scaffold

### 3. Reactive Intermediates & Instability - 15 features

Critical safety and reactivity markers:

- `epoxide_present` - Strained ring, electrophilic
- `aziridine_present` - Strained N-ring, toxic
- `diazo_compound_present` - Carbene precursor, explosive
- `azide_present` - Click chemistry, explosive
- `peroxide_present` - Oxidizing, explosive
- `hydroperoxide_present` - Reactive autoxidation product
- `isocyanate_present` - Electrophilic, toxic
- `isothiocyanate_present` - Electrophilic, pungent
- `sulfonyl_chloride_present` - Highly reactive
- `acid_anhydride_present` - Acylating agent
- `alpha_beta_unsaturated_carbonyl_present` - Michael acceptor
- `alpha_beta_unsaturated_ester_present` - Conjugate addition
- `alpha_beta_unsaturated_nitrile_present` - Electrophilic
- `alpha_halo_carbonyl_present` - Alkylating agent
- `strained_ring_present` - 3-4 membered rings

### 4. Physicochemical Properties (ADME) - 10 features

Drug-likeness and bioavailability predictors:

#### Lipinski Rules
- `lipinski_hbd_compliant` - ≤5 H-bond donors
- `lipinski_hba_compliant` - ≤10 H-bond acceptors
- `lipinski_mw_compliant` - MW < 500 Da
- `lipinski_logp_compliant` - logP < 5

#### Veber Rules
- `veber_rotb_compliant` - ≤10 rotatable bonds
- `veber_tpsa_compliant` - TPSA ≤ 140 Ų

#### Ionization
- `ionizable_basic_group_present` - Basic groups (pKa 7-11)
- `ionizable_acidic_group_present` - Acidic groups (pKa 2-10)
- `quaternary_ammonium_present` - Permanently charged
- `zwitterion_potential` - Both acidic and basic groups

### 5. Medicinal Chemistry Features - 13 features

Drug design and optimization markers:

- `aromatic_ring_count` - Lipophilicity indicator
- `fsp3_high` - >50% sp3 carbons, 3D character
- `fsp3_low` - <30% sp3, flat/aromatic
- `fluorine_present` - Metabolic blocker
- `trifluoromethyl_present` - CF3, lipophilic, stable
- `difluoromethyl_present` - CHF2, H-bond donor mimic
- `sulfone_present` - Polar, stable, HBA
- `sulfonamide_present` - Acid bioisostere
- `urea_present` - HBD/HBA, kinase motif
- `carbamate_present` - Ester-amide hybrid

#### PAINS Alerts
- `pains_aldehyde_alert` - Reactive aldehydes
- `pains_catechol_alert` - Redox-active ortho-diphenol
- `pains_quinone_alert` - Redox-active quinone

### 6. Stereochemistry & Topology - 6 features

3D structure and chirality indicators:

- `alkene_e_isomer_present` - E/trans configuration
- `alkene_z_isomer_present` - Z/cis configuration
- `atropisomer_potential` - Axial chirality from hindered rotation
- `spiro_center_present` - Spirocyclic complexity
- `bridgehead_nitrogen_present` - Conformationally restricted
- `quaternary_carbon_present` - All-carbon substituted center

### 7. Redox & Photochemistry - 8 features

Oxidation/reduction and photoreactivity markers:

#### Oxidation-Prone Sites
- `benzylic_ch_present` - Weak C-H BDE, metabolic hotspot
- `allylic_ch_present` - Weak BDE, radical site
- `propargylic_ch_present` - Acidic C-H
- `alpha_amino_ch_present` - SET oxidation site
- `alpha_oxy_ch_present` - Oxidation-prone

#### Photoactivity
- `extended_conjugation_present` - Chromophore, UV-active
- `aromatic_ketone_present` - Photoactive, Norrish reactions
- `naphthoquinone_present` - Redox-active, ROS generator

### 8. Additional Organometallics - 8 features

Expanded coupling partner detection:

- `aryl_boron_present` - Suzuki aryl donor
- `vinyl_boron_present` - Suzuki vinyl donor
- `alkyl_boron_present` - Activated coupling
- `aryl_stannane_present` - Stille aryl donor
- `vinyl_stannane_present` - Stille vinyl donor (stereochemistry preserved)
- `aryl_silane_present` - Hiyama coupling
- `vinyl_silane_present` - Hiyama/Peterson
- `organocuprate_potential` - Soft nucleophile

### 9. Sulfur & Phosphorus - 8 features

Expanded S/P functionality:

- `sulfide_present` - Thioether, oxidizable
- `sulfoxide_present` - Chiral, polar, Pummerer substrate
- `sulfonic_acid_present` - Strong acid, very polar
- `sulfinic_acid_present` - Reducing agent, SO2 source
- `disulfide_present` - Redox-active, reductive cleavage
- `phosphate_present` - Biologically relevant
- `phosphonate_present` - Phosphate bioisostere
- `phosphonium_salt_present` - Wittig precursor

## New Derived Shortcuts (21 added)

Logical combinations for complex queries:

### Protection Strategy Queries
- `protected_alcohol_present` - Any alcohol protecting group
- `protected_amine_present` - Any amine protecting group
- `protected_carbonyl_present` - Any carbonyl protecting group
- `protected_acid_present` - Any carboxylic acid protecting group

### Heterocycle Classification
- `five_membered_heterocycle_present` - Any 5-membered heterocycle
- `six_membered_heterocycle_present` - Any 6-membered heterocycle
- `fused_heterocycle_present` - Any fused heterocycle system

### Drug-Likeness
- `lipinski_compliant` - Passes all 4 Lipinski rules
- `veber_compliant` - Passes both Veber rules
- `drug_like` - Lipinski + Veber + basic PAINS filters

### Reactivity Patterns
- `michael_acceptor_present` - Any Michael acceptor system
- `electrophilic_warhead_present` - Covalent inhibitor motifs
- `metabolic_soft_spot_present` - Oxidation-prone sites

### Functional Motifs
- `fluorinated_motif_present` - Any fluorine-containing group
- `cross_coupling_electrophile` - Suitable Pd-coupling electrophile
- `cross_coupling_nucleophile` - Suitable Pd-coupling nucleophile

### Stability Indicators
- `oxidation_prone` - Contains oxidizable groups
- `reduction_prone` - Contains reducible groups
- `moisture_sensitive` - Requires anhydrous conditions
- `air_sensitive` - Requires inert atmosphere
- `explosive_risk` - Potentially explosive groups

## Usage Examples

### Basic Feature Detection

```python
from chemtools.featurizers.calculable import detect_all_features, detect_feature

# Detect all features for a molecule
features = detect_all_features("c1ccc(Br)cc1")
# Returns: {'sp2_bromide_present': True, 'aryl_halide_present': True, ...}

# Check specific feature
has_protecting_group = detect_feature("CC(=O)Oc1ccccc1", "acetate_ester_present")
# Returns: True
```

### Drug-Likeness Assessment

```python
features = detect_all_features("CC(C)Cc1ccc(cc1)C(C)C(=O)O")  # Ibuprofen

drug_like = features.get('drug_like', False)
lipinski = features.get('lipinski_compliant', False)
pains_alerts = any([
    features.get('pains_aldehyde_alert', False),
    features.get('pains_catechol_alert', False),
    features.get('pains_quinone_alert', False)
])
```

### Protection Strategy Planning

```python
features = detect_all_features("Nc1ccc(O)cc1")  # 4-Aminophenol

needs_protection = {
    'amine': not features.get('protected_amine_present', False),
    'phenol': not features.get('protected_alcohol_present', False)
}

# Check what protecting groups are available
available_pg = {
    'Boc': True,  # Can add
    'Cbz': True,
    'Fmoc': True,
    'Acetate': True
}
```

### Cross-Coupling Reaction Planning

```python
# Check electrophile
elec_features = detect_all_features("Brc1ccc(F)cc1")
is_suitable_electrophile = elec_features.get('cross_coupling_electrophile', False)
# True (has aryl bromide)

# Check nucleophile
nuc_features = detect_all_features("c1ccc(B(O)O)cc1")
is_suitable_nucleophile = nuc_features.get('cross_coupling_nucleophile', False)
# True (has aryl boron)

# Check for potential issues
ortho_hindered = elec_features.get('ortho_substitution_present', False)
```

### Safety Assessment

```python
features = detect_all_features("CC(=O)N=[N+]=[N-]")  # Acetyl azide

safety_flags = {
    'explosive': features.get('explosive_risk', False),
    'moisture_sensitive': features.get('moisture_sensitive', False),
    'air_sensitive': features.get('air_sensitive', False)
}

# All return True - handle with extreme caution!
```

## LLM Reasoning Integration

Each feature includes a `why` field explaining its chemical significance, enabling LLM-based reasoning:

```python
from chemtools.featurizers.calculable import get_feature_spec

spec = get_feature_spec()

# Get explanation for a feature
for feature in spec['features']:
    if feature['token'] == 'pains_catechol_alert':
        print(feature['why'])
        # "PAINS alert: catechol (ortho-diphenol) is redox-active, metal chelator"
```

### Example LLM Prompt Integration

```python
features = detect_all_features(smiles)
active_features = [k for k, v in features.items() if v]

prompt = f"""
Analyze this molecule: {smiles}

Detected features:
{chr(10).join('- ' + feat for feat in active_features)}

Provide:
1. Synthesis strategy considering protecting groups needed
2. Potential reactivity issues
3. Cross-coupling reaction suitability
4. Drug-likeness assessment
5. Safety considerations
"""
```

## Category-Based Queries

Features are organized into 17 categories for systematic analysis:

```python
spec = get_feature_spec()

# Get all protecting group features
pg_features = [f for f in spec['features'] if f.get('category') == 'protecting_groups']

# Get all heterocycle features
hetero_features = [f for f in spec['features'] if f.get('category') == 'heterocycles']

# Get all ADME-relevant features
adme_features = [f for f in spec['features'] if f.get('category') == 'physicochemical']
```

## Testing

Run comprehensive tests:

```bash
pytest tests/test_calculable_features.py -v
```

## Version History

- **v3.0-comprehensive** (2025-11-04): Expanded from 107 to 244+ tokens
  - Added 115 new direct detection features
  - Added 21 new derived shortcuts
  - Organized into 17 categories
  - Enhanced LLM reasoning support

- **v2.2** (2025-11-02): Original implementation
  - 91 base features
  - 17 derived shortcuts
  - Basic coverage for cross-coupling reactions

## Future Enhancements

Potential areas for further expansion:

1. **Specific Reaction Detection**: Named reaction patterns (Wittig, Grubbs, etc.)
2. **Isotope Labeling**: Detection of D, ¹³C, ¹⁵N labels
3. **Metal Complexes**: Coordination chemistry features
4. **Polymer Chemistry**: Repeat units, end groups
5. **Carbohydrate Features**: Glycosidic bonds, anomeric centers
6. **Peptide Features**: Amino acid residues, backbone patterns
7. **Natural Product Features**: Polyketide, terpene, alkaloid motifs
8. **Fragmentability Scores**: Mass spec fragmentation predictors
9. **Synthetic Accessibility**: Retrosynthetic complexity scores
10. **Green Chemistry Metrics**: Atom economy, E-factor proxies

## Integration with Existing Workflows

The expanded feature set integrates seamlessly with existing chemtools modules:

- **Router**: Family detection enhanced with protecting group awareness
- **Recommend**: Condition selection improved with safety/stability markers
- **Constraints**: Filtering enhanced with drug-likeness and PAINS alerts
- **Explain**: Reasoning enriched with comprehensive feature explanations
- **LLM Tools**: Chemistry agents can leverage detailed feature descriptions

## References

- Greene's Protective Groups in Organic Synthesis
- Lipinski's Rule of Five (2001)
- Veber's Rules (2002)
- PAINS filters (Baell & Holloway, 2010)
- Cross-Coupling Reactions (Hartwig, 2010)

## Contact & Contributions

For questions, bug reports, or feature requests, please open an issue in the repository.

---

**Note**: This comprehensive expansion maintains backward compatibility with existing code while adding extensive new capabilities for organic chemistry reasoning and planning.
