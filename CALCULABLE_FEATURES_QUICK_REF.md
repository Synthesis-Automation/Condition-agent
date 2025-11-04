# Calculable Features v3.0 - Quick Reference

## 🚀 Quick Start

```python
from chemtools.featurizers.calculable import detect_all_features, detect_feature

# Detect all features
features = detect_all_features("c1ccc(Br)cc1")

# Check specific feature
has_pg = detect_feature("CC(C)(C)OC(=O)Nc1ccccc1", "boc_present")
```

## 📋 Feature Categories (244 total tokens)

| Category | Count | Examples |
|----------|-------|----------|
| **Protecting Groups** | 27 | Boc, Cbz, Fmoc, TBS, benzyl, acetate, phthalimide |
| **Heterocycles** | 20 | Furan, thiophene, quinoline, morpholine, triazole |
| **Reactive Intermediates** | 15 | Epoxides, azides, peroxides, Michael acceptors |
| **ADME Properties** | 10 | Lipinski rules, Veber rules, ionizable groups |
| **Medicinal Chemistry** | 13 | Fsp3, fluorinated motifs, PAINS alerts, ureas |
| **Stereochemistry** | 6 | E/Z alkenes, atropisomers, spiro centers |
| **Redox/Photo** | 8 | Benzylic/allylic CH, chromophores, quinones |
| **Organometallics** | 8 | Aryl/vinyl/alkyl boron, stannanes, silanes |
| **S/P Functionality** | 8 | Sulfides, sulfoxides, phosphates, disulfides |
| **Halides** | 20+ | Aryl/vinyl/alkyl halides, triflates, tosylates |
| **Amines** | 10+ | Primary, secondary, tertiary, anilines |
| **Carbonyls** | 15+ | Ketones, aldehydes, esters, amides, acids |
| **Derived Features** | 38 | Drug-like, explosive risk, cross-coupling ready |

## 🎯 Common Use Cases

### 1. Drug-Likeness Assessment
```python
features = detect_all_features(smiles)
checks = {
    'lipinski': features.get('lipinski_compliant'),  # Rule of 5
    'veber': features.get('veber_compliant'),        # Additional rules
    'drug_like': features.get('drug_like'),          # Combined + PAINS
    'pains_alerts': any([
        features.get('pains_aldehyde_alert'),
        features.get('pains_catechol_alert'),
        features.get('pains_quinone_alert')
    ])
}
```

### 2. Protection Strategy Planning
```python
features = detect_all_features("Nc1ccc(O)cc1")  # 4-Aminophenol
strategy = {
    'amine_protected': features.get('protected_amine_present'),
    'alcohol_protected': features.get('protected_alcohol_present'),
    'available_pg_types': {
        'amine': ['Boc', 'Cbz', 'Fmoc', 'Ts', 'Ns'],
        'alcohol': ['TBS', 'TBDPS', 'Bn', 'Ac', 'Bz']
    }
}
```

### 3. Cross-Coupling Reaction Planning
```python
# Check electrophile
elec = detect_all_features("Brc1ccccc1")
suitable_elec = elec.get('cross_coupling_electrophile')  # True

# Check nucleophile
nuc = detect_all_features("c1ccc(B(O)O)cc1")
suitable_nuc = nuc.get('cross_coupling_nucleophile')  # True

# Check compatibility
issues = {
    'ortho_hindered': elec.get('ortho_substitution_present'),
    'pyridine_poison': nuc.get('pyridine_poison_risk'),
    'moisture_sensitive': nuc.get('moisture_sensitive')
}
```

### 4. Safety Assessment
```python
features = detect_all_features("CCN=[N+]=[N-]")  # Ethyl azide
safety = {
    'explosive_risk': features.get('explosive_risk'),      # True!
    'moisture_sensitive': features.get('moisture_sensitive'),
    'air_sensitive': features.get('air_sensitive'),
    'specific_hazards': [
        'azide_present',
        'diazo_compound_present',
        'peroxide_present'
    ]
}
```

### 5. Metabolic Hotspot Identification
```python
features = detect_all_features("CC(C)NCC(COc1ccccc1)O")  # Propranolol
metabolism = {
    'soft_spots': features.get('metabolic_soft_spot_present'),
    'oxidation_sites': features.get('oxidation_prone'),
    'specific': {
        'benzylic': features.get('benzylic_ch_present'),
        'allylic': features.get('allylic_ch_present'),
        'alpha_amino': features.get('alpha_amino_ch_present'),
        'alpha_oxy': features.get('alpha_oxy_ch_present')
    }
}
```

## 📊 Key Derived Features

### Protection Queries
- `protected_alcohol_present` - Any OH protecting group
- `protected_amine_present` - Any NH protecting group
- `protected_carbonyl_present` - Any C=O protecting group
- `protected_acid_present` - Any COOH protecting group (esters)

### Heterocycle Queries
- `five_membered_heterocycle_present`
- `six_membered_heterocycle_present`
- `fused_heterocycle_present`

### Drug-Likeness Queries
- `lipinski_compliant` - All 4 Lipinski rules
- `veber_compliant` - Both Veber rules
- `drug_like` - Lipinski + Veber + basic PAINS

### Reactivity Queries
- `michael_acceptor_present` - Any α,β-unsaturated electrophile
- `electrophilic_warhead_present` - Covalent inhibitor motifs
- `cross_coupling_electrophile` - Suitable for Pd coupling
- `cross_coupling_nucleophile` - Suitable for Pd coupling

### Safety Queries
- `explosive_risk` - Azides, diazo, peroxides, nitro groups
- `moisture_sensitive` - Acyl halides, anhydrides, organometallics
- `air_sensitive` - Phosphines, thiols, aldehydes, Grignards
- `oxidation_prone` - Benzylic/allylic CH, sulfides, thiols
- `reduction_prone` - Nitro, nitrile, carbonyl, unsaturation

### Stability Queries
- `metabolic_soft_spot_present` - Oxidation-prone metabolic sites
- `fluorinated_motif_present` - Metabolic blocking with F

## 🔍 Feature Detection Methods

Features are detected via:
1. **SMARTS patterns** - Direct substructure matching
2. **Heuristics** - Calculated properties (MW, logP, TPSA, etc.)
3. **Derived logic** - Boolean combinations of other features

All features include:
- `token` - Feature name
- `type` - bool, int, or float
- `scope` - Usually "global"
- `detect` - Detection method (SMARTS or heuristic)
- `why` - Chemical reasoning explanation (for LLMs)
- `category` - Organizational category (optional)

## 📈 Coverage Statistics

- **244 total tokens** (206 base + 38 derived)
- **126% expansion** from v2.2 (107 → 244)
- **17 categories** for systematic organization
- **20+ chemistry domains** covered
- **LLM-ready** with comprehensive explanations

## 🛠️ Advanced Usage

### Get Feature Specification
```python
from chemtools.featurizers.calculable import get_feature_spec

spec = get_feature_spec()
version = spec['version']
all_features = spec['features']
derived_features = spec['derived_shortcuts']
categories = spec['schema_notes']['categories']
```

### Filter by Category
```python
spec = get_feature_spec()
pg_features = [f for f in spec['features'] 
               if f.get('category') == 'protecting_groups']
```

### Get Feature Explanations
```python
for feature in spec['features']:
    if feature['token'] == 'boc_present':
        print(feature['why'])
        # "tert-Butoxycarbonyl (Boc) protecting group on nitrogen"
```

## 📚 Related Documentation

- **CALCULABLE_FEATURES_V3_EXPANSION.md** - Comprehensive guide (400+ lines)
- **CALCULABLE_FEATURES_V3_SUMMARY.md** - Executive summary
- **chemtools/featurizers/calculable.py** - Implementation
- **tests/test_calculable_features.py** - Test suite

## ⚠️ Important Notes

1. **RDKit Required**: Most features require RDKit for SMARTS matching
2. **Graceful Fallback**: Returns False/0 if RDKit unavailable
3. **Caching**: Compiled SMARTS patterns are cached for performance
4. **Boolean Logic**: Use derived features for complex queries
5. **Case Sensitivity**: SMARTS patterns are case-sensitive

## 🎓 Example Workflows

### Complete Molecule Analysis
```python
def analyze_molecule(smiles):
    features = detect_all_features(smiles)
    
    return {
        'drug_likeness': {
            'lipinski': features.get('lipinski_compliant'),
            'veber': features.get('veber_compliant'),
            'pains': any([features.get(f'pains_{x}_alert') 
                         for x in ['aldehyde', 'catechol', 'quinone']])
        },
        'protection_needs': {
            'amine': not features.get('protected_amine_present'),
            'alcohol': not features.get('protected_alcohol_present'),
            'carbonyl': not features.get('protected_carbonyl_present')
        },
        'reactivity': {
            'michael_acceptor': features.get('michael_acceptor_present'),
            'oxidation_sites': features.get('oxidation_prone'),
            'reduction_sites': features.get('reduction_prone')
        },
        'safety': {
            'explosive': features.get('explosive_risk'),
            'moisture': features.get('moisture_sensitive'),
            'air': features.get('air_sensitive')
        },
        'coupling': {
            'electrophile': features.get('cross_coupling_electrophile'),
            'nucleophile': features.get('cross_coupling_nucleophile')
        }
    }
```

### LLM Integration
```python
def generate_synthesis_prompt(smiles):
    features = detect_all_features(smiles)
    spec = get_feature_spec()
    
    active = [k for k, v in features.items() if v]
    explanations = {}
    
    for feat in spec['features'] + spec['derived_shortcuts']:
        if feat['token'] in active:
            explanations[feat['token']] = feat.get('why', '')
    
    prompt = f"""
    Molecule: {smiles}
    
    Detected Features:
    {chr(10).join(f"- {k}: {explanations.get(k, '')}" for k in active[:20])}
    
    Provide:
    1. Retrosynthetic analysis
    2. Protecting group strategy
    3. Safety considerations
    4. Recommended reaction conditions
    """
    
    return prompt
```

---

**Version**: 2025-11-04.v3.0-comprehensive  
**Quick Ref**: For full details see CALCULABLE_FEATURES_V3_EXPANSION.md
