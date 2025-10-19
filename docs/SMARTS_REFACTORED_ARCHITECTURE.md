# Chemistry-Aware SMARTS Generation - Refactored Architecture

## 🎯 Design Principle: Separation of Concerns

You're absolutely right! The chemistry intelligence should be in **reusable utilities**, not buried in application code.

### **The Right Way to Organize**

```
chemtools/                      # Reusable chemistry libraries
├── util/                       # Core chemistry utilities (REUSABLE)
│   ├── functional_groups.py   # ✅ Already exists (80+ FG detection)
│   ├── rdkit_helpers.py        # ✅ Already exists (RDKit wrappers)
│   ├── substrate_classifier.py # 🆕 NEW: Classify substrates
│   └── smarts_builders.py      # 🆕 NEW: Build context-aware SMARTS
│
├── protocol/                   # Protocol-specific applications
│   └── smarts_generator_cli.py # Uses util/* for chemistry logic
│
├── featurizers/                # ML feature extraction
│   └── *.py                    # Could also benefit from substrate_classifier
│
└── recommend/                  # Recommendation engine
    └── *.py                    # Could use smarts_builders for matching
```

---

## 📦 Module Breakdown

### **Module 1: `chemtools/util/substrate_classifier.py`** (NEW - REUSABLE)

**Purpose**: Classify chemical substrates by type and context

**Exports**:
```python
class SubstrateClassifier:
    """Classify chemical substrates into meaningful categories"""
    
    def classify(self, mol_or_smiles) -> SubstrateInfo:
        """Complete substrate classification"""
        
    def get_carbon_types(self, mol) -> Dict[int, str]:
        """Map atom index → 'sp3'|'sp2'|'sp'|'aromatic'"""
        
    def find_special_positions(self, mol) -> SpecialPositions:
        """Identify benzylic, allylic, propargylic positions"""
        
    def analyze_heteroatom_environment(self, atom, mol) -> Environment:
        """Analyze chemical environment of heteroatom"""

# Data classes
@dataclass
class SubstrateInfo:
    substrate_class: str  # 'alkyl_halide', 'aryl_halide', 'aniline', etc.
    carbon_types: Dict[int, str]  # atom_idx → hybridization
    functional_groups: List[str]  # detected functional groups
    special_positions: SpecialPositions
    reactive_centers: List[int]  # atom indices

@dataclass
class SpecialPositions:
    benzylic: List[int]       # [CH2] next to aromatic
    allylic: List[int]        # [CH2] next to C=C
    propargylic: List[int]    # [CH2] next to C≡C
    alpha_to_carbonyl: List[int]
    ortho_to_heteroatom: List[int]
```

**Who uses this?**
- ✅ SMARTS generator (protocol scope definition)
- ✅ Featurizers (ML feature extraction)
- ✅ Recommendation engine (substrate matching)
- ✅ Reaction type detector (understanding reactants)

---

### **Module 2: `chemtools/util/smarts_builders.py`** (NEW - REUSABLE)

**Purpose**: Build context-aware SMARTS patterns for functional groups

**Exports**:
```python
class SmartsBuilder:
    """Build chemically meaningful SMARTS patterns"""
    
    def build_for_substrate(self, substrate_info: SubstrateInfo) -> str:
        """Generate SMARTS for entire substrate"""
        
    def build_halide_smarts(self, halide_atom, mol, context: SubstrateInfo) -> str:
        """Context-aware halide SMARTS"""
        
    def build_amine_smarts(self, n_atom, mol, context: SubstrateInfo) -> str:
        """Context-aware amine SMARTS"""
        
    def build_carbonyl_smarts(self, c_atom, mol, context: SubstrateInfo) -> str:
        """Context-aware carbonyl SMARTS"""
    
    def build_alcohol_smarts(self, o_atom, mol, context: SubstrateInfo) -> str:
        """Context-aware alcohol SMARTS"""
    
    def build_boron_smarts(self, b_atom, mol, context: SubstrateInfo) -> str:
        """Context-aware boron SMARTS"""

class SmartsPatternMatcher:
    """Match molecules against SMARTS patterns with explanations"""
    
    def match(self, mol, smarts: str) -> MatchResult:
        """Match molecule against SMARTS"""
        
    def explain_match(self, mol, smarts: str) -> str:
        """Explain why molecule matches/doesn't match"""
        
    def find_matching_atoms(self, mol, smarts: str) -> List[int]:
        """Return atom indices that match pattern"""
```

**Who uses this?**
- ✅ SMARTS generator CLI (main user)
- ✅ Protocol database queries (find compatible protocols)
- ✅ Substrate validation (check if substrate matches scope)
- ✅ Test frameworks (validate patterns)

---

### **Module 3: `chemtools/util/functional_groups.py`** (ALREADY EXISTS ✅)

**Current state**: Already has 80+ functional group SMARTS patterns!

**Enhancement needed**: Add classification helper
```python
# ADD THIS to existing file:

def classify_functional_group_context(mol, fg_name: str) -> str:
    """Given a functional group, classify its chemical context
    
    Examples:
        - 'amine' → 'aniline' | 'aliphatic_primary' | 'aliphatic_secondary'
        - 'halide' → 'aryl' | 'alkyl' | 'vinyl' | 'benzylic' | 'allylic'
        - 'alcohol' → 'phenol' | 'benzylic' | 'allylic' | 'aliphatic'
    """
```

---

### **Module 4: `chemtools/protocol/smarts_generator_cli.py`** (APPLICATION LAYER)

**Purpose**: CLI tool for protocol scope definition (uses util modules)

**Responsibilities**:
- ✅ Parse command-line arguments
- ✅ Parse reaction SMILES
- ✅ Call `SubstrateClassifier` to understand chemistry
- ✅ Call `SmartsBuilder` to generate patterns
- ✅ Suggest guard patterns based on classification
- ✅ Generate visualizations
- ✅ Format output as JSON for protocol database

**Should NOT contain**:
- ❌ Chemistry classification logic → `substrate_classifier.py`
- ❌ SMARTS pattern building → `smarts_builders.py`
- ❌ Functional group detection → `functional_groups.py`

**Should be thin wrapper**:
```python
class SmartsGenerator:
    def __init__(self, reaction_smiles: str):
        self.reaction_smiles = reaction_smiles
        self.classifier = SubstrateClassifier()  # from util
        self.builder = SmartsBuilder()            # from util
    
    def generate_core_smarts(self) -> str:
        # 1. Parse reaction
        reactants, products = self._parse_reaction()
        
        # 2. Classify substrates (REUSABLE)
        reactant_info = self.classifier.classify(reactants[0])
        product_info = self.classifier.classify(products[0])
        
        # 3. Build SMARTS (REUSABLE)
        reactant_smarts = self.builder.build_for_substrate(reactant_info)
        product_smarts = self.builder.build_for_substrate(product_info)
        
        return f"{reactant_smarts}>>{product_smarts}"
    
    def suggest_guard_patterns(self) -> List[str]:
        # Use classifier results to suggest guards
        reactant_info = self.classifier.classify(self.reactants_smiles)
        return self._generate_guards_from_classification(reactant_info)
```

---

## 🎯 Benefits of This Architecture

### **1. Reusability**
Other modules can use the same chemistry intelligence:

```python
# In featurizers/custom_featurizer.py
from chemtools.util.substrate_classifier import SubstrateClassifier

classifier = SubstrateClassifier()
info = classifier.classify(molecule)
if info.substrate_class == 'benzylic_halide':
    # Add special features for benzylic substrates
    ...

# In recommend/substrate_validator.py
from chemtools.util.smarts_builders import SmartsPatternMatcher

matcher = SmartsPatternMatcher()
if matcher.match(substrate, protocol.applicability.core):
    # Substrate is compatible with protocol
    ...
```

### **2. Testability**
Each module can be tested independently:

```python
# tests/test_substrate_classifier.py
def test_classify_primary_alkyl_iodide():
    classifier = SubstrateClassifier()
    info = classifier.classify("CCCCI")
    assert info.substrate_class == 'primary_alkyl_halide'
    assert info.special_positions.benzylic == []
    assert info.special_positions.allylic == []

# tests/test_smarts_builders.py  
def test_build_halide_smarts():
    builder = SmartsBuilder()
    mol = Chem.MolFromSmiles("CCCCI")
    iodine = mol.GetAtomWithIdx(4)
    smarts = builder.build_halide_smarts(iodine, mol, context)
    assert smarts == "[CX4;H2,H3]-[I]"
```

### **3. Maintainability**
Chemistry logic is centralized, not scattered:

```
❌ BAD (current):
  - Chemistry logic mixed with CLI code
  - Hard to find and update
  - Can't reuse in other modules

✅ GOOD (refactored):
  - Chemistry logic in chemtools/util/*
  - Clear separation of concerns
  - Easy to update and extend
```

### **4. Documentation**
Clear module boundaries with specific purposes:

```
chemtools/util/substrate_classifier.py  → "What is this molecule?"
chemtools/util/smarts_builders.py       → "How do I match this chemistry?"
chemtools/util/functional_groups.py     → "What functional groups are present?"
chemtools/protocol/smarts_generator_cli.py → "Generate protocol scope patterns"
```

---

## 📋 Implementation Plan (Refactored)

### **Phase 1: Create Reusable Utilities** (Week 1)

#### **Day 1-2: SubstrateClassifier**
**File**: `chemtools/util/substrate_classifier.py`

```python
# Core classification logic
- classify() - main entry point
- get_carbon_types() - sp3, sp2, sp, aromatic
- find_special_positions() - benzylic, allylic, propargylic
- analyze_heteroatom_environment() - what's around N, O, halide, etc.
- _is_benzylic(), _is_allylic(), _is_propargylic() - helper methods

# Uses existing:
- functional_groups.detect_all()
- rdkit_helpers.parse_smiles()
```

**Tests**: `tests/test_substrate_classifier.py` (50+ test cases)

#### **Day 3-4: SmartsBuilder**
**File**: `chemtools/util/smarts_builders.py`

```python
# SMARTS generation for each functional group family
- build_halide_smarts() - alkyl, aryl, vinyl, benzylic, allylic
- build_amine_smarts() - aniline, aliphatic, amide distinction
- build_carbonyl_smarts() - acid, ester, amide, ketone, aldehyde
- build_alcohol_smarts() - phenol, benzylic, aliphatic
- build_boron_smarts() - boronic acid, Bpin, BF3K
- build_for_substrate() - orchestrates builders

# Pattern matching utilities
- SmartsPatternMatcher.match()
- SmartsPatternMatcher.explain_match()
```

**Tests**: `tests/test_smarts_builders.py` (50+ test cases)

#### **Day 5: Enhance functional_groups.py**
**File**: `chemtools/util/functional_groups.py` (add to existing)

```python
# ADD classification context helper
def classify_functional_group_context(mol, fg_name: str) -> str:
    """Classify the chemical context of a functional group"""
```

**Tests**: Add to `tests/test_functional_groups.py`

---

### **Phase 2: Refactor CLI Application** (Week 2)

#### **Day 6-7: Refactor smarts_generator_cli.py**
**File**: `chemtools/protocol/smarts_generator_cli.py`

```python
# REMOVE chemistry logic (move to util/)
# KEEP only:
- CLI argument parsing
- Reaction SMILES parsing
- Visualization generation
- JSON output formatting
- Guard pattern suggestions (using classifier results)

# NOW USES:
from chemtools.util.substrate_classifier import SubstrateClassifier
from chemtools.util.smarts_builders import SmartsBuilder
```

**Tests**: Update `tests/test_smarts_generator.py`

---

### **Phase 3: Enable Reuse** (Week 3)

#### **Day 8: Add to featurizers**
Show how other modules can use the new utilities:

```python
# chemtools/featurizers/custom_cn_coupling.py
from chemtools.util.substrate_classifier import SubstrateClassifier

def extract_substrate_features(mol):
    classifier = SubstrateClassifier()
    info = classifier.classify(mol)
    
    features = {
        'is_benzylic': len(info.special_positions.benzylic) > 0,
        'is_allylic': len(info.special_positions.allylic) > 0,
        'substrate_class': info.substrate_class,
        # ... more features
    }
    return features
```

#### **Day 9: Add to recommendation engine**
```python
# chemtools/recommend/substrate_validator.py
from chemtools.util.smarts_builders import SmartsPatternMatcher

def validate_substrate_for_protocol(substrate_smiles, protocol):
    matcher = SmartsPatternMatcher()
    
    # Check if substrate matches protocol scope
    if not matcher.match(substrate_smiles, protocol.applicability.core):
        return False, matcher.explain_match(substrate_smiles, protocol.applicability.core)
    
    # Check guard patterns (forbidden)
    for guard in protocol.applicability.guards_forbid:
        if matcher.match(substrate_smiles, guard):
            return False, f"Substrate matches forbidden pattern: {guard}"
    
    return True, "Substrate is compatible"
```

#### **Day 10: Documentation & Examples**
- `docs/SUBSTRATE_CLASSIFIER_GUIDE.md`
- `docs/SMARTS_BUILDERS_GUIDE.md`
- `docs/REUSING_CHEMISTRY_UTILS.md`

---

## 🎓 Example Usage Across Modules

### **Use Case 1: Protocol Scope Definition (CLI)**
```python
# chemtools/protocol/smarts_generator_cli.py
from chemtools.util.substrate_classifier import SubstrateClassifier
from chemtools.util.smarts_builders import SmartsBuilder

classifier = SubstrateClassifier()
builder = SmartsBuilder()

info = classifier.classify("CCCCCCI")
smarts = builder.build_for_substrate(info)
# Result: [CX4;H2,H3]-[I]
```

### **Use Case 2: ML Feature Extraction**
```python
# chemtools/featurizers/advanced_featurizer.py
from chemtools.util.substrate_classifier import SubstrateClassifier

classifier = SubstrateClassifier()
info = classifier.classify(substrate_smiles)

features = {
    'has_benzylic_position': len(info.special_positions.benzylic) > 0,
    'has_allylic_position': len(info.special_positions.allylic) > 0,
    'carbon_sp3_count': sum(1 for t in info.carbon_types.values() if t == 'sp3'),
    'substrate_class': info.substrate_class,
}
```

### **Use Case 3: Protocol Recommendation**
```python
# chemtools/recommend/protocol_matcher.py
from chemtools.util.smarts_builders import SmartsPatternMatcher

matcher = SmartsPatternMatcher()

compatible_protocols = []
for protocol in protocol_database:
    if matcher.match(substrate, protocol.applicability.core):
        compatible_protocols.append(protocol)
```

### **Use Case 4: Reaction Type Detection**
```python
# chemtools/reaction_type_detector.py
from chemtools.util.substrate_classifier import SubstrateClassifier

classifier = SubstrateClassifier()
reactant_info = classifier.classify(reactants[0])

if reactant_info.substrate_class in ['aryl_halide', 'aryl_triflate']:
    if 'amine' in reactant_info.functional_groups:
        reaction_type = 'Ullmann CN Coupling'
```

---

## ✅ Summary

### **What Goes Where**

| Component | Location | Purpose | Reusable? |
|-----------|----------|---------|-----------|
| **Substrate Classification** | `chemtools/util/substrate_classifier.py` | Classify molecules by chemistry | ✅ YES |
| **SMARTS Pattern Building** | `chemtools/util/smarts_builders.py` | Build context-aware SMARTS | ✅ YES |
| **Functional Group Detection** | `chemtools/util/functional_groups.py` | Detect 80+ functional groups | ✅ ALREADY EXISTS |
| **RDKit Helpers** | `chemtools/util/rdkit_helpers.py` | RDKit wrappers | ✅ ALREADY EXISTS |
| **CLI Tool** | `chemtools/protocol/smarts_generator_cli.py` | Protocol scope definition | ❌ Application-specific |

### **Benefits**

1. ✅ **Reusable** - Other modules can use the same chemistry intelligence
2. ✅ **Testable** - Each utility module tested independently
3. ✅ **Maintainable** - Chemistry logic centralized
4. ✅ **Discoverable** - Clear module boundaries
5. ✅ **Extensible** - Easy to add new substrate types or SMARTS builders

### **Next Steps**

**Start with**: `chemtools/util/substrate_classifier.py`
- Most reusable
- Needed by both SMARTS builder and other modules
- Clear, well-defined scope
- Can be tested immediately

**Then**: `chemtools/util/smarts_builders.py`
- Uses SubstrateClassifier results
- Provides pattern generation for all modules

**Finally**: Refactor CLI to use the new utilities

---

**This is the right architecture!** Agree? Should I start implementing `SubstrateClassifier`? 🚀
