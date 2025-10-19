# Chemistry-Aware SMARTS Pattern Generation - Comprehensive Plan

## 🎯 Problem Statement

The current SMARTS generator is **too simplistic** - it doesn't understand chemical context. We need it to recognize:

### Current Issues
1. **No substrate class recognition**: Can't distinguish alkyl vs aryl vs vinyl
2. **No functional group context**: Doesn't recognize amides, esters, alcohols, etc.
3. **No chemical logic**: Treats all C-I bonds the same (alkyl, aryl, benzylic should be different)
4. **No reaction context**: Pattern should reflect the actual chemistry happening

### What We Need
Generate patterns that capture **chemical meaning**:

```
If reaction uses:           Pattern should specify:
------------------------------------------------------------
Primary alkyl iodide    →   [CX4;H2,H3]-I  (sp³, not aromatic)
Aryl bromide           →   c-Br           (aromatic carbon)
Benzylic chloride      →   [CH2]-[c]-Cl   (next to aromatic)
Allylic iodide         →   [CH2]-C=C-I    (next to double bond)
Aniline                →   c-[NH2]        (aromatic amine)
Aliphatic amine        →   [CX4]-[NH2]    (sp³ carbon)
Amide                  →   N-C(=O)        (carbonyl, not amine)
Ester                  →   O-C(=O)        (carbonyl ester, not acid)
```

---

## 🔬 Proposed Architecture

### Phase 1: Substrate Classification System

Create a **ChemicalContextAnalyzer** that identifies substrate types:

```python
class SubstrateAnalyzer:
    """Analyzes substrates to understand chemical context"""
    
    def analyze_molecule(self, mol) -> SubstrateContext:
        """Return comprehensive substrate classification"""
        return {
            'carbon_type': 'sp3' | 'sp2' | 'sp' | 'aromatic',
            'functional_groups': ['amine', 'alcohol', 'halide', ...],
            'special_positions': {
                'benzylic': [atom_idx, ...],
                'allylic': [atom_idx, ...],
                'propargylic': [atom_idx, ...],
            },
            'hybridization_map': {atom_idx: 'sp3|sp2|sp|aromatic'},
            'heteroatom_environment': {...}
        }
```

**Leverage Existing Code:**
- `chemtools/util/functional_groups.py` - 80+ functional group SMARTS patterns
- `chemtools/features/role/smarts.py` - Reactive center detection
- `chemtools/featurizers/ullmann.py` - Electrophile classification
- `chemtools/features/role/role_feats/` - Aryl halide, amine features

### Phase 2: Context-Aware Pattern Generator

Upgrade `_mol_to_generic_smarts()` to be **chemistry-aware**:

```python
def _generate_context_aware_smarts(self, mol, role='reactant') -> str:
    """Generate SMARTS that captures chemical essence"""
    
    # 1. Identify all functional groups
    fg_detector = FunctionalGroupDetector()
    functional_groups = fg_detector.detect_all(mol)
    
    # 2. Classify carbon types (sp3, sp2, aromatic)
    carbon_types = self._classify_carbons(mol)
    
    # 3. Find reactive centers
    reactive_centers = self._find_reactive_centers(mol, functional_groups)
    
    # 4. Analyze chemical environment
    for center_idx in reactive_centers:
        atom = mol.GetAtomWithIdx(center_idx)
        
        # Determine specificity level
        if self._is_halide(atom):
            return self._generate_halide_smarts(atom, mol)
        elif self._is_amine(atom):
            return self._generate_amine_smarts(atom, mol)
        elif self._is_carbonyl(atom):
            return self._generate_carbonyl_smarts(atom, mol)
        # ... more functional group handlers
    
    return self._generate_generic_smarts(mol)
```

### Phase 3: Functional Group-Specific SMARTS Generators

Create specialized generators for each functional group family:

#### **Halide Handler**
```python
def _generate_halide_smarts(self, halide_atom, mol) -> str:
    """Context-aware halide SMARTS"""
    X = halide_atom.GetSymbol()  # Cl, Br, I, F
    
    # Find carbon neighbor
    carbon = self._get_carbon_neighbor(halide_atom)
    if not carbon:
        return f"[{X}]"
    
    # Check carbon environment
    if carbon.GetIsAromatic():
        # Aryl halide: c-X
        return f"c-[{X}]"
    
    elif self._is_benzylic(carbon, mol):
        # Benzylic: [CH2]-c-X or similar
        h_count = carbon.GetTotalNumHs()
        return f"[C;H{h_count}]-[c]-[{X}]"
    
    elif self._is_allylic(carbon, mol):
        # Allylic: [CH2]-C=C or C=C-[CH2]-X
        return f"[C;H2]-[C]=[C]-[{X}]"
    
    elif self._is_propargylic(carbon, mol):
        # Propargylic: C#C-[CH2]-X
        return f"[C]#[C]-[C;H2]-[{X}]"
    
    else:
        # Aliphatic (sp3)
        h_count = carbon.GetTotalNumHs()
        if h_count >= 2:
            return f"[CX4;H2,H3]-[{X}]"  # Primary
        elif h_count == 1:
            return f"[CX4;H1]-[{X}]"      # Secondary
        else:
            return f"[CX4;H0]-[{X}]"      # Tertiary
```

#### **Amine Handler**
```python
def _generate_amine_smarts(self, n_atom, mol) -> str:
    """Context-aware amine SMARTS"""
    
    # Check if amide (exclude from amine)
    if self._is_amide(n_atom, mol):
        return "[N]-[C](=O)"
    
    # Check if aniline
    if self._is_aniline(n_atom, mol):
        h_count = n_atom.GetTotalNumHs()
        if h_count == 2:
            return "c-[NH2]"         # Primary aniline
        elif h_count == 1:
            return "c-[NH1]"         # Secondary aniline
        else:
            return "c-[N]"           # Tertiary aniline
    
    # Aliphatic amine
    h_count = n_atom.GetTotalNumHs()
    if h_count == 2:
        return "[CX4]-[NH2]"         # Primary aliphatic
    elif h_count == 1:
        return "[NX3;H1]"            # Secondary
    else:
        return "[NX3;H0]"            # Tertiary
```

#### **Carbonyl Handler**
```python
def _generate_carbonyl_smarts(self, c_atom, mol) -> str:
    """Context-aware carbonyl SMARTS"""
    
    neighbors = c_atom.GetNeighbors()
    has_o_double = any(n.GetSymbol() == 'O' and 
                       mol.GetBondBetweenAtoms(c_atom.GetIdx(), n.GetIdx()).GetBondType() == 2
                       for n in neighbors)
    
    if not has_o_double:
        return "[C]"
    
    # Check what's attached to C=O
    other_neighbors = [n for n in neighbors if n.GetSymbol() != 'O']
    
    # Carboxylic acid: C(=O)OH
    if any(n.GetSymbol() == 'O' and n.GetTotalNumHs() == 1 for n in neighbors):
        return "[C](=O)[OH]"
    
    # Ester: C(=O)O[C]
    if any(n.GetSymbol() == 'O' and n.GetTotalNumHs() == 0 for n in neighbors):
        return "[C](=O)[O]-[C]"
    
    # Amide: C(=O)N
    if any(n.GetSymbol() == 'N' for n in neighbors):
        return "[C](=O)-[N]"
    
    # Acyl chloride: C(=O)Cl
    if any(n.GetSymbol() in ['Cl', 'Br', 'I', 'F'] for n in neighbors):
        X = next(n.GetSymbol() for n in neighbors if n.GetSymbol() in ['Cl', 'Br', 'I', 'F'])
        return f"[C](=O)-[{X}]"
    
    # Aldehyde or ketone
    if c_atom.GetTotalNumHs() == 1:
        return "[C;H1](=O)"          # Aldehyde
    else:
        return "[C](=O)"             # Ketone
```

### Phase 4: Intelligent Guard Pattern Suggestions

Make guard patterns **context-aware** too:

```python
def suggest_context_aware_guards(self, substrate_analysis) -> List[str]:
    """Suggest guards based on actual substrate type"""
    
    guards = []
    
    # If primary alkyl halide used, exclude secondary/tertiary/benzylic/allylic
    if substrate_analysis['type'] == 'primary_alkyl_halide':
        X = substrate_analysis['halogen']
        guards.extend([
            f"[CX4;H1]-[{X}]  # Exclude secondary",
            f"[CX4;H0]-[{X}]  # Exclude tertiary",
            f"[CH2;$([CH2][c])]-[{X}]  # Exclude benzylic",
            f"[CH2;$([CH2]C=C)]-[{X}]  # Exclude allylic",
        ])
    
    # If aryl halide used, exclude vinyl/alkyl
    elif substrate_analysis['type'] == 'aryl_halide':
        X = substrate_analysis['halogen']
        guards.extend([
            f"[CX4]-[{X}]  # Exclude aliphatic halides",
            f"[CX3]=[CX3]-[{X}]  # Exclude vinyl halides",
        ])
    
    # If aniline used, exclude aliphatic amines
    elif substrate_analysis['type'] == 'aniline':
        guards.extend([
            "[CX4]-[NH2]  # Exclude aliphatic primary amines",
            "[CX4]-[NH1]  # Exclude aliphatic secondary amines",
        ])
    
    # If amide used, exclude simple amines
    elif substrate_analysis['type'] == 'amide':
        guards.extend([
            "[NH2;!$(N-C(=O))]  # Exclude non-amide amines",
        ])
    
    return guards
```

---

## 🛠️ Implementation Plan

### **Step 1: Create Substrate Analysis Module** ✅ Easy
**File:** `chemtools/protocol/substrate_analyzer.py`

- Integrate with existing `functional_groups.py`
- Add carbon hybridization detection
- Add special position detection (benzylic, allylic, propargylic)
- Add heteroatom environment analysis

**Time:** 2-3 hours

### **Step 2: Create Functional Group Handlers** ⚙️ Medium
**File:** `chemtools/protocol/smarts_context_generators.py`

- Halide handler (alkyl, aryl, vinyl, benzylic, allylic)
- Amine handler (aniline, aliphatic, amide)
- Carbonyl handler (acid, ester, amide, ketone, aldehyde)
- Alcohol/phenol handler
- Boron handler (boronic acid, boronate ester, BF3K)
- Misc handlers (thiol, sulfonic acid, etc.)

**Time:** 4-6 hours

### **Step 3: Upgrade SmartsGenerator** 🔧 Medium
**File:** `chemtools/protocol/smarts_generator_cli.py`

- Replace `_mol_to_generic_smarts()` with `_generate_context_aware_smarts()`
- Integrate substrate analyzer
- Route to appropriate functional group handler
- Update guard pattern suggestions

**Time:** 3-4 hours

### **Step 4: Comprehensive Testing** 🧪 Critical
**File:** `tests/test_smarts_context_aware.py`

Test cases covering:
- ✅ Primary alkyl iodide → `[CX4;H2,H3]-I`
- ✅ Aryl bromide → `c-Br`
- ✅ Benzylic chloride → `[CH2;$([CH2][c])]-Cl`
- ✅ Aniline → `c-[NH2]`
- ✅ Amide → `[N]-C(=O)`
- ✅ Boronic acid → `[B](O)(O)`
- ✅ Phenol → `c-[OH]`
- ✅ Allylic alcohol → `C=C-[CH2]-[OH]`

**Time:** 3-4 hours

### **Step 5: Documentation & Examples** 📚
**Files:** 
- `docs/SMARTS_CHEMISTRY_AWARE_GUIDE.md`
- `docs/SMARTS_CONTEXT_EXAMPLES.md`

**Time:** 2 hours

---

## 📊 Expected Outcomes

### Before (Current Simple Approach)
```python
# Input: CCCCCCCCI >> CCCCCCCCB
# Output: [C;H2,H3]-[I]>>[C;H2,H3]-[B]
# Problem: Too generic, doesn't specify sp3 carbon
```

### After (Chemistry-Aware Approach)
```python
# Input: CCCCCCCCI >> CCCCCCCCB
# Analysis: Primary alkyl iodide (sp3 carbon, not aromatic, not benzylic)
# Output: [CX4;H2,H3]-[I]>>[CX4;H2,H3]-[B]
# Guards: [CX4;H1]-I, [CX4;H0]-I, [CH2;$([CH2][c])]-I, [CH2;$([CH2]C=C)]-I
```

```python
# Input: c1ccccc1Br >> c1ccccc1B(O)O
# Analysis: Aryl bromide (aromatic carbon)
# Output: c-[Br]>>c-[B](O)(O)
# Guards: [CX4]-Br, [CX3]=[CX3]-Br  # Exclude alkyl & vinyl
```

```python
# Input: c1ccccc1N >> c1ccccc1NC(=O)C
# Analysis: Aniline (aromatic amine, not amide)
# Output: c-[NH2;!$(N-C(=O))]>>c-[N]-[C](=O)
# Guards: [CX4]-[NH2]  # Exclude aliphatic amines
```

---

## 🎓 Key Principles

### 1. **Substrate Class First**
Always determine the **class** before generating pattern:
- Alkyl vs Aryl vs Vinyl
- Primary vs Secondary vs Tertiary
- Aliphatic vs Aromatic vs Heteroaromatic

### 2. **Functional Group Context**
Understand the **role** of heteroatoms:
- Amine vs Amide vs Aniline (all have N, different chemistry)
- Alcohol vs Ether vs Ester (all have O, different chemistry)
- Alkyl halide vs Aryl halide vs Acyl halide

### 3. **Special Positions Matter**
Recognize **activated positions**:
- Benzylic (next to aromatic) - more reactive
- Allylic (next to C=C) - more reactive
- Propargylic (next to C≡C) - special chemistry

### 4. **Exclude by Default**
When pattern specifies primary alkyl:
- Explicitly exclude secondary, tertiary
- Explicitly exclude benzylic, allylic
- Explicitly exclude aromatic

---

## 🚀 Migration Strategy

### Phase A: Add Without Breaking (Week 1)
1. Create new modules (`substrate_analyzer.py`, `smarts_context_generators.py`)
2. Add tests for new functionality
3. Keep old `_mol_to_generic_smarts()` as fallback

### Phase B: Integrate (Week 2)
1. Add `--chemistry-aware` flag to CLI
2. Route to new logic when flag enabled
3. Compare old vs new output in tests
4. Document differences

### Phase C: Make Default (Week 3)
1. Switch default to chemistry-aware mode
2. Keep `--simple-mode` flag for old behavior
3. Update all documentation
4. Announce change to users

---

## ✅ Success Criteria

1. **Correct Classification**: 95%+ accuracy on substrate type detection
2. **Meaningful Patterns**: Patterns specify chemical context, not just atoms
3. **Accurate Guards**: Guard patterns correctly exclude incompatible substrates
4. **Test Coverage**: 100 test cases covering major functional groups
5. **User Satisfaction**: Patterns work out-of-box for real protocols

---

## 📝 Example Test Cases

```python
def test_primary_alkyl_iodide():
    gen = SmartsGenerator("CCCCCCI>>CCCCCB")
    pattern = gen.generate_context_aware()
    assert pattern['core'] == "[CX4;H2,H3]-[I]>>[CX4;H2,H3]-[B]"
    assert "[CX4;H1]-[I]" in pattern['guards_forbid']  # No secondary
    
def test_aryl_bromide():
    gen = SmartsGenerator("c1ccccc1Br>>c1ccccc1B(O)O")
    pattern = gen.generate_context_aware()
    assert pattern['core'] == "c-[Br]>>c-[B](O)(O)"
    assert "[CX4]-[Br]" in pattern['guards_forbid']  # No alkyl
    
def test_aniline():
    gen = SmartsGenerator("c1ccccc1N>>c1ccccc1NC(=O)C")
    pattern = gen.generate_context_aware()
    assert "c-[NH2" in pattern['core']
    assert "[CX4]-[NH2]" in pattern['guards_forbid']  # No aliphatic
    
def test_amide_not_amine():
    gen = SmartsGenerator("CC(=O)N>>CC(=O)NCC")
    pattern = gen.generate_context_aware()
    assert "[N]-[C](=O)" in pattern['core']
    assert "[NH2;!$(N-C(=O))]" in pattern['guards_forbid']  # Exclude simple amines
```

---

## 🎯 Priority Order

1. **Halides** (most common) - alkyl, aryl, vinyl, benzylic, allylic
2. **Amines** (very common) - aniline, aliphatic, amide distinction
3. **Carbonyls** (complex) - acid, ester, amide, ketone, aldehyde
4. **Alcohols** (common) - aliphatic, benzylic, allylic
5. **Boron** (growing importance) - boronic acid, Bpin, BF3K
6. **Others** (as needed) - thiols, sulfonates, phosphines

---

## 💡 This Is The Right Approach Because:

1. ✅ **Leverages existing code** - `functional_groups.py`, featurizers
2. ✅ **Modular design** - each functional group is a separate handler
3. ✅ **Testable** - clear expected outputs for each substrate type
4. ✅ **Extensible** - easy to add new functional group handlers
5. ✅ **Backward compatible** - can keep simple mode as fallback
6. ✅ **Chemistry-correct** - patterns reflect actual chemical meaning

---

**Ready to implement?** We start with Step 1: Substrate Analysis Module! 🚀
