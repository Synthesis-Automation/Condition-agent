# Available Reaction Datasets for Analysis

## Dataset Overview

| Dataset | Reactions | Size (KB) | Reaction Type | Status |
|---------|-----------|-----------|---------------|--------|
| **Amide Formation** | 41,427 | 40,795 | Amide coupling | 🔥 **LARGEST** |
| **Suzuki Combined** | 50,215 | 41,984 | C-C coupling | 🔥 **LARGEST** |
| **Ullmann** | 5,552 | 4,151 | C-N/C-O coupling (Cu) | ✅ Medium |
| **Buchwald** | 1,343 | 1,085 | C-N coupling (Pd) | ✅ **DONE** |
| **Ni Amination** | 1,131 | 798 | C-N coupling (Ni) | ✅ Small |

---

## Priority Recommendations

### 🥇 **Option 1: Suzuki Coupling (50,215 reactions)**

**Why this is excellent:**
- **Massive dataset** - 50k reactions = robust ML training
- **Well-studied reaction** - C-C bond formation via Suzuki-Miyaura
- **Rich condition variety** - bases, solvents, catalysts
- **High yields** - generally reliable reactions (good for ML)

**Data Structure:**
```json
{
  "reaction_type": "Suzuki",
  "condition_core": "Base: Et3N",
  "reagents": [...],
  "solvents": [...],
  "conditions": {"temperature_c": null, "time_h": null, "yield_pct": 84.0},
  "smiles": {"reactants": "...", "products": "..."}
}
```

**Challenges:**
- Very large dataset (may need sampling or incremental training)
- Missing temperature/time data in many reactions
- Boron reagent handling (boronic acids, esters, MIDA boronates)

---

### 🥈 **Option 2: Amide Formation (41,427 reactions)**

**Why this is excellent:**
- **Massive dataset** - 41k reactions
- **Critical reaction** - amide bond = most common in pharmaceuticals
- **Diverse conditions** - coupling reagents (HATU, EDC, HBTU, PyBOP), bases, solvents
- **Varied substrates** - carboxylic acids + amines (primary, secondary, anilines)

**Data Structure:**
```json
{
  "reaction_type": "Amide formation",
  "condition_core": "Base: DIPEA",
  "reagents": [
    {"name": "HATU", "cas": "148893-10-1", "role": "COUPLING_REAGENT"},
    {"name": "N,N-Diisopropylethylamine", "cas": "7087-68-5", "role": "BASE"}
  ],
  "conditions": {"yield_pct": 89.0}
}
```

**Unique Features:**
- **Rule-based system already exists** (`chemtools/recommend.py` has amide-specific logic!)
- Can validate ML against existing expert rules
- Diverse coupling reagent vocabulary

**Challenges:**
- Coupling reagent selection is key (HATU vs EDC vs HBTU vs DCC)
- Protecting group complexity
- Racemization concerns for chiral acids

---

### 🥉 **Option 3: Ullmann C-N Coupling (5,552 reactions)**

**Why this is good:**
- **Medium-sized** - manageable for quick iteration
- **Complementary to Buchwald** - both are C-N, but Cu vs Pd
- **Simpler catalysts** - Cu salts, fewer ligand variations
- **Different substrates** - often heteroaryls, nucleophiles

**Data Structure:**
```json
{
  "reaction_type": "Ullman",
  "condition_core": "Cu",
  "reagents": [{"name": "Cesium carbonate", "cas": "534-17-8", "role": "BASE"}],
  "conditions": {"yield_pct": 75.0}
}
```

**Comparison with Buchwald:**
| Aspect | Buchwald (Pd) | Ullmann (Cu) |
|--------|---------------|--------------|
| Catalyst | Pd + ligands | Cu salts |
| Conditions | 80-120°C | 80-150°C |
| Substrates | Aryl halides | Aryl halides |
| Ligands | Complex (XPhos, etc.) | Simple (none or phenanthroline) |

---

### 🔬 **Option 4: Ni-Catalyzed Amination (1,131 reactions)**

**Why this is interesting:**
- **Emerging chemistry** - Ni catalysis is modern alternative to Pd
- **Different mechanism** - Ni(0)/Ni(II) vs Pd(0)/Pd(II)
- **Cost advantage** - Ni is cheaper than Pd
- **Manageable size** - quick to train and test

**Challenges:**
- Small dataset (may have limited diversity)
- Less mature than Buchwald/Ullmann

---

## 📊 Recommended Analysis Workflow

### **Phase 1: Quick Win - Ullmann (5,552 reactions)**

**Why start here:**
- ✅ Similar size to Buchwald (5.5k vs 1.3k)
- ✅ Can reuse C-N coupling infrastructure
- ✅ Validation framework already built
- ✅ 2-3 days to complete

**Tasks:**
1. Extract Ullmann reactions
2. Train DRFP model (reuse Buchwald code)
3. Test and verify (reuse verification script)
4. Compare Pd vs Cu recommendations

---

### **Phase 2: Big Impact - Amide Formation (41k reactions)**

**Why this is high priority:**
- ✅ Largest pharmaceutical relevance
- ✅ Can validate against existing rule-based system
- ✅ Rich coupling reagent vocabulary
- ✅ 1-2 weeks to complete

**Tasks:**
1. Analyze coupling reagent distribution
2. Feature engineering for acid/amine classes
3. Train multi-output model (yield + coupling reagent + base + solvent)
4. Validate against `chemtools/recommend.py` amide rules

---

### **Phase 3: Massive Scale - Suzuki (50k reactions)**

**Why save for later:**
- ⚠️ Very large (50k reactions)
- ⚠️ May need sampling strategy
- ⚠️ Boron chemistry adds complexity
- ⚠️ 2-3 weeks to complete

**Tasks:**
1. Sample representative subset (5k-10k reactions)
2. Handle boron reagent types (boronic acid vs ester vs MIDA)
3. Multi-condition prediction
4. Incremental learning approach

---

## 🎯 My Recommendation: Start with Ullmann

**Why Ullmann is the perfect next step:**

1. **Quick validation** - Use existing verification framework
2. **Complementary chemistry** - C-N coupling like Buchwald, but Cu instead of Pd
3. **Manageable scope** - 5.5k reactions = 4x Buchwald size
4. **Interesting comparison** - Can we find patterns where Cu > Pd or vice versa?

**Expected Results:**
- Train Ullmann DRFP model in 1-2 hours
- Test on sample reactions
- Verify with rule-based system
- Generate insights on Cu vs Pd catalyst selection

**Timeline:** 2-3 days to complete end-to-end

---

## 🚀 Next Steps

**Option A: Ullmann (Recommended)**
```bash
# 1. Create Ullmann test script (similar to Buchwald)
python scripts/train_ullmann_drfp.py

# 2. Test on sample reactions
python scripts/test_ullmann_reactions.py

# 3. Verify predictions
python scripts/verify_ullmann_ml_with_rules.py
```

**Option B: Amide Formation (High Impact)**
```bash
# 1. Analyze coupling reagent distribution
python scripts/analyze_amide_reagents.py

# 2. Train amide DRFP model
python scripts/train_amide_drfp.py

# 3. Validate against existing rules
python scripts/verify_amide_ml_with_rules.py
```

**Option C: Suzuki (Big Data)**
```bash
# 1. Sample representative subset
python scripts/sample_suzuki_dataset.py --n=10000

# 2. Train on sampled data
python scripts/train_suzuki_drfp.py

# 3. Test and verify
python scripts/test_suzuki_reactions.py
```

---

## 💡 Key Questions to Answer

### For Ullmann:
1. **Cu vs Pd:** When does Cu outperform Pd for C-N coupling?
2. **Ligand-free:** Can simple Cu salts match Pd/ligand performance?
3. **Substrate scope:** Which electrophiles prefer Cu?

### For Amide Formation:
1. **Reagent selection:** Can ML predict optimal coupling reagent (HATU vs EDC)?
2. **Rule validation:** Does ML agree with expert amide formation rules?
3. **Racemization:** Can ML predict conditions that avoid racemization?

### For Suzuki:
1. **Boron reagent:** Does boronic acid vs ester vs MIDA matter?
2. **Base selection:** Can ML predict optimal base (K2CO3 vs Cs2CO3)?
3. **Scale:** Does model performance improve with 50k reactions vs 5k?

---

## 📁 Files Ready to Create

1. **Ullmann:** `scripts/train_ullmann_drfp.py`, `scripts/test_ullmann_reactions.py`
2. **Amide:** `scripts/train_amide_drfp.py`, `scripts/analyze_amide_reagents.py`
3. **Suzuki:** `scripts/sample_suzuki_dataset.py`, `scripts/train_suzuki_drfp.py`

---

## ✅ Conclusion

**I recommend starting with Ullmann** for quick validation, then moving to **Amide Formation** for high pharmaceutical impact, and finally **Suzuki** for big data exploration.

**Which dataset would you like to work on first?**
