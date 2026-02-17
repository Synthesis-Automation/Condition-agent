# Phase 1: Experiments Bridge - REVISED PLAN

## Critical Discoveries from HTERecommender Testing

### Test Results Summary

```
TEST: Suzuki Coupling (Brc1ccccc1 + B(O)(O)c1ccccc1 >> biphenyl)

WITH Product SMILES:
✅ Detected reaction type: "Suzuki_miyaura" (0.95 confidence)
✅ Reacted motifs: ['Ar-B(OH)2', 'Ar-Br']
✅ Formed motifs: ['Ar-Ar']
✅ Recommendations: 5 conditions found
✅ Reaction key: |Ar-B(OH)2|Ar-Br -> Ar-Ar | events: LGDisp+C-C

WITHOUT Product SMILES:
❌ Predicted reaction type: None
❌ Recommendations: 0 conditions found
❌ Cannot match database
```

### **CRITICAL REQUIREMENT: Product SMILES is MANDATORY**

The HTERecommender **requires product SMILES** to:
1. Call `featurize_reaction()` internally
2. Detect reaction type automatically (0.95 confidence for Suzuki)
3. Identify reacted_motifs and formed_motifs for matching
4. Generate reaction_key for database lookup

**Without product**: Zero recommendations, no reaction type detection.

---

## HTERecommender Internal Workflow (Discovered)

```
Input: reactant_a_smiles, reactant_b_smiles, product_smiles

Step 1: Construct full reaction SMILES
    reaction_smiles = f"{reactant_a}.{reactant_b}>>{product}"

Step 2: Call featurize_reaction() [KEY STEP]
    → Detects reaction_type: "Suzuki_miyaura"
    → Extracts reacted_motifs: ["Ar-Br", "Ar-B(OH)2"]
    → Extracts formed_motifs: ["Ar-Ar"]
    → Generates reaction_key: "|Ar-B(OH)2|Ar-Br -> Ar-Ar|..."

Step 3: Match against HTE database
    → Primary: reaction_key matching
    → Secondary: reaction_type + motif matching
    → Tertiary: motif-only matching (fallback)

Step 4: Rank and return top_k recommendations
    → Z-score based ranking
    → Success rate, yield statistics
```

**Key Insight**: The recommender doesn't just match SMILES strings - it uses the **taxonomy system** to extract motifs and reaction signatures, then matches those against the database.

---

## Taxonomy System Overview

### Motif Detection (V2 Taxonomy)

```python
featurize_molecule("Brc1ccccc1")
# Returns: motifs = [{"id": "Ar-Br", ...}]

featurize_molecule("B(O)(O)c1ccccc1")
# Returns: motifs = [{"id": "Ar-B(OH)2", ...}]
```

**Common motifs in HTE database:**
- `Ar-Br`, `Ar-Cl`, `Ar-I` (aryl halides)
- `Ar-B(OH)2`, `Ar-B(OR)2` (boronic acids/esters)
- `Ar-NH2`, `Ar-NHR`, `Ar-OH` (aryl nucleophiles)
- `Alkyl-NHR`, `R2CH-NH2` (alkyl amines)
- `Ar-H` (C-H activation substrates)

### Reaction Type Detection

```python
featurize_reaction("Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1")
# Returns: reaction_type = "Suzuki_miyaura" (confidence: 0.95)
```

**Database reaction types** (41 total):
- `Suzuki_miyaura` (not "Suzuki-Miyaura"!)
- `C_N_Coupling`
- `Amide_formation`
- `Snar`
- `C_H_arylation`
- `Buchwald_Hartwig`
- etc.

### Reaction Key Format

```
|Ar-B(OH)2|Ar-Br -> Ar-Ar | events: LGDisp+C-C | mech: oa_based_coupling
└─ reactants ─┘  └─ products ─┘  └─ reaction events ─┘  └─ mechanism ─┘
```

The reaction_key encodes:
- Reacted motifs (consumed in reaction)
- Formed motifs (created in reaction)
- Reaction events (bond changes, leaving group displacement)
- Mechanism classification

This enables hierarchical matching:
1. **Best**: Exact reaction_key match
2. **Good**: reaction_type + motif match
3. **Fallback**: Motif-only match

---

## Revised Bridge Implementation

### Gap Analysis

#### What Analysis Currently Provides

From `analyze_reaction_smiles()` result:
```python
{
    "input": {
        "rxn_smiles_clean": "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        # ✅ Has full reaction SMILES with product!
    },
    "quick_glance": {  # Tier 2
        "reaction_types": ["Suzuki-Miyaura cross-coupling"],  # Natural language
        ...
    },
    "interpretation": {  # Tier 3
        "overall_class": "cross_coupling",
        ...
    }
}
```

✅ **Good news**: `rxn_smiles_clean` contains full reaction SMILES including product!

#### What Bridge Needs to Do

1. **Parse rxn_smiles_clean** → Extract reactant_a, reactant_b, product
2. **Pass to HTERecommender** → Let it handle taxonomy/motif detection internally
3. **NO manual reaction type mapping needed!** → The recommender detects it automatically

**Simplified Bridge Logic:**

```python
def recommend_conditions(analysis_result):
    # Step 1: Parse SMILES
    rxn_smiles = analysis_result['input']['rxn_smiles_clean']
    parts = rxn_smiles.split('>>')
    reactants_str = parts[0]
    product_str = parts[1]

    reactants = reactants_str.split('.')
    reactant_a = reactants[0]
    reactant_b = reactants[1] if len(reactants) > 1 else None
    product = product_str.split('.')[0]  # First product

    # Step 2: Call HTERecommender (it handles everything internally!)
    return recommender.recommend(
        reactant_a_smiles=reactant_a,
        reactant_b_smiles=reactant_b,
        product_smiles=product,  # CRITICAL!
        top_k=10,
        source_group="experiments"
        # NO reaction_type_filter needed - it auto-detects!
    )
```

---

## Simplified Implementation Plan

### MAJOR SIMPLIFICATION

The original plan was too complex because I didn't understand the taxonomy system.

**OLD (incorrect) assumptions:**
- ❌ Need to manually classify reactants into motifs
- ❌ Need to map Tier 2 reaction types to database names
- ❌ Need complex fallback logic for matching

**NEW (correct) understanding:**
- ✅ HTERecommender handles ALL taxonomy/motif detection internally
- ✅ Just pass reactant_a, reactant_b, product SMILES
- ✅ Recommender auto-detects reaction type with 0.95 confidence
- ✅ Simple bridge: parse SMILES → call recommender → display

---

## Day 1: Create Bridge Module

**File**: `reaction_agent/condition_bridge.py`

```python
"""
Minimal bridge connecting reaction analysis to HTE experiments.

The HTERecommender handles all taxonomy/motif detection internally.
Bridge just parses SMILES and passes to recommender.
"""

from typing import Dict, Any, Optional, Tuple
from pathlib import Path
from chemtools.recommend.recommender import HTERecommender, HTERecommendationResult


class ConditionBridge:
    """
    Simple bridge to HTE experiments database.

    Leverages HTERecommender's internal taxonomy system.
    """

    def __init__(self, hte_db_path: str = "data/HTE_db"):
        self.hte_db_path = Path(hte_db_path)
        self.recommender = HTERecommender(str(self.hte_db_path))

    def extract_smiles(
        self,
        analysis_result: Dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Extract reactant and product SMILES from analysis.

        Args:
            analysis_result: From analyze_reaction_smiles()

        Returns:
            (reactant_a, reactant_b, product) SMILES
        """
        input_data = analysis_result.get('input', {})
        rxn_smiles_clean = input_data.get('rxn_smiles_clean', '')

        if not rxn_smiles_clean or '>>' not in rxn_smiles_clean:
            return None, None, None

        # Split reaction SMILES
        parts = rxn_smiles_clean.split('>>')
        reactants_str = parts[0] if len(parts) > 0 else ''
        products_str = parts[1] if len(parts) > 1 else ''

        # Parse reactants (split by '.')
        reactants = [r.strip() for r in reactants_str.split('.') if r.strip()]
        reactant_a = reactants[0] if len(reactants) > 0 else None
        reactant_b = reactants[1] if len(reactants) > 1 else None

        # Parse products (use first product)
        products = [p.strip() for p in products_str.split('.') if p.strip()]
        product = products[0] if len(products) > 0 else None

        return reactant_a, reactant_b, product

    def recommend_conditions(
        self,
        analysis_result: Dict[str, Any],
        top_k: int = 10,
        min_experiments: int = 2,
        reaction_type_filter: Optional[str] = None
    ) -> HTERecommendationResult:
        """
        Get condition recommendations from HTE experiments.

        The HTERecommender internally:
        - Calls featurize_reaction() to detect reaction type
        - Extracts reacted/formed motifs using taxonomy
        - Matches against database using reaction_key + motifs

        Args:
            analysis_result: From analyze_reaction_smiles()
            top_k: Number of recommendations
            min_experiments: Minimum experiments per condition
            reaction_type_filter: Optional override (usually not needed)

        Returns:
            HTERecommendationResult with ranked recommendations
        """
        # Extract SMILES
        reactant_a, reactant_b, product = self.extract_smiles(analysis_result)

        if not reactant_a:
            raise ValueError("Could not extract reactant SMILES from analysis")

        if not product:
            raise ValueError(
                "Product SMILES required for recommendation. "
                "Recommender needs full reaction for taxonomy detection."
            )

        # Call HTERecommender (handles all taxonomy internally)
        return self.recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,  # CRITICAL for reaction type detection
            top_k=top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,  # Usually None (auto-detect)
            source_group="experiments"  # Phase 1: experiments only
        )

    def analyze_and_recommend(
        self,
        rxn_smiles: str,
        analyzer: 'ReactionSMILESAnalyzer',
        top_k: int = 10,
        validate: bool = True
    ) -> Dict[str, Any]:
        """
        End-to-end: Analyze → Recommend.

        Args:
            rxn_smiles: Reaction SMILES
            analyzer: ReactionSMILESAnalyzer instance
            top_k: Number of recommendations
            validate: Enable Tier 4 validation

        Returns:
            Dict with analysis + recommendations + metadata
        """
        import time

        start_time = time.time()

        # Analyze reaction
        analysis_result = analyzer.analyze(rxn_smiles, validate=validate)
        analysis_time = time.time() - start_time

        # Get recommendations
        rec_start = time.time()
        recommendations = self.recommend_conditions(analysis_result, top_k=top_k)
        rec_time = time.time() - rec_start

        return {
            "analysis": analysis_result,
            "recommendations": recommendations,
            "metadata": {
                "analysis_time_s": analysis_time,
                "recommendation_time_s": rec_time,
                "total_time_s": time.time() - start_time,
                "detected_reaction_type": recommendations.predicted_reaction_type,
                "reaction_type_confidence": recommendations.reaction_type_confidence,
                "reacted_motifs": recommendations.reacted_motifs,
                "formed_motifs": recommendations.formed_motifs
            }
        }


__all__ = ['ConditionBridge']
```

**Key Simplifications:**
- ❌ No `REACTION_TYPE_MAP` needed (auto-detected)
- ❌ No `map_reaction_type()` method needed
- ❌ No complex fallback logic needed
- ✅ Just parse SMILES and call recommender
- ✅ Recommender handles all taxonomy internally

---

## Day 1 Afternoon: CLI Integration

Same as before, but even simpler since no reaction_type_filter needed.

**Add to cli.py:**

```python
# Around line 655 - Add flags
parser.add_argument(
    '--recommend',
    action='store_true',
    help='Recommend conditions from HTE experiments'
)
parser.add_argument(
    '--top-conditions',
    type=int,
    default=10,
    help='Number of recommendations (default: 10)'
)

# Display function (same as before)
def display_recommendations(result: HTERecommendationResult, top_k: int = 10):
    """Display recommendations from HTE database."""
    # ... same implementation as before ...

# Integrate into analyze_reaction_interactive()
if recommend:
    from reaction_agent.condition_bridge import ConditionBridge

    print_header("CONDITION RECOMMENDATIONS")

    # Show what was detected
    print(f"Detected reaction type: {result.recommendations.predicted_reaction_type}")
    print(f"Confidence: {result.recommendations.reaction_type_confidence:.2f}")
    print(f"Reacted motifs: {', '.join(result.recommendations.reacted_motifs)}")
    print(f"Formed motifs: {', '.join(result.recommendations.formed_motifs)}")
    print()

    display_recommendations(result.recommendations, top_conditions)
```

---

## Day 2: Testing

```bash
# Test 1: Simple Suzuki (should auto-detect)
python -m reaction_agent.cli \
  --reaction "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" \
  --recommend \
  --validate

# Expected:
# - Detected reaction type: Suzuki_miyaura (0.95 confidence)
# - Reacted motifs: Ar-Br, Ar-B(OH)2
# - Formed motifs: Ar-Ar
# - 5+ recommendations with Pd catalysts

# Test 2: Buchwald-Hartwig
python -m reaction_agent.cli \
  --reaction "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1" \
  --recommend

# Expected:
# - Detected reaction type: C_N_Coupling or Buchwald_Hartwig
# - Recommendations with Pd catalysts + ligands

# Test 3: Amide formation
python -m reaction_agent.cli \
  --reaction "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1" \
  --recommend

# Expected:
# - Detected reaction type: Amide_formation
# - Recommendations with coupling reagents (EDC, HATU, etc.)
```

---

## Summary of Key Changes from Original Plan

### What Changed

**BEFORE** (over-engineered):
1. Manual motif detection
2. Reaction type mapping table
3. Complex fallback logic
4. Multiple matching strategies

**AFTER** (simplified):
1. ✅ Parse SMILES (trivial)
2. ✅ Pass to HTERecommender (handles everything)
3. ✅ Display results
4. ✅ Done!

### Why So Much Simpler?

The HTERecommender **already has a sophisticated taxonomy system** built-in:
- `featurize_molecule()` for motif detection
- `featurize_reaction()` for reaction type detection
- Hierarchical matching (reaction_key → motifs → fallback)
- Z-score based ranking

**We just need to pass it SMILES strings!**

---

## Success Criteria (Unchanged)

✅ Can run: `python -m reaction_agent.cli --reaction "SMILES" --recommend`
✅ Shows top 10 conditions from experiments
✅ Displays catalyst, ligand, base, solvent, performance
✅ Works for common reactions (Suzuki, Buchwald-Hartwig, amide)
✅ Total time < 30 seconds

---

## Risk Assessment

### Risk 1: Product SMILES Missing
**Likelihood**: Low (analysis always provides clean reaction SMILES)
**Impact**: High (no recommendations without product)
**Mitigation**:
- Validate that rxn_smiles_clean has '>>' and products
- Clear error message if missing

### Risk 2: Multiple Products
**Likelihood**: Medium (some reactions have multiple products)
**Impact**: Low (we just use first product)
**Mitigation**:
- Document that first product is used
- Could enhance later to try all products

### Risk 3: Reaction Type Not in Database
**Likelihood**: Medium (41 reaction types, but chemistry is vast)
**Impact**: Medium (may get generic recommendations)
**Mitigation**:
- Recommender falls back to motif-only matching
- Still gets recommendations, just less specific

### Risk 4: Reactant Order Ambiguity
**Likelihood**: Low (order doesn't matter much)
**Impact**: Low (taxonomy handles swapped reactants)
**Mitigation**:
- Recommender tries both orders internally
- Will match regardless

---

## Files Summary

**NEW**:
- `reaction_agent/condition_bridge.py` (~150 lines, much simpler!)
- `test_recommender_workflow.py` (testing script)

**MODIFIED**:
- `reaction_agent/cli.py` (~100 lines added)

**TOTAL**: ~250 lines of simple code

**Timeline**: Still 2-3 days, but much less complex

---

## Next Steps After Phase 1

**Phase 1B**: Add literature and rules
- Change `source_group="experiments"` to None (all sources)
- Aggregate from multiple sources

**Phase 2**: Structured analysis (electrophile/nucleophile profiling)
- Optional enhancement
- Current auto-detection already works well

**Phase 3**: Intelligent filtering
- Filter by mechanism (SN2 → polar aprotic)
- Boost compatible conditions

This revised plan is **much simpler** and **more robust** because it leverages the existing taxonomy infrastructure rather than trying to recreate it.
