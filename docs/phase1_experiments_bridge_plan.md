# Phase 1: Bridge to Experiments-Based Condition Recommendation

## Simplified Plan - MVP Implementation

**Goal**: Connect reaction analysis to HTE experiments database only (simplest case)

**Timeline**: 2-3 days

---

## Current Situation Analysis

### What We Have

#### 1. **Reaction Analysis Output Structure**
```python
{
    "schema_version": "reaction_analysis.v1",
    "input": {
        "rxn_smiles_raw": "...",
        "rxn_smiles_clean": "reactants>>products",
        "warnings": [...]
    },
    "auto_interpretation": {  # Tier 1
        "interpretation": {...},
        "report": "..."
    },
    "quick_glance": {  # Tier 2 (DeepSeek-v3.2)
        "reaction_types": ["Suzuki-Miyaura", "deprotection", ...],
        "all_transformations": [...],  # If thorough=True
        "complexity": "simple|moderate|complex|tandem",
        "summary": "...",
        "confidence": 0.85,
        "success": True
    },
    "interpretation": {  # Tier 3
        "overall_class": "cross_coupling",
        "events": [...],
        "confidence": 0.80
    },
    "validation": {  # Tier 4 (optional)
        "gate": {"status": "pass|warning|fail"}
    },
    "metadata": {...}
}
```

#### 2. **HTERecommender API** (in `chemtools/recommend/recommender.py`)
```python
def recommend(
    self,
    reactant_a_smiles: str,          # REQUIRED
    reactant_b_smiles: Optional[str] = None,
    product_smiles: Optional[str] = None,
    top_k: int = 10,
    min_experiments: int = 2,
    reaction_type_filter: Optional[str] = None,  # e.g., "Suzuki_Miyaura_coupling"
    source_group: Optional[str] = None,  # Filter: "experiments" | "literature" | "rules"
    ...
) -> HTERecommendationResult
```

#### 3. **HTE Experiments Database** (`data/HTE_db/experiments/HTE_canonical.csv`)
```csv
reaction_type,yield,z_score,reactant_1,reactant_2,catalyst,ligand,base,solvent,additive
C_H_arylation,48.01,0.42,Ar-H,Ar-Br,XantPhos Pd(allyl)Cl,XantPhos,K2CO3,"DMAc, water",KOPiv
Suzuki_Miyaura_coupling,82.5,1.56,Ar-Br,Ar-B(OH)2,Pd(OAc)2,SPhos,K3PO4,dioxane,
...
```

### What's Missing (The Gaps)

#### Gap 1: Reactant SMILES Extraction
**Problem**: Analysis provides `rxn_smiles_clean: "A.B.C>>D"` but HTERecommender needs separate `reactant_a_smiles` and `reactant_b_smiles`

**Solution**: Parse and split reactants:
```python
# From: "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
# Extract:
#   reactant_a_smiles = "Brc1ccccc1"
#   reactant_b_smiles = "B(O)(O)c1ccccc1"
#   product_smiles = "c1ccc(-c2ccccc2)cc1"
```

**Challenge**: How to decide which is A vs B?
- For 2 reactants: arbitrary order is OK (recommender handles both)
- For 3+ reactants: may need heuristic (largest first, or use Tier 2 to identify roles)

#### Gap 2: Reaction Type Name Mapping
**Problem**: Tier 2 reports "Suzuki-Miyaura" but HTE database uses "Suzuki_Miyaura_coupling"

**Solution**: Create mapping table:
```python
REACTION_TYPE_MAP = {
    # Tier 2 name → HTE database name
    "Suzuki-Miyaura": "Suzuki_Miyaura_coupling",
    "Suzuki": "Suzuki_Miyaura_coupling",
    "Buchwald-Hartwig": "Buchwald_Hartwig_amination",
    "Buchwald-Hartwig C-N coupling": "Buchwald_Hartwig_amination",
    "Amide formation": "Amide_formation",
    "C-H arylation": "C_H_arylation",
    "Chan-Lam": "ChanLam_CN",
    # ... expand as needed
}
```

**Challenge**: Tier 2 may return multiple types for tandem reactions
- Strategy: Try first type, if no results try others

#### Gap 3: Multiple Reaction Types (Tandem/Complex)
**Problem**: Complex reactions like "Suzuki + THP deprotection" have multiple types

**Solution**: Prioritize types:
1. Try most specific type first (e.g., "Suzuki_Miyaura_coupling")
2. If no experiments found, try secondary types
3. If still no results, fall back to reactant-only matching (no type filter)

#### Gap 4: Display Format
**Problem**: CLI needs to show condition recommendations clearly

**Solution**: Add new display section:
```
================================================================================
  CONDITION RECOMMENDATIONS (from 127 experiments)
================================================================================

Rank 1: [Score: 1.56 | Confidence: 95% | 23 experiments]
  Catalyst: Pd(OAc)2
  Ligand: SPhos
  Base: K3PO4
  Solvent: dioxane
  Performance: 82% success rate, 78% avg yield

Rank 2: [Score: 1.42 | Confidence: 92% | 18 experiments]
  ...
```

---

## Implementation Plan

### Step 1: Create Bridge Module (Day 1 Morning)

**File**: `reaction_agent/condition_bridge.py`

**Components**:

```python
"""
Bridge module connecting reaction analysis to HTE condition recommendation.

Workflow:
    Reaction SMILES → analyze_reaction_smiles() → AnalysisResult
                                                       ↓
    AnalysisResult → extract_smiles() → (reactants, product)
                                                       ↓
    (reactants, product) → map_reaction_type() → reaction_type_filter
                                                       ↓
    HTERecommender.recommend() → Condition recommendations
"""

from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
from chemtools.recommend.recommender import HTERecommender, HTERecommendationResult

# Reaction type mapping: Tier 2 names → HTE database names
REACTION_TYPE_MAP = {
    "Suzuki-Miyaura": "Suzuki_Miyaura_coupling",
    "Suzuki-Miyaura cross-coupling": "Suzuki_Miyaura_coupling",
    "Suzuki": "Suzuki_Miyaura_coupling",
    "Buchwald-Hartwig": "Buchwald_Hartwig_amination",
    "Buchwald-Hartwig C-N coupling": "Buchwald_Hartwig_amination",
    "Buchwald-Hartwig amination": "Buchwald_Hartwig_amination",
    "Amide formation": "Amide_formation",
    "Amidation": "Amide_formation",
    "C-H arylation": "C_H_arylation",
    "Chan-Lam coupling": "ChanLam_CN",
    "Chan-Lam": "ChanLam_CN",
    "Sonogashira": "Sonogashira",
    "Sonogashira coupling": "Sonogashira",
    "Heck": "Heck",
    "Heck coupling": "Heck",
    "Reductive amination": "Reductive_amination",
    # Add more as discovered
}


class ConditionBridge:
    """
    Bridge between reaction analysis and HTE experiments database.

    Simplified version - only connects to experiments (not literature/rules).
    """

    def __init__(self, hte_db_path: str = "data/HTE_db"):
        """
        Initialize bridge with HTE database.

        Args:
            hte_db_path: Path to HTE database directory
        """
        self.hte_db_path = Path(hte_db_path)
        self.recommender = HTERecommender(str(self.hte_db_path))

    def extract_smiles(
        self,
        analysis_result: Dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Extract reactant and product SMILES from analysis result.

        Args:
            analysis_result: Output from analyze_reaction_smiles()

        Returns:
            Tuple of (reactant_a_smiles, reactant_b_smiles, product_smiles)
        """
        input_data = analysis_result.get('input', {})
        rxn_smiles_clean = input_data.get('rxn_smiles_clean', '')

        if not rxn_smiles_clean or '>>' not in rxn_smiles_clean:
            return None, None, None

        # Parse reaction SMILES
        parts = rxn_smiles_clean.split('>>')
        reactants_str = parts[0] if len(parts) > 0 else ''
        products_str = parts[1] if len(parts) > 1 else ''

        # Split reactants by '.'
        reactants = [r.strip() for r in reactants_str.split('.') if r.strip()]

        # Extract first two reactants (assume A and B)
        # For simplicity, order doesn't matter - HTERecommender handles this
        reactant_a = reactants[0] if len(reactants) > 0 else None
        reactant_b = reactants[1] if len(reactants) > 1 else None

        # Extract first product (main product)
        products = [p.strip() for p in products_str.split('.') if p.strip()]
        product = products[0] if len(products) > 0 else None

        return reactant_a, reactant_b, product

    def map_reaction_type(
        self,
        analysis_result: Dict[str, Any]
    ) -> Optional[str]:
        """
        Map reaction type from analysis to HTE database name.

        Priority:
        1. Tier 2 quick_glance reaction_types
        2. Tier 3 interpretation overall_class
        3. None (fall back to reactant-only matching)

        Args:
            analysis_result: Output from analyze_reaction_smiles()

        Returns:
            HTE database reaction type name, or None if not found
        """
        # Try Tier 2 first (most reliable)
        quick_glance = analysis_result.get('quick_glance', {})
        if quick_glance.get('success'):
            reaction_types = quick_glance.get('reaction_types', [])
            for rxn_type in reaction_types:
                mapped = REACTION_TYPE_MAP.get(rxn_type)
                if mapped:
                    return mapped

        # Try Tier 3 as fallback
        interpretation = analysis_result.get('interpretation', {})
        overall_class = interpretation.get('overall_class')
        if overall_class:
            # Try to map overall_class (e.g., "cross_coupling" → might match several types)
            # For now, just check if it's in the map
            mapped = REACTION_TYPE_MAP.get(overall_class)
            if mapped:
                return mapped

        # No mapping found - will use reactant-only matching
        return None

    def recommend_conditions(
        self,
        analysis_result: Dict[str, Any],
        top_k: int = 10,
        min_experiments: int = 2
    ) -> HTERecommendationResult:
        """
        Get condition recommendations based on analysis result.

        Args:
            analysis_result: Output from analyze_reaction_smiles()
            top_k: Number of recommendations to return
            min_experiments: Minimum experiments for a condition

        Returns:
            HTERecommendationResult with ranked recommendations
        """
        # Step 1: Extract SMILES
        reactant_a, reactant_b, product = self.extract_smiles(analysis_result)

        if not reactant_a:
            raise ValueError("Could not extract reactant SMILES from analysis")

        # Step 2: Map reaction type
        reaction_type = self.map_reaction_type(analysis_result)

        # Step 3: Call HTERecommender (experiments only)
        return self.recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,
            top_k=top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type,
            source_group="experiments"  # Only experiments for Phase 1
        )

    def analyze_and_recommend(
        self,
        rxn_smiles: str,
        analyzer: 'ReactionSMILESAnalyzer',
        top_k: int = 10,
        validate: bool = True
    ) -> Dict[str, Any]:
        """
        End-to-end: Analyze reaction → Recommend conditions.

        Args:
            rxn_smiles: Reaction SMILES
            analyzer: ReactionSMILESAnalyzer instance
            top_k: Number of recommendations
            validate: Enable Tier 4 validation

        Returns:
            Dict with:
            - analysis: Full reaction analysis (Tiers 1-4)
            - recommendations: HTERecommendationResult
            - metadata: Timing info
        """
        import time

        start_time = time.time()

        # Step 1: Analyze reaction
        analysis_result = analyzer.analyze(rxn_smiles, validate=validate)
        analysis_time = time.time() - start_time

        # Step 2: Get recommendations
        rec_start = time.time()
        recommendations = self.recommend_conditions(
            analysis_result,
            top_k=top_k
        )
        rec_time = time.time() - rec_start

        return {
            "analysis": analysis_result,
            "recommendations": recommendations,
            "metadata": {
                "analysis_time_s": analysis_time,
                "recommendation_time_s": rec_time,
                "total_time_s": time.time() - start_time
            }
        }


# Export public API
__all__ = ['ConditionBridge', 'REACTION_TYPE_MAP']
```

**Testing approach**:
```python
# Quick test (add to module docstring)
if __name__ == "__main__":
    from llmtools.clients import LLMClient
    from reaction_agent import ReactionSMILESAnalyzer

    client = LLMClient(provider="openai", model="gpt-4o-mini")
    analyzer = ReactionSMILESAnalyzer(client)
    bridge = ConditionBridge()

    # Test Suzuki coupling
    rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
    result = bridge.analyze_and_recommend(rxn, analyzer, top_k=5)

    print(f"Analysis: {result['analysis']['quick_glance']['summary']}")
    print(f"Recommendations: {len(result['recommendations'].recommendations)}")
    for i, rec in enumerate(result['recommendations'].recommendations[:3], 1):
        print(f"{i}. {rec.catalyst} + {rec.ligand} (score: {rec.avg_z_score:.2f})")
```

---

### Step 2: CLI Integration (Day 1 Afternoon)

**Modify**: `reaction_agent/cli.py`

**Changes**:

#### 2.1 Add CLI Flag (around line 655)
```python
parser.add_argument(
    '--recommend',
    action='store_true',
    help='Recommend reaction conditions from HTE experiments database'
)
parser.add_argument(
    '--top-conditions',
    type=int,
    default=10,
    help='Number of condition recommendations to show (default: 10)'
)
```

#### 2.2 Add Display Function (around line 100, after existing display functions)
```python
def display_recommendations(
    recommendations: HTERecommendationResult,
    top_k: int = 10
):
    """
    Display condition recommendations from HTE database.

    Args:
        recommendations: HTERecommendationResult object
        top_k: Number of recommendations to display
    """
    if not recommendations.recommendations:
        print(f"{Colors.YELLOW}No condition recommendations found{Colors.END}")
        print(f"Try using different reactants or removing reaction_type filter")
        return

    # Header
    total_expts = recommendations.total_matching_experiments
    print(f"\nFound {len(recommendations.recommendations)} condition sets from {total_expts} experiments")
    print(f"Showing top {min(top_k, len(recommendations.recommendations))}:\n")

    # Display each recommendation
    for i, rec in enumerate(recommendations.recommendations[:top_k], 1):
        # Color code by score
        if rec.avg_z_score > 1.0:
            score_color = Colors.GREEN
        elif rec.avg_z_score > 0.0:
            score_color = Colors.YELLOW
        else:
            score_color = Colors.RED

        # Header
        print(f"{Colors.BOLD}Rank {i}:{Colors.END}")
        print(f"  {score_color}Score: {rec.avg_z_score:.2f}{Colors.END} | "
              f"Confidence: {rec.confidence_score:.0f}% | "
              f"Experiments: {rec.num_experiments}")

        # Conditions
        print(f"\n  {Colors.BOLD}Conditions:{Colors.END}")
        if rec.catalyst:
            print(f"    Catalyst: {rec.catalyst}")
        if rec.ligand:
            print(f"    Ligand: {rec.ligand}")
        if rec.base:
            print(f"    Base: {rec.base}")
        if rec.solvent:
            print(f"    Solvent: {rec.solvent}")
        if rec.additive:
            print(f"    Additive: {rec.additive}")

        # Performance
        print(f"\n  {Colors.BOLD}Performance:{Colors.END}")
        print(f"    Success rate: {rec.success_rate:.1f}%")
        print(f"    Avg yield: {rec.avg_yield:.1f}%")
        print(f"    Median yield: {rec.median_yield:.1f}%")

        print()  # Blank line between recommendations
```

#### 2.3 Integrate into Main Flow (in main() function, around line 800)
```python
# After analysis completes in analyze_reaction_interactive()
def analyze_reaction_interactive(
    analyzer: ReactionSMILESAnalyzer,
    rxn_smiles: str,
    save_output: Optional[Path] = None,
    mode: str = "auto",
    validate: bool = False,
    retry_config: Optional['RetryConfig'] = None,
    recommend: bool = False,  # NEW
    top_conditions: int = 10  # NEW
) -> Dict[str, Any]:
    """Analyze reaction with optional condition recommendation."""

    # ... existing analysis code ...

    # NEW: Add recommendation section
    if recommend:
        from reaction_agent.condition_bridge import ConditionBridge

        print_header("CONDITION RECOMMENDATIONS")

        try:
            bridge = ConditionBridge()
            recommendations = bridge.recommend_conditions(
                result,
                top_k=top_conditions
            )
            display_recommendations(recommendations, top_conditions)

            # Optionally add to result
            result['condition_recommendations'] = {
                'recommendations': recommendations,
                'source': 'experiments',
                'top_k': top_conditions
            }

        except Exception as e:
            print(f"{Colors.RED}✗ Recommendation failed: {e}{Colors.END}")
            import traceback
            traceback.print_exc()

    return result
```

#### 2.4 Wire Up CLI Args (in main() around line 800)
```python
elif args.reaction:
    # Single reaction mode
    # ... existing retry config code ...

    analyze_reaction_interactive(
        analyzer,
        args.reaction,
        save_output=args.output,
        mode=effective_mode,
        validate=args.validate,
        retry_config=retry_config,
        recommend=args.recommend,  # NEW
        top_conditions=args.top_conditions  # NEW
    )
```

---

### Step 3: Testing & Refinement (Day 2)

#### 3.1 Manual Testing
```bash
# Test 1: Simple Suzuki coupling (should find matches)
python -m reaction_agent.cli \
  --reaction "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" \
  --validate \
  --recommend \
  --top-conditions 5

# Expected: Should see Suzuki conditions (Pd catalysts, SPhos/XPhos ligands, etc.)

# Test 2: Complex tandem reaction
python -m reaction_agent.cli \
  --reaction "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C.Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1" \
  --validate \
  --recommend

# Expected: Should find Suzuki conditions (deprotection not in experiments DB)

# Test 3: Unknown reaction type (fallback to reactants)
python -m reaction_agent.cli \
  --reaction "CC(C)C[C@H](NC(=O)OC(C)(C)C)C(=O)O.NCc1ccccc1>>CC(C)C[C@H](NC(=O)OC(C)(C)C)C(=O)NCc1ccccc1" \
  --validate \
  --recommend

# Expected: Should find amide formation conditions
```

#### 3.2 Expand Reaction Type Mapping
Based on testing, add more entries to `REACTION_TYPE_MAP`:
- Check `data/HTE_db/experiments/HTE_canonical.csv` for all reaction_types
- Test common reaction names from Tier 2
- Add aliases (e.g., "Pd-catalyzed cross-coupling" → "Suzuki_Miyaura_coupling")

#### 3.3 Handle Edge Cases
- **No recommendations found**: Show helpful message with suggestions
- **Parsing errors**: Validate SMILES extraction, show warning
- **Multiple products**: Currently uses first product, may need to handle better

---

### Step 4: Documentation (Day 2 Afternoon)

#### 4.1 Update README
Add section:
```markdown
## Condition Recommendation

Get experimental condition recommendations based on reaction analysis:

```bash
python -m reaction_agent.cli \
  --reaction "YOUR_REACTION_SMILES" \
  --validate \
  --recommend \
  --top-conditions 10
```

The system will:
1. Analyze the reaction (Tiers 1-4)
2. Match against HTE experiments database
3. Recommend top conditions by Z-score
4. Show catalyst, ligand, base, solvent, performance stats

**Supported reaction types**: Suzuki-Miyaura, Buchwald-Hartwig, Amide formation, C-H arylation, etc.
```

#### 4.2 Create Usage Example
Add to `docs/examples/`:
```python
"""
Example: End-to-end reaction analysis + condition recommendation
"""

from llmtools.clients import LLMClient
from reaction_agent import ReactionSMILESAnalyzer
from reaction_agent.condition_bridge import ConditionBridge

# Initialize
client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)
bridge = ConditionBridge()

# Analyze + Recommend
rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
result = bridge.analyze_and_recommend(rxn, analyzer, validate=True)

# Display
print(f"Reaction type: {result['analysis']['quick_glance']['reaction_types']}")
print(f"Found {len(result['recommendations'].recommendations)} conditions")

for i, rec in enumerate(result['recommendations'].recommendations[:3], 1):
    print(f"\n{i}. Score: {rec.avg_z_score:.2f}")
    print(f"   Catalyst: {rec.catalyst}, Ligand: {rec.ligand}")
    print(f"   Success rate: {rec.success_rate:.1f}%, Avg yield: {rec.avg_yield:.1f}%")
```

---

## Success Criteria

✅ **Phase 1 Complete When:**

1. ✅ Can run: `python -m reaction_agent.cli --reaction "SMILES" --recommend`
2. ✅ Shows top 10 conditions from experiments database
3. ✅ Displays catalyst, ligand, base, solvent, performance metrics
4. ✅ Works for common reactions (Suzuki, Buchwald-Hartwig, amide)
5. ✅ Handles edge cases gracefully (no matches, parsing errors)
6. ✅ Total time < 30 seconds (analysis + recommendation)
7. ✅ Documentation updated with usage examples

---

## Deferred to Later Phases

**Not included in Phase 1 (Experiments-only MVP):**
- ❌ Literature conditions (Phase 1B)
- ❌ Rule-based conditions (Phase 1B)
- ❌ Multi-source aggregation (Phase 1B)
- ❌ Structured reaction understanding (Phase 2)
- ❌ Intelligent filtering by mechanism (Phase 3)
- ❌ User feedback tracking (Phase 4)
- ❌ Web interface (Phase 5)

**These will be added incrementally after Phase 1 is proven working.**

---

## Analysis Module Improvements Needed

### Issue 1: Reactant Role Identification

**Current**: Analysis splits reactants by '.' but doesn't identify which is electrophile vs nucleophile

**Fix Needed**: Add to Tier 2 prompt to identify reactant roles
```python
# In quick_reaction_glance comprehensive mode:
# Request: "Identify reactant roles: electrophile, nucleophile, catalyst, base, etc."
# Add to output:
{
    "reactant_roles": {
        "reactant_1": "aryl_halide",  # First SMILES component
        "reactant_2": "boronic_acid",  # Second SMILES component
        "reactant_3": "base"  # Third (spectator)
    }
}
```

**Priority**: LOW (not critical for Phase 1 - order doesn't matter much)

### Issue 2: Reaction Type Standardization

**Current**: Tier 2 returns natural language names (e.g., "Suzuki-Miyaura cross-coupling")

**Fix Needed**: Add standardized names to output
```python
# In quick_reaction_glance output:
{
    "reaction_types": ["Suzuki-Miyaura cross-coupling"],
    "standardized_types": ["Suzuki_Miyaura_coupling"],  # HTE database format
}
```

**Priority**: MEDIUM (nice to have, but mapping table works for now)

### Issue 3: Confidence in Reaction Type

**Current**: Single confidence for entire analysis

**Fix Needed**: Per-reaction-type confidence
```python
{
    "reaction_types": [
        {"name": "Suzuki-Miyaura", "confidence": 0.95},
        {"name": "deprotection", "confidence": 0.85}
    ]
}
```

**Priority**: LOW (can use overall confidence for now)

---

## Risk Mitigation

### Risk 1: No Matches Found
**Likelihood**: Medium (for unusual reactions)
**Impact**: High (user gets no recommendations)
**Mitigation**:
- Fall back to reactant-only matching (no reaction_type filter)
- Show message: "No exact matches found, showing similar reactions"
- Log unmapped reaction types for future expansion

### Risk 2: Incorrect Reaction Type Detection
**Likelihood**: Low (Tier 2 is accurate)
**Impact**: Medium (wrong conditions recommended)
**Mitigation**:
- Always include Tier 2 summary in output
- User can verify reaction type is correct
- Show match quality in recommendations

### Risk 3: Multiple Reaction Types (Tandem)
**Likelihood**: Medium (complex reactions common)
**Impact**: Medium (unclear which conditions to show)
**Mitigation**:
- Try each reaction type in order
- Aggregate results from multiple types
- Label which type each condition came from

---

## Next Steps After Phase 1

**Phase 1B: Add Literature & Rules**
- Extend `source_group` to include "literature" and "rules"
- Aggregate recommendations from multiple sources
- Weight by source reliability

**Phase 2: Structured Analysis**
- Implement electrophile/nucleophile profiling
- Add selectivity risk assessment
- Use structure to filter conditions

**Phase 3: Intelligent Filtering**
- Filter incompatible conditions (e.g., no protic solvent for SN2)
- Boost compatible conditions (e.g., polar aprotic for SN2)
- Add mechanism-specific warnings

---

## Summary

**Phase 1 Deliverable**: Working end-to-end CLI for reaction analysis + experiments-based condition recommendation

**Timeline**: 2-3 days
- Day 1 AM: Bridge module implementation
- Day 1 PM: CLI integration
- Day 2 AM: Testing and refinement
- Day 2 PM: Documentation

**Key Files**:
- NEW: `reaction_agent/condition_bridge.py` (~200 lines)
- MODIFY: `reaction_agent/cli.py` (~100 lines added)

**Testing Command**:
```bash
python -m reaction_agent.cli \
  --reaction "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1" \
  --validate \
  --recommend \
  --top-conditions 5
```

This simplified Phase 1 focuses on the minimal viable connection: analysis → experiments → display. Later phases will add sophistication (structured analysis, intelligent filtering, multi-source aggregation).
