# Condition Recommendation Implementation Plan

## Executive Summary

**Goal**: Connect the newly-built reaction analysis system (Tiers 1-4 + Retry) to the existing HTE condition recommendation system, implementing the structured reaction understanding framework from to_do.md.

**Current State:**

- ✅ **Reaction Analysis** (Phase 1): Three-tier analysis + validation + retry system complete
- ✅ **HTE Recommendation System** (Production): Sophisticated recommender with 100+ datasets, rule-based templates, motif matching
- ⚠ **Gap**: These two systems are disconnected - no bridge between analysis → recommendation

**Vision** (from to_do.md):
Decompose reactions into orthogonal aspects:

1. Net transformation (bond changes)
2. Electrophile analysis (carbon + leaving group)
3. Nucleophile/basicity profile
4. Driving forces & thermodynamics
5. Mechanistic class
6. Selectivity map
7. **Conditions as constraints** ← Critical for recommendation
8. Side reactions & failure modes

---

## Phase 1: Bridge Reaction Analysis → Condition Recommendation

### Objective

Connect reaction analysis output to HTERecommender input, creating end-to-end workflow.

### 1.1 New Module: `reaction_agent/condition_bridge.py`

**Purpose**: Map reaction analysis results to HTERecommender API

```python
"""
Bridge module connecting reaction analysis to condition recommendation.

Workflow:
    Reaction SMILES → analyze_reaction_smiles() → AnalysisResult
                                                        ↓
    AnalysisResult → extract_recommendation_inputs() → RecommendationQuery
                                                        ↓
    RecommendationQuery → HTERecommender.recommend() → Condition recommendations
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass
from chemtools.recommend.recommender import HTERecommender, HTERecommendationResult


@dataclass
class RecommendationQuery:
    """
    Structured input for condition recommendation.

    Maps reaction analysis output to HTERecommender input format.
    """
    # Core reaction SMILES
    reactants_smiles: str
    products_smiles: str
    reaction_smiles: str

    # Reaction classification (from Tier 2/3)
    reaction_types: List[str]  # From quick_glance or interpretation
    primary_reaction_type: str  # Most confident

    # Reactant types (from Tier 1)
    reactant_a_smiles: Optional[str] = None
    reactant_b_smiles: Optional[str] = None
    reactant_c_smiles: Optional[str] = None

    # Mechanistic insights (from Tier 3)
    overall_class: Optional[str] = None  # cross_coupling, substitution, etc.
    mechanistic_events: Optional[List[Dict]] = None

    # Context from validation (Tier 4)
    validation_warnings: Optional[List[str]] = None
    atom_balance: Optional[Dict] = None

    # Search strategy
    top_k: int = 10
    include_rules: bool = True
    include_literature: bool = True
    include_experiments: bool = True


class RecommendationBridge:
    """
    Bridge between reaction analysis and condition recommendation.

    Key responsibilities:
    1. Extract recommendation inputs from analysis results
    2. Map terminology (Tier 2/3 names → HTE database names)
    3. Handle ambiguous/multiple reaction types
    4. Aggregate recommendations from multiple strategies
    """

    def __init__(self, hte_db_path: str = "data/HTE_db"):
        """
        Initialize recommendation bridge.

        Args:
            hte_db_path: Path to HTE database directory
        """
        self.recommender = HTERecommender(hte_db_path)

        # Reaction type mapping: Tier 2/3 names → HTE database names
        self.reaction_type_map = {
            "Suzuki-Miyaura cross-coupling": "Suzuki_Miyaura_coupling",
            "Suzuki-Miyaura": "Suzuki_Miyaura_coupling",
            "Suzuki": "Suzuki_Miyaura_coupling",
            "Buchwald-Hartwig C-N coupling": "Buchwald_Hartwig_amination",
            "Buchwald-Hartwig": "Buchwald_Hartwig_amination",
            "Amide formation": "Amide_formation",
            "Amidation": "Amide_formation",
            "Sonogashira coupling": "Sonogashira",
            "Heck coupling": "Heck",
            "Chan-Lam coupling": "ChanLam_CN",
            "Reductive amination": "Reductive_amination",
            # ... expand this mapping
        }

    def extract_recommendation_inputs(
        self,
        analysis_result: Dict[str, Any]
    ) -> RecommendationQuery:
        """
        Extract recommendation query from analysis result.

        Args:
            analysis_result: Output from analyze_reaction_smiles()

        Returns:
            RecommendationQuery ready for HTERecommender
        """
        # Extract input SMILES
        input_data = analysis_result.get('input', {})
        rxn_smiles_clean = input_data.get('rxn_smiles_clean', '')

        # Parse reactants and products
        parts = rxn_smiles_clean.split('>>')
        reactants_smiles = parts[0] if len(parts) > 0 else ''
        products_smiles = parts[1] if len(parts) > 1 else ''

        # Split reactants (handle multiple reactants)
        reactant_smiles_list = [s.strip() for s in reactants_smiles.split('.') if s.strip()]
        reactant_a = reactant_smiles_list[0] if len(reactant_smiles_list) > 0 else None
        reactant_b = reactant_smiles_list[1] if len(reactant_smiles_list) > 1 else None
        reactant_c = reactant_smiles_list[2] if len(reactant_smiles_list) > 2 else None

        # Extract reaction types from Tier 2 (most reliable)
        quick_glance = analysis_result.get('quick_glance', {})
        reaction_types_raw = quick_glance.get('reaction_types', [])

        # Map to HTE database names
        reaction_types = [
            self.reaction_type_map.get(rt, rt)
            for rt in reaction_types_raw
        ]

        # Primary reaction type (first one, highest confidence)
        primary_reaction_type = reaction_types[0] if reaction_types else None

        # Mechanistic insights from Tier 3
        interpretation = analysis_result.get('interpretation', {})
        overall_class = interpretation.get('overall_class')
        mechanistic_events = interpretation.get('events', [])

        # Validation context from Tier 4
        validation = analysis_result.get('validation', {})
        rdkit = validation.get('rdkit', {})
        validation_warnings = rdkit.get('warnings', [])
        atom_balance = rdkit.get('atom_balance', {})

        return RecommendationQuery(
            reactants_smiles=reactants_smiles,
            products_smiles=products_smiles,
            reaction_smiles=rxn_smiles_clean,
            reaction_types=reaction_types,
            primary_reaction_type=primary_reaction_type,
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            reactant_c_smiles=reactant_c,
            overall_class=overall_class,
            mechanistic_events=mechanistic_events,
            validation_warnings=validation_warnings,
            atom_balance=atom_balance
        )

    def recommend_conditions(
        self,
        query: RecommendationQuery,
        **kwargs
    ) -> HTERecommendationResult:
        """
        Get condition recommendations for reaction query.

        Args:
            query: RecommendationQuery from reaction analysis
            **kwargs: Additional HTERecommender.recommend() parameters

        Returns:
            HTERecommendationResult with ranked condition recommendations
        """
        # Build source filter
        sources = []
        if query.include_rules:
            sources.append("rules")
        if query.include_literature:
            sources.append("literature")
        if query.include_experiments:
            sources.append("experiments")

        # Call HTERecommender with reaction type filter
        return self.recommender.recommend(
            reactant_a_smiles=query.reactant_a_smiles,
            reactant_b_smiles=query.reactant_b_smiles or query.reactant_c_smiles,  # Fallback to C if no B
            product_smiles=query.products_smiles,
            top_k=query.top_k,
            reaction_type_filter=query.primary_reaction_type,
            source_groups=sources,
            **kwargs
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
            top_k: Number of condition recommendations
            validate: Enable Tier 4 validation

        Returns:
            Dict with:
            - analysis: Reaction analysis result (Tiers 1-4)
            - query: Extracted recommendation query
            - recommendations: HTERecommendationResult
            - metadata: Timing, costs, etc.
        """
        import time

        start_time = time.time()

        # Step 1: Analyze reaction
        analysis_result = analyzer.analyze(rxn_smiles, validate=validate)
        analysis_time = time.time() - start_time

        # Step 2: Extract recommendation inputs
        query = self.extract_recommendation_inputs(analysis_result)

        # Step 3: Get condition recommendations
        recommendation_start = time.time()
        recommendations = self.recommend_conditions(query, top_k=top_k)
        recommendation_time = time.time() - recommendation_start

        return {
            "analysis": analysis_result,
            "query": query,
            "recommendations": recommendations,
            "metadata": {
                "analysis_time_s": analysis_time,
                "recommendation_time_s": recommendation_time,
                "total_time_s": time.time() - start_time,
                "analysis_cost": analysis_result.get('metadata', {}).get('total_tokens', 0) * 0.000001,  # Rough estimate
                "top_k": top_k
            }
        }


# Export public API
__all__ = [
    'RecommendationQuery',
    'RecommendationBridge'
]
```

### 1.2 CLI Integration

Add new command to `reaction_agent/cli.py`:

```python
# New command: --recommend
parser.add_argument(
    '--recommend',
    action='store_true',
    help='Recommend reaction conditions after analysis'
)
parser.add_argument(
    '--top-k',
    type=int,
    default=10,
    help='Number of condition recommendations (default: 10)'
)
parser.add_argument(
    '--include-rules',
    action='store_true',
    default=True,
    help='Include rule-based conditions'
)
parser.add_argument(
    '--include-literature',
    action='store_true',
    default=True,
    help='Include literature conditions'
)
parser.add_argument(
    '--include-experiments',
    action='store_true',
    default=True,
    help='Include HTE experimental conditions'
)

# In main() after analysis:
if args.recommend:
    from reaction_agent.condition_bridge import RecommendationBridge

    bridge = RecommendationBridge()
    full_result = bridge.analyze_and_recommend(
        rxn_smiles=args.reaction,
        analyzer=analyzer,
        top_k=args.top_k,
        validate=args.validate
    )

    # Display recommendations
    print_header("CONDITION RECOMMENDATIONS")
    display_recommendations(full_result['recommendations'], args.top_k)
```

### 1.3 Display Function for Recommendations

```python
def display_recommendations(
    result: HTERecommendationResult,
    top_k: int = 10
):
    """Display condition recommendations in CLI."""

    if not result.recommendations:
        print(f"{Colors.YELLOW}No recommendations found{Colors.END}")
        return

    print(f"\nFound {len(result.recommendations)} condition sets")
    print(f"Showing top {min(top_k, len(result.recommendations))}:\n")

    for i, rec in enumerate(result.recommendations[:top_k], 1):
        # Ranking info
        quality_color = Colors.GREEN if rec.avg_z_score > 1.0 else Colors.YELLOW if rec.avg_z_score > 0.0 else Colors.RED
        print(f"\n{Colors.BOLD}Rank {i}:{Colors.END}")
        print(f"  {quality_color}Score: {rec.avg_z_score:.2f}{Colors.END} | "
              f"Confidence: {rec.confidence_score:.0f}% | "
              f"Experiments: {rec.num_experiments}")

        # Core conditions
        print(f"\n  {Colors.BOLD}Conditions:{Colors.END}")
        if rec.catalyst:
            print(f"    Catalyst: {rec.catalyst}")
        if rec.ligand:
            print(f"    Ligand: {rec.ligand}")
        if rec.base:
            print(f"    Base: {rec.base}")
        if rec.solvent:
            print(f"    Solvent: {rec.solvent}")
        if rec.secondary_solvent:
            print(f"    Co-solvent: {rec.secondary_solvent}")
        if rec.additive:
            print(f"    Additive: {rec.additive}")
        if rec.coupling_reagent:
            print(f"    Coupling Reagent: {rec.coupling_reagent}")

        # Performance statistics
        print(f"\n  {Colors.BOLD}Performance:{Colors.END}")
        print(f"    Success rate: {rec.success_rate:.0f}%")
        print(f"    Avg yield: {rec.avg_yield:.1f}%")
        print(f"    Median yield: {rec.median_yield:.1f}%")

        # Match quality
        if hasattr(rec, 'match_score'):
            match_color = Colors.GREEN if rec.match_score > 0.8 else Colors.YELLOW
            print(f"    {match_color}Match quality: {rec.match_score:.2f}{Colors.END}")
```

---

## Phase 2: Implement Structured Reaction Understanding (to_do.md)

### Objective

Build deterministic analysis modules for electrophile, nucleophile, mechanism, and selectivity analysis.

### 2.1 New Module: `reaction_agent/structure_analysis.py`

**Purpose**: Implement the to_do.md framework using RDKit + existing chemtools

```python
"""
Structured reaction understanding framework.

Implements to_do.md decomposition:
1. Net transformation (bond changes)
2. Electrophile analysis
3. Nucleophile/basicity profile
4. Driving forces & thermodynamics
5. Mechanistic class
6. Selectivity map
7. Conditions as constraints
8. Side reactions & failure modes
"""

from typing import Dict, Any, List, Optional, Set, Tuple
from dataclasses import dataclass
from chemtools.util.rdkit_helpers import parse_smiles
from chemtools.featurizers.molecule import detect_reaction_motifs
from rdkit import Chem


@dataclass
class ElectrophileProfile:
    """
    Electrophile center analysis (to_do.md #2).

    Characterizes carbon + leaving group:
    - Hybridization (sp, sp2, sp3)
    - Degree (methyl, 1°, 2°, 3°)
    - Stabilization (benzylic, allylic, resonance)
    - Electronic activation (EWG presence)
    - Leaving group quality
    - Steric hindrance
    """
    center_atom_idx: int
    element: str  # Usually 'C'
    hybridization: str  # 'sp', 'sp2', 'sp3'
    degree: str  # 'methyl', '1°', '2°', '3°'

    # Stabilization
    is_benzylic: bool
    is_allylic: bool
    is_propargylic: bool
    has_resonance: bool

    # Electronic effects
    ewg_neighbors: List[str]  # EWGs within 2 bonds
    is_activated: bool

    # Leaving group
    leaving_group: Optional[str]  # 'Br', 'I', 'OTf', 'OTs', etc.
    lg_quality_score: float  # 0-1, OTf=1.0, F=0.1

    # Sterics
    neopentyl_penalty: bool
    backside_access: str  # 'good', 'hindered', 'blocked'

    # Mechanism hints
    sn1_favorable: bool
    sn2_favorable: bool
    e1_favorable: bool
    e2_favorable: bool
    snar_favorable: bool
    oxidative_addition_favorable: bool


@dataclass
class NucleophileProfile:
    """
    Nucleophile/base analysis (to_do.md #3).

    Characterizes attacking species:
    - Hard/soft character
    - Charge and aggregation
    - Steric bulk
    - Ambident nucleophiles
    - Basicity vs nucleophilicity
    """
    atom_idx: int
    element: str  # 'N', 'O', 'S', 'C', etc.

    # Character
    hardness: str  # 'hard', 'soft', 'borderline'
    charge: int  # -1, 0, +1

    # Sterics
    steric_bulk: str  # 'small', 'medium', 'large'
    eliminates_preferred: bool  # tBuO- prefers E2

    # Basicity
    pka_proxy: Optional[float]
    is_strong_base: bool  # pKa of conjugate acid
    is_weak_nucleophile: bool  # DBU, etc.

    # Ambident
    is_ambident: bool  # Enolate, cyanide, nitrite, etc.
    nucleophilic_sites: List[int]  # Atom indices

    # Role prediction
    role_scores: Dict[str, float]  # nucleophile, base, reductant, ligand


@dataclass
class MechanisticClass:
    """
    Minimal viable mechanism (to_do.md #5).

    Assigns reaction to mechanistic family:
    - Polar 2e pathways (SN1, SN2, E1, E2, additions, acyl substitution)
    - Aromatic substitution (SNAr, benzyne, SEAr)
    - Radical 1e pathways (SRN1, Giese, HAT, radical-polar crossover)
    - Organometallic cycles (OA/TM/RE, MI, β-H elim)
    """
    primary_class: str  # 'sn2', 'snar', 'oxidative_addition', etc.
    confidence: float  # 0-1

    # Diagnostic evidence
    evidence: List[str]  # Why this mechanism was chosen

    # Competing mechanisms
    alternatives: List[Tuple[str, float]]  # (mechanism, likelihood)

    # Key diagnostic answers
    is_aryl_vinyl_electrophile: bool
    has_strong_base: bool
    has_beta_h: bool
    has_redox_reagents: bool
    has_metal_catalyst: bool


@dataclass
class SelectivityRisks:
    """
    Selectivity and competing pathways (to_do.md #6, #8).

    Identifies potential issues:
    - Chemoselectivity: which functional group reacts
    - Regioselectivity: Markovnikov, SNAr site, etc.
    - Stereoselectivity: inversion, retention, E/Z
    - Competing pathways: substitution vs elimination
    - Side reactions: overalkylation, rearrangements
    """
    # Selectivity concerns
    chemoselectivity_issues: List[str]
    regioselectivity_issues: List[str]
    stereoselectivity_issues: List[str]

    # Competing pathways
    substitution_vs_elimination: Optional[float]  # Ratio estimate
    overreaction_risk: bool
    rearrangement_risk: bool

    # Side reactions
    side_reactions: List[str]  # Hydrolysis, polymerization, etc.
    catalyst_poisoning_risk: bool


@dataclass
class StructuredReactionAnalysis:
    """
    Complete structured reaction understanding.

    Integrates all aspects from to_do.md framework.
    """
    # Reaction identity
    reaction_smiles: str
    reactants_smiles: str
    products_smiles: str

    # Core transformation
    bond_changes: List[Dict]  # Formed/broken bonds
    functional_group_changes: List[Tuple[str, str]]  # (from, to)
    redox_changes: Optional[str]  # Oxidation/reduction

    # Profiles
    electrophiles: List[ElectrophileProfile]
    nucleophiles: List[NucleophileProfile]

    # Mechanism
    mechanistic_class: MechanisticClass

    # Selectivity
    selectivity_risks: SelectivityRisks

    # Thermodynamics (simple heuristics)
    driving_forces: List[str]  # LG departure, strong bonds, aromaticity, etc.
    is_thermodynamically_favorable: bool

    # Condition constraints (from to_do.md #7)
    solvent_requirements: List[str]  # Polar aprotic, protic, etc.
    temperature_range: Optional[Tuple[float, float]]  # (min, max) °C
    catalyst_requirements: List[str]  # Metal, acid, base, etc.
    atmosphere_requirements: List[str]  # N2, Air, etc.

    # Integration with existing analysis
    tier1_auto_interpretation: Optional[Dict]
    tier2_quick_glance: Optional[Dict]
    tier3_llm_interpretation: Optional[Dict]


class StructuredReactionAnalyzer:
    """
    Deterministic structured reaction analyzer.

    Implements to_do.md framework using:
    - RDKit for molecular structure analysis
    - chemtools for motif detection and taxonomy
    - Rule-based heuristics for mechanism assignment
    """

    def __init__(self):
        """Initialize analyzer with rule databases."""
        self._load_leaving_group_table()
        self._load_nucleophile_table()
        self._load_mechanism_rules()

    def analyze(
        self,
        rxn_smiles: str,
        existing_analysis: Optional[Dict[str, Any]] = None
    ) -> StructuredReactionAnalysis:
        """
        Perform structured reaction analysis.

        Args:
            rxn_smiles: Reaction SMILES
            existing_analysis: Optional result from analyze_reaction_smiles()

        Returns:
            StructuredReactionAnalysis with complete decomposition
        """
        # Parse reaction
        parts = rxn_smiles.split('>>')
        reactants_smiles = parts[0]
        products_smiles = parts[1]

        # Integrate existing analysis if provided
        tier1 = existing_analysis.get('auto_interpretation') if existing_analysis else None
        tier2 = existing_analysis.get('quick_glance') if existing_analysis else None
        tier3 = existing_analysis.get('interpretation') if existing_analysis else None

        # Step 1: Identify bond changes (net transformation)
        bond_changes = self._identify_bond_changes(rxn_smiles)
        fg_changes = self._identify_functional_group_changes(bond_changes)
        redox = self._identify_redox_changes(reactants_smiles, products_smiles)

        # Step 2: Analyze electrophile centers
        electrophiles = self._analyze_electrophiles(
            reactants_smiles,
            bond_changes,
            tier2_hints=tier2
        )

        # Step 3: Analyze nucleophiles
        nucleophiles = self._analyze_nucleophiles(
            reactants_smiles,
            bond_changes,
            tier2_hints=tier2
        )

        # Step 4: Assign mechanistic class
        mechanism = self._assign_mechanism(
            electrophiles,
            nucleophiles,
            bond_changes,
            tier3_hints=tier3
        )

        # Step 5: Assess selectivity risks
        selectivity = self._assess_selectivity_risks(
            electrophiles,
            nucleophiles,
            mechanism
        )

        # Step 6: Identify driving forces
        driving_forces = self._identify_driving_forces(
            bond_changes,
            fg_changes,
            reactants_smiles,
            products_smiles
        )

        # Step 7: Infer condition requirements
        solvent_req = self._infer_solvent_requirements(mechanism)
        temp_range = self._infer_temperature_range(mechanism, electrophiles)
        catalyst_req = self._infer_catalyst_requirements(mechanism)
        atmosphere_req = self._infer_atmosphere_requirements(mechanism)

        return StructuredReactionAnalysis(
            reaction_smiles=rxn_smiles,
            reactants_smiles=reactants_smiles,
            products_smiles=products_smiles,
            bond_changes=bond_changes,
            functional_group_changes=fg_changes,
            redox_changes=redox,
            electrophiles=electrophiles,
            nucleophiles=nucleophiles,
            mechanistic_class=mechanism,
            selectivity_risks=selectivity,
            driving_forces=driving_forces,
            is_thermodynamically_favorable=len(driving_forces) > 0,
            solvent_requirements=solvent_req,
            temperature_range=temp_range,
            catalyst_requirements=catalyst_req,
            atmosphere_requirements=atmosphere_req,
            tier1_auto_interpretation=tier1,
            tier2_quick_glance=tier2,
            tier3_llm_interpretation=tier3
        )

    # Private methods implementing each analysis step
    def _identify_bond_changes(self, rxn_smiles: str) -> List[Dict]:
        """Identify formed/broken bonds using reaction differencing."""
        # Use existing chemtools.reaction_center or rxnmapper
        pass

    def _analyze_electrophiles(
        self,
        reactants_smiles: str,
        bond_changes: List[Dict],
        tier2_hints: Optional[Dict] = None
    ) -> List[ElectrophileProfile]:
        """
        Analyze electrophile centers (to_do.md #2).

        For each carbon with leaving group:
        1. Determine hybridization
        2. Count substitution degree
        3. Check for stabilization (benzylic, allylic, resonance)
        4. Detect EWG activation
        5. Rate leaving group quality
        6. Assess sterics
        7. Predict mechanism preferences
        """
        # Implementation using RDKit
        pass

    def _assign_mechanism(
        self,
        electrophiles: List[ElectrophileProfile],
        nucleophiles: List[NucleophileProfile],
        bond_changes: List[Dict],
        tier3_hints: Optional[Dict] = None
    ) -> MechanisticClass:
        """
        Assign mechanistic class using decision tree.

        Diagnostic questions (to_do.md #5):
        1. Is electrophile aryl/vinyl? → SNAr or metal-catalyzed
        2. Strong base + β-H? → E2
        3. Redox reagents/light? → Radical
        4. etc.
        """
        # Rule-based mechanism assignment
        pass

    # ... more helper methods
```

### 2.2 Integration Point

Modify `reaction_agent/agent.py` to optionally run structured analysis:

```python
def analyze(
    self,
    rxn_smiles: str,
    mode: str = "auto",
    validate: bool = False,
    retry_config: Optional['RetryConfig'] = None,
    structured_analysis: bool = False  # NEW
) -> Dict[str, Any]:
    """
    Analyze a reaction SMILES string.

    Args:
        structured_analysis: Enable deterministic structured analysis (electrophile, nucleophile, etc.)
    """
    # ... existing analysis code ...

    # Optional: Add structured analysis
    if structured_analysis:
        from reaction_agent.structure_analysis import StructuredReactionAnalyzer

        analyzer = StructuredReactionAnalyzer()
        structured = analyzer.analyze(rxn_smiles, existing_analysis=result)
        result['structured_analysis'] = structured

    return result
```

---

## Phase 3: Intelligent Condition Filtering Based on Structured Analysis

### Objective

Use structured analysis to filter/rank condition recommendations intelligently.

### 3.1 Enhanced Bridge: `condition_bridge.py` (Enhancement)

Add filtering logic based on structured analysis:

```python
class RecommendationBridge:
    # ... existing methods ...

    def recommend_with_structure(
        self,
        query: RecommendationQuery,
        structured_analysis: StructuredReactionAnalysis,
        **kwargs
    ) -> HTERecommendationResult:
        """
        Enhanced recommendation using structured analysis.

        Uses electrophile/nucleophile profiles to:
        1. Filter incompatible conditions
        2. Boost compatible conditions
        3. Add warnings for risky conditions
        """
        # Get base recommendations
        recommendations = self.recommend_conditions(query, **kwargs)

        # Apply filters based on mechanism
        mechanism = structured_analysis.mechanistic_class

        if mechanism.primary_class == 'sn2':
            # Boost polar aprotic solvents
            recommendations = self._boost_solvents(
                recommendations,
                ['DMF', 'DMSO', 'DMA', 'NMP'],
                boost_factor=1.5
            )

        elif mechanism.primary_class == 'snar':
            # Require activated aryl electrophile
            if not any(e.is_activated for e in structured_analysis.electrophiles):
                recommendations = self._add_warning(
                    recommendations,
                    "SNAr requires electron-withdrawing activation - consider using metal-catalyzed coupling instead"
                )

        elif mechanism.primary_class == 'oxidative_addition':
            # Require Pd/Ni catalyst
            recommendations = self._filter_by_catalyst(
                recommendations,
                ['Pd', 'Ni'],
                required=True
            )

        # Apply selectivity-based filters
        if structured_analysis.selectivity_risks.substitution_vs_elimination:
            ratio = structured_analysis.selectivity_risks.substitution_vs_elimination
            if ratio > 0.5:  # E2 likely
                recommendations = self._add_warning(
                    recommendations,
                    f"High elimination risk ({ratio:.0%}) - consider bulky nucleophile or lower temperature"
                )

        return recommendations
```

---

## Phase 4: Learning from User Feedback

### Objective

Track which recommendations users select and use feedback to improve ranking.

### 4.1 Feedback Module: `reaction_agent/feedback.py`

```python
"""
User feedback tracking for recommendation improvement.

Stores:
- Which conditions were recommended
- Which conditions were selected by users
- Which conditions worked (if reported)
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional


class FeedbackTracker:
    """
    Track user selections and outcomes.

    Storage format: JSONL in results/feedback/
    """

    def __init__(self, feedback_dir: Path = Path("results/feedback")):
        self.feedback_dir = feedback_dir
        self.feedback_dir.mkdir(parents=True, exist_ok=True)

    def record_recommendation(
        self,
        reaction_smiles: str,
        query: 'RecommendationQuery',
        recommendations: 'HTERecommendationResult',
        user_selection: Optional[int] = None,
        outcome: Optional[str] = None,
        notes: Optional[str] = None
    ):
        """
        Record a recommendation event.

        Args:
            reaction_smiles: Input reaction
            query: Recommendation query
            recommendations: Returned recommendations
            user_selection: Index of selected condition (0-indexed)
            outcome: "success", "failed", "partial", etc.
            notes: User notes
        """
        record = {
            "timestamp": datetime.now().isoformat(),
            "reaction_smiles": reaction_smiles,
            "query": {
                "reaction_types": query.reaction types,
                "primary_type": query.primary_reaction_type,
                "top_k": query.top_k
            },
            "recommendations": [
                {
                    "rank": i + 1,
                    "catalyst": rec.catalyst,
                    "ligand": rec.ligand,
                    "base": rec.base,
                    "solvent": rec.solvent,
                    "avg_z_score": rec.avg_z_score,
                    "confidence_score": rec.confidence_score
                }
                for i, rec in enumerate(recommendations.recommendations[:query.top_k])
            ],
            "user_selection": user_selection,
            "outcome": outcome,
            "notes": notes
        }

        # Append to JSONL file
        feedback_file = self.feedback_dir / "feedback.jsonl"
        with open(feedback_file, 'a') as f:
            f.write(json.dumps(record) + '\n')
```

---

## Phase 5: Web Interface (Optional)

### Objective

Provide web UI for reaction analysis + condition recommendation.

### 5.1 Simple Flask App: `app/web.py`

```python
from flask import Flask, render_template, request, jsonify
from reaction_agent import ReactionSMILESAnalyzer
from reaction_agent.condition_bridge import RecommendationBridge
from llmtools.clients import LLMClient

app = Flask(__name__)

# Initialize on startup
client = LLMClient(provider="openai", model="gpt-4o-mini")
analyzer = ReactionSMILESAnalyzer(client)
bridge = RecommendationBridge()

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    """Analyze reaction + recommend conditions."""
    data = request.json
    rxn_smiles = data.get('reaction_smiles')
    top_k = data.get('top_k', 10)
    validate = data.get('validate', True)

    # Run analysis + recommendation
    result = bridge.analyze_and_recommend(
        rxn_smiles=rxn_smiles,
        analyzer=analyzer,
        top_k=top_k,
        validate=validate
    )

    # Format for frontend
    return jsonify({
        "reaction_analysis": {
            "tier1": result['analysis'].get('auto_interpretation'),
            "tier2": result['analysis'].get('quick_glance'),
            "tier3": result['analysis'].get('interpretation'),
            "validation": result['analysis'].get('validation')
        },
        "recommendations": [
            {
                "rank": i + 1,
                "catalyst": rec.catalyst,
                "ligand": rec.ligand,
                "base": rec.base,
                "solvent": rec.solvent,
                "score": rec.avg_z_score,
                "confidence": rec.confidence_score,
                "success_rate": rec.success_rate,
                "avg_yield": rec.avg_yield
            }
            for i, rec in enumerate(result['recommendations'].recommendations)
        ],
        "metadata": result['metadata']
    })

if __name__ == '__main__':
    app.run(debug=True, port=5000)
```

---

## Implementation Timeline

| Phase | Effort | Dependencies | Priority |
|-------|--------|--------------|----------|
| **Phase 1: Bridge** | 2-3 days | Existing systems | **HIGH** |
| **Phase 2: Structured Analysis** | 1-2 weeks | RDKit, chemtools | **MEDIUM** |
| **Phase 3: Intelligent Filtering** | 3-5 days | Phase 1, 2 | **MEDIUM** |
| **Phase 4: Feedback** | 2-3 days | Phase 1 | **LOW** |
| **Phase 5: Web UI** | 1 week | Phase 1 | **LOW** |

---

## Quick Start: Minimum Viable Product (MVP)

**Goal**: Get end-to-end workflow working in 3 days

### Day 1: Basic Bridge

- Create `condition_bridge.py` with `RecommendationBridge` class
- Implement `extract_recommendation_inputs()`
- Implement `recommend_conditions()` (thin wrapper around HTERecommender)
- Test: Reaction SMILES → Analysis → Recommendation (manual)

### Day 2: CLI Integration

- Add `--recommend` flag to CLI
- Implement `display_recommendations()` function
- Test: `python -m reaction_agent.cli --reaction "SMILES" --recommend --validate`

### Day 3: Polish & Test

- Add reaction type mapping table
- Handle edge cases (no recommendations, parsing errors)
- Test on 10 diverse reactions (Suzuki, amide, Buchwald-Hartwig, etc.)
- Document usage

**Deliverable**: Working end-to-end CLI for reaction analysis + condition recommendation

---

## Success Metrics

1. **Coverage**: % of analyzed reactions that get ≥5 condition recommendations
2. **Relevance**: User satisfaction with recommended conditions (need feedback loop)
3. **Speed**: Total time < 30 seconds (analysis + recommendation)
4. **Cost**: < $0.01 per reaction (analysis + recommendation)
5. **Accuracy**: Recommended conditions match literature for known reactions

---

## Future Enhancements (Beyond Phase 5)

1. **Active Learning**: Train ML model on feedback data to re-rank recommendations
2. **Condition Optimization**: Multi-objective optimization (yield, cost, safety, sustainability)
3. **Protocol Generation**: Automatically generate full experimental protocols
4. **Batch Planning**: Design screening sets for HTE plates
5. **Integration with ELN**: Export to electronic lab notebooks
6. **Real-time Literature Search**: Query Reaxys/SciFinder APIs for latest conditions
7. **Retrosynthesis Integration**: Recommend conditions for each step in synthesis plan

---

This implementation plan provides a clear path from the current state (reaction analysis) to a complete condition recommendation system, leveraging the existing HTE infrastructure while adding the structured understanding framework from to_do.md.
