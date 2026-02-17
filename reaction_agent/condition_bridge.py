"""
Condition recommendation bridge connecting reaction analysis to HTE experiments.

This module provides a simple interface between the reaction analysis system
(Tiers 1-4) and the HTE-based condition recommendation system.

The HTERecommender handles all taxonomy/motif detection internally via:
- featurize_molecule() for reactant motif detection
- featurize_reaction() for reaction type and event detection
- Hierarchical matching (reaction_key → motifs → fallback)

The bridge simply:
1. Parses reaction SMILES from analysis output
2. Passes to HTERecommender
3. Returns recommendations

Key Requirement:
    Product SMILES is MANDATORY for recommendation. The HTERecommender
    requires full reaction SMILES to detect reaction type and motifs
    via the taxonomy system.
"""

from typing import Dict, Any, Optional, Tuple
from pathlib import Path
from chemtools.recommend.recommender import HTERecommender, HTERecommendationResult


class ConditionBridge:
    """
    Bridge between reaction analysis and HTE experiments database.

    Phase 1 Implementation:
    - Connects to experiments database only (source_group="experiments")
    - Leverages HTERecommender's internal taxonomy system
    - Supports basic recommendation workflow

    Example:
        >>> from llmtools.clients import LLMClient
        >>> from reaction_agent import ReactionSMILESAnalyzer
        >>> from reaction_agent.condition_bridge import ConditionBridge
        >>>
        >>> client = LLMClient(provider="openai", model="gpt-4o-mini")
        >>> analyzer = ReactionSMILESAnalyzer(client)
        >>> bridge = ConditionBridge()
        >>>
        >>> rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
        >>> result = bridge.analyze_and_recommend(rxn, analyzer, validate=True)
        >>>
        >>> print(f"Reaction type: {result['metadata']['detected_reaction_type']}")
        >>> print(f"Found {len(result['recommendations'].recommendations)} conditions")
    """

    def __init__(self, hte_db_path: str = "data/HTE_db"):
        """
        Initialize bridge with HTE database.

        Args:
            hte_db_path: Path to HTE database directory (default: "data/HTE_db")
        """
        self.hte_db_path = Path(hte_db_path)
        self.recommender = HTERecommender(str(self.hte_db_path))

    def extract_smiles(
        self,
        analysis_result: Dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str], Optional[str]]:
        """
        Extract reactant and product SMILES from analysis result.

        Parses rxn_smiles_clean from analysis to get separate reactants and product.
        For multi-reactant reactions, uses first two reactants (arbitrary order is OK,
        HTERecommender handles this).

        Args:
            analysis_result: Output from analyze_reaction_smiles()

        Returns:
            Tuple of (reactant_a_smiles, reactant_b_smiles, product_smiles)
            Returns (None, None, None) if parsing fails

        Example:
            >>> result = analyzer.analyze("Brc1ccccc1.B(O)(O)c1ccccc1>>biphenyl")
            >>> reactant_a, reactant_b, product = bridge.extract_smiles(result)
            >>> # reactant_a = "Brc1ccccc1"
            >>> # reactant_b = "B(O)(O)c1ccccc1"
            >>> # product = "c1ccc(-c2ccccc2)cc1"
        """
        input_data = analysis_result.get('input', {})
        rxn_smiles_clean = input_data.get('rxn_smiles_clean', '')

        if not rxn_smiles_clean or '>>' not in rxn_smiles_clean:
            return None, None, None

        # Split reaction SMILES: "reactants>>products"
        parts = rxn_smiles_clean.split('>>')
        reactants_str = parts[0] if len(parts) > 0 else ''
        products_str = parts[1] if len(parts) > 1 else ''

        # Parse reactants (split by '.')
        reactants = [r.strip() for r in reactants_str.split('.') if r.strip()]
        reactant_a = reactants[0] if len(reactants) > 0 else None
        reactant_b = reactants[1] if len(reactants) > 1 else None

        # Parse products (use first product as main product)
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
        Get condition recommendations from HTE experiments database.

        The HTERecommender internally:
        1. Constructs full reaction SMILES from reactants + product
        2. Calls featurize_reaction() to detect reaction type and motifs
        3. Extracts reacted_motifs, formed_motifs, reaction_key
        4. Matches against HTE database using taxonomy-based matching
        5. Ranks by Z-score and returns top_k recommendations

        Args:
            analysis_result: Output from analyze_reaction_smiles()
            top_k: Number of recommendations to return (default: 10)
            min_experiments: Minimum experiments required for a condition (default: 2)
            reaction_type_filter: Optional reaction type override (usually not needed,
                                 recommender auto-detects with high confidence)

        Returns:
            HTERecommendationResult with:
            - recommendations: List of ranked condition sets
            - predicted_reaction_type: Auto-detected reaction type
            - reaction_type_confidence: Confidence score (0.0-1.0)
            - reacted_motifs: Consumed motifs (e.g., ["Ar-Br", "Ar-B(OH)2"])
            - formed_motifs: Created motifs (e.g., ["Ar-Ar"])
            - total_matching_experiments: Total experiments in database

        Raises:
            ValueError: If reactant or product SMILES cannot be extracted

        Example:
            >>> result = analyzer.analyze("Brc1ccccc1.B(O)(O)c1ccccc1>>biphenyl")
            >>> recommendations = bridge.recommend_conditions(result, top_k=5)
            >>> print(f"Reaction: {recommendations.predicted_reaction_type}")
            >>> print(f"Confidence: {recommendations.reaction_type_confidence:.2f}")
            >>> for rec in recommendations.recommendations[:3]:
            ...     print(f"  {rec.catalyst} + {rec.ligand} (score: {rec.avg_z_score:.2f})")
        """
        # Extract SMILES components
        reactant_a, reactant_b, product = self.extract_smiles(analysis_result)

        if not reactant_a:
            raise ValueError(
                "Could not extract reactant SMILES from analysis result. "
                "Expected 'input.rxn_smiles_clean' with format 'A.B>>C'"
            )

        if not product:
            raise ValueError(
                "Product SMILES required for condition recommendation. "
                "The HTERecommender needs full reaction SMILES (with product) "
                "to detect reaction type and motifs via the taxonomy system."
            )

        # Call HTERecommender
        # It handles all taxonomy/motif detection internally via featurize_reaction()
        return self.recommender.recommend(
            reactant_a_smiles=reactant_a,
            reactant_b_smiles=reactant_b,
            product_smiles=product,  # CRITICAL for taxonomy-based detection
            top_k=top_k,
            min_experiments=min_experiments,
            reaction_type_filter=reaction_type_filter,
            source_group="experiments"  # Phase 1: experiments only
        )

    def analyze_and_recommend(
        self,
        rxn_smiles: str,
        analyzer: 'ReactionSMILESAnalyzer',
        top_k: int = 10,
        validate: bool = True,
        mode: str = "auto"
    ) -> Dict[str, Any]:
        """
        End-to-end workflow: Analyze reaction → Recommend conditions.

        Combines reaction analysis (Tiers 1-4) with condition recommendation
        in a single call.

        Args:
            rxn_smiles: Reaction SMILES string (reactants>>products)
            analyzer: ReactionSMILESAnalyzer instance
            top_k: Number of condition recommendations (default: 10)
            validate: Enable Tier 4 RDKit validation (default: True)
            mode: Analysis mode - "auto", "tier2", "tier3" (default: "auto")

        Returns:
            Dict with:
            - analysis: Full reaction analysis result (Tiers 1-4)
            - recommendations: HTERecommendationResult object
            - metadata: Dict with timing and detection info
                - analysis_time_s: Time for analysis
                - recommendation_time_s: Time for recommendation
                - total_time_s: Total time
                - detected_reaction_type: Auto-detected reaction type
                - reaction_type_confidence: Confidence score
                - reacted_motifs: Consumed motifs tuple
                - formed_motifs: Created motifs tuple

        Raises:
            ValueError: If SMILES parsing or analysis fails

        Example:
            >>> from llmtools.clients import LLMClient
            >>> from reaction_agent import ReactionSMILESAnalyzer
            >>> from reaction_agent.condition_bridge import ConditionBridge
            >>>
            >>> client = LLMClient(provider="openai", model="gpt-4o-mini")
            >>> analyzer = ReactionSMILESAnalyzer(client)
            >>> bridge = ConditionBridge()
            >>>
            >>> rxn = "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"
            >>> result = bridge.analyze_and_recommend(rxn, analyzer, top_k=5)
            >>>
            >>> print(f"Analysis: {result['analysis']['quick_glance']['summary']}")
            >>> print(f"Reaction type: {result['metadata']['detected_reaction_type']}")
            >>> print(f"Confidence: {result['metadata']['reaction_type_confidence']:.2f}")
            >>> print(f"Found {len(result['recommendations'].recommendations)} conditions")
        """
        import time

        start_time = time.time()

        # Step 1: Analyze reaction (Tiers 1-4)
        analysis_result = analyzer.analyze(rxn_smiles, mode=mode, validate=validate)
        analysis_time = time.time() - start_time

        # Step 2: Get condition recommendations
        rec_start = time.time()
        recommendations = self.recommend_conditions(analysis_result, top_k=top_k)
        rec_time = time.time() - rec_start

        # Step 3: Package results
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
                "formed_motifs": recommendations.formed_motifs,
                "total_experiments": recommendations.total_matching_experiments
            }
        }


__all__ = ['ConditionBridge']
