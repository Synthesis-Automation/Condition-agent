#!/usr/bin/env python
"""
Quantitative Validation Framework for Reaction Analysis

Provides multiple independent metrics to assess model/workflow reliability
beyond self-reported LLM confidence.

Key Validation Approaches:
1. Deterministic metrics (atom mapping quality)
2. Cross-model consistency
3. Specificity scores (named reactions, detailed events)
4. Literature comparison (optional)
5. Ensemble confidence
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional
import json
import numpy as np
from collections import Counter

project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from llmtools import LLMClient
from reaction_agent import ReactionSMILESAnalyzer, analyze_deterministic


class QuantitativeValidator:
    """Quantitative validation of reaction analysis quality."""

    def __init__(self):
        self.validation_results = {}

    def compute_deterministic_quality_score(self, result: Dict) -> float:
        """
        Compute quality score from deterministic analysis (independent of LLM).

        This is the most reliable metric as it's based on cheminformatics tools,
        not LLM self-reporting.

        Returns: 0.0-1.0 score based on mapping quality and bond changes
        """
        tool_facts = result.get('tool_facts', {})

        # 1. Atom mapping confidence (from rxnmapper, independent)
        mapping_qc = tool_facts.get('mapping_qc', {})
        mapping_conf = mapping_qc.get('confidence', 0.0)
        mapping_ok = mapping_qc.get('ok', False)

        # 2. Bond change detection quality
        bond_changes = tool_facts.get('bond_changes', [])
        num_bond_changes = len(bond_changes)

        # 3. Reaction center identification
        reaction_center = tool_facts.get('reaction_center_atoms', [])
        has_reaction_center = len(reaction_center) > 0

        # 4. Parse warnings (fewer is better)
        input_data = result.get('input', {})
        parse_warnings = len(input_data.get('parse_warnings', []))

        # Compute composite score
        score = 0.0

        # Mapping confidence is most important (50% weight)
        if mapping_ok:
            score += 0.5 * mapping_conf

        # Bond changes detected (25% weight)
        if num_bond_changes > 0:
            score += 0.25 * min(num_bond_changes / 5.0, 1.0)

        # Reaction center identified (15% weight)
        if has_reaction_center:
            score += 0.15

        # No parse errors (10% weight)
        if parse_warnings == 0:
            score += 0.10

        return min(score, 1.0)

    def compute_specificity_score(self, result: Dict) -> float:
        """
        Compute specificity score based on detail level.

        More specific = more reliable (named reactions, multiple events, etc.)

        Returns: 0.0-1.0 score
        """
        interp = result.get('interpretation', {})

        score = 0.0

        # 1. Named reaction tags (30% weight)
        tags = interp.get('tags', [])
        named_reactions = ['Suzuki', 'Buchwald', 'Sonogashira', 'SNAr', 'Wittig',
                          'Grignard', 'Heck', 'Stille', 'Negishi', 'Ullmann']

        has_named = any(any(nr.lower() in tag.lower() for nr in named_reactions)
                       for tag in tags)

        if has_named:
            score += 0.30
        elif len(tags) > 0:
            score += 0.15  # Generic tags still better than none

        # 2. Event detection (30% weight)
        events = interp.get('events', [])
        num_events = len(events)

        if num_events > 0:
            score += 0.30 * min(num_events / 3.0, 1.0)

        # 3. Mechanism detail (20% weight)
        mechanism = interp.get('mechanism_summary', [])
        num_steps = len(mechanism)

        if num_steps > 0:
            score += 0.20 * min(num_steps / 3.0, 1.0)

        # 4. Role identification (20% weight)
        roles = interp.get('roles', {})
        filled_roles = sum(1 for v in roles.values() if v and v != "N/A")

        if filled_roles > 0:
            score += 0.20 * min(filled_roles / 3.0, 1.0)

        return min(score, 1.0)

    def compute_cross_model_consistency(self, results: List[Dict]) -> Dict[str, float]:
        """
        Compute consistency across multiple models.

        If multiple models agree, more likely to be correct.

        Args:
            results: List of analysis results from different models

        Returns:
            Dict with consistency metrics
        """
        if len(results) < 2:
            return {'consistency_score': 0.0, 'note': 'Need 2+ models'}

        # Extract key features from each model
        classes = [r.get('interpretation', {}).get('overall_class', 'unknown')
                  for r in results]

        tags_lists = [set(r.get('interpretation', {}).get('tags', []))
                     for r in results]

        confidences = [r.get('interpretation', {}).get('confidence', 0.0)
                      for r in results]

        num_events = [len(r.get('interpretation', {}).get('events', []))
                     for r in results]

        # 1. Class agreement
        class_counts = Counter(classes)
        most_common_class = class_counts.most_common(1)[0]
        class_agreement = most_common_class[1] / len(classes)

        # 2. Tag overlap (Jaccard similarity)
        if len(tags_lists) >= 2:
            all_tags = set().union(*tags_lists)
            if len(all_tags) > 0:
                intersection = set.intersection(*tags_lists)
                tag_overlap = len(intersection) / len(all_tags)
            else:
                tag_overlap = 0.0
        else:
            tag_overlap = 0.0

        # 3. Confidence variance (lower is better)
        conf_std = np.std(confidences) if len(confidences) > 1 else 0.0
        conf_consistency = 1.0 - min(conf_std * 2, 1.0)  # Scale to 0-1

        # 4. Event count consistency
        event_std = np.std(num_events) if len(num_events) > 1 else 0.0
        event_consistency = 1.0 - min(event_std / 2.0, 1.0)

        # Overall consistency score
        consistency_score = (
            0.40 * class_agreement +
            0.30 * tag_overlap +
            0.20 * conf_consistency +
            0.10 * event_consistency
        )

        return {
            'consistency_score': consistency_score,
            'class_agreement': class_agreement,
            'tag_overlap': tag_overlap,
            'confidence_std': conf_std,
            'event_count_std': event_std,
            'models_tested': len(results),
            'consensus_class': most_common_class[0],
            'confidence_range': (min(confidences), max(confidences))
        }

    def compute_warning_score(self, result: Dict) -> float:
        """
        Compute score based on warnings (fewer = better).

        Returns: 0.0-1.0 (1.0 = no warnings)
        """
        interp = result.get('interpretation', {})
        warnings = interp.get('warnings', [])

        # Categorize warnings by severity
        critical_warnings = ['mapping_failed', 'parse_failed']
        moderate_warnings = ['mapping_low_confidence', 'uncertain_bond_changes']

        critical_count = sum(1 for w in warnings
                           if any(cw in str(w).lower() for cw in critical_warnings))
        moderate_count = sum(1 for w in warnings
                           if any(mw in str(w).lower() for mw in moderate_warnings))
        other_count = len(warnings) - critical_count - moderate_count

        # Penalty for warnings
        penalty = (critical_count * 0.4 +
                  moderate_count * 0.2 +
                  other_count * 0.1)

        return max(1.0 - penalty, 0.0)

    def compute_ensemble_confidence(self, results: List[Dict]) -> Dict[str, float]:
        """
        Compute ensemble confidence from multiple models.

        Uses weighted voting based on model reliability and agreement.

        Args:
            results: List of analysis results from different models

        Returns:
            Dict with ensemble metrics
        """
        if len(results) < 2:
            return {
                'ensemble_confidence': results[0].get('interpretation', {}).get('confidence', 0.0),
                'note': 'Single model only'
            }

        # Model reliability weights (based on our testing)
        model_weights = {
            'o3': 1.0,
            'o3-mini': 0.95,
            'gpt-4o': 0.90,
            'gpt-5.2': 0.85,
            'gpt-4o-mini': 0.80,
        }

        weighted_confidences = []
        weights = []

        for result in results:
            model = result.get('metadata', {}).get('model', '')
            confidence = result.get('interpretation', {}).get('confidence', 0.0)

            # Find matching weight
            weight = 0.75  # Default
            for model_key, model_weight in model_weights.items():
                if model_key in model.lower():
                    weight = model_weight
                    break

            weighted_confidences.append(confidence * weight)
            weights.append(weight)

        # Ensemble confidence (weighted average)
        if sum(weights) > 0:
            ensemble_conf = sum(weighted_confidences) / sum(weights)
        else:
            ensemble_conf = np.mean([r.get('interpretation', {}).get('confidence', 0.0)
                                    for r in results])

        # Confidence variance (lower is better)
        confidences = [r.get('interpretation', {}).get('confidence', 0.0) for r in results]
        conf_variance = np.var(confidences)
        conf_std = np.std(confidences)

        # Adjust ensemble confidence based on agreement
        agreement_factor = 1.0 - min(conf_std, 0.5)  # Penalize disagreement
        adjusted_ensemble = ensemble_conf * agreement_factor

        return {
            'ensemble_confidence': ensemble_conf,
            'adjusted_ensemble': adjusted_ensemble,
            'confidence_std': conf_std,
            'confidence_variance': conf_variance,
            'agreement_factor': agreement_factor,
            'individual_confidences': confidences
        }

    def compute_comprehensive_score(self,
                                   result: Dict,
                                   cross_model_results: Optional[List[Dict]] = None) -> Dict:
        """
        Compute comprehensive quality score combining all metrics.

        This is the most reliable overall assessment.

        Args:
            result: Single analysis result
            cross_model_results: Optional list of results from multiple models

        Returns:
            Dict with all quality metrics and overall score
        """
        # Individual scores
        deterministic_score = self.compute_deterministic_quality_score(result)
        specificity_score = self.compute_specificity_score(result)
        warning_score = self.compute_warning_score(result)

        # LLM self-reported confidence (for comparison)
        llm_confidence = result.get('interpretation', {}).get('confidence', 0.0)

        # Cross-model metrics (if available)
        if cross_model_results and len(cross_model_results) > 1:
            consistency = self.compute_cross_model_consistency(cross_model_results)
            ensemble = self.compute_ensemble_confidence(cross_model_results)

            consistency_score = consistency['consistency_score']
            ensemble_confidence = ensemble['adjusted_ensemble']
        else:
            consistency_score = None
            ensemble_confidence = None
            consistency = None
            ensemble = None

        # Compute overall validated score
        # Weight deterministic most heavily (it's independent/objective)
        if consistency_score is not None:
            overall_score = (
                0.35 * deterministic_score +      # Most important
                0.25 * consistency_score +         # Cross-validation
                0.20 * specificity_score +         # Detail level
                0.10 * warning_score +             # Quality flags
                0.10 * llm_confidence              # LLM opinion (least trusted)
            )
        else:
            # Single model only
            overall_score = (
                0.45 * deterministic_score +
                0.30 * specificity_score +
                0.15 * warning_score +
                0.10 * llm_confidence
            )

        # Reliability assessment
        if deterministic_score >= 0.8 and overall_score >= 0.7:
            reliability = "HIGH"
            recommendation = "Results are reliable for production use"
        elif deterministic_score >= 0.6 and overall_score >= 0.5:
            reliability = "MEDIUM"
            recommendation = "Results are reasonable but validate important details"
        elif deterministic_score >= 0.4 and overall_score >= 0.3:
            reliability = "LOW"
            recommendation = "Results are uncertain, manual review required"
        else:
            reliability = "VERY_LOW"
            recommendation = "Results unreliable, trust with caution or reanalyze"

        return {
            'overall_score': overall_score,
            'reliability': reliability,
            'recommendation': recommendation,
            'individual_scores': {
                'deterministic_quality': deterministic_score,
                'specificity': specificity_score,
                'warning_score': warning_score,
                'llm_confidence': llm_confidence,
                'consistency': consistency_score,
                'ensemble_confidence': ensemble_confidence
            },
            'detailed_metrics': {
                'deterministic': {
                    'mapping_confidence': result.get('tool_facts', {}).get('mapping_qc', {}).get('confidence', 0.0),
                    'mapping_ok': result.get('tool_facts', {}).get('mapping_qc', {}).get('ok', False),
                    'num_bond_changes': len(result.get('tool_facts', {}).get('bond_changes', [])),
                    'reaction_center_size': len(result.get('tool_facts', {}).get('reaction_center_atoms', []))
                },
                'specificity': {
                    'has_named_reactions': len(result.get('interpretation', {}).get('tags', [])) > 0,
                    'num_events': len(result.get('interpretation', {}).get('events', [])),
                    'num_mechanism_steps': len(result.get('interpretation', {}).get('mechanism_summary', []))
                },
                'warnings': {
                    'num_warnings': len(result.get('interpretation', {}).get('warnings', [])),
                    'warning_list': result.get('interpretation', {}).get('warnings', [])
                },
                'cross_model': consistency,
                'ensemble': ensemble
            }
        }


def validate_reaction(rxn_smiles: str,
                     models: List[str] = ['gpt-4o-mini', 'gpt-4o', 'o3-mini'],
                     max_tokens: int = 3000) -> Dict:
    """
    Perform comprehensive validation of a reaction analysis.

    Tests multiple models and computes quantitative quality metrics.

    Args:
        rxn_smiles: Reaction SMILES string
        models: List of models to test (for cross-validation)
        max_tokens: Max tokens per model

    Returns:
        Dict with validation results and recommendations
    """
    print(f"\n{'='*80}")
    print(f"  Quantitative Validation")
    print(f"{'='*80}\n")
    print(f"Reaction: {rxn_smiles[:60]}...")
    print(f"Testing {len(models)} models for cross-validation\n")

    validator = QuantitativeValidator()
    results = []

    # Test each model
    for model in models:
        print(f"Testing {model}...", end=' ', flush=True)

        try:
            client = LLMClient(provider="openai", model=model)
            analyzer = ReactionSMILESAnalyzer(client, max_tokens=max_tokens)

            result = analyzer.analyze(rxn_smiles)
            results.append(result)

            conf = result.get('interpretation', {}).get('confidence', 0.0)
            print(f"✓ (confidence: {conf:.3f})")

        except Exception as e:
            print(f"✗ Failed: {e}")

    if not results:
        return {"error": "All models failed"}

    # Compute comprehensive validation
    print(f"\nComputing validation metrics...")

    # Use first result as primary
    primary_result = results[0]

    # Compute comprehensive score
    comprehensive = validator.compute_comprehensive_score(
        primary_result,
        cross_model_results=results if len(results) > 1 else None
    )

    # Print results
    print(f"\n{'='*80}")
    print(f"  VALIDATION RESULTS")
    print(f"{'='*80}\n")

    print(f"Overall Validated Score: {comprehensive['overall_score']:.3f}")
    print(f"Reliability: {comprehensive['reliability']}")
    print(f"Recommendation: {comprehensive['recommendation']}\n")

    print(f"Individual Scores:")
    for metric, score in comprehensive['individual_scores'].items():
        if score is not None:
            print(f"  {metric:25s}: {score:.3f}")

    print(f"\nKey Metrics:")
    det = comprehensive['detailed_metrics']['deterministic']
    print(f"  Mapping Confidence: {det['mapping_confidence']:.3f} ({'OK' if det['mapping_ok'] else 'Failed'})")
    print(f"  Bond Changes:       {det['num_bond_changes']}")
    print(f"  Reaction Center:    {det['reaction_center_size']} atoms")

    spec = comprehensive['detailed_metrics']['specificity']
    print(f"  Named Reactions:    {spec['has_named_reactions']}")
    print(f"  Events Detected:    {spec['num_events']}")
    print(f"  Mechanism Steps:    {spec['num_mechanism_steps']}")

    if comprehensive['detailed_metrics']['cross_model']:
        cm = comprehensive['detailed_metrics']['cross_model']
        print(f"\nCross-Model Validation:")
        print(f"  Consistency Score:  {cm['consistency_score']:.3f}")
        print(f"  Class Agreement:    {cm['class_agreement']:.3f}")
        print(f"  Consensus Class:    {cm['consensus_class']}")
        print(f"  Confidence Range:   {cm['confidence_range'][0]:.3f} - {cm['confidence_range'][1]:.3f}")

    print(f"\n{'='*80}\n")

    return comprehensive


def compare_to_literature(result: Dict, known_class: str, known_mechanism: List[str]) -> Dict:
    """
    Compare analysis results to known literature values.

    This provides ground truth validation if you have known correct answers.

    Args:
        result: Analysis result
        known_class: Known correct reaction class
        known_mechanism: List of known mechanistic steps

    Returns:
        Dict with comparison metrics
    """
    interp = result.get('interpretation', {})

    # Class match
    predicted_class = interp.get('overall_class', '')
    class_match = predicted_class.lower() == known_class.lower()

    # Mechanism similarity (simple keyword matching)
    predicted_mechanism = interp.get('mechanism_summary', [])

    # Count matching keywords
    known_keywords = set(' '.join(known_mechanism).lower().split())
    predicted_keywords = set(' '.join(predicted_mechanism).lower().split())

    if len(known_keywords.union(predicted_keywords)) > 0:
        mechanism_similarity = len(known_keywords.intersection(predicted_keywords)) / \
                             len(known_keywords.union(predicted_keywords))
    else:
        mechanism_similarity = 0.0

    # Overall literature agreement
    lit_score = (0.6 * (1.0 if class_match else 0.0) +
                0.4 * mechanism_similarity)

    return {
        'literature_agreement': lit_score,
        'class_match': class_match,
        'mechanism_similarity': mechanism_similarity,
        'predicted_class': predicted_class,
        'known_class': known_class
    }


if __name__ == "__main__":
    # Example usage
    if not os.getenv("OPENAI_API_KEY"):
        print("Set OPENAI_API_KEY to run validation")
        sys.exit(1)

    # Test reaction
    rxn = "Clc1nc2ccccc2s1.Cn1ccnc1>>CN1C=C[N+](C2=NC3=CC=CC=C3S2)=C1"

    # Run validation
    validation = validate_reaction(rxn, models=['gpt-4o-mini', 'gpt-4o', 'o3-mini'])

    # Save results
    output_dir = Path("reaction_agent/results")
    output_dir.mkdir(parents=True, exist_ok=True)

    with open(output_dir / "validation_example.json", 'w') as f:
        json.dump(validation, f, indent=2)

    print(f"Validation results saved to: {output_dir}/validation_example.json")
