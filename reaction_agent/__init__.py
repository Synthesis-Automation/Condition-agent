"""
Reaction SMILES Analysis Agent - POC

A self-contained module for analyzing reaction SMILES using:
- Deterministic cheminformatics (RDKit + rxnmapper)
- LLM interpretation (GPT-4 family)

Design: docs/reaction_smiles_analysis_agent_simple_v1.md

Components:
    core: Deterministic analysis (cleaning, mapping, bond changes)
    prompts: LLM prompt templates
    agent: LLM integration and orchestration

Usage:
    from llmtools.clients import LLMClient
    from reaction_agent import ReactionSMILESAnalyzer

    client = LLMClient(provider="openai", model="gpt-4o-mini")
    analyzer = ReactionSMILESAnalyzer(client)

    result = analyzer.analyze("Brc1ccccc1.Nc1ccccc1>>c1ccc(Nc2ccccc2)cc1")
    print(result["interpretation"]["overall_class"])
"""

__version__ = "0.1.0"

from reaction_agent.smiles_pipeline import (
    ReactionPipeline,
    PipelineResult,
    NormalizationResult,
    FeaturizationResult,
    LLMFallbackResult,
)
from reaction_agent.pipeline_eval import (
    QualityConfig,
    QualityReport,
    QualityEvaluator,
)
from reaction_agent.taxonomy_prompts import TaxonomyContext
from reaction_agent.reactivity_profile import (
    ReactivityProfile,
    ElectrophileProfile,
    NucleophileProfile,
    MechanismAnalysis,
    TransformationRecord,
)
from reaction_agent.reasoning_agent import ReactionReasoningAgent
from reaction_agent.core import (
    clean_reaction_smiles,
    map_reaction,
    extract_bond_changes,
    analyze_reaction_smiles as analyze_deterministic,
    CleanReport,
    MappingReport,
    BondChange,
    BondChangeReport,
    SPECTATORS,
)
from reaction_agent.prompts import get_template, REACTION_SMILES_ANALYSIS
from reaction_agent.agent import (
    ReactionSMILESAnalyzer,
    analyze_reaction_smiles,
)

__all__ = [
    # Pipeline
    "ReactionPipeline",
    "PipelineResult",
    "NormalizationResult",
    "FeaturizationResult",
    "LLMFallbackResult",
    "QualityConfig",
    "QualityReport",
    "QualityEvaluator",
    "TaxonomyContext",
    # Reasoning agent
    "ReactionReasoningAgent",
    "ReactivityProfile",
    "ElectrophileProfile",
    "NucleophileProfile",
    "MechanismAnalysis",
    "TransformationRecord",
    # Core deterministic functions
    "clean_reaction_smiles",
    "map_reaction",
    "extract_bond_changes",
    "analyze_deterministic",
    # Data classes
    "CleanReport",
    "MappingReport",
    "BondChange",
    "BondChangeReport",
    # Constants
    "SPECTATORS",
    # Prompts
    "get_template",
    "REACTION_SMILES_ANALYSIS",
    # Agent
    "ReactionSMILESAnalyzer",
    "analyze_reaction_smiles",
    # Version
    "__version__",
]
