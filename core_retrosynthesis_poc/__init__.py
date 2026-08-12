"""Reaction-core-derived one-step retrosynthesis proof of concept."""

from .compiler import CompilationResult, compile_core_templates
from .condition_ranking import (
    ConditionRankedRetrosynthesisCandidate,
    RetrosynthesisConditionEvidence,
    rank_retrosynthesis_candidates_with_conditions,
    recommend_retrosynthesis_conditions,
)
from .comparison import run_comparison, split_by_reference
from .ensemble import EnsembleCandidate, disconnect_ensemble
from .coverage_audit import audit_operator_library_coverage
from .full_scale import (
    FullScaleBuildConfig,
    build_full_scale_operator_library,
    compile_operator_shard,
    merge_operator_shards,
)
from .generic_compiler import (
    GenericCompilationResult,
    GenericReactionIdentity,
    analyze_generic_reaction,
    compile_generic_templates,
)
from .generic_library import (
    build_generic_library,
    load_generic_library,
    save_generic_library,
)
from .generic_models import (
    GenericCoreTemplate,
    GenericDisconnectionCandidate,
    GenericGraphOperator,
    GenericRetrievalIndex,
    GenericSearchDiagnostics,
    GenericHandleCompletionGroup,
    GenericTemplateLibrary,
)
from .generic_search import (
    disconnect_generic_target,
    disconnect_generic_target_detailed,
    disconnect_operator_ladder,
    rank_operator_site_diverse,
)
from .html_report import render_comparison_html, write_comparison_html
from .library import build_library, load_library, save_library
from .models import (
    CenterReactivityContext,
    CoreDisconnectionCandidate,
    CoreLibraryBuildReport,
    CoreTemplate,
    CoreTemplateLibrary,
    CoreTemplatePrecedent,
    TemplateContext,
)
from .multistep import (
    MultistepRetrosynthesisResult,
    MultistepRetrosynthesisRoute,
    MultistepSearchDiagnostics,
    RetrosynthesisRouteStep,
    StartingMaterialAssessment,
    plan_multistep_routes,
)
from .operator_benchmark import (
    load_operator_rows,
    run_operator_coverage_benchmark,
)
from .ranking_policy import (
    RetrosynthesisRankingPolicy,
    load_retrosynthesis_ranking_policy,
)
from .search import disconnect_target
from .selectivity_poc import (
    ChoiceModelTrainingReport,
    ConditionalEditChoiceModel,
    FunctionalGroupCompetitionOutcome,
    FunctionalGroupCompetitionWarning,
    RankedReactionOutcome,
    ReactionChoiceSet,
    ReactionOutcomeCandidate,
    SelectivityAssessment,
    build_reaction_choice_set,
    build_reaction_choice_set_from_record,
    condition_tokens_from_mapping,
    condition_tokens_from_recipe,
    detect_functional_group_competition,
)
from .sources import (
    LIBRARY_MODES,
    iter_library_rows,
    resolve_library_mode,
    source_shard_files,
)

__all__ = [
    "CenterReactivityContext",
    "CompilationResult",
    "ChoiceModelTrainingReport",
    "ConditionalEditChoiceModel",
    "ConditionRankedRetrosynthesisCandidate",
    "CoreDisconnectionCandidate",
    "CoreLibraryBuildReport",
    "CoreTemplate",
    "CoreTemplateLibrary",
    "CoreTemplatePrecedent",
    "EnsembleCandidate",
    "FullScaleBuildConfig",
    "FunctionalGroupCompetitionOutcome",
    "FunctionalGroupCompetitionWarning",
    "GenericCompilationResult",
    "GenericCoreTemplate",
    "GenericDisconnectionCandidate",
    "GenericGraphOperator",
    "GenericRetrievalIndex",
    "GenericReactionIdentity",
    "GenericSearchDiagnostics",
    "GenericHandleCompletionGroup",
    "GenericTemplateLibrary",
    "LIBRARY_MODES",
    "MultistepRetrosynthesisResult",
    "MultistepRetrosynthesisRoute",
    "MultistepSearchDiagnostics",
    "RankedReactionOutcome",
    "ReactionChoiceSet",
    "ReactionOutcomeCandidate",
    "RetrosynthesisConditionEvidence",
    "RetrosynthesisRouteStep",
    "RetrosynthesisRankingPolicy",
    "SelectivityAssessment",
    "StartingMaterialAssessment",
    "TemplateContext",
    "analyze_generic_reaction",
    "audit_operator_library_coverage",
    "build_full_scale_operator_library",
    "build_library",
    "build_reaction_choice_set",
    "build_reaction_choice_set_from_record",
    "build_generic_library",
    "compile_generic_templates",
    "compile_core_templates",
    "compile_operator_shard",
    "condition_tokens_from_mapping",
    "condition_tokens_from_recipe",
    "disconnect_target",
    "detect_functional_group_competition",
    "disconnect_ensemble",
    "disconnect_generic_target",
    "disconnect_generic_target_detailed",
    "disconnect_operator_ladder",
    "load_library",
    "load_generic_library",
    "load_operator_rows",
    "load_retrosynthesis_ranking_policy",
    "iter_library_rows",
    "merge_operator_shards",
    "plan_multistep_routes",
    "rank_retrosynthesis_candidates_with_conditions",
    "recommend_retrosynthesis_conditions",
    "rank_operator_site_diverse",
    "run_comparison",
    "run_operator_coverage_benchmark",
    "render_comparison_html",
    "resolve_library_mode",
    "save_library",
    "save_generic_library",
    "split_by_reference",
    "source_shard_files",
    "write_comparison_html",
]
