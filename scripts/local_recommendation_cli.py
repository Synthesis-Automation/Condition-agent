"""
Local Recommendation Tester
===========================

This script mirrors `interactive_recommendation_cli.py` but calls the
ChemTools APIs directly rather than going through the FastAPI server.
It prompts for a reaction SMILES and reaction family, runs rule-based,
ML, fusion, and protocol recommendation pipelines, saves their JSON results,
and prints a compact summary to the console.

Usage:
    python scripts/local_recommendation_cli.py

Requirements:
    - ChemTools library dependencies installed (same as the FastAPI app)
    - Optional SCDB JSON database for rule-based matching
    - Protocol index built (python -m chemtools.protocol.cli build)
"""

from __future__ import annotations

import argparse
import io
import os
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chemtools import chem, output_formatter, recommend

try:
    from chemtools.protocol import ProtocolRecommender
    HAS_PROTOCOL = True
except ImportError:
    HAS_PROTOCOL = False
    ProtocolRecommender = None

try:
    from llmtools.clients import LLMClient
    from llmtools.recommendation_llm import (
        synthesize_recommendations_llm,
        convert_llm_synthesis_to_standard_format,
    )
    HAS_LLM = True
except ImportError:
    HAS_LLM = False
    LLMClient = None
    synthesize_recommendations_llm = None
    convert_llm_synthesis_to_standard_format = None

try:
    from chemtools.rule_scdb_matcher.loader import SchemeConditionDBError
except Exception:  # pragma: no cover - fallback path
    class SchemeConditionDBError(RuntimeError):
        """Fallback when SCDB loader is unavailable."""
        pass

try:
    from scripts.recommendation_cli_utils import (
        DEFAULT_SCDB_PATH,
        FUSION_VARIANTS_DEFAULT,
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_fusion,
        summarize_ml,
        summarize_rule,
        summarize_protocol,
        summarize_llm_synthesis,
    )
except ModuleNotFoundError:
    sys.path.append(str(HERE))
    from recommendation_cli_utils import (
        DEFAULT_SCDB_PATH,
        FUSION_VARIANTS_DEFAULT,
        K_DEFAULT,
        LIMIT_DEFAULT,
        choose_reaction_type,
        prompt_smiles,
        save_json,
        slugify_label,
        summarize_fusion,
        summarize_ml,
        summarize_rule,
        summarize_protocol,
        summarize_llm_synthesis,
    )


if sys.platform == "win32":
    # Ensure UTF-8 output when running in Windows terminals.
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")


_ENV_SCDB_DEFAULT = (
    os.environ.get("SCDB_MATCH_DB_PATH", "cn_coupling_pd_db.json").strip()
    or "cn_coupling_pd_db.json"
)


def _resolve_rule_db(candidate: Optional[str]) -> Optional[str]:
    """Determine the best available rule DB path."""
    if candidate:
        return candidate
    if Path(DEFAULT_SCDB_PATH).exists():
        return DEFAULT_SCDB_PATH
    if Path(_ENV_SCDB_DEFAULT).exists():
        return _ENV_SCDB_DEFAULT
    return None


def local_rule_based_match(reaction: str, db_path: Optional[str]) -> Dict[str, Any]:
    """Replicate the /match endpoint using in-process ChemTools calls."""
    target_db = _resolve_rule_db(db_path) or _ENV_SCDB_DEFAULT
    start_time = time.perf_counter()

    try:
        database = chem.rules.load_database(target_db, cache=True)
        result = chem.rules.match(database, reaction)
        payload = result.to_json_dict()
        processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)
        database_label = Path(target_db).name if target_db else "SchemeConditionDB"
        return output_formatter.format_rule_match_result(
            reaction_smiles=reaction,
            match_result=payload,
            requested_type=None,
            database_name=database_label,
            processing_time_ms=processing_time_ms,
        )
    except FileNotFoundError as exc:
        return {
            "error": f"Database file not found: {target_db}",
            "exception": repr(exc),
        }
    except SchemeConditionDBError as exc:
        return {"error": f"Rule-based match failed: {exc}"}
    except Exception as exc:  # pragma: no cover - unexpected path
        return {"error": f"Unexpected error during rule-based match: {exc}"}


def _format_conditions_output(
    raw_data: Dict[str, Any],
    reaction: str,
    reaction_type: Optional[str],
    limit: int,
    elapsed_ms: float,
) -> Dict[str, Any]:
    """Apply the same formatting as /api/v1/recommend/conditions."""
    detection = raw_data.get("detection", {})
    # Try multiple keys for detected type
    detected_type = (
        detection.get("family") 
        or detection.get("detected_reaction_type")
        or detection.get("reaction_type")
        or "Unknown"
    )
    confidence = detection.get("confidence", 0.0)

    recommendations_data: List[Dict[str, Any]] = []
    for rec in raw_data.get("recommendations", [])[:limit]:
        summary = rec.get("summary", {})
        chemicals = rec.get("chemicals", [])
        conditions_info = rec.get("conditions", {})

        reagents = []
        for chemical in chemicals:
            reagents.append(
                {
                    "id": chemical.get("uid", chemical.get("cas")),
                    "role": chemical.get("role", "reagent"),
                    "name": chemical.get("name"),
                    "abbreviation": None,
                    "cas": chemical.get("cas"),
                    "smiles": chemical.get("smiles"),
                    "equivalents": None,
                }
            )

        conditions: Dict[str, Any] = {}
        if conditions_info.get("temperature") is not None:
            conditions["temperature"] = {
                "value": conditions_info["temperature"],
                "unit": "°C",
            }
        if conditions_info.get("time") is not None:
            conditions["time"] = {
                "value": conditions_info["time"],
                "unit": "hours",
            }
        if conditions_info.get("atmosphere"):
            conditions["atmosphere"] = conditions_info["atmosphere"]

        recommendations_data.append(
            {
                "rank": rec.get("rank", len(recommendations_data) + 1),
                "confidence": (
                    summary.get("confidence", 0.0) / 100.0
                    if summary.get("confidence")
                    else 0.0
                ),
                "reagents": reagents,
                "conditions": conditions,
                "precedent_count": (
                    summary.get("support", {}).get("count")
                    if isinstance(summary.get("support"), dict)
                    else summary.get("support", 0)
                ),
            }
        )

    formatted_output = output_formatter.format_ml_output(
        reaction_smiles=reaction,
        requested_type=reaction_type,
        detected_type=detected_type,
        detection_confidence=confidence,
        recommendations_data=recommendations_data,
        processing_time_ms=elapsed_ms,
    )
    
    # Preserve precedents_used section from raw_data if available
    if "precedents_used" in raw_data:
        formatted_output["precedents_used"] = raw_data["precedents_used"]
    
    return formatted_output


def local_ml_recommendation(
    reaction: str,
    reaction_type: Optional[str],
    k_value: int,
    limit: int,
    rerank_strategy: str = 'rule',
    filter_unknown_reagents: bool = False,
) -> Dict[str, Any]:
    """Replicate the /api/v1/recommend/conditions endpoint locally."""
    start_time = time.perf_counter()
    try:
        raw_data = chem.recommend.conditions(
            reaction=reaction,
            reaction_type=reaction_type,
            k=k_value,
            limit=limit,
            relax={},
            constraints={},
            rerank_strategy=rerank_strategy,
            filter_unknown_reagents=filter_unknown_reagents,
        )
    except Exception as exc:
        return {"error": f"Local ML recommendation failed: {exc}"}

    elapsed_ms = (time.perf_counter() - start_time) * 1000
    try:
        return _format_conditions_output(raw_data, reaction, reaction_type, limit, elapsed_ms)
    except Exception as exc:  # pragma: no cover - formatting failure
        return {"error": f"Failed to format ML recommendation: {exc}"}


def local_fusion_recommendation(
    reaction: str,
    k_value: int,
    max_variants: int,
) -> Dict[str, Any]:
    """
    DEPRECATED: Fusion recommendation has been replaced.
    
    This function now redirects to the standard recommendation with rule-based reranking.
    Use local_ml_recommendation() with rerank_strategy='rule' instead.
    """
    import warnings
    warnings.warn(
        "Fusion recommendation is deprecated. Use local_ml_recommendation() "
        "with rerank_strategy='rule' instead.",
        DeprecationWarning,
        stacklevel=2
    )
    
    start_time = time.perf_counter()
    try:
        # Redirect to standard recommendation with rule reranking
        result = recommend.recommend_from_reaction(
            reaction=reaction,
            k=k_value,
            rerank_strategy='rule',  # Use rule reranking instead of fusion
            filter_unknown_reagents=False,
            max_variants=max_variants,
            relax={},
            constraint_rules={},
        )
    except Exception as exc:
        return {"error": f"Local fusion recommendation failed: {exc}"}

    processing_time_ms = round((time.perf_counter() - start_time) * 1000, 2)

    # Add deprecation notice
    result['_deprecated'] = {
        'message': 'Fusion is deprecated. Use rerank_strategy="rule" instead.',
        'migration': 'Use --rerank rule instead of --strategy fusion'
    }

    return output_formatter.format_fusion_output(
        reaction_smiles=reaction,
        requested_type=None,
        fusion_result=result,
        processing_time_ms=processing_time_ms,
    )


def local_protocol_recommendation(
    reaction: str,
    k_value: int,
    tags: Optional[List[str]] = None,
    reaction_family: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run protocol-based recommendation locally.
    
    Uses DRFP similarity to find matching experimental protocols.
    Returns standard JSON format (same as ML/Rule modes).
    
    Args:
        reaction: Reaction SMILES
        k_value: Number of recommendations to return
        tags: Optional list of tags to filter by
        reaction_family: Optional reaction family to filter by
    
    Returns:
        Standard format result with meta, input, detection, recommended_conditions
    """
    if not HAS_PROTOCOL:
        return {
            "error": "Protocol recommendation not available. Install with: pip install drfp",
            "meta": {
                "model": "Protocol-DRFP",
                "status": "error"
            }
        }
    
    try:
        # Initialize recommender (loads index)
        recommender = ProtocolRecommender()
        
        # Get recommendations with details (includes extracted conditions)
        result = recommender.recommend_with_details(
            reaction_smiles=reaction,
            k=k_value,
            reaction_family=reaction_family,
            tags=tags,
            use_standard_format=True  # Ensure standard format
        )
        
        return result
        
    except FileNotFoundError as exc:
        return {
            "error": f"Protocol index not found. Run: python -m chemtools.protocol.cli build",
            "details": str(exc),
            "meta": {
                "model": "Protocol-DRFP",
                "status": "error"
            }
        }
    except Exception as exc:
        return {
            "error": f"Protocol recommendation failed: {exc}",
            "meta": {
                "model": "Protocol-DRFP",
                "status": "error"
            }
        }


def local_llm_synthesis(
    reaction: str,
    ml_result: Optional[Dict[str, Any]] = None,
    rule_result: Optional[Dict[str, Any]] = None,
    protocol_result: Optional[Dict[str, Any]] = None,
    constraints: Optional[Dict[str, Any]] = None,
    llm_provider: str = "aliyun",
    llm_model: str = "deepseek-v3.2-exp",
    prompt_version: str = "v2",
    requested_type: Optional[str] = None,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Run LLM-enhanced multi-source synthesis locally.
    
    Combines ML, Rule, and Protocol recommendations using LLM intelligence
    to provide chemistry-aware synthesis with explanations, warnings, and backups.
    
    Args:
        reaction: Reaction SMILES
        ml_result: ML recommendation result (from local_ml_recommendation)
        rule_result: Rule-based result (from local_rule_based_match)
        protocol_result: Protocol result (from local_protocol_recommendation)
        constraints: Optional user constraints (scale, cost, air_sensitivity, etc.)
        llm_provider: LLM provider ("aliyun" or "openai")
        llm_model: LLM model name (default: "deepseek-v3.2-exp")
        prompt_version: Prompt version ("v1" or "v2", default "v2" for optimized)
        requested_type: Requested reaction type for standard format output
    
    Returns:
        Tuple of (analysis_result, standard_format_result):
            - analysis_result: LLM synthesis analysis (original format)
            - standard_format_result: Standard format for robotic execution
    """
    error_result = {
        "error": "LLM synthesis not available. Install with: pip install dashscope openai",
        "meta": {
            "model": "LLM-Synthesis",
            "status": "error"
        }
    }
    
    if not HAS_LLM:
        return error_result, error_result
    
    try:
        # Initialize LLM client
        llm_client = LLMClient(provider=llm_provider, model=llm_model)
        
        # Run synthesis
        start_time = time.perf_counter()
        analysis_result = synthesize_recommendations_llm(
            reaction_smiles=reaction,
            ml_results=ml_result,
            rule_results=rule_result,
            protocol_results=protocol_result,
            constraints=constraints,
            llm_client=llm_client,
            prompt_version=prompt_version,
        )
        elapsed_ms = (time.perf_counter() - start_time) * 1000
        
        # Add processing time
        if analysis_result.get('status') == 'success':
            if 'llm_metadata' not in analysis_result:
                analysis_result['llm_metadata'] = {}
            analysis_result['llm_metadata']['processing_time_ms'] = elapsed_ms
        
        # Convert to standard format for robotic execution
        standard_result = convert_llm_synthesis_to_standard_format(
            reaction_smiles=reaction,
            synthesis_result=analysis_result,
            requested_type=requested_type,
            processing_time_ms=elapsed_ms,
        )
        
        return analysis_result, standard_result
        
    except Exception as exc:
        error_result = {
            "error": f"LLM synthesis failed: {exc}",
            "meta": {
                "model": f"LLM-Synthesis-{llm_model}",
                "status": "error"
            }
        }
        return error_result, error_result


def main() -> None:
    """Main entry point with optional command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Local Recommendation Tester - Test ChemTools recommendation pipelines",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode (prompts for input):
  python scripts/local_recommendation_cli.py
  
  # Provide reaction and type via command line:
  python scripts/local_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --family Buchwald_CN
  
  # Auto-detect reaction type:
  python scripts/local_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1"
  
  # Specify output directory:
  python scripts/local_recommendation_cli.py --rxn "Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1" --save-dir my_results
        """
    )
    
    parser.add_argument(
        "--rxn", "--reaction",
        type=str,
        default=None,
        help="Reaction SMILES (reactants>>products). If not provided, will prompt interactively."
    )
    
    parser.add_argument(
        "--family", "--type",
        type=str,
        default=None,
        choices=[None, "Suzuki", "Suzuki_CC", "C_N_Coupling_Cu", "Ullmann_CN", 
                 "C_N_Coupling_Pd", "Buchwald_CN", "C_N_Coupling_Ni", "Amide_formation"],
        help="Reaction family/type. If not provided, will prompt interactively or auto-detect."
    )
    
    parser.add_argument(
        "--k",
        type=int,
        default=K_DEFAULT,
        help=f"Number of precedents to retrieve for ML (default: {K_DEFAULT})"
    )
    
    parser.add_argument(
        "--limit",
        type=int,
        default=LIMIT_DEFAULT,
        help=f"Number of ML recommendations to return (default: {LIMIT_DEFAULT})"
    )
    
    parser.add_argument(
        "--fusion-variants",
        type=int,
        default=FUSION_VARIANTS_DEFAULT,
        help=f"Number of fusion recommendation variants (default: {FUSION_VARIANTS_DEFAULT})"
    )
    
    parser.add_argument(
        "--save-dir",
        type=str,
        default="results",
        help="Directory to save output JSON files (default: results)"
    )
    
    parser.add_argument(
        "--strategy",
        type=str,
        default="all",
        choices=["all", "rule", "ml", "fusion", "protocol", "llm"],
        help="Which recommendation strategy to run (default: all). "
             "Use 'llm' for multi-source LLM synthesis. "
             "NOTE: 'fusion' is deprecated, use --rerank rule instead."
    )
    
    parser.add_argument(
        "--rerank",
        type=str,
        default="rule",
        choices=["none", "rule", "analytics"],
        help="Reranking strategy for ML recommendations: 'none' (similarity only), "
             "'rule' (boost by rule match), 'analytics' (boost by popularity). Default: rule"
    )
    
    parser.add_argument(
        "--filter-unknown",
        action="store_true",
        help="Filter out precedents with unknown base/solvent reagents (not in database)"
    )
    
    parser.add_argument(
        "--protocol-tags",
        type=str,
        default=None,
        help="Comma-separated tags to filter protocols (e.g., 'suzuki,palladium')"
    )
    
    # LLM synthesis options
    parser.add_argument(
        "--llm-provider",
        type=str,
        default="aliyun",
        choices=["aliyun", "openai"],
        help="LLM provider for multi-source synthesis (default: aliyun)"
    )
    
    parser.add_argument(
        "--llm-model",
        type=str,
        default="deepseek-v3.2-exp",
        help="LLM model for synthesis (default: deepseek-v3.2-exp)"
    )
    
    parser.add_argument(
        "--llm-prompt-version",
        type=str,
        default="v2",
        choices=["v1", "v2"],
        help="LLM prompt version: v1 (original) or v2 (optimized, 30%% faster). Default: v2"
    )
    
    parser.add_argument(
        "--constraints",
        type=str,
        default=None,
        help="User constraints as JSON string (e.g., '{\"scale\": \"multigram\", \"cost\": \"low\"}')"
    )
    
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print JSON output files"
    )
    
    args = parser.parse_args()
    
    print("Local Recommendation Test")
    print("-------------------------")
    
    # Get reaction SMILES - from args or prompt
    if args.rxn:
        reaction = args.rxn.strip()
        print(f"Reaction SMILES: {reaction}")
    else:
        reaction = prompt_smiles()
    
    # Get reaction type - from args or prompt
    if args.family:
        reaction_type = args.family
        selected_label = args.family
        print(f"Reaction type: {selected_label}")
    else:
        selected_label, reaction_type = choose_reaction_type()
        print(f"Selected reaction type: {selected_label}")

    k_value = args.k
    limit_value = args.limit
    fusion_variants = args.fusion_variants

    db_path = DEFAULT_SCDB_PATH if Path(DEFAULT_SCDB_PATH).exists() else None
    if db_path:
        print(f"Using rule database: {db_path}")
    else:
        print("Rule database: default resolver (environment or ChemTools defaults)")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    type_label = reaction_type or "auto"
    label = slugify_label(type_label)

    print("\nRunning local pipelines...\n")

    # Update output directory from args
    output_dir = Path(args.save_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a custom save function with the specified directory
    def save_to_dir(data: Dict[str, Any], filename: str) -> Path:
        """Save JSON data to the specified output directory."""
        output_path = output_dir / filename
        import json
        output_path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
        return output_path

    # Run selected strategies
    run_rule = args.strategy in ["all", "rule", "llm"]
    run_ml = args.strategy in ["all", "ml", "llm"]
    run_fusion = args.strategy in ["all", "fusion"]
    run_protocol = args.strategy in ["all", "protocol", "llm"]
    run_llm = args.strategy in ["all", "llm"]
    
    rule_result = None
    ml_result = None
    fusion_result = None
    protocol_result = None
    llm_result = None
    
    if run_rule:
        rule_result = local_rule_based_match(reaction, db_path)
        rule_file = save_to_dir(rule_result, f"{timestamp}_{label}_rule_local.json")
    
    if run_ml:
        ml_result = local_ml_recommendation(
            reaction, 
            reaction_type, 
            k_value, 
            limit_value,
            rerank_strategy=args.rerank,
            filter_unknown_reagents=args.filter_unknown
        )
        ml_file = save_to_dir(ml_result, f"{timestamp}_{label}_ml_local.json")
    
    if run_fusion:
        fusion_result = local_fusion_recommendation(reaction, k_value, fusion_variants)
        fusion_file = save_to_dir(fusion_result, f"{timestamp}_{label}_fusion_local.json")
    
    if run_protocol:
        # Parse protocol tags if provided
        protocol_tags = None
        if args.protocol_tags:
            protocol_tags = [tag.strip() for tag in args.protocol_tags.split(',')]
        
        # Protocol recommendation doesn't strictly need reaction_family
        # Only pass it if it might help with filtering
        protocol_result = local_protocol_recommendation(
            reaction=reaction,
            k_value=k_value,
            tags=protocol_tags,
            reaction_family=None  # Let protocol module auto-detect from similarity
        )
        protocol_file = save_to_dir(protocol_result, f"{timestamp}_{label}_protocol_local.json")
    
    # Initialize LLM result variables
    llm_analysis_result = None
    llm_standard_result = None
    llm_analysis_file = None
    llm_standard_file = None
    
    if run_llm:
        # Parse constraints if provided
        constraints = None
        if args.constraints:
            import json
            try:
                constraints = json.loads(args.constraints)
            except json.JSONDecodeError as e:
                print(f"Warning: Failed to parse constraints JSON: {e}")
                constraints = None
        
        # Run LLM synthesis - returns both analysis and standard format
        llm_analysis_result, llm_standard_result = local_llm_synthesis(
            reaction=reaction,
            ml_result=ml_result,
            rule_result=rule_result,
            protocol_result=protocol_result,
            constraints=constraints,
            llm_provider=args.llm_provider,
            llm_model=args.llm_model,
            prompt_version=args.llm_prompt_version,
            requested_type=args.family,
        )
        
        # Save both formats
        llm_analysis_file = save_to_dir(llm_analysis_result, f"{timestamp}_{label}_llm_analysis.json")
        llm_standard_file = save_to_dir(llm_standard_result, f"{timestamp}_{label}_llm_local.json")

    print("Summary\n-------")
    
    if rule_result:
        summarize_rule(rule_result)
        print()
    
    if ml_result:
        summarize_ml(ml_result)
        print()
    
    if fusion_result:
        summarize_fusion(fusion_result)
        print()
    
    if protocol_result:
        summarize_protocol(protocol_result)
        print()
    
    if llm_analysis_result:
        summarize_llm_synthesis(llm_analysis_result)

    print("\nSaved outputs:")
    if run_rule:
        print(f"  Rule JSON:            {rule_file}")
    if run_ml:
        print(f"  ML JSON:              {ml_file}")
    if run_fusion:
        print(f"  Fusion JSON:          {fusion_file}")
    if run_protocol:
        print(f"  Protocol JSON:        {protocol_file}")
    if run_llm:
        print(f"  LLM Analysis JSON:    {llm_analysis_file}")
        print(f"  LLM Standard JSON:    {llm_standard_file}  (for robotic execution)")
    
    print("\nDone.")


if __name__ == "__main__":
    main()
