"""
LangGraph agent for ChemTools featurization, analysis, and reagent registry access.
"""

import json
import os
import re
import ast
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from dotenv import load_dotenv
from langchain_core.messages import AIMessage, BaseMessage, HumanMessage, SystemMessage
from langchain_openai import ChatOpenAI
from pydantic import BaseModel, Field

from .chemtools_wrapper import CHEMTOOLS_TOOLS

try:
    from langgraph.prebuilt import create_agent as _langgraph_agent_factory
    _LANGGRAPH_AGENT_FACTORY_NAME = "create_agent"
except (ImportError, AttributeError):
    try:
        from langgraph.prebuilt import create_react_agent as _langgraph_agent_factory  # type: ignore[attr-defined]
        _LANGGRAPH_AGENT_FACTORY_NAME = "create_react_agent"
    except (ImportError, AttributeError):
        from langchain.agents import create_agent as _langgraph_agent_factory  # type: ignore[attr-defined]
        _LANGGRAPH_AGENT_FACTORY_NAME = "create_react_agent"

LANGGRAPH_AGENT_FACTORY: Callable[..., object] = _langgraph_agent_factory

load_dotenv()

DEFAULT_ALIYUN_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"
_REACTION_SMILES_RE = re.compile(
    r"([A-Za-z0-9@+\-\[\]\(\)=#%/.]+>>[A-Za-z0-9@+\-\[\]\(\)=#%/.]*)"
)
_RUNTIME_POLICY_PREFIX = "Runtime policy (installed by "
_RUNTIME_POLICY_MARKER = (
    "Runtime policy (installed by chem_assistant.chemtools_agent)"
)


def get_llm_client(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0,
) -> ChatOpenAI:
    """
    Get LLM client using configured provider and model.

    Environment variables:
        LLM_PROVIDER: "openai" or "aliyun" (default: "openai")
        LLM_MODEL: Model name (default: "gpt-4o")
        OPENAI_API_KEY: OpenAI API key
        OPENAI_BASE_URL: OpenAI base URL (optional)
        ALIYUN_API_KEY: Aliyun API key
        ALIYUN_BASE_URL: Aliyun base URL (optional)
    """
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "openai")
    if model is None:
        model = os.getenv("LLM_MODEL", "gpt-4o")

    if provider == "aliyun":
        api_key = os.getenv("ALIYUN_API_KEY")
        base_url = os.getenv("ALIYUN_BASE_URL", DEFAULT_ALIYUN_BASE_URL)
        if not api_key:
            raise RuntimeError("ALIYUN_API_KEY environment variable not set")
    elif provider == "openai":
        api_key = os.getenv("OPENAI_API_KEY")
        base_url = os.getenv("OPENAI_BASE_URL", DEFAULT_OPENAI_BASE_URL)
        if not api_key:
            raise RuntimeError("OPENAI_API_KEY environment variable not set")
    else:
        raise ValueError(f"Unsupported provider '{provider}'. Use 'openai' or 'aliyun'.")

    return ChatOpenAI(
        model=model,
        api_key=api_key,
        base_url=base_url,
        temperature=temperature,
    )


CHEMISTRY_SYSTEM_PROMPT = """You are ChemBot, an expert chemistry assistant focused on featurization and analysis.

You have access to the following tools (featurization/analysis/reagent registry):

- analysis_normalize_smiles: normalize SMILES strings
- analysis_normalize_reaction: normalize reaction SMILES
- analysis_analyze_reaction: normalize + taxonomy matches + family detection
- analysis_classify_reactant_smiles: best reactant match
- analysis_classify_reactant_category: reactant category id
- analysis_classify_reactant_group: reactant group label
- analysis_classify_reactant_batch: batch reactant classification
- analysis_get_reactant_category_matches: all reactant categories matched
- analysis_get_all_reactant_matches: all SMARTS matches
- analysis_normalize_reactant_identifier: alias to canonical reactant id
- analysis_normalize_reaction_type: alias to canonical reaction id
- analysis_resolve_reaction_family: alias to canonical family id
- analysis_canonical_family_label: canonical family resolution helper
- analysis_classify_reactants_with_context: context-aware reactant roles
- analysis_reactant_summary: summary of context-aware reactant roles

- detection_detect_reaction_types: reaction detection from reaction SMILES
- detection_detect_reaction_types_from_smiles: reaction detection from reactant/product lists
- detection_detect_motif_ids_from_smiles: motif id detection from reactants

- unified_featurize_molecule: unified molecule bundle (motifs, nearby groups, RDKit props)
- unified_featurize_reaction: unified reaction bundle (detection, roles, aggregates)
- motif_featurize_molecule: motif-based steric/electronic analysis
- motif_featurize_reaction: motif-based reaction featurization
- reaction_pair_featurize_pair: electrophile/nucleophile pair features
- reaction_pair_featurize_flat: flat pair features
- calculable_detect_all_features: calculable feature flags
- calculable_get_reactant_type_features: reactant type features
- calculable_classify_reactant_smiles: calculable-based reactant classification
- hte_recommend_conditions: HTE-based condition recommendations
- hte_dataset_summary: dataset-level condition summaries (filtered by reaction/ reactant types)
- hte_database_stats: HTE database summary statistics
- reagent_lookup: lookup reagent by name, abbreviation, or CAS
- reagent_list_roles: list available reagent roles (with counts)
- reagent_list_by_role: list reagents for a role (supports limit)
- reagent_list_by_family: list reagents for a role/family
- rag_search: retrieve curated knowledge base snippets (RAG)
- kb_recommend_conditions: extract summary condition tables from the knowledge base

Workflow (stepwise, always follow):
1) Identify the task type and required inputs.
2) Call the minimum necessary tools. Do not guess tool outputs.
3) Extract evidence from tool outputs (fields + metrics + provenance).
4) Respond with evidence first, then interpretations tied to the evidence.
5) If inputs are missing, ask for them explicitly and stop.

Tool selection rubric:
- Molecule question -> unified_featurize_molecule first. Add motif_featurize_molecule or calculable_* only if asked.
- Reaction question -> unified_featurize_reaction first.
  - Also run analysis_analyze_reaction for taxonomy context.
  - If the user wants reactant roles, add analysis_classify_reactants_with_context or analysis_reactant_summary.
- Electrophile/nucleophile pair -> reaction_pair_featurize_pair (or reaction_pair_featurize_flat for compact output).
- Only use detection_* tools when the user asks for reaction typing without full featurization.
- HTE data or condition screening -> hte_recommend_conditions (use reaction_smiles when available).
- HTE dataset insights (e.g., top conditions for aryl chlorides) -> hte_dataset_summary.
- If the user wants literature-only HTE results, set source_group="literature" in hte_recommend_conditions.
- HTE database questions -> hte_database_stats.
- Reagent lookup or validation -> reagent_lookup (optionally set role or include_all).
- Reagent inventory -> reagent_list_roles or reagent_list_by_role.
- Reagent family queries -> reagent_list_by_family.
- If the user asks for specific protocols, rules, or literature-style guidance, call rag_search and cite the snippets.
- If the user asks for condition summaries or top conditions, call kb_recommend_conditions.

HTE explanation guidance:
- If asked to explain why conditions were chosen, cite the HTE metrics (avg_z_score, match_score, num_experiments).
- If no recommendations were returned, explain that no HTE matches exist and cite diagnostics
  (reactant types/categories, matched motifs, database coverage) plus next steps.

Consistency checks for reactions:
- Compare unified_featurize_reaction.reaction.reaction_type with analysis_analyze_reaction.family.canonical_id.
- If they disagree or confidence is low, report both and state uncertainty explicitly.

Response templates:
- Molecule:
  - Evidence: list tool fields + values used (cite tool name)
  - Input: <smiles>
  - Highlights: motifs, nearby groups, key RDKit props
  - Workflow: render workflow.steps from the tool output in order with data (no paraphrasing)
  - Notes: any errors or assumptions
- Reaction:
  - Evidence: list tool fields + values used (cite tool name)
  - Input: <reaction_smiles>
  - Type: unified reaction_type + analysis family (and confidence)
  - Reactants: key roles or categories
  - Aggregates: max steric/electronic summary
  - Notes: disagreements or edge cases
- Pair:
  - Evidence: list tool fields + values used (cite tool name)
  - Input: electrophile + nucleophile
  - Pair features: LG, nuc_class, sterics
  - Flags: any key tags or calculable signals
- Reagent:
  - Evidence: list tool fields + values used (cite tool name)
  - Query: <name or CAS>
  - Result: matched role, name, CAS, and key fields
  - Notes: multiple matches or missing entries

Condition recommendations (HTE) output should be structured JSON by default.
If the user explicitly asks for a human-readable summary, return a readable summary instead of JSON.
Use this Pydantic response format and populate only from tool outputs when returning JSON:

class ConditionEvidence(BaseModel):
    source: str  # "hte_recommend_conditions"
    reaction_smiles: Optional[str]
    reactant_a_smiles: Optional[str]
    reactant_b_smiles: Optional[str]
    product_smiles: Optional[str]
    predicted_reaction_type: Optional[str]
    reaction_type_confidence: float
    total_matching_experiments: int
    database_coverage: float
    is_fallback_match: bool
    matched_motifs: Optional[List[str]]

class ConditionRecommendationItem(BaseModel):
    rank: int
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str]
    additive: Optional[str]
    coupling_reagent: Optional[str]
    success_rate: float
    avg_yield: float
    median_yield: float
    num_experiments: int
    avg_z_score: float
    confidence_score: float
    match_score: float
    reaction_type: Optional[str]
    reaction_category: Optional[str]
    reaction_id: Optional[str]
    reactant_types: List[str]
    z_score_range: List[float]

class KBConditionItem(BaseModel):
    rank: int
    context: str
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str]
    additive: Optional[str]
    coupling_reagent: Optional[str]
    temperature_c: Optional[str]
    time_h: Optional[str]
    notes: Optional[str]
    reaction_id: Optional[str]
    reaction_type: Optional[str]
    score: float
    source_doc: Optional[str]
    source_path: Optional[str]
    tags: List[str]
    extras: Dict[str, str]

class ConditionRecommendationResponse(BaseModel):
    query: str
    evidence: ConditionEvidence
    recommendations: List[ConditionRecommendationItem]
    knowledge_base: List["KBConditionItem"]
    warnings: List[str]

Guidelines:
- Use unified_featurize_molecule/reaction as the primary entry points.
- Use analysis tools for normalization and taxonomy-level reasoning.
- Use reaction_pair tools when the user provides electrophile/nucleophile pairs.
- Keep answers concise and focused on the analysis output.
- When workflow.steps is available, output it with the actual data (or a clearly labeled, truncated subset). Do not narrate or invent details.
"""


class ConditionEvidence(BaseModel):
    """Structured evidence for HTE condition recommendations."""

    source: str = Field(default="hte_recommend_conditions")
    reaction_smiles: Optional[str] = None
    reactant_a_smiles: Optional[str] = None
    reactant_b_smiles: Optional[str] = None
    product_smiles: Optional[str] = None
    predicted_reaction_type: Optional[str] = None
    reaction_type_confidence: float = 0.0
    total_matching_experiments: int = 0
    database_coverage: float = 0.0
    is_fallback_match: bool = False
    matched_motifs: Optional[List[str]] = None


class ConditionRecommendationItem(BaseModel):
    """Structured single-condition recommendation."""

    rank: int
    catalyst: str
    ligand: str
    base: str
    solvent: str
    secondary_solvent: Optional[str] = None
    additive: Optional[str] = None
    coupling_reagent: Optional[str] = None
    success_rate: float = 0.0
    avg_yield: float = 0.0
    median_yield: float = 0.0
    num_experiments: int = 0
    avg_z_score: float = 0.0
    confidence_score: float = 0.0
    match_score: float = 0.0
    reaction_type: Optional[str] = None
    reaction_category: Optional[str] = None
    reaction_id: Optional[str] = None
    reactant_types: List[str] = Field(default_factory=list)
    z_score_range: List[float] = Field(default_factory=list)


class KBConditionItem(BaseModel):
    """Knowledge base condition summary row."""

    rank: int
    context: str = ""
    catalyst: str = ""
    ligand: str = ""
    base: str = ""
    solvent: str = ""
    secondary_solvent: Optional[str] = None
    additive: Optional[str] = None
    coupling_reagent: Optional[str] = None
    temperature_c: Optional[str] = None
    time_h: Optional[str] = None
    notes: Optional[str] = None
    reaction_id: Optional[str] = None
    reaction_type: Optional[str] = None
    score: float = 0.0
    source_doc: Optional[str] = None
    source_path: Optional[str] = None
    tags: List[str] = Field(default_factory=list)
    extras: Dict[str, str] = Field(default_factory=dict)


class ConditionRecommendationResponse(BaseModel):
    """Structured response for HTE condition recommendations."""

    query: str
    evidence: ConditionEvidence
    recommendations: List[ConditionRecommendationItem] = Field(default_factory=list)
    knowledge_base: List[KBConditionItem] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)


def _safe_json_loads(value: str) -> Optional[Any]:
    if not isinstance(value, str):
        return None
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return None


def _coerce_mapping(value: Any) -> Optional[Dict[str, Any]]:
    if isinstance(value, dict):
        return value
    if not isinstance(value, str):
        return None
    text = value.strip()
    if not text:
        return None

    parsed = _safe_json_loads(text)
    if isinstance(parsed, dict):
        return parsed

    try:
        literal = ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return None
    return literal if isinstance(literal, dict) else None


def _coerce_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


def _wants_human_readable(query: str) -> bool:
    text = (query or "").lower()
    triggers = [
        "easier to read",
        "readable",
        "human",
        "explain",
        "why",
        "summary",
        "format",
        "table",
    ]
    return any(trigger in text for trigger in triggers)


def _wants_json(query: str) -> bool:
    text = (query or "").lower()
    return "json" in text or "raw" in text


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _registered_tool_names() -> List[str]:
    names = {
        str(getattr(tool, "name", "")).strip()
        for tool in CHEMTOOLS_TOOLS
        if str(getattr(tool, "name", "")).strip()
    }
    return sorted(names)


@lru_cache(maxsize=1)
def _runtime_capabilities() -> Dict[str, bool]:
    root = Path(__file__).resolve().parents[1]
    hte_ready = (root / "data" / "HTE_db").exists()
    kb_ready = (root / "data" / "knowledge_base").exists()
    taxonomy_ready = False
    try:
        from chemtools.taxonomy.reaction_catalog import load_reaction_catalog

        rxn_defs, _ = load_reaction_catalog()
        taxonomy_ready = bool(rxn_defs)
    except Exception:
        taxonomy_ready = False
    return {
        "hte_ready": hte_ready,
        "kb_ready": kb_ready,
        "taxonomy_ready": taxonomy_ready,
    }


def _append_runtime_tool_policy(system_prompt: str) -> str:
    prompt = str(system_prompt or "").strip()
    if not prompt:
        return prompt
    if _RUNTIME_POLICY_PREFIX in prompt:
        return prompt

    tool_names = _registered_tool_names()
    capabilities = _runtime_capabilities()
    has_hte_tool = "hte_recommend_conditions" in tool_names
    has_kb_tools = any(
        name in tool_names for name in ("rag_search", "kb_recommend_conditions")
    )

    policy_lines = [
        _RUNTIME_POLICY_MARKER,
        f"- Registered tools: {len(tool_names)}",
        "- Always call the minimum number of tools needed for the user request.",
        "- Start with one primary tool, then add follow-up tools only when evidence is missing.",
        "- If required inputs are missing, ask for them explicitly instead of guessing.",
    ]
    if has_hte_tool and capabilities["hte_ready"]:
        policy_lines.append(
            "- HTE data is available; for condition screening prefer hte_recommend_conditions first."
        )
    elif has_hte_tool:
        policy_lines.append(
            "- hte_recommend_conditions exists but HTE data is unavailable; require explicit hte_db_path before use."
        )
    else:
        policy_lines.append(
            "- HTE tool is unavailable; avoid promising condition recommendations."
        )
    if has_kb_tools and capabilities["kb_ready"]:
        policy_lines.append(
            "- Knowledge base is available; use rag_search/kb_recommend_conditions only for protocol, rules, or literature requests."
        )
    elif has_kb_tools:
        policy_lines.append(
            "- KB tools exist but knowledge_base data is unavailable; avoid KB retrieval until a valid root is provided."
        )
    else:
        policy_lines.append(
            "- KB retrieval tools are unavailable; avoid KB-specific claims."
        )
    if capabilities["taxonomy_ready"]:
        policy_lines.append(
            "- Taxonomy is available; prefer taxonomy-consistent reaction naming and filter values."
        )
    else:
        policy_lines.append(
            "- Taxonomy appears unavailable; report uncertainty for reaction-family conclusions."
        )

    return f"{prompt}\n\n" + "\n".join(policy_lines)


def _extract_tool_payload(
    messages: List[BaseMessage],
    tool_name: str,
) -> Optional[Dict[str, Any]]:
    call_map: Dict[str, str] = {}
    for message in messages:
        if isinstance(message, AIMessage):
            for call in getattr(message, "tool_calls", []) or []:
                call_payload = _coerce_mapping(call)
                if not call_payload:
                    continue
                call_id = call_payload.get("id") or call_payload.get("tool_call_id")
                name = call_payload.get("name")
                if call_id and name:
                    call_map[str(call_id)] = name

    for message in reversed(messages):
        tool_call_id = getattr(message, "tool_call_id", None)
        if not tool_call_id:
            continue
        if call_map.get(str(tool_call_id)) != tool_name:
            continue
        artifact = getattr(message, "artifact", None)
        artifact_payload = _coerce_mapping(artifact)
        if artifact_payload is not None:
            return artifact_payload
        content = getattr(message, "content", None)
        content_payload = _coerce_mapping(content)
        if content_payload is not None:
            return content_payload
        if isinstance(content, list):
            if len(content) == 1 and isinstance(content[0], dict):
                return content[0]
            joined = "".join(str(item) for item in content)
            payload = _coerce_mapping(joined)
            if payload is not None:
                return payload
        if isinstance(content, str):
            payload = _coerce_mapping(content)
            if payload is not None:
                return payload
    return None


def _build_condition_response(
    query: str,
    tool_payload: Dict[str, Any],
    kb_payload: Optional[Dict[str, Any]] = None,
) -> ConditionRecommendationResponse:
    tool_payload = tool_payload if isinstance(tool_payload, dict) else {}

    warnings: List[str] = []
    success = bool(tool_payload.get("success", True))
    if not success:
        error = tool_payload.get("error") or {}
        detail = error.get("details") if isinstance(error, dict) else None
        warnings.append("hte_recommend_conditions failed.")
        if detail:
            warnings.append(str(detail))
        elif isinstance(error, str) and error.strip():
            warnings.append(error.strip())

    input_block_raw = tool_payload.get("input")
    diagnostics_raw = tool_payload.get("diagnostics")
    input_block = input_block_raw if isinstance(input_block_raw, dict) else {}
    diagnostics = diagnostics_raw if isinstance(diagnostics_raw, dict) else {}
    evidence = ConditionEvidence(
        reaction_smiles=input_block.get("reaction_smiles"),
        reactant_a_smiles=input_block.get("reactant_a_smiles")
        or tool_payload.get("reactant_a_smiles"),
        reactant_b_smiles=input_block.get("reactant_b_smiles")
        or tool_payload.get("reactant_b_smiles"),
        product_smiles=input_block.get("product_smiles")
        or tool_payload.get("product_smiles"),
        predicted_reaction_type=tool_payload.get("predicted_reaction_type")
        or diagnostics.get("predicted_reaction_type"),
        reaction_type_confidence=_to_float(
            tool_payload.get("reaction_type_confidence")
            or diagnostics.get("reaction_type_confidence"),
            0.0,
        ),
        total_matching_experiments=_to_int(
            tool_payload.get("total_matching_experiments")
            or diagnostics.get("total_matching_experiments"),
            0,
        ),
        database_coverage=_to_float(
            tool_payload.get("database_coverage") or diagnostics.get("database_coverage"),
            0.0,
        ),
        is_fallback_match=bool(
            tool_payload.get("is_fallback_match", diagnostics.get("is_fallback_match", False))
        ),
        matched_motifs=[str(v) for v in _coerce_list(tool_payload.get("matched_motifs"))]
        or None,
    )

    if evidence.is_fallback_match:
        warnings.append("Fallback match used; similarity may be weak.")

    recommendations: List[ConditionRecommendationItem] = []
    for idx, rec in enumerate(_coerce_list(tool_payload.get("recommendations")), start=1):
        if not isinstance(rec, dict):
            continue
        recommendations.append(
            ConditionRecommendationItem(
                rank=idx,
                catalyst=str(rec.get("catalyst") or ""),
                ligand=str(rec.get("ligand") or ""),
                base=str(rec.get("base") or ""),
                solvent=str(rec.get("solvent") or ""),
                secondary_solvent=rec.get("secondary_solvent") or None,
                additive=rec.get("additive") or None,
                coupling_reagent=rec.get("coupling_reagent") or None,
                success_rate=_to_float(rec.get("success_rate"), 0.0),
                avg_yield=_to_float(rec.get("avg_yield"), 0.0),
                median_yield=_to_float(rec.get("median_yield"), 0.0),
                num_experiments=_to_int(rec.get("num_experiments"), 0),
                avg_z_score=_to_float(rec.get("avg_z_score"), 0.0),
                confidence_score=_to_float(rec.get("confidence_score"), 0.0),
                match_score=_to_float(rec.get("match_score"), 0.0),
                reaction_type=rec.get("reaction_type"),
                reaction_category=rec.get("reaction_category"),
                reaction_id=rec.get("reaction_id"),
                reactant_types=[str(v) for v in _coerce_list(rec.get("reactant_types"))],
                z_score_range=[
                    _to_float(v, 0.0) for v in _coerce_list(rec.get("z_score_range"))
                ],
            )
        )

    if not recommendations:
        warnings.append("No condition recommendations returned.")
        reactant_a_type = diagnostics.get("reactant_a_type") or ""
        reactant_b_type = diagnostics.get("reactant_b_type") or ""
        reactant_a_category = diagnostics.get("reactant_a_category") or ""
        reactant_b_category = diagnostics.get("reactant_b_category") or ""
        if reactant_a_type or reactant_b_type:
            warnings.append(
                "No HTE matches for reactant types: "
                f"{reactant_a_type or 'unknown'} / {reactant_b_type or 'unknown'}."
            )
        elif reactant_a_category or reactant_b_category:
            warnings.append(
                "No HTE matches for reactant categories: "
                f"{reactant_a_category or 'unknown'} / {reactant_b_category or 'unknown'}."
            )
        else:
            warnings.append(
                "Reactant motifs could not be detected; verify SMILES or provide neutral forms."
            )

    knowledge_base: List[KBConditionItem] = []
    if kb_payload:
        kb_payload_map = (
            kb_payload
            if isinstance(kb_payload, dict)
            else _coerce_mapping(kb_payload)
        ) or {}
        kb_success = bool(kb_payload_map.get("success", True))
        if not kb_success:
            warnings.append("kb_recommend_conditions failed.")
        for item in _coerce_list(kb_payload_map.get("results")):
            if not isinstance(item, dict):
                continue
            condition_raw = item.get("condition")
            extras_raw = item.get("extras")
            condition = condition_raw if isinstance(condition_raw, dict) else {}
            extras = extras_raw if isinstance(extras_raw, dict) else {}
            knowledge_base.append(
                KBConditionItem(
                    rank=_to_int(item.get("rank"), len(knowledge_base) + 1),
                    context=str(item.get("context") or ""),
                    catalyst=str(condition.get("catalyst") or ""),
                    ligand=str(condition.get("ligand") or ""),
                    base=str(condition.get("base") or ""),
                    solvent=str(condition.get("solvent") or ""),
                    secondary_solvent=condition.get("secondary_solvent") or None,
                    additive=condition.get("additive") or None,
                    coupling_reagent=condition.get("coupling_reagent") or None,
                    temperature_c=condition.get("temperature_c") or None,
                    time_h=condition.get("time_h") or None,
                    notes=condition.get("notes") or None,
                    reaction_id=condition.get("reaction_id") or None,
                    reaction_type=condition.get("reaction_type") or None,
                    score=_to_float(item.get("score"), 0.0),
                    source_doc=str(item.get("doc_id") or ""),
                    source_path=str(item.get("path") or ""),
                    tags=[str(tag) for tag in _coerce_list(item.get("tags"))],
                    extras={str(k): str(v) for k, v in extras.items()},
                )
            )

    return ConditionRecommendationResponse(
        query=query,
        evidence=evidence,
        recommendations=recommendations,
        knowledge_base=knowledge_base,
        warnings=warnings,
    )


def _split_condition_items(value: Optional[str]) -> List[str]:
    if not value:
        return []
    parts = re.split(r"[,/;]+", value)
    return [part.strip() for part in parts if part and part.strip()]


def _format_reagent_detail(name: str, role: str) -> str:
    from chemtools import reagent as reagent_tools

    if not name:
        return "not specified"

    details: List[str] = []
    for item in _split_condition_items(name):
        info = reagent_tools.enrich_reagent_info(item, role)
        if not info.get("found"):
            details.append(f"{item} (not found)")
            continue
        role_payload = (info.get("roles") or {}).get(role, {})
        families = role_payload.get("families") or []
        family_id = families[0] if families else None
        tag = role_payload.get("tag") or None
        pieces = [info.get("name") or item]
        cas = info.get("cas")
        if cas:
            pieces.append(f"CAS {cas}")
        if family_id:
            pieces.append(f"family {family_id}")
        if tag:
            pieces.append(f"tag {tag}")
        details.append(" | ".join(pieces))
    return "; ".join(details) if details else "not specified"


def _format_condition_summary(response: ConditionRecommendationResponse) -> str:
    evidence = response.evidence
    lines = [
        "**HTE Conditions**",
        f"Reaction: {evidence.reaction_smiles or 'unknown'}",
        "Evidence:",
        f"- predicted reaction type: {evidence.predicted_reaction_type or 'unknown'} (confidence {evidence.reaction_type_confidence:.2f})",
        f"- matched motifs: {', '.join(evidence.matched_motifs or []) or 'none'}",
        f"- total matching experiments: {evidence.total_matching_experiments}",
        f"- database coverage: {evidence.database_coverage:.2f}%",
        f"- fallback match: {'yes' if evidence.is_fallback_match else 'no'}",
    ]
    if response.warnings:
        lines.append("Warnings:")
        for warning in response.warnings:
            lines.append(f"- {warning}")

    if not response.recommendations:
        lines.append("No condition recommendations available.")
        if not response.knowledge_base:
            return "\n".join(lines)
        lines.append("")
        lines.append("**Knowledge Base Summaries**")
        for item in response.knowledge_base:
            lines.append(f"- {item.context or 'summary'}")
            lines.append(
                "  "
                f"Catalyst={item.catalyst or 'n/a'}; "
                f"Ligand={item.ligand or 'n/a'}; "
                f"Base={item.base or 'n/a'}; "
                f"Solvent={item.solvent or 'n/a'}"
            )
            if item.notes:
                lines.append(f"  Notes: {item.notes}")
        return "\n".join(lines)

    lines.append("Conditions:")
    for rec in response.recommendations:
        lines.append(f"Condition {rec.rank}:")
        lines.append(f"Catalyst: {rec.catalyst or 'not specified'}")
        lines.append(f"Ligand: {rec.ligand or 'not specified'}")
        lines.append(f"Base: {rec.base or 'not specified'}")
        lines.append(f"Solvent: {rec.solvent or 'not specified'}")
        lines.append(f"Secondary solvent: {rec.secondary_solvent or 'not specified'}")
        lines.append(f"Additive: {rec.additive or 'not specified'}")
        lines.append(f"Coupling reagent: {rec.coupling_reagent or 'not specified'}")
        lines.append(
            "Metrics: "
            f"avg_z_score={rec.avg_z_score:.2f}, "
            f"success_rate={rec.success_rate:.1f}%, "
            f"avg_yield={rec.avg_yield:.1f}%, "
            f"median_yield={rec.median_yield:.1f}%, "
            f"num_experiments={rec.num_experiments}, "
            f"match_score={rec.match_score:.2f}, "
            f"confidence_score={rec.confidence_score:.2f}"
        )
        lines.append("Reagent lookup:")
        lines.append(f"- catalyst: {_format_reagent_detail(rec.catalyst, 'metal_catalyst')}")
        lines.append(f"- ligand: {_format_reagent_detail(rec.ligand, 'ligand')}")
        lines.append(f"- base: {_format_reagent_detail(rec.base, 'base')}")
        lines.append(f"- solvent: {_format_reagent_detail(rec.solvent, 'solvent')}")
        lines.append(f"- additive: {_format_reagent_detail(rec.additive or '', 'additive')}")
        lines.append(f"- coupling reagent: {_format_reagent_detail(rec.coupling_reagent or '', 'other_reagent')}")

    if response.knowledge_base:
        lines.append("")
        lines.append("**Knowledge Base Summaries**")
        for item in response.knowledge_base:
            lines.append(f"- {item.context or 'summary'}")
            lines.append(
                "  "
                f"Catalyst={item.catalyst or 'n/a'}; "
                f"Ligand={item.ligand or 'n/a'}; "
                f"Base={item.base or 'n/a'}; "
                f"Solvent={item.solvent or 'n/a'}"
            )
            if item.temperature_c or item.time_h:
                lines.append(
                    "  "
                    f"Temp={item.temperature_c or 'n/a'}; "
                    f"Time={item.time_h or 'n/a'}"
                )
            if item.notes:
                lines.append(f"  Notes: {item.notes}")

    return "\n".join(lines)


class ChemToolsAgent:
    """LangGraph agent for chemistry analysis using ChemTools featurizers."""

    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0,
        system_prompt: Optional[str] = None,
        verbose: bool = False,
        auto_check: bool = True,
        runtime_tool_policy: bool = True,
    ):
        self.llm = get_llm_client(provider, model, temperature)
        base_prompt = system_prompt or CHEMISTRY_SYSTEM_PROMPT
        self.system_prompt = (
            _append_runtime_tool_policy(base_prompt)
            if runtime_tool_policy
            else base_prompt
        )
        self.verbose = verbose
        self.auto_check = auto_check
        self.runtime_tool_policy = runtime_tool_policy
        self.agent_factory_name = _LANGGRAPH_AGENT_FACTORY_NAME
        self.tools = CHEMTOOLS_TOOLS
        self.agent = LANGGRAPH_AGENT_FACTORY(
            self.llm,
            self.tools,
            prompt=self.system_prompt,
        )

    def run(
        self,
        query: str,
        history: Optional[List[BaseMessage]] = None,
        recursion_limit: int = 15,
    ) -> str:
        try:
            messages = list(history or [])
            preflight = self._maybe_run_preflight(query)
            if preflight:
                messages.append(SystemMessage(content=preflight))
            messages.append(HumanMessage(content=query))

            result = self.agent.invoke(
                {"messages": messages},
                config={"recursion_limit": recursion_limit},
            )
            result_payload = (
                result
                if isinstance(result, dict)
                else _coerce_mapping(result)
            ) or {}
            result_messages = result_payload.get("messages")
            if not isinstance(result_messages, list):
                result_messages = []
            extracted_payload = _extract_tool_payload(
                result_messages, "hte_recommend_conditions"
            )
            structured_payload = (
                extracted_payload
                if isinstance(extracted_payload, dict)
                else _coerce_mapping(extracted_payload)
            )
            if structured_payload is not None:
                kb_payload = None
                try:
                    from chem_assistant.chemtools_wrapper import kb_recommend_conditions

                    diagnostics_raw = structured_payload.get("diagnostics")
                    diagnostics = (
                        diagnostics_raw
                        if isinstance(diagnostics_raw, dict)
                        else {}
                    )
                    predicted = (
                        structured_payload.get("predicted_reaction_type")
                        or diagnostics.get("predicted_reaction_type")
                    )
                    kb_query = query
                    predicted_text = str(predicted).strip() if predicted is not None else ""
                    if predicted_text and predicted_text.lower() not in (query or "").lower():
                        kb_query = f"{query} {predicted_text}"
                    kb_payload = kb_recommend_conditions.invoke({"query": kb_query})
                except Exception:
                    kb_payload = None
                structured_response = _build_condition_response(
                    query, structured_payload, kb_payload
                )
                if _wants_human_readable(query) and not _wants_json(query):
                    return _format_condition_summary(structured_response)
                return structured_response.model_dump_json(indent=2)
            if not result_messages:
                return str(result)
            final_message = result_messages[-1]
            if isinstance(final_message, AIMessage):
                return final_message.content
            final_content = getattr(final_message, "content", final_message)
            return str(final_content)
        except Exception as exc:
            if self.verbose:
                print(f"Agent error: {exc}")
            return f"Error: {exc}\n\nPlease rephrase your question or provide more details."

    @staticmethod
    def _extract_reaction_smiles(text: str) -> Optional[str]:
        match = _REACTION_SMILES_RE.search(text or "")
        if match:
            return match.group(1).strip()
        return None

    def _maybe_run_preflight(self, query: str) -> Optional[str]:
        if not self.auto_check:
            return None
        reaction_smiles = self._extract_reaction_smiles(query)
        if not reaction_smiles:
            return None

        from chem_assistant.chemtools_wrapper import (
            analysis_analyze_reaction,
            unified_featurize_reaction,
        )

        unified_payload = None
        analysis_payload = None
        unified_error = None
        analysis_error = None

        try:
            unified_result_raw = unified_featurize_reaction.invoke(
                {"reaction_smiles": reaction_smiles}
            )
            unified_result = (
                unified_result_raw
                if isinstance(unified_result_raw, dict)
                else _coerce_mapping(unified_result_raw)
            ) or {}
            if unified_result.get("success"):
                unified_payload = unified_result
            else:
                unified_error = unified_result.get("error") or (
                    f"unexpected payload: {type(unified_result_raw).__name__}"
                )
        except Exception as exc:
            unified_error = str(exc)

        try:
            analysis_result_raw = analysis_analyze_reaction.invoke(
                {"reaction_smiles": reaction_smiles}
            )
            analysis_result = (
                analysis_result_raw
                if isinstance(analysis_result_raw, dict)
                else _coerce_mapping(analysis_result_raw)
            ) or {}
            if analysis_result.get("success"):
                analysis_payload = analysis_result
            else:
                analysis_error = analysis_result.get("error") or (
                    f"unexpected payload: {type(analysis_result_raw).__name__}"
                )
        except Exception as exc:
            analysis_error = str(exc)

        unified_type = None
        unified_conf = None
        if unified_payload:
            reaction_raw = unified_payload.get("reaction")
            reaction = reaction_raw if isinstance(reaction_raw, dict) else {}
            reaction_type_raw = reaction.get("reaction_type")
            reaction_type = (
                reaction_type_raw
                if isinstance(reaction_type_raw, dict)
                else {}
            )
            unified_type = reaction_type.get("reaction_type") or reaction_type.get("name")
            unified_conf = reaction_type.get("confidence")

        analysis_family = None
        if analysis_payload:
            family_raw = analysis_payload.get("family")
            family = family_raw if isinstance(family_raw, dict) else {}
            analysis_family = family.get("canonical_id")

        agree = (
            unified_type is not None
            and analysis_family is not None
            and unified_type == analysis_family
        )

        lines = [
            "Preflight reaction cross-check (auto):",
            f"- reaction_smiles: {reaction_smiles}",
            f"- unified.reaction_type: {unified_type or 'unknown'}"
            + (f" (confidence {unified_conf})" if unified_conf is not None else ""),
            f"- analysis.family.canonical_id: {analysis_family or 'unknown'}",
            f"- agree: {'yes' if agree else 'no'}",
        ]
        if unified_error:
            lines.append(f"- unified error: {unified_error}")
        if analysis_error:
            lines.append(f"- analysis error: {analysis_error}")

        return "\n".join(lines)

    def chat(
        self,
        query: str,
        history: List[BaseMessage],
        recursion_limit: int = 15,
    ) -> tuple[str, List[BaseMessage]]:
        response = self.run(query, history, recursion_limit)
        updated_history = history + [HumanMessage(content=query), AIMessage(content=response)]
        return response, updated_history


def create_agent(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0,
    **kwargs: Any,
) -> ChemToolsAgent:
    """Convenience function to create a ChemToolsAgent."""
    return ChemToolsAgent(
        provider=provider,
        model=model,
        temperature=temperature,
        **kwargs,
    )


def quick_query(
    query: str,
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0,
    **kwargs: Any,
) -> str:
    """Quick one-shot query without maintaining state."""
    agent = ChemToolsAgent(
        provider=provider,
        model=model,
        temperature=temperature,
        **kwargs,
    )
    return agent.run(query)


if __name__ == "__main__":
    print("=" * 70)
    print("ChemTools Agent Test")
    print("=" * 70)
    try:
        agent = ChemToolsAgent(verbose=True)
        test_queries = [
            "Normalize this SMILES: c1ccccc1",
            "Analyze this reaction: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "Featurize molecule: c1ccccc1O",
        ]
        for i, query in enumerate(test_queries, 1):
            print(f"\n{'=' * 70}")
            print(f"Test {i}: {query}")
            print("=" * 70)
            response = agent.run(query)
            print(response)
        print(f"\n{'=' * 70}")
        print("All tests completed")
        print("=" * 70)
    except Exception as exc:
        print(f"Error: {exc}")
        print("Make sure your API keys are set:")
        print("  - OPENAI_API_KEY or ALIYUN_API_KEY")
        print("  - Optionally set LLM_PROVIDER and LLM_MODEL")
        raise SystemExit(1)
