"""
LangGraph agent for ChemTools featurization, analysis, and reagent registry access.
"""

import json
import os
import re
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
- hte_database_stats: HTE database summary statistics
- reagent_lookup: lookup reagent by name, abbreviation, or CAS
- reagent_list_roles: list available reagent roles (with counts)
- reagent_list_by_role: list reagents for a role (supports limit)
- reagent_list_by_family: list reagents for a role/family
- rag_search: retrieve curated knowledge base snippets (RAG)

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
- HTE database questions -> hte_database_stats.
- Reagent lookup or validation -> reagent_lookup (optionally set role or include_all).
- Reagent inventory -> reagent_list_roles or reagent_list_by_role.
- Reagent family queries -> reagent_list_by_family.
- If the user asks for specific protocols, rules, or literature-style guidance, call rag_search and cite the snippets.

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

Condition recommendations (HTE) output MUST be structured JSON only (no extra text).
Use this Pydantic response format and populate only from tool outputs:

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

class ConditionRecommendationResponse(BaseModel):
    query: str
    evidence: ConditionEvidence
    recommendations: List[ConditionRecommendationItem]
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


class ConditionRecommendationResponse(BaseModel):
    """Structured response for HTE condition recommendations."""

    query: str
    evidence: ConditionEvidence
    recommendations: List[ConditionRecommendationItem] = Field(default_factory=list)
    warnings: List[str] = Field(default_factory=list)


def _safe_json_loads(value: str) -> Optional[Dict[str, Any]]:
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return None


def _coerce_list(value: Any) -> List[Any]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]


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


def _extract_tool_payload(
    messages: List[BaseMessage],
    tool_name: str,
) -> Optional[Dict[str, Any]]:
    call_map: Dict[str, str] = {}
    for message in messages:
        if isinstance(message, AIMessage):
            for call in getattr(message, "tool_calls", []) or []:
                call_id = call.get("id") or call.get("tool_call_id")
                name = call.get("name")
                if call_id and name:
                    call_map[str(call_id)] = name

    for message in reversed(messages):
        tool_call_id = getattr(message, "tool_call_id", None)
        if not tool_call_id:
            continue
        if call_map.get(str(tool_call_id)) != tool_name:
            continue
        artifact = getattr(message, "artifact", None)
        if isinstance(artifact, dict):
            return artifact
        content = getattr(message, "content", None)
        if isinstance(content, dict):
            return content
        if isinstance(content, list):
            if len(content) == 1 and isinstance(content[0], dict):
                return content[0]
            joined = "".join(str(item) for item in content)
            payload = _safe_json_loads(joined)
            if payload is not None:
                return payload
        if isinstance(content, str):
            payload = _safe_json_loads(content)
            if payload is not None:
                return payload
    return None


def _build_condition_response(
    query: str,
    tool_payload: Dict[str, Any],
) -> ConditionRecommendationResponse:
    warnings: List[str] = []
    success = bool(tool_payload.get("success", True))
    if not success:
        error = tool_payload.get("error") or {}
        detail = error.get("details") if isinstance(error, dict) else None
        warnings.append("hte_recommend_conditions failed.")
        if detail:
            warnings.append(str(detail))

    input_block = tool_payload.get("input") or {}
    evidence = ConditionEvidence(
        reaction_smiles=input_block.get("reaction_smiles"),
        reactant_a_smiles=input_block.get("reactant_a_smiles")
        or tool_payload.get("reactant_a_smiles"),
        reactant_b_smiles=input_block.get("reactant_b_smiles")
        or tool_payload.get("reactant_b_smiles"),
        product_smiles=input_block.get("product_smiles")
        or tool_payload.get("product_smiles"),
        predicted_reaction_type=tool_payload.get("predicted_reaction_type"),
        reaction_type_confidence=_to_float(
            tool_payload.get("reaction_type_confidence"), 0.0
        ),
        total_matching_experiments=_to_int(
            tool_payload.get("total_matching_experiments"), 0
        ),
        database_coverage=_to_float(tool_payload.get("database_coverage"), 0.0),
        is_fallback_match=bool(tool_payload.get("is_fallback_match", False)),
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

    return ConditionRecommendationResponse(
        query=query,
        evidence=evidence,
        recommendations=recommendations,
        warnings=warnings,
    )


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
    ):
        self.llm = get_llm_client(provider, model, temperature)
        self.system_prompt = system_prompt or CHEMISTRY_SYSTEM_PROMPT
        self.verbose = verbose
        self.auto_check = auto_check
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
            structured_payload = _extract_tool_payload(
                result.get("messages", []), "hte_recommend_conditions"
            )
            if structured_payload is not None:
                structured_response = _build_condition_response(query, structured_payload)
                return structured_response.model_dump_json(indent=2)
            final_message = result["messages"][-1]
            if isinstance(final_message, AIMessage):
                return final_message.content
            return str(final_message.content)
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
            unified_result = unified_featurize_reaction.invoke(
                {"reaction_smiles": reaction_smiles}
            )
            if unified_result.get("success"):
                unified_payload = unified_result
            else:
                unified_error = unified_result.get("error")
        except Exception as exc:
            unified_error = str(exc)

        try:
            analysis_result = analysis_analyze_reaction.invoke(
                {"reaction_smiles": reaction_smiles}
            )
            if analysis_result.get("success"):
                analysis_payload = analysis_result
            else:
                analysis_error = analysis_result.get("error")
        except Exception as exc:
            analysis_error = str(exc)

        unified_type = None
        unified_conf = None
        if unified_payload:
            reaction = unified_payload.get("reaction") or {}
            reaction_type = reaction.get("reaction_type") or {}
            unified_type = reaction_type.get("reaction_type") or reaction_type.get("name")
            unified_conf = reaction_type.get("confidence")

        analysis_family = None
        if analysis_payload:
            family = analysis_payload.get("family") or {}
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
