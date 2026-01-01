"""
LangGraph agent for ChemTools.

This module provides an agent that can use ChemTools functions to analyze
reactions, recommend conditions, and answer chemistry questions.

The agent can:
    - Normalize SMILES and reaction SMILES
    - Detect reaction families and types
    - Recommend optimal reaction conditions
    - Search for similar precedent reactions
    - Look up reagent information
    - Analyze reactants and functional groups

Usage:
    from lang_chain.chemtools_agent import ChemToolsAgent

    agent = ChemToolsAgent()
    result = agent.run("What conditions should I use for this reaction: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
    print(result)
"""

import json
import os
import re
from typing import Any, Callable, Dict, List, Optional, Set, Union

from langchain_openai import ChatOpenAI
from langchain_core.messages import AIMessage, BaseMessage, HumanMessage, SystemMessage
from langchain_core.tools import BaseTool, StructuredTool
from dotenv import load_dotenv

from .chemtools_wrapper import CHEMTOOLS_TOOLS
from .planner import ReactionInput, auto_conditions
from .constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    format_constraints_for_prompt,
    merge_specs,
)

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


# Load environment variables
load_dotenv()


# ============================================================================
# LLM Configuration
# ============================================================================

DEFAULT_ALIYUN_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"

_REACTION_SMILES_RE = re.compile(
    r"([A-Za-z0-9@+\-\[\]\(\)=#%/.]+>>[A-Za-z0-9@+\-\[\]\(\)=#%/.]*)"
)


def get_llm_client(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0
) -> ChatOpenAI:
    """
    Get LLM client using configured provider and model.
    
    Args:
        provider: LLM provider ("openai" or "aliyun"). If None, reads from env.
        model: Model name. If None, reads from env or uses default.
        temperature: Sampling temperature (0.0-1.0)
    
    Returns:
        Configured ChatOpenAI client
    
    Environment variables:
        LLM_PROVIDER: "openai" or "aliyun" (default: "openai")
        LLM_MODEL: Model name (default: "gpt-4o")
        OPENAI_API_KEY: OpenAI API key
        OPENAI_BASE_URL: OpenAI base URL (optional)
        ALIYUN_API_KEY: Aliyun API key
        ALIYUN_BASE_URL: Aliyun base URL (optional)
    """
    # Determine provider from parameter or environment
    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "openai")
    
    # Determine model from parameter or environment
    if model is None:
        model = os.getenv("LLM_MODEL", "gpt-4o")
    
    # Get API credentials and base URL
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
        temperature=temperature
    )


# ============================================================================
# System Prompts
# ============================================================================

CHEMISTRY_SYSTEM_PROMPT = """You are ChemBot, an expert chemistry assistant with access to powerful ChemTools for analyzing reactions and recommending conditions.

You have access to the following tools:

1. **normalize_smiles_tool**: Canonicalize SMILES strings
2. **normalize_reaction_tool**: Canonicalize reaction SMILES
3. **detect_reaction_family_tool**: Identify reaction type (Suzuki, Buchwald, etc.)
4. **classify_reactant_tool**: Classify reactant types (aryl halide, amine, etc.)
5. **get_functional_groups_tool**: Detect functional groups in molecules
6. **calculable_features_tool**: Evaluate curated calculable feature library for a molecule
7. **molpipeline_featurize_tool**: Generate molecular features (optional MolPipeline vectors)
8. **analyze_bond_changes_tool**: Analyze bond breaking/formation in reactions
9. **reaction_similarity_tool**: Compare two reactions using DRFP similarity
10. **recommend_conditions_tool**: Get ML-based condition recommendations
11. **rule_based_conditions_tool**: Deterministic rule-engine guidance with feature-level reasoning
12. **search_precedents_tool**: Find similar precedent reactions
13. **protocol_recommendation_tool**: Retrieve full experimental protocols from literature
14. **unified_recommender_tool**: Unified protocol + rule search (DRFP similarity)
15. **reaction_dataset_analytics_tool**: Summarize reaction dataset composition, yields, and popular reagents
16. **find_reagent_tool**: Look up reagent information from database
17. **reagent_database_analytics_tool**: Summarize reagent database composition and completeness
18. **list_supported_cores_tool**: List catalyst cores observed in similar precedents
19. **list_all_families_tool**: List all reaction families available in the dataset
20. **add_reagent_tool**: Add or dry-run reagent entries in the taxonomy registry
21. **rule_builder_autofill_tool**: Convert protocol text into draft rule databases (LLM-assisted + deterministic validation)
22. **hte_recommend_tool**: Get data-driven condition recommendations from 66K+ HTE experiments with catalyst/reaction filtering
23. **hte_analytics_tool**: Explore HTE database reactant pairs, catalysts, and success patterns
24. **hte_conditions_tool**: Get detailed experimental conditions for specific substrate combinations
25. **hte_screening_set_tool**: Generate diverse condition sets for HTE screening plates (up to 24 conditions for parallel testing)

**How to help users:**

For reaction condition recommendations:
1. First normalize the reaction SMILES if needed
2. Detect the reaction family to understand the reaction type
3. Choose the appropriate recommendation tool:
   - **hte_screening_set_tool**: For generating DIVERSE condition sets for HTE screening plates (12-24 conditions for parallel testing). USE THIS when users want to "screen", "test multiple conditions", "set up a plate", or need a "group of conditions". This is the PRIMARY HTE use case.
   - **hte_recommend_tool**: For data-driven top-K recommendations from 66K+ HTE experiments (fast, <100ms, statistical evidence with Z-scores). USE THIS when users explicitly mention "HTE system" or want experimental data-backed conditions for a SINGLE best condition or small set.
   - **unified_recommender_tool**: Default best-effort protocol + rule search when users want practical conditions and/or automation-ready output.
   - **protocol_recommendation_tool**: For full literature protocols (stepwise procedures) when users ask for "protocol", "procedure", or "paper conditions".
   - **recommend_conditions_tool**: For ML-based predictions with diverse precedent coverage
   - **rule_based_conditions_tool**: For deterministic guidance with feature-level reasoning
4. Include catalyst and reagent constraints via the tool parameters when the user mentions preferences (use constraint_text / allow_metals / exclude_metals / prefer_metals / search_all_families).
5. Call list_supported_cores_tool when you need to understand which catalyst cores exist before setting constraints.
6. Optionally search for precedents to provide context
7. Explain the recommendations in clear, practical terms

Before calling any tool, always pause to analyze the user's intent:
- If the request clearly requires ChemTools data (e.g., reaction conditions, reagent lookup, functional groups, reagent database updates), proceed with the relevant tool workflow.
- If the request falls outside ChemTools coverage (e.g., high-level synthetic planning, theory explanation, unrelated topics), respond using your own knowledge without invoking tools, but still provide a thoughtful, chemistry-aware answer.

For reagent database updates:
1. Collect whatever structured data the user already supplied (CAS, role, name, synonyms, family).
2. If some fields are missing, still run add_reagent_tool with `dry_run=True` and `auto_resolve=True` so the resolver can fill gaps. Share the previewed entry with the user.
3. When the user explicitly approves the write, rerun add_reagent_tool with `dry_run=False` (and optional `taxonomy_dir`) to persist the entry, then relay the status and destination path.
4. If the user declines to persist, keep the dry-run preview only and offer further adjustments.

For dataset analytics:
1. Use reaction_dataset_analytics_tool to answer questions about dataset-wide statistics (e.g., counts of solvents, top reagents, yield distributions)
2. When users ask about reagent inventory totals, prefer reagent_database_analytics_tool

For reagent questions:
1. Use find_reagent_tool to look up properties and roles
2. Use reagent_database_analytics_tool for aggregate counts (e.g., solvent totals, missing data)
3. Provide CAS numbers, structures, and usage information

For reactant/molecule analysis:
1. Normalize SMILES first
2. Use classify_reactant_tool or get_functional_groups_tool
3. Explain the structural features clearly

For HTE (High-Throughput Experimentation) queries:
1. **hte_recommend_tool**: Primary tool for HTE condition recommendations
   - Automatically extracts reactants from reaction SMILES format
   - Filters by catalyst metal (Cu, Pd, Ni, etc.) and reaction type
   - Returns conditions ranked by Z-score (primary metric) and confidence
   - Fast queries (<100ms) with statistical evidence
   - USE THIS when users say "use HTE system" or "HTE data"
2. **hte_analytics_tool**: Explore database coverage and patterns
   - List available reactant pairs for specific catalysts/reactions
   - Analyze catalyst usage and success rates
   - Understand metal distribution (Pd vs Cu vs Ni)
3. **hte_conditions_tool**: Get detailed conditions for specific substrate pairs
   - Returns full experimental details (temperature, time, concentration)
   - Best for when users need complete protocols

**Important guidelines:**
- Always normalize SMILES before analysis
- Tools return JSON - parse it and present results clearly
- For reactions, use the format: "reactants>>products"
- Explain recommendations with practical context
- When tools return errors, explain the issue and suggest fixes
- Be concise but informative
- Focus on actionable insights

Remember: You're helping chemists design better experiments!
"""

_CONSTRAINT_TOOL_FIELDS: Set[str] = {
    "allow_metals",
    "exclude_metals",
    "prefer_metals",
    "search_all_families",
    "constraint_rules",
}


def _schema_field_names(args_schema: object) -> Set[str]:
    """Return pydantic model field names (v1/v2 compatible)."""
    fields = getattr(args_schema, "model_fields", None)
    if isinstance(fields, dict):
        return set(fields.keys())
    legacy_fields = getattr(args_schema, "__fields__", None)
    if isinstance(legacy_fields, dict):
        return set(legacy_fields.keys())
    return set()


# ============================================================================
# ChemTools Agent Class
# ============================================================================

class ChemToolsAgent:
    """
    LangGraph agent for chemistry tasks using ChemTools.

    This agent combines LLM reasoning with deterministic ChemTools functions.
    It prefers LangGraph's `create_agent` helper (auto-selecting tool-calling
    vs. ReAct) and falls back to `create_react_agent` on older releases.

    Example:
        >>> agent = ChemToolsAgent()
        >>> result = agent.run("Recommend conditions for Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1")
        >>> print(result)

        >>> # With conversation history
        >>> history = []
        >>> result = agent.run("What is the CAS number for Cs2CO3?", history=history)
        >>> result = agent.run("What role does it play?", history=history)
    """
    
    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0,
        system_prompt: Optional[str] = None,
        verbose: bool = False,
        proactive: bool = False,
        proactive_top_k: int = 5,
        proactive_max_protocols: int = 3,
        proactive_build_protocols: bool = True,
        session_constraints: Optional[Union[ConstraintSpec, dict, str]] = None,
    ):
        """
        Initialize ChemTools agent.
        
        Args:
            provider: LLM provider ("openai" or "aliyun")
            model: Model name
            temperature: Sampling temperature (0.0-1.0)
            system_prompt: Custom system prompt (default: CHEMISTRY_SYSTEM_PROMPT)
            verbose: Print debug information
            proactive: If True, precompute deterministic auto-conditions for reactions
            proactive_top_k: Number of DRFP protocols to retrieve in preflight
            proactive_max_protocols: Number of protocols to format in preflight
            proactive_build_protocols: Whether to build protocol additions in preflight
        """
        self.llm = get_llm_client(provider, model, temperature)
        self.system_prompt = system_prompt or CHEMISTRY_SYSTEM_PROMPT
        self.verbose = verbose
        self.proactive = proactive
        self.proactive_top_k = proactive_top_k
        self.proactive_max_protocols = proactive_max_protocols
        self.proactive_build_protocols = proactive_build_protocols
        self.session_constraints = self._coerce_constraints(session_constraints)
        self.agent_factory_name = _LANGGRAPH_AGENT_FACTORY_NAME
        self._active_constraints = ConstraintSpec()
        
        # Create agent graph with tools (prefers LangGraph create_agent when available)
        self.tools = self._wrap_tools_with_constraint_defaults(CHEMTOOLS_TOOLS)
        self.agent = LANGGRAPH_AGENT_FACTORY(
            self.llm,
            self.tools,
            prompt=self.system_prompt
        )
    
    def run(
        self,
        query: str,
        history: Optional[List[BaseMessage]] = None,
        recursion_limit: int = 15,
        constraints: Optional[Union[ConstraintSpec, dict, str]] = None
    ) -> str:
        """
        Run the agent on a query.
        
        Args:
            query: User question or task
            history: Conversation history (list of messages)
            recursion_limit: Maximum reasoning steps
            constraints: Optional one-off constraint specification for this call
        
        Returns:
            Agent's response text
        
        Example:
            >>> agent = ChemToolsAgent()
            >>> response = agent.run("What conditions for Suzuki coupling?")
        """
        try:
            messages = list(history or [])

            # Merge session and call-specific constraints.
            call_spec = self._coerce_constraints(constraints)
            active_spec = merge_specs(self.session_constraints, call_spec)

            system_messages: List[BaseMessage] = []
            constraint_prompt = format_constraints_for_prompt(active_spec)
            if constraint_prompt:
                system_messages.append(SystemMessage(content=constraint_prompt))

            preflight = self._maybe_run_proactive_preflight(query, active_spec)
            if preflight:
                system_messages.append(SystemMessage(content=preflight))

            messages = system_messages + messages
            messages.append(HumanMessage(content=query))
        
            # Invoke agent with query and history
            self._active_constraints = active_spec
            result = self.agent.invoke(
                {"messages": messages},
                config={"recursion_limit": recursion_limit}
            )
            
            # Extract final AI message
            final_message = result["messages"][-1]
            
            if isinstance(final_message, AIMessage):
                return final_message.content
            else:
                return str(final_message.content)
                
        except Exception as e:
            error_msg = f"Error: {str(e)}\n\nPlease rephrase your question or provide more details."
            if self.verbose:
                print(f"Agent error: {e}")
            return error_msg
        finally:
            self._active_constraints = ConstraintSpec()
    
    def chat(
        self,
        query: str,
        history: List[BaseMessage],
        recursion_limit: int = 15,
        constraints: Optional[Union[ConstraintSpec, dict, str]] = None
    ) -> tuple[str, List[BaseMessage]]:
        """
        Chat with the agent and maintain conversation history.
        
        Args:
            query: User question
            history: Current conversation history
            recursion_limit: Maximum reasoning steps
            constraints: Optional one-off constraint specification for this turn
        
        Returns:
            Tuple of (response, updated_history)
        
        Example:
            >>> agent = ChemToolsAgent()
            >>> history = []
            >>> response, history = agent.chat("Normalize c1ccccc1", history)
            >>> response, history = agent.chat("What functional groups?", history)
        """
        response = self.run(query, history, recursion_limit, constraints=constraints)
        
        # Update history
        updated_history = history + [
            HumanMessage(content=query),
            AIMessage(content=response)
        ]
        
        return response, updated_history

    # ======================================================================
    # Constraint helpers
    # ======================================================================

    def _coerce_constraints(
        self,
        value: Optional[Union[ConstraintSpec, dict, str]]
    ) -> ConstraintSpec:
        """Normalize value into a ConstraintSpec instance."""
        if isinstance(value, ConstraintSpec):
            return value
        if isinstance(value, dict):
            return build_constraint_spec(
                text=value.get("constraint_text"),
                allow_metals=value.get("allow_metals"),
                exclude_metals=value.get("exclude_metals"),
                prefer_metals=value.get("prefer_metals"),
                search_all_families=value.get("search_all_families"),
                base_constraint_rules=value.get("constraint_rules"),
            )
        if isinstance(value, str):
            return build_constraint_spec(text=value)
        return ConstraintSpec()

    def update_session_constraints(self, **kwargs) -> ConstraintSpec:
        """
        Merge new parameters into the session-level constraint defaults.
        
        Returns the updated constraint specification.
        """
        spec = build_constraint_spec(
            text=kwargs.get("constraint_text"),
            allow_metals=kwargs.get("allow_metals"),
            exclude_metals=kwargs.get("exclude_metals"),
            prefer_metals=kwargs.get("prefer_metals"),
            search_all_families=kwargs.get("search_all_families"),
            base_constraint_rules=kwargs.get("constraint_rules"),
        )
        self.session_constraints = merge_specs(self.session_constraints, spec)
        return self.session_constraints

    def clear_session_constraints(self) -> None:
        """Remove all stored session-level constraints."""
        self.session_constraints = ConstraintSpec()

    # ======================================================================
    # Proactive preflight helpers
    # ======================================================================

    @staticmethod
    def _extract_reaction_smiles(text: str) -> Optional[str]:
        """Extract a reaction SMILES from free-form text."""
        match = _REACTION_SMILES_RE.search(text)
        if not match:
            return None
        return match.group(1)

    @staticmethod
    def _model_to_dict(obj: Any) -> Dict[str, Any]:
        """Convert a pydantic model to a plain dict (v1/v2 compatible)."""
        if hasattr(obj, "model_dump"):
            return obj.model_dump()  # type: ignore[no-any-return]
        if hasattr(obj, "dict"):
            return obj.dict()  # type: ignore[no-any-return]
        return dict(obj) if isinstance(obj, dict) else {"value": obj}

    def _maybe_run_proactive_preflight(
        self,
        query: str,
        constraints: ConstraintSpec,
    ) -> Optional[str]:
        """Run deterministic auto-conditions and return a summary string."""
        if not self.proactive:
            return None

        reaction_smiles = self._extract_reaction_smiles(query)
        if not reaction_smiles:
            return None

        try:
            result = auto_conditions(
                ReactionInput(reaction_smiles=reaction_smiles),
                constraints=constraints,
                top_k_protocols=self.proactive_top_k,
                build_protocols=self.proactive_build_protocols,
                max_protocols=self.proactive_max_protocols,
            )
        except Exception as exc:
            if self.verbose:
                print(f"Proactive preflight failed: {exc}")
            return None

        summary: Dict[str, Any] = {
            "reaction_smiles": reaction_smiles,
            "family": self._model_to_dict(result.family),
            "plan": self._model_to_dict(result.plan),
            "counts": {
                "rule_candidates": len(result.rule_candidates),
                "protocol_candidates": len(result.protocol_candidates),
            },
            "hte_summary": self._model_to_dict(result.hte_summary) if result.hte_summary else None,
            "top_protocols": [],
        }
        for proto in result.protocols:
            summary["top_protocols"].append(
                {
                    "candidate_id": proto.candidate_id,
                    "source": proto.source,
                    "additions": proto.additions[:6],
                    "notes": proto.notes,
                }
            )

        payload = json.dumps(summary, indent=2, sort_keys=True)
        return (
            "Deterministic preflight summary (auto_conditions):\n"
            f"{payload}\n"
            "Use this to guide tool selection and explanations."
        )

    # ======================================================================
    # Tool wrapping (constraint defaults)
    # ======================================================================

    def _wrap_tools_with_constraint_defaults(self, tools: List[BaseTool]) -> List[BaseTool]:
        """Wrap tools so constraint fields default to the active ConstraintSpec."""
        wrapped: List[BaseTool] = []
        for tool_obj in tools:
            args_schema = getattr(tool_obj, "args_schema", None)
            if not args_schema:
                wrapped.append(tool_obj)
                continue

            schema_fields = _schema_field_names(args_schema)
            if not (schema_fields & _CONSTRAINT_TOOL_FIELDS):
                wrapped.append(tool_obj)
                continue

            wrapped.append(self._make_constrained_tool(tool_obj, schema_fields, args_schema))
        return wrapped

    def _make_constrained_tool(
        self,
        tool_obj: BaseTool,
        schema_fields: Set[str],
        args_schema: object,
    ) -> BaseTool:
        """Create a new tool object that injects defaults from active constraints."""

        def _call(**kwargs: Any) -> Any:
            merged_payload = self._merge_constraints_into_payload(dict(kwargs), schema_fields)
            return tool_obj.invoke(merged_payload)

        return StructuredTool.from_function(
            func=_call,
            name=tool_obj.name,
            description=tool_obj.description,
            args_schema=args_schema,  # type: ignore[arg-type]
            return_direct=getattr(tool_obj, "return_direct", False),
            infer_schema=False,
        )

    def _merge_constraints_into_payload(
        self,
        payload: Dict[str, Any],
        schema_fields: Set[str],
    ) -> Dict[str, Any]:
        """Merge the active ConstraintSpec into a tool payload."""
        spec = self._active_constraints
        if not isinstance(spec, ConstraintSpec):
            return payload

        merged = dict(payload)

        if "allow_metals" in schema_fields and spec.allow_metals:
            existing_allow = merged.get("allow_metals")
            existing_allow_set = (
                set(existing_allow or []) if isinstance(existing_allow, (list, tuple, set)) else set()
            )
            merged_allow = existing_allow_set & spec.allow_metals if existing_allow_set else set(spec.allow_metals)
            merged["allow_metals"] = sorted(merged_allow or spec.allow_metals)

        if "exclude_metals" in schema_fields and spec.exclude_metals:
            existing_exclude = merged.get("exclude_metals")
            existing_exclude_set = (
                set(existing_exclude or []) if isinstance(existing_exclude, (list, tuple, set)) else set()
            )
            merged["exclude_metals"] = sorted(existing_exclude_set | spec.exclude_metals)

        if "prefer_metals" in schema_fields and spec.prefer_metals:
            existing_prefer = merged.get("prefer_metals")
            existing_prefer_set = (
                set(existing_prefer or []) if isinstance(existing_prefer, (list, tuple, set)) else set()
            )
            merged["prefer_metals"] = sorted(existing_prefer_set | spec.prefer_metals)

        if "search_all_families" in schema_fields and spec.search_all_families:
            merged["search_all_families"] = True

        if "constraint_rules" in schema_fields and spec.constraint_rules:
            existing_rules = merged.get("constraint_rules")
            merged_rules: Dict[str, Any] = dict(existing_rules) if isinstance(existing_rules, dict) else {}
            merged_rules.update(spec.constraint_rules)
            merged["constraint_rules"] = merged_rules

        # Keep allow/exclude consistent when both fields exist.
        if "allow_metals" in schema_fields and "exclude_metals" in schema_fields:
            allow_values = merged.get("allow_metals")
            exclude_values = merged.get("exclude_metals")
            allow_set = set(allow_values or []) if isinstance(allow_values, (list, tuple, set)) else set()
            exclude_set = set(exclude_values or []) if isinstance(exclude_values, (list, tuple, set)) else set()
            if allow_set and exclude_set:
                merged["exclude_metals"] = sorted(exclude_set - allow_set)

        return merged


# ============================================================================
# Convenience Functions
# ============================================================================

def create_agent(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    proactive: bool = False,
    proactive_top_k: int = 5,
    proactive_max_protocols: int = 3,
    proactive_build_protocols: bool = True,
    session_constraints: Optional[Union[ConstraintSpec, dict, str]] = None,
    **kwargs
) -> ChemToolsAgent:
    """
    Convenience function to create a ChemTools agent.
    
    Args:
        provider: LLM provider ("openai" or "aliyun")
        model: Model name
        session_constraints: Optional default constraints applied to each query
        **kwargs: Additional arguments for ChemToolsAgent
    
    Returns:
        Configured ChemToolsAgent instance
    """
    return ChemToolsAgent(
        provider=provider,
        model=model,
        proactive=proactive,
        proactive_top_k=proactive_top_k,
        proactive_max_protocols=proactive_max_protocols,
        proactive_build_protocols=proactive_build_protocols,
        session_constraints=session_constraints,
        **kwargs
    )


def quick_query(
    query: str,
    constraints: Optional[Union[ConstraintSpec, dict, str]] = None,
    proactive: bool = False,
    proactive_top_k: int = 5,
    proactive_max_protocols: int = 3,
    proactive_build_protocols: bool = True,
    **kwargs
) -> str:
    """
    Quick one-shot query without maintaining state.
    
    Args:
        query: Question or task
        constraints: Optional one-off constraint specification
        **kwargs: Arguments for ChemToolsAgent
    
    Returns:
        Agent response
    
    Example:
        >>> from lang_chain.chemtools_agent import quick_query
        >>> result = quick_query("Recommend conditions for Suzuki coupling")
    """
    agent = ChemToolsAgent(
        proactive=proactive,
        proactive_top_k=proactive_top_k,
        proactive_max_protocols=proactive_max_protocols,
        proactive_build_protocols=proactive_build_protocols,
        **kwargs,
    )
    return agent.run(query, constraints=constraints)


# ============================================================================
# Main (for testing)
# ============================================================================

if __name__ == "__main__":
    import sys
    
    # Simple test
    print("=" * 70)
    print("ChemTools Agent Test")
    print("=" * 70)
    
    try:
        agent = ChemToolsAgent(verbose=True)
        
        # Test queries
        test_queries = [
            "Normalize this SMILES: c1ccccc1",
            "What reaction family is this: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            "Look up Cs2CO3 in the reagent database",
        ]
        
        for i, query in enumerate(test_queries, 1):
            print(f"\n{'=' * 70}")
            print(f"Test {i}: {query}")
            print('=' * 70)
            
            response = agent.run(query)
            print(response)
        
        print(f"\n{'=' * 70}")
        print("✅ All tests completed!")
        print('=' * 70)
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        print("\nMake sure your API keys are set:")
        print("  - OPENAI_API_KEY or ALIYUN_API_KEY")
        print("  - Optionally set LLM_PROVIDER and LLM_MODEL")
        sys.exit(1)
