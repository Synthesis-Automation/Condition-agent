"""
ReactionReasoningAgent — "Claude Code for Organic Reactions"

An agentic chemistry reasoner that uses RDKit-based tools + LLM reasoning to
build a comprehensive ReactivityProfile for any reaction. The agent:

1. Inspects molecules (functional groups, descriptors, electronics, sterics)
2. Analyzes the net transformation (bond changes, FG diff)
3. Characterizes electrophile/nucleophile
4. Determines mechanism class with evidence
5. Identifies selectivity risks and condition implications
6. Maps to taxonomy (reaction type + motif labels)

Uses LangGraph for the tool-calling loop, following the same pattern as
chem_assistant/chemtools_agent.py but specialized for reactivity analysis.
"""

from __future__ import annotations

import json
import logging
import os
import re
import time
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)

# System prompt — bond-first reasoning with explicit disambiguation rules
REASONING_SYSTEM_PROMPT = """You are an expert organic chemistry reasoning agent.
You have RDKit-based tools for molecular inspection and a reaction taxonomy database.

Your goal: build a complete ReactivityProfile and identify the CORRECT reaction type.

## CORE PRINCIPLE — BOND-FORMED FIRST

The single most reliable signal for reaction type is the KEY BOND FORMED in the product.
Always determine this first. Never select a reaction type that contradicts it.

  C-C bond formed  →  cross-coupling family (Suzuki, Negishi, Stille, Heck, Sonogashira…)
  C-N bond formed  →  C_N_Coupling, Chan_Lam, Reductive_amination, Amide_formation…
  C-O bond formed  →  C_O_Coupling, Esterification, Acylation_ester…
  C-S bond formed  →  C_S_Coupling
  C-B bond formed  →  Miyaura_borylation
  C-halogen formed →  Halogenation_aromatic, Azide_coupling, Cyanation_coupling…
  no new C-X bond  →  oxidation, reduction, elimination, addition, cycloaddition…

## NUCLEOPHILE → REACTION FAMILY LOOKUP

After identifying the key bond, confirm by checking the nucleophile type:

  Nucleophile contains -B(OH)2 / -Bpin / -BF3K   →  C-C formed  →  Suzuki_miyaura family
  Nucleophile contains -NH2 / -NHR / -NR2 / amide →  C-N formed  →  C_N_Coupling
  Nucleophile contains -OH / -OR (alcohol/alkoxide) →  C-O formed  →  C_O_Coupling
  Nucleophile contains -SH / -SR (thiol/sulfide)   →  C-S formed  →  C_S_Coupling
  Nucleophile contains -ZnX (organozinc)            →  C-C formed  →  Negishi
  Nucleophile contains -SnR3 (organotin)            →  C-C formed  →  Stille
  Nucleophile contains -MgX (Grignard)              →  C-C formed  →  Kumada
  Nucleophile is terminal alkyne (-C≡CH)            →  C-C formed  →  Sonogashira
  Nucleophile is alkene (Heck partner)              →  C-C formed  →  Heck

If bond formed contradicts the nucleophile class above, flag a warning and
re-examine — one of the assignments is wrong.

## ANALYSIS STEPS

### 1. INSPECT MOLECULES
- Call inspect_functional_groups() on EACH reactant and the product separately
- Call compute_molecular_descriptors() on each molecule
- List which functional groups are present in each reactant

### 2. NET TRANSFORMATION  ← MOST CRITICAL STEP
- Call compare_reactant_product() → examine groups_removed and groups_added
- Call analyze_bond_changes()     → examine bonds_broken and bonds_formed
- Explicitly state: "KEY BOND FORMED: [C-C / C-N / C-O / C-S / other]"
- This determination drives ALL subsequent reasoning — do not skip it

### 3. ELECTROPHILE ANALYSIS
- Which reactant has the leaving group (halide, pseudohalide, or activated group)?
- Hybridization of electrophilic carbon: sp2 (aryl/vinyl) or sp3 (alkyl)?
- Call analyze_electronics() on the electrophilic reactant
- Call analyze_steric() on the electrophilic reactant

### 4. NUCLEOPHILE CLASSIFICATION
- Which reactant provides the electron pair / new bond partner?
- Identify the attacking atom (N, C, O, S, B…)
- Apply the NUCLEOPHILE → REACTION FAMILY table above
- Confirm this matches the KEY BOND FORMED from Step 2

### 5. REACTION PATTERN TYPE
- Call get_reaction_pattern_types() with no argument to get all 12 definitions
- Select ONE pattern type consistent with Step 2 + Step 4:
  coupling_substitution | condensation | addition | cycloaddition | cyclization |
  annulation | electrophilic_aromatic_substitution | nucleophilic_aromatic_substitution |
  oxidation | reduction | olefination | elimination

### 6. MECHANISM CLASS
- If sp2 aryl electrophile + nucleophile → likely oxidative addition / reductive elimination cycle
- If sp2 + activated ring + nucleophile displacement → call check_snar_feasibility()
- If sp3 electrophile → SN2/SN1/E2 depending on steric/base
- If C=O + nucleophile → acyl substitution or addition

### 7. SELECTIVITY & RISKS
- Competing pathways (substitution vs elimination, homo-coupling, over-reaction)
- Regioselectivity, chemoselectivity concerns

### 8. CONDITION IMPLICATIONS
- Catalyst type implied by mechanism (Pd, Cu, Ni…)
- Base, solvent, temperature requirements

### 9. TAXONOMY MAPPING
SEARCH STRATEGY — always include KEY BOND TYPE + nucleophile class in query:
  Good:  "C-N coupling amine aryl halide Pd"   (for C-N formed, amine nucleophile)
  Bad:   "Pd coupling aryl halide"              (ignores nucleophile → Suzuki bias)

Procedure:
a) Call search_reaction_types("<bond_type> <nucleophile_class> <optional catalyst>")
b) Call search_motifs() for each reactant based on FG analysis
c) If a broad motif returned, call get_motif_hierarchy() to find the specific label
d) Select the best-matching reaction_type ID

CROSS-VALIDATION (mandatory before outputting):
- Does the selected reaction_type match the KEY BOND FORMED? If not → reject and re-search
- Call get_reaction_pattern_types(reaction_type_id="<your_ID>") to verify
  leaving_groups match what you observed in bond analysis
- If leaving_groups don't match → the reaction type is wrong, search again

ONLY use taxonomy identifiers that appear in search_reaction_types() results.

## COMMON MISTAKES TO AVOID

× Selecting Suzuki_miyaura when C-N bond is formed (Suzuki always forms C-C)
× Selecting Suzuki_miyaura when no boron species is present in reactants
× Picking a reaction type based solely on the electrophile (Ar-I/Br/Cl appear in
  many different reaction families — the NUCLEOPHILE distinguishes them)
× Ignoring the bond change analysis output

## OUTPUT FORMAT

After completing all steps, output ONLY a valid JSON object.
No markdown fences, no explanatory text before or after. Start with { end with }.

{
  "reaction_type": "taxonomy_ID or null",
  "reaction_type_confidence": 0.0-1.0,
  "reaction_pattern_type": "coupling_substitution|condensation|addition|cycloaddition|cyclization|annulation|electrophilic_aromatic_substitution|nucleophilic_aromatic_substitution|oxidation|reduction|olefination|elimination",
  "reacted_motifs": ["motif1", "motif2"],
  "formed_motifs": ["motif1"],
  "taxonomy_reasoning": "Bond formed: C-N. Nucleophile: secondary amide (N-H). Electrophile: aryl iodide. → C_N_Coupling confirmed by leaving_groups match.",

  "electrophile": {
    "center_atom": "description",
    "hybridization": "sp2_aryl|sp2_vinyl|sp3|carbonyl",
    "leaving_group": "I|Br|Cl|OTf|...",
    "leaving_group_quality": "excellent|good|moderate|poor",
    "steric_score": 0.0,
    "steric_class": "unhindered|moderate|hindered",
    "electronic_score": 5.0,
    "electronic_class": "electron_poor|neutral|electron_rich",
    "activation_tags": []
  },

  "nucleophile": {
    "identity": "secondary amide / boronic acid / primary amine / ...",
    "attacking_atom": "N|C|O|S|B",
    "hardsoft": "soft|hard|borderline",
    "is_also_base": false,
    "steric_bulk": "small|moderate|bulky",
    "functional_groups": []
  },

  "mechanism": {
    "primary_class": "oa_reductive_elimination|snar|sn2|sn1|acyl_sub|...",
    "evidence": ["C-N bond formed per analyze_bond_changes", "N-H present in reactant 2"],
    "confidence": 0.0-1.0,
    "alternative_mechanisms": [],
    "requires_catalyst": true,
    "likely_catalyst_metals": ["Pd"]
  },

  "transformation": {
    "bonds_broken": ["C-I (sp2, aryl)", "N-H (amide)"],
    "bonds_formed": ["C-N (sp2-N, aryl-amide)"],
    "key_bond_type": "C-N",
    "fg_removed": ["aryl_iodide", "secondary_amide_NH"],
    "fg_formed": ["N-aryl_amide"],
    "redox_change": "neutral"
  },

  "selectivity_risks": [],
  "competing_pathways": [],
  "driving_forces": [],
  "condition_implications": {},
  "is_tandem": false,
  "confidence": 0.0-1.0,
  "warnings": []
}
"""


def _get_llm_client(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    temperature: float = 0,
):
    """Get a LangChain ChatOpenAI client.

    Model recommendations for the reasoning agent (multi-step tool-calling):
      - openai:  "o4-mini" or "o3"  (best reasoning; o4-mini is cost-effective)
                 "gpt-4.1"          (good balance, supports tool use well)
                 "gpt-4o"           (acceptable but may misclassify edge cases)
      - aliyun:  "qwen-max" or provider-specific reasoning model

    Set LLM_MODEL env var or pass model= to override. Default: o4-mini.
    """
    from langchain_openai import ChatOpenAI

    if provider is None:
        provider = os.getenv("LLM_PROVIDER", "openai")
    if model is None:
        model = os.getenv("LLM_MODEL", "o4-mini")

    if provider == "aliyun":
        api_key = os.getenv("ALIYUN_API_KEY")
        base_url = os.getenv(
            "ALIYUN_BASE_URL",
            "https://dashscope.aliyuncs.com/compatible-mode/v1",
        )
    elif provider == "openai":
        api_key = os.getenv("OPENAI_API_KEY")
        base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")
    else:
        raise ValueError(f"Unsupported provider '{provider}'. Use 'openai' or 'aliyun'.")

    if not api_key:
        raise RuntimeError(f"{provider.upper()}_API_KEY environment variable not set")

    # OpenAI o-series reasoning models (o1, o3, o4-mini, etc.) do not accept
    # a temperature parameter — they only support the default value of 1.
    _reasoning_model_prefixes = ("o1", "o3", "o4")
    is_reasoning_model = any(
        model.startswith(p) for p in _reasoning_model_prefixes
    )

    kwargs: Dict[str, Any] = {
        "model": model,
        "api_key": api_key,
        "base_url": base_url,
    }
    if not is_reasoning_model:
        kwargs["temperature"] = temperature

    return ChatOpenAI(**kwargs)


def _get_agent_factory():
    """Import the LangGraph agent factory (handles API changes across versions)."""
    try:
        from langgraph.prebuilt import create_agent as factory
        return factory
    except (ImportError, AttributeError):
        pass
    try:
        from langgraph.prebuilt import create_react_agent as factory
        return factory
    except (ImportError, AttributeError):
        pass
    from langchain.agents import create_react_agent as factory
    return factory


class ReactionReasoningAgent:
    """
    Agentic chemistry reasoner — uses RDKit tools + LLM to build
    a comprehensive ReactivityProfile for any reaction.

    Like Claude Code, but for organic reactions: the agent decides what to
    investigate based on what it sees, calling tools iteratively.
    """

    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        temperature: float = 0,
        max_iterations: int = 8,
        system_prompt: Optional[str] = None,
        verbose: bool = False,
    ):
        """
        Args:
            provider: LLM provider ("openai" or "aliyun")
            model: Model name (e.g. "gpt-4o", "gpt-5.2")
            temperature: LLM temperature (0 for deterministic)
            max_iterations: Maximum tool-calling turns
            system_prompt: Override the default reasoning prompt
            verbose: Print debug info during execution
        """
        self.provider = provider
        self.model_name = model or os.getenv("LLM_MODEL", "o4-mini")
        self.max_iterations = max_iterations
        self.verbose = verbose

        # Build LLM + tools + agent
        self.llm = _get_llm_client(provider, model, temperature)

        from reaction_agent.reasoning_tools import REASONING_TOOLS
        self.tools = REASONING_TOOLS

        factory = _get_agent_factory()
        self.agent = factory(
            self.llm,
            self.tools,
            prompt=system_prompt or REASONING_SYSTEM_PROMPT,
        )

    def analyze(self, reaction_smiles: str) -> "ReactivityProfile":
        """
        Run the reasoning agent on a reaction SMILES.

        The agent will:
        1. Call tools to inspect molecules
        2. Reason about mechanism and reactivity
        3. Map to taxonomy identifiers
        4. Return a rich ReactivityProfile

        Args:
            reaction_smiles: Reaction SMILES (reactants>>products)

        Returns:
            ReactivityProfile with comprehensive analysis
        """
        from reaction_agent.reactivity_profile import ReactivityProfile

        start_time = time.time()

        # Parse input
        if ">>" not in reaction_smiles:
            profile = ReactivityProfile(
                reaction_smiles=reaction_smiles,
                warnings=["Invalid reaction SMILES: missing '>>'"],
            )
            return profile

        reactants_str, products_str = reaction_smiles.split(">>", 1)
        reactant_smiles = [s.strip() for s in reactants_str.split(".") if s.strip()]
        product_smiles = [s.strip() for s in products_str.split(".") if s.strip()]

        # Build the query message
        query = (
            f"Analyze this reaction and build a complete reactivity profile:\n\n"
            f"Reaction SMILES: {reaction_smiles}\n"
            f"Reactants: {reactant_smiles}\n"
            f"Products: {product_smiles}\n\n"
            f"Follow the systematic checklist in your instructions. "
            f"Call the tools, reason about the chemistry, and return the JSON profile."
        )

        # Run the LangGraph agent
        try:
            from langchain_core.messages import HumanMessage, AIMessage
            messages = [HumanMessage(content=query)]

            result = self.agent.invoke(
                {"messages": messages},
                config={"recursion_limit": self.max_iterations * 2 + 2},
            )

            # Extract result
            result_payload = result if isinstance(result, dict) else {}
            result_messages = result_payload.get("messages", [])
            if not isinstance(result_messages, list):
                result_messages = []

            # Extract tool call records for the reasoning chain
            reasoning_chain = []
            tools_called = []
            for msg in result_messages:
                msg_type = type(msg).__name__
                if msg_type == "ToolMessage":
                    tool_name = getattr(msg, "name", "unknown")
                    tools_called.append(tool_name)
                    reasoning_chain.append({
                        "tool": tool_name,
                        "output_preview": str(getattr(msg, "content", ""))[:200],
                    })

            # Get the final AI message content
            final_content = ""
            for msg in reversed(result_messages):
                if isinstance(msg, AIMessage) and msg.content:
                    final_content = msg.content
                    break

            if not final_content:
                final_content = str(result)

            # Parse JSON from the final response
            profile = self._parse_profile(
                final_content,
                reaction_smiles=reaction_smiles,
                reactant_smiles=reactant_smiles,
                product_smiles=product_smiles,
                reasoning_chain=reasoning_chain,
                tools_called=tools_called,
            )

            elapsed = time.time() - start_time
            profile.llm_model = self.model_name
            profile.total_tool_calls = len(tools_called)
            if self.verbose:
                logger.info(
                    f"ReasoningAgent completed: {len(tools_called)} tool calls, "
                    f"{elapsed:.1f}s, model={self.model_name}"
                )

            return profile

        except Exception as exc:
            elapsed = time.time() - start_time
            logger.error(f"ReasoningAgent failed after {elapsed:.1f}s: {exc}")
            return ReactivityProfile(
                reaction_smiles=reaction_smiles,
                reactant_smiles=reactant_smiles,
                product_smiles=product_smiles,
                llm_model=self.model_name,
                warnings=[f"Agent execution failed: {exc}"],
            )

    def _parse_profile(
        self,
        raw_content: str,
        reaction_smiles: str,
        reactant_smiles: List[str],
        product_smiles: List[str],
        reasoning_chain: List[Dict[str, Any]],
        tools_called: List[str],
    ) -> "ReactivityProfile":
        """Parse the agent's JSON output into a ReactivityProfile."""
        from reaction_agent.reactivity_profile import (
            ReactivityProfile,
            ElectrophileProfile,
            NucleophileProfile,
            MechanismAnalysis,
            TransformationRecord,
        )
        from reaction_agent.taxonomy_prompts import TaxonomyContext

        warnings: List[str] = []

        # Try to extract JSON from the response
        text = raw_content.strip()
        text = re.sub(r'^```(?:json)?\s*\n', '', text, flags=re.MULTILINE)
        text = re.sub(r'\n```\s*$', '', text, flags=re.MULTILINE)
        text = text.strip()

        # Find JSON object in the text
        json_match = re.search(r'\{[\s\S]*\}', text)
        if json_match:
            text = json_match.group(0)

        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            warnings.append(f"Could not parse agent JSON output: {e}")
            data = {}

        # Validate taxonomy IDs
        taxonomy = TaxonomyContext()
        reaction_type = data.get("reaction_type")
        if reaction_type and not taxonomy.validate_reaction_type(reaction_type):
            warnings.append(f"Invalid reaction_type '{reaction_type}' — not in taxonomy")
            reaction_type = None

        raw_reacted = data.get("reacted_motifs", [])
        if not isinstance(raw_reacted, list):
            raw_reacted = []
        valid_reacted, inv_reacted = taxonomy.filter_motifs(raw_reacted)
        if inv_reacted:
            warnings.append(f"Invalid reacted_motifs discarded: {inv_reacted}")

        raw_formed = data.get("formed_motifs", [])
        if not isinstance(raw_formed, list):
            raw_formed = []
        valid_formed, inv_formed = taxonomy.filter_motifs(raw_formed)
        if inv_formed:
            warnings.append(f"Invalid formed_motifs discarded: {inv_formed}")

        # Build sub-profiles
        elec_data = data.get("electrophile", {})
        nuc_data = data.get("nucleophile", {})
        mech_data = data.get("mechanism", {})
        trans_data = data.get("transformation", {})

        electrophile = ElectrophileProfile(
            center_atom=str(elec_data.get("center_atom", "")),
            hybridization=str(elec_data.get("hybridization", "")),
            leaving_group=str(elec_data.get("leaving_group", "")),
            leaving_group_quality=str(elec_data.get("leaving_group_quality", "")),
            steric_score=float(elec_data.get("steric_score", 0.0)),
            steric_class=str(elec_data.get("steric_class", "")),
            electronic_score=float(elec_data.get("electronic_score", 5.0)),
            electronic_class=str(elec_data.get("electronic_class", "")),
            activation_tags=elec_data.get("activation_tags", []),
        )

        nucleophile = NucleophileProfile(
            identity=str(nuc_data.get("identity", "")),
            attacking_atom=str(nuc_data.get("attacking_atom", "")),
            hardsoft=str(nuc_data.get("hardsoft", "")),
            is_also_base=bool(nuc_data.get("is_also_base", False)),
            steric_bulk=str(nuc_data.get("steric_bulk", "")),
            functional_groups=nuc_data.get("functional_groups", []),
        )

        mechanism = MechanismAnalysis(
            primary_class=str(mech_data.get("primary_class", "")),
            evidence=mech_data.get("evidence", []),
            confidence=float(mech_data.get("confidence", 0.0)),
            alternative_mechanisms=mech_data.get("alternative_mechanisms", []),
            requires_catalyst=bool(mech_data.get("requires_catalyst", False)),
            likely_catalyst_metals=mech_data.get("likely_catalyst_metals", []),
        )

        # SNAr feasibility if present
        snar = data.get("snar_feasibility")
        if snar:
            mechanism.snar_feasibility = snar

        transformation = TransformationRecord(
            bonds_broken=trans_data.get("bonds_broken", []),
            bonds_formed=trans_data.get("bonds_formed", []),
            key_bond_type=str(trans_data.get("key_bond_type", "")),
            fg_removed=trans_data.get("fg_removed", []),
            fg_formed=trans_data.get("fg_formed", []),
            redox_change=str(trans_data.get("redox_change", "")),
        )

        # Validate reaction_pattern_type against known values
        valid_pattern_types = {
            "coupling_substitution", "condensation", "addition", "cycloaddition",
            "cyclization", "annulation", "electrophilic_aromatic_substitution",
            "nucleophilic_aromatic_substitution", "oxidation", "reduction",
            "olefination", "elimination",
        }
        reaction_pattern_type = str(data.get("reaction_pattern_type", ""))
        if reaction_pattern_type and reaction_pattern_type not in valid_pattern_types:
            warnings.append(
                f"Unknown reaction_pattern_type '{reaction_pattern_type}' — cleared"
            )
            reaction_pattern_type = ""

        return ReactivityProfile(
            reaction_smiles=reaction_smiles,
            reactant_smiles=reactant_smiles,
            product_smiles=product_smiles,
            transformation=transformation,
            electrophile=electrophile,
            nucleophile=nucleophile,
            mechanism=mechanism,
            selectivity_risks=data.get("selectivity_risks", []),
            competing_pathways=data.get("competing_pathways", []),
            driving_forces=data.get("driving_forces", []),
            condition_implications=data.get("condition_implications", {}),
            functional_group_inventory=data.get("functional_group_inventory", {}),
            molecular_descriptors=data.get("molecular_descriptors", {}),
            reaction_type=reaction_type,
            reaction_type_confidence=float(data.get("reaction_type_confidence", 0.0)),
            reaction_pattern_type=reaction_pattern_type,
            reacted_motifs=valid_reacted,
            formed_motifs=valid_formed,
            taxonomy_reasoning=str(data.get("taxonomy_reasoning", "")),
            is_tandem=bool(data.get("is_tandem", False)),
            tandem_steps=data.get("tandem_steps", []),
            reasoning_chain=reasoning_chain,
            tools_called=tools_called,
            confidence=float(data.get("confidence", 0.0)),
            warnings=warnings + data.get("warnings", []),
        )


__all__ = [
    "ReactionReasoningAgent",
    "REASONING_SYSTEM_PROMPT",
]
