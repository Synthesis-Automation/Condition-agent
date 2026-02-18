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
from typing import Any, Dict, List, Optional, Tuple

from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Message / JSON extraction helpers
# ---------------------------------------------------------------------------

def _get_message_text(msg: Any) -> str:
    """Extract plain text from a LangChain message.

    Handles:
    - content as str (standard OpenAI)
    - content as list of dicts (multi-modal / some providers)
    - empty content on tool-call-only messages
    """
    content = getattr(msg, "content", "")
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts: List[str] = []
        for item in content:
            if isinstance(item, dict):
                t = item.get("text") or item.get("content", "")
                if t:
                    parts.append(str(t))
            elif isinstance(item, str):
                parts.append(item)
        return "\n".join(parts)
    return ""


def _extract_json_from_text(text: str) -> Optional[Dict[str, Any]]:
    """Robustly extract a JSON object from model output text.

    Handles:
    - Plain JSON response (ideal case)
    - JSON wrapped in ```json ... ``` or ``` ... ``` fences
    - JSON preceded/followed by explanatory text
    - Control characters embedded in JSON (MiniMax issue)
    """
    if not text:
        return None

    # Strip ASCII control characters (fixes MiniMax invalid-control-char issue)
    # Keep \t, \n, \r but remove other control chars (\x00-\x08, \x0b-\x0c, \x0e-\x1f)
    cleaned = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]", "", text)

    # Strategy 1: JSON inside markdown code fence
    fence_match = re.search(r"```(?:json)?\s*\n?([\s\S]*?)\n?```", cleaned)
    if fence_match:
        candidate = fence_match.group(1).strip()
        if candidate:
            try:
                return json.loads(candidate)
            except json.JSONDecodeError:
                pass

    # Strategy 2: All {...} blocks — try last (most recent) first
    brace_spans = list(re.finditer(r"\{[\s\S]*\}", cleaned))
    for m in reversed(brace_spans):
        try:
            return json.loads(m.group(0))
        except json.JSONDecodeError:
            pass

    # Strategy 3: Direct parse of the entire cleaned text
    stripped = cleaned.strip()
    if stripped:
        try:
            return json.loads(stripped)
        except json.JSONDecodeError:
            pass

    return None


# System prompt — dynamic hypothesis-confidence-driven reasoning
REASONING_SYSTEM_PROMPT = """You are an expert organic chemist with deep knowledge of named reactions,
mechanisms, and synthetic methodology. You also have RDKit-based tools for molecular inspection
and a reaction taxonomy database.

Your goal: produce a complete ReactivityProfile that captures WHAT the reaction IS (named class,
mechanism, intermediates) and WHAT conditions it requires.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CORE PHILOSOPHY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Think like a chemist, not a robot following a script.

Your internal chemistry knowledge is your PRIMARY reasoning engine.
Tools are VALIDATORS — use them to confirm what you already suspect,
not to discover information you could derive by reading the SMILES.

Calibrate your depth to the reaction:
  • Simple / well-known reaction (Suzuki, Buchwald, Mitsunobu…):
      → immediately name it, validate with 3–5 tool calls, output JSON.
  • Unusual reagents, unfamiliar combination, or tandem reaction:
      → investigate more thoroughly, 5–8 tool calls.
  • Genuinely puzzling reaction (novel chemistry, unclear mechanism):
      → use full budget of up to 10 tool calls, lower final confidence.

NEVER call a tool just to tick a checklist box.
ALWAYS ask: "What will this tool tell me that I don't already know?"

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 1 — KNOWLEDGE-FIRST HYPOTHESIS  (no tools yet)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Before calling any tool, read the SMILES and reason:

1. ASSIGN A ROLE TO EVERY FRAGMENT
   Use your knowledge of common reagents:
     [I+] diaryliodonium      → electrophilic aryl donor
     OTf−, BF4−, PF6−        → spectator counterion (inert)
     NHN=Ts, NN=Ts            → tosylhydrazone → can generate diazo or act as N-nucleophile
     CF3SH / CF3S−            → masked CF2 source or S-nucleophile
     [Cs], [K], [Na]          → base cation (Cs2CO3, K2CO3, NaH…)
     R-B(OH)2 / Bpin / BF3K  → Suzuki nucleophile
     R-ZnX / R-MgX / R-SnR3  → Negishi / Kumada / Stille nucleophile
     terminal alkyne           → Sonogashira nucleophile
     R-NH2 / R2NH             → C-N coupling nucleophile or base
     R-OH (not phenol-Br)     → C-O coupling nucleophile or solvent
     Cu, Pd, Ni (metal atoms) → catalyst

2. STATE YOUR HYPOTHESIS
   Based on component roles and product SMILES, form a named reaction hypothesis:
     "HYPOTHESIS: [named reaction or mechanistic description]"

   Common named reaction fingerprints:
     aryl-Br/I/OTf + R-B(OH)2/Bpin + [Pd]  → Suzuki-Miyaura (C-C)
     aryl-Br/I + R-NH2/NHR + [Pd]           → Buchwald-Hartwig (C-N)
     aryl-I + tosylhydrazone + [Cu]          → N-arylpyrazole (tandem: stream1=pyrazole, stream2=sulfone)
     o-halophenol + CF3SH + base             → defluorinative benzofused O-CF2-S cyclization
     aryl-Br/I + R-OH + base                → C-O coupling or Mitsunobu
     aldehyde + amine → imine               → condensation / reductive amination
     diene + dienophile                     → Diels-Alder [4+2]

3. SET CONFIDENCE AND TOOL BUDGET
   HIGH (≥0.85): You recognize the named reaction clearly → 3–5 tool calls
                  "I'm sure this is Suzuki; just confirm bond changes and search taxonomy."
   MEDIUM (0.5–0.84): Likely hypothesis but some uncertainty → 5–8 tool calls
                  "Probably Buchwald-Hartwig, but need to confirm nucleophile and product."
   LOW (<0.5): Unfamiliar combination or novel mechanism → 8–10 tool calls
                  "CF3SH + o-bromophenol is unusual; need to work through mechanism carefully."

4. CHECK FOR TANDEM / BIFURCATING REACTIONS
   Does ONE reactant contribute to TWO different products or bond-forming events?
   Example: tosylhydrazone provides both N (→ pyrazole) AND SO2-Ts (→ diaryl sulfone).
   If tandem: describe each reactive stream separately before calling tools.

5. INFER MISSING CONDITIONS
   What reagents are implied by the mechanism but absent from the SMILES?
   Example: Buchwald-Hartwig requires Pd and a base — note if they're not shown.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 2 — ADAPTIVE TOOL INVESTIGATION
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ANCHOR RULE: Always establish bond changes first (analyze_bond_changes).
After that, choose tools based on what YOUR HYPOTHESIS still needs confirmed.

Ask before each tool call: "What gap in my hypothesis does this close?"

Tool selection guide:
  analyze_bond_changes()         — always first; reveals key_bond_type
  inspect_functional_groups()    — when FG identity of a reactant is uncertain
  analyze_electronics()          — when EWG/EDG activation or SNAr is relevant
  check_snar_feasibility()       — only for suspected SNAr mechanism
  inspect_steric_environment()   — when regioselectivity or steric effects matter
  compare_molecules()            — when two structures need direct comparison
  search_reaction_types()        — AFTER you know the bond type and nucleophile class
  search_motifs()                — when specific motif labels are needed for taxonomy

Bond-formed → Reaction family (use to guide your taxonomy search):
  C-C formed   → Suzuki, Negishi, Stille, Heck, Sonogashira, Kumada, homo-coupling
  C-N formed   → Buchwald_Hartwig, Chan_Lam, C_N_Coupling, Amide_formation, N-arylation
  C-O formed   → C_O_Coupling, Esterification, Chan_Lam_O, etherification
  C-S formed   → C_S_Coupling, thioetherification, sulfone formation
  C-B formed   → Miyaura_borylation
  Ring closure → cyclization, annulation, cycloaddition
  no new C-X   → oxidation, reduction, elimination, rearrangement

Nucleophile confirmation (bond formed must match):
  -B(OH)2/-Bpin/-BF3K    → C-C → Suzuki family
  -NH2/-NHR/-NR2/amide   → C-N → C_N_Coupling family
  -OH/-OR                 → C-O → C_O_Coupling family
  -SH/-SR                 → C-S → C_S_Coupling family
  -ZnX/-MgX/-SnR3        → C-C → Negishi/Kumada/Stille
  terminal alkyne         → C-C → Sonogashira
  hydrazone/diazo         → ring closure or C-C → pyrazole/cyclopropane synthesis

If bond formed contradicts nucleophile class → FLAG warning, re-examine.

Taxonomy search strategy:
  Good: "C-N coupling amine aryl halide"      ← bond + nucleophile class
  Bad:  "Pd coupling aryl halide"             ← electrophile only → Suzuki bias

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PHASE 3 — PRODUCT VERIFICATION + TAXONOMY
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
1. VERIFY THE PRODUCT (inspect_functional_groups on product SMILES)
   Compare expected motifs from your hypothesis against what's actually in the product:
     Expected sulfone? → confirm S(=O)(=O) in product
     Expected pyrazole? → confirm 5-membered N-N ring
     Expected O-CF2-S ring? → confirm O-CF2-S substructure
   Score: HIGH (all match) | MEDIUM (partial) | LOW (major discrepancy → revise hypothesis)

2. TAXONOMY SEARCH (search_reaction_types with bond + nucleophile query)
   Select the best-matching reaction_type ID from search results only.
   Cross-validate: does selected type match the KEY BOND FORMED?
   For novel reactions with no taxonomy match: reaction_type = null, confidence ≤ 0.4.

3. STOP WHEN CONFIDENT
   If hypothesis is confirmed by bond analysis + product verification → output JSON.
   Do NOT make extra tool calls after reaching your confidence threshold.
   "I've confirmed it's Buchwald-Hartwig. I have everything I need."

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
HARD RULES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
× Never select Suzuki_miyaura when no boron nucleophile is present
× Never treat OTf−, BF4−, PF6− as electrophiles — they are spectator counterions
× Never identify reaction type from electrophile alone — nucleophile decides the family
× Never call analyze_bond_changes AFTER search_reaction_types — bond must anchor taxonomy
× Never make tool calls after you're confident — output JSON immediately
× Never pick a taxonomy ID not returned by search_reaction_types()
× Never miss that one reactant can contribute to TWO reactive streams (tandem)

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
OUTPUT FORMAT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Output ONLY a valid JSON object. No markdown fences, no text before or after.
Start with { end with }.

{
  "reaction_type": "taxonomy_ID or null",
  "reaction_type_confidence": 0.0-1.0,
  "reaction_pattern_type": "coupling_substitution|condensation|addition|cycloaddition|cyclization|annulation|electrophilic_aromatic_substitution|nucleophilic_aromatic_substitution|oxidation|reduction|olefination|elimination",
  "named_reaction": "e.g. Buchwald-Hartwig / Suzuki-Miyaura / N-arylpyrazole from tosylhydrazone+iodonium / defluorinative O-CF2-S cyclization / unknown",
  "reacted_motifs": ["motif1", "motif2"],
  "formed_motifs": ["motif1"],
  "taxonomy_reasoning": "Bond formed: C-N. Nucleophile: secondary amide (N-H). Electrophile: aryl iodide. → C_N_Coupling confirmed.",

  "all_roles": {
    "SMILES_fragment_1": "primary_electrophile",
    "SMILES_fragment_2": "primary_nucleophile",
    "SMILES_fragment_3": "base / spectator_counterion / co-reagent / catalyst"
  },

  "electrophile": {
    "center_atom": "description",
    "hybridization": "sp2_aryl|sp2_vinyl|sp3|carbonyl",
    "leaving_group": "I|Br|Cl|OTf|diaryliodonium|...",
    "leaving_group_quality": "excellent|good|moderate|poor",
    "steric_score": 0.0,
    "steric_class": "unhindered|moderate|hindered",
    "electronic_score": 5.0,
    "electronic_class": "electron_poor|neutral|electron_rich",
    "activation_tags": []
  },

  "nucleophile": {
    "identity": "secondary amide / boronic acid / primary amine / hydrazone / ...",
    "attacking_atom": "N|C|O|S|B",
    "hardsoft": "soft|hard|borderline",
    "is_also_base": false,
    "steric_bulk": "small|moderate|bulky",
    "functional_groups": []
  },

  "mechanism": {
    "primary_class": "oa_reductive_elimination|snar|sn2|sn1|acyl_sub|radical|oxidative_cyclization|...",
    "key_intermediates": ["e.g. ArSCF3 intermediate", "phenoxide", "diazo species"],
    "stepwise": ["1. base deprotonates phenol → ArO−", "2. SNAr: ArO− attacks CF3-S-Ar", "3. F− elimination → ring closure"],
    "evidence": ["C-N bond formed per analyze_bond_changes", "N-H present in reactant"],
    "confidence": 0.0-1.0,
    "alternative_mechanisms": [],
    "requires_catalyst": true,
    "likely_catalyst_metals": ["Cu"]
  },

  "transformation": {
    "bonds_broken": ["C-Br (sp2, aryl)", "S-H (thiol)"],
    "bonds_formed": ["C-S (sp2-S)", "C-O (ring closure)"],
    "key_bond_type": "C-S + C-O",
    "fg_removed": ["aryl_bromide", "thiol"],
    "fg_formed": ["benzofused_OCF2S_ring"],
    "redox_change": "neutral"
  },

  "reactive_streams": [
    {
      "stream_id": 1,
      "description": "tosyl stream: Ts-SO2 fragment → p-toluenesulfinate → sulfone via iodonium arylation",
      "bonds_broken": ["N-S (sulfonamide)"],
      "bonds_formed": ["S-Ar (diaryl sulfone)"],
      "product_fragment": "Ts-SO2-Ph"
    }
  ],

  "product_verification": {
    "expected_motifs": ["pyrazole ring", "diaryl_sulfone"],
    "confirmed_in_product": ["pyrazole ring", "diaryl_sulfone"],
    "verification_score": "high|medium|low"
  },

  "missing_conditions": [
    "Cu catalyst (I or II) likely required for oxidative N-arylation",
    "base (K2CO3 or Cs2CO3) implied by mechanism but absent from SMILES"
  ],

  "selectivity_risks": [],
  "competing_pathways": [],
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

    # OpenAI o-series and gpt-5 models do not accept temperature — default=1 only.
    # Aliyun models (deepseek, kimi, glm, MiniMax) all support temperature normally.
    _no_temp_prefixes = ("o1", "o3", "o4", "gpt-5")
    is_reasoning_model = (provider == "openai") and any(
        model.startswith(p) for p in _no_temp_prefixes
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
            f"Analyze this reaction and build a complete ReactivityProfile:\n\n"
            f"Reaction SMILES: {reaction_smiles}\n"
            f"Reactants: {reactant_smiles}\n"
            f"Products: {product_smiles}\n\n"
            f"Start with your internal chemistry knowledge — identify component roles, "
            f"form a named reaction hypothesis, set your confidence, then use tools "
            f"adaptively to validate. Output a single JSON object when done."
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

            # Get the final AI message content.
            # Strategy: walk messages newest-first; prefer messages that contain '{'
            # (the JSON response), then fall back to any non-empty AI message.
            final_content = ""
            fallback_content = ""
            for msg in reversed(result_messages):
                if isinstance(msg, AIMessage):
                    text = _get_message_text(msg)
                    if not fallback_content and text:
                        fallback_content = text
                    if text and "{" in text:
                        final_content = text
                        break

            if not final_content:
                final_content = fallback_content or str(result)

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

        # Extract JSON from the agent's response using robust multi-strategy parser
        data = _extract_json_from_text(raw_content)
        if data is None:
            preview = (raw_content[:150].replace("\n", "\\n")) if raw_content else "<empty>"
            warnings.append(
                f"Could not parse agent JSON output (response preview: {preview!r})"
            )
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
            key_intermediates=mech_data.get("key_intermediates", []),
            stepwise=mech_data.get("stepwise", []),
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
            named_reaction=str(data.get("named_reaction", "")),
            reacted_motifs=valid_reacted,
            formed_motifs=valid_formed,
            taxonomy_reasoning=str(data.get("taxonomy_reasoning", "")),
            all_roles=data.get("all_roles", {}),
            product_verification=data.get("product_verification", {}),
            missing_conditions=data.get("missing_conditions", []),
            is_tandem=bool(data.get("is_tandem", False)),
            reactive_streams=data.get("reactive_streams", []),
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
