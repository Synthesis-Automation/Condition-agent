"""
Chemistry-Specific LLM Agents
==============================

Intelligent agents that combine LLM capabilities with chemtools functions
for advanced chemistry tasks:

- ChemistryAgent: General-purpose chemistry assistant
- ConditionOptimizer: Optimize reaction conditions
- RetrosynthesisPlanner: Plan synthetic routes
- ReactionAnalyzer: Analyze and explain reactions
- LiteratureAgent: Search and summarize literature

Agents can:
- Use chemtools functions for deterministic calculations
- Call LLMs for reasoning and generation
- Chain multiple steps together
- Validate outputs against chemistry rules
"""

from typing import Any, Dict, List, Optional, Union

from llmtools.clients import LLMClient, LLMResponse
from llmtools.prompts import get_template


class ChemistryAgent:
    """
    General-purpose chemistry agent combining LLM and chemtools.
    
    This agent can:
    - Recommend reaction conditions
    - Explain mechanisms
    - Suggest alternatives
    - Troubleshoot problems
    - Generate protocols
    
    Examples:
        from llmtools import LLMClient, ChemistryAgent
        
        client = LLMClient(provider="openai", model="gpt-4o")
        agent = ChemistryAgent(client)
        
        # Recommend conditions
        result = agent.suggest_conditions(
            reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
            reaction_type="Buchwald_CN"
        )
        
        # Explain mechanism
        explanation = agent.explain_mechanism(
            reaction="...",
            reaction_type="Suzuki"
        )
    """
    
    def __init__(
        self,
        client: LLMClient,
        use_chemtools: bool = True,
        verbose: bool = False,
    ):
        """
        Initialize chemistry agent.
        
        Args:
            client: LLM client for reasoning
            use_chemtools: Whether to integrate chemtools functions
            verbose: Print intermediate steps
        """
        self.client = client
        self.use_chemtools = use_chemtools
        self.verbose = verbose
        
        # Import chemtools if available and requested
        if use_chemtools:
            try:
                from chemtools import chem
                self.chem = chem
                self._has_chemtools = True
            except ImportError:
                if verbose:
                    print("Warning: chemtools not available, running in LLM-only mode")
                self._has_chemtools = False
        else:
            self._has_chemtools = False
    
    def suggest_conditions(
        self,
        reaction: str,
        reaction_type: Optional[str] = None,
        context: Optional[str] = None,
        temperature: float = 0.7,
    ) -> Dict[str, Any]:
        """
        Suggest reaction conditions using LLM + chemtools precedents.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family (optional, will auto-detect)
            context: Additional context or constraints
            temperature: LLM sampling temperature
            
        Returns:
            Dict with conditions, rationale, and metadata
        """
        # Step 1: Get precedents from chemtools if available
        precedents_info = ""
        if self._has_chemtools:
            try:
                # Get ML-based precedents
                ml_result = self.chem.recommend.conditions(
                    reaction=reaction,
                    reaction_type=reaction_type,
                    k=5,
                    limit=3,
                )
                
                # Extract condition patterns
                if "recommended_conditions" in ml_result:
                    precedents_info = "\n\n**Precedents from Database**:\n"
                    for i, cond in enumerate(ml_result["recommended_conditions"][:3], 1):
                        precedents_info += f"\n{i}. Similarity: {cond.get('confidence', 0):.2f}\n"
                        if "conditions" in cond:
                            c = cond["conditions"]
                            precedents_info += f"   - Catalyst: {c.get('catalyst', 'N/A')}\n"
                            precedents_info += f"   - Ligand: {c.get('ligand', 'N/A')}\n"
                            precedents_info += f"   - Base: {c.get('base', 'N/A')}\n"
                            precedents_info += f"   - Solvent: {c.get('solvent', 'N/A')}\n"
                            if c.get('temperature'):
                                precedents_info += f"   - Temp: {c['temperature']}\n"
            except Exception as e:
                if self.verbose:
                    print(f"Could not fetch precedents: {e}")
        
        # Step 2: Use LLM to generate recommendations
        prompt_template = get_template("condition_recommendation")
        prompt = prompt_template.format(
            reaction_smiles=reaction,
            reaction_type=reaction_type or "Unknown (please identify)",
            context=(context or "None") + precedents_info,
        )
        
        if self.verbose:
            print("Generating LLM recommendations...")
        
        response = self.client.chat(
            prompt=prompt,
            temperature=temperature,
        )
        
        return {
            "reaction": reaction,
            "reaction_type": reaction_type,
            "llm_response": response.content,
            "precedents_used": bool(precedents_info),
            "model": response.model,
            "provider": response.provider,
            "tokens": response.total_tokens,
            "latency_ms": response.latency_ms,
        }
    
    def explain_mechanism(
        self,
        reaction: str,
        reaction_type: str,
        reagents: Optional[str] = None,
        detail_level: str = "Detailed with all intermediates",
    ) -> Dict[str, Any]:
        """
        Explain reaction mechanism.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            reagents: Catalyst and reagents used
            detail_level: Level of detail in explanation
            
        Returns:
            Dict with mechanism explanation and metadata
        """
        prompt_template = get_template("mechanism")
        prompt = prompt_template.format(
            reaction_smiles=reaction,
            reaction_type=reaction_type,
            reagents=reagents or "Not specified",
            detail_level=detail_level,
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.5)
        
        return {
            "reaction": reaction,
            "reaction_type": reaction_type,
            "mechanism_explanation": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }
    
    def troubleshoot_reaction(
        self,
        reaction: str,
        reaction_type: str,
        problem: str,
        current_conditions: Optional[str] = None,
        observed_result: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Help troubleshoot a problematic reaction.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            problem: Description of the problem
            current_conditions: Current reaction conditions
            observed_result: What was observed
            
        Returns:
            Dict with troubleshooting advice
        """
        prompt_template = get_template("troubleshooting")
        prompt = prompt_template.format(
            problem_description=problem,
            reaction_smiles=reaction,
            reaction_type=reaction_type,
            current_conditions=current_conditions or "Not specified",
            observed_result=observed_result or "Poor yield or no conversion",
            expected_result="High yield of desired product",
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.7)
        
        return {
            "reaction": reaction,
            "problem": problem,
            "troubleshooting_advice": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }
    
    def generate_protocol(
        self,
        reaction: str,
        reaction_type: str,
        conditions: Optional[Dict[str, Any]] = None,
        scale: str = "1 mmol",
    ) -> Dict[str, Any]:
        """
        Generate detailed experimental protocol.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            conditions: Recommended conditions dict
            scale: Reaction scale
            
        Returns:
            Dict with detailed protocol
        """
        # Format conditions if provided
        conditions_str = "Not specified"
        if conditions:
            conditions_str = "\n".join([f"- {k}: {v}" for k, v in conditions.items()])
        
        prompt_template = get_template("protocol")
        prompt = prompt_template.format(
            reaction_smiles=reaction,
            reaction_type=reaction_type,
            conditions=conditions_str,
            scale=scale,
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.3)
        
        return {
            "reaction": reaction,
            "reaction_type": reaction_type,
            "scale": scale,
            "protocol": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }
    
    def search_literature(
        self,
        reaction: str,
        reaction_type: str,
        focus: str = "General conditions and scope",
        year_range: str = "Last 10 years",
    ) -> Dict[str, Any]:
        """
        Search and summarize relevant literature.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            focus: What to focus on
            year_range: Time range for search
            
        Returns:
            Dict with literature summary
        """
        prompt_template = get_template("literature")
        prompt = prompt_template.format(
            reaction_smiles=reaction,
            reaction_type=reaction_type,
            focus=focus,
            year_range=year_range,
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.5)
        
        return {
            "reaction": reaction,
            "reaction_type": reaction_type,
            "literature_summary": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }
    
    def assess_safety(
        self,
        reaction: str,
        reagents: str,
        solvents: str,
        scale: str = "Laboratory (mg-g scale)",
    ) -> Dict[str, Any]:
        """
        Assess safety considerations.
        
        Args:
            reaction: Reaction SMILES
            reagents: Reagents used
            solvents: Solvents used
            scale: Reaction scale
            
        Returns:
            Dict with safety assessment
        """
        prompt_template = get_template("safety")
        prompt = prompt_template.format(
            reaction_smiles=reaction,
            reagents=reagents,
            solvents=solvents,
            scale=scale,
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.3)
        
        return {
            "reaction": reaction,
            "safety_assessment": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }
    
    def get_usage_summary(self) -> Dict[str, Any]:
        """Get LLM usage summary."""
        return self.client.get_usage_summary()


class ConditionOptimizer(ChemistryAgent):
    """
    Specialized agent for optimizing reaction conditions.
    
    Uses iterative LLM reasoning + chemtools validation.
    """
    
    def optimize(
        self,
        reaction: str,
        reaction_type: str,
        constraints: Optional[Dict[str, Any]] = None,
        max_iterations: int = 3,
    ) -> Dict[str, Any]:
        """
        Iteratively optimize reaction conditions.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            constraints: Optimization constraints
            max_iterations: Maximum optimization rounds
            
        Returns:
            Dict with optimized conditions and reasoning chain
        """
        # TODO: Implement iterative optimization loop
        # 1. Get initial recommendations from LLM + precedents
        # 2. Validate against constraints
        # 3. Refine based on feedback
        # 4. Repeat until converged or max iterations
        
        return self.suggest_conditions(reaction, reaction_type)


class RetrosynthesisPlanner(ChemistryAgent):
    """
    Specialized agent for retrosynthetic planning.
    
    Combines LLM reasoning with chemtools validation.
    """
    
    def plan_route(
        self,
        target_smiles: str,
        max_steps: int = 5,
        availability: str = "commercially available",
    ) -> Dict[str, Any]:
        """
        Plan retrosynthetic route.
        
        Args:
            target_smiles: Target molecule SMILES
            max_steps: Maximum synthetic steps
            availability: Starting material availability
            
        Returns:
            Dict with retrosynthetic tree and forward route
        """
        prompt_template = get_template("retrosynthesis")
        prompt = prompt_template.format(
            target_smiles=target_smiles,
            max_steps=max_steps,
            availability=availability,
        )
        
        response = self.client.chat(prompt=prompt, temperature=0.7)
        
        return {
            "target": target_smiles,
            "retrosynthesis": response.content,
            "model": response.model,
            "tokens": response.total_tokens,
        }


class ReactionAnalyzer(ChemistryAgent):
    """
    Specialized agent for analyzing reactions.
    
    Provides detailed mechanism, literature, and safety info.
    """
    
    def analyze(
        self,
        reaction: str,
        reaction_type: str,
        include_mechanism: bool = True,
        include_literature: bool = True,
        include_safety: bool = True,
    ) -> Dict[str, Any]:
        """
        Comprehensive reaction analysis.
        
        Args:
            reaction: Reaction SMILES
            reaction_type: Reaction family
            include_mechanism: Include mechanism explanation
            include_literature: Include literature search
            include_safety: Include safety assessment
            
        Returns:
            Dict with comprehensive analysis
        """
        analysis = {
            "reaction": reaction,
            "reaction_type": reaction_type,
        }
        
        if include_mechanism:
            mech = self.explain_mechanism(reaction, reaction_type)
            analysis["mechanism"] = mech["mechanism_explanation"]
        
        if include_literature:
            lit = self.search_literature(reaction, reaction_type)
            analysis["literature"] = lit["literature_summary"]
        
        if include_safety:
            safety = self.assess_safety(reaction, "TBD", "TBD")
            analysis["safety"] = safety["safety_assessment"]
        
        return analysis
