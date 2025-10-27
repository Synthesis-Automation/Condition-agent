"""
LangGraph ReAct agent for ChemTools.

This module provides a ReAct (Reasoning + Acting) agent that can use
ChemTools functions to analyze reactions, recommend conditions, and
answer chemistry questions.

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

import os
from typing import List, Optional
from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, AIMessage, BaseMessage
from langgraph.prebuilt import create_react_agent
from dotenv import load_dotenv

from .chemtools_wrapper import CHEMTOOLS_TOOLS


# Load environment variables
load_dotenv()


# ============================================================================
# LLM Configuration
# ============================================================================

DEFAULT_ALIYUN_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"


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
6. **recommend_conditions_tool**: Get ML-based condition recommendations
7. **search_precedents_tool**: Find similar precedent reactions
8. **find_reagent_tool**: Look up reagent information from database

**How to help users:**

For reaction condition recommendations:
1. First normalize the reaction SMILES if needed
2. Detect the reaction family to understand the reaction type
3. Use recommend_conditions_tool to get comprehensive recommendations
4. Optionally search for precedents to provide context
5. Explain the recommendations in clear, practical terms

For reagent questions:
1. Use find_reagent_tool to look up properties and roles
2. Provide CAS numbers, structures, and usage information

For reactant/molecule analysis:
1. Normalize SMILES first
2. Use classify_reactant_tool or get_functional_groups_tool
3. Explain the structural features clearly

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


# ============================================================================
# ChemTools Agent Class
# ============================================================================

class ChemToolsAgent:
    """
    ReAct agent for chemistry tasks using ChemTools.
    
    This agent combines LLM reasoning with deterministic ChemTools functions
    to provide accurate chemistry analysis and recommendations.
    
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
        verbose: bool = False
    ):
        """
        Initialize ChemTools agent.
        
        Args:
            provider: LLM provider ("openai" or "aliyun")
            model: Model name
            temperature: Sampling temperature (0.0-1.0)
            system_prompt: Custom system prompt (default: CHEMISTRY_SYSTEM_PROMPT)
            verbose: Print debug information
        """
        self.llm = get_llm_client(provider, model, temperature)
        self.system_prompt = system_prompt or CHEMISTRY_SYSTEM_PROMPT
        self.verbose = verbose
        
        # Create ReAct agent with tools
        self.agent = create_react_agent(
            self.llm,
            CHEMTOOLS_TOOLS,
            prompt=self.system_prompt
        )
    
    def run(
        self,
        query: str,
        history: Optional[List[BaseMessage]] = None,
        recursion_limit: int = 15
    ) -> str:
        """
        Run the agent on a query.
        
        Args:
            query: User question or task
            history: Conversation history (list of messages)
            recursion_limit: Maximum reasoning steps
        
        Returns:
            Agent's response text
        
        Example:
            >>> agent = ChemToolsAgent()
            >>> response = agent.run("What conditions for Suzuki coupling?")
        """
        if history is None:
            history = []
        
        try:
            # Invoke agent with query and history
            result = self.agent.invoke(
                {"messages": history + [HumanMessage(content=query)]},
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
    
    def chat(
        self,
        query: str,
        history: List[BaseMessage],
        recursion_limit: int = 15
    ) -> tuple[str, List[BaseMessage]]:
        """
        Chat with the agent and maintain conversation history.
        
        Args:
            query: User question
            history: Current conversation history
            recursion_limit: Maximum reasoning steps
        
        Returns:
            Tuple of (response, updated_history)
        
        Example:
            >>> agent = ChemToolsAgent()
            >>> history = []
            >>> response, history = agent.chat("Normalize c1ccccc1", history)
            >>> response, history = agent.chat("What functional groups?", history)
        """
        response = self.run(query, history, recursion_limit)
        
        # Update history
        updated_history = history + [
            HumanMessage(content=query),
            AIMessage(content=response)
        ]
        
        return response, updated_history


# ============================================================================
# Convenience Functions
# ============================================================================

def create_agent(
    provider: Optional[str] = None,
    model: Optional[str] = None,
    **kwargs
) -> ChemToolsAgent:
    """
    Convenience function to create a ChemTools agent.
    
    Args:
        provider: LLM provider ("openai" or "aliyun")
        model: Model name
        **kwargs: Additional arguments for ChemToolsAgent
    
    Returns:
        Configured ChemToolsAgent instance
    """
    return ChemToolsAgent(provider=provider, model=model, **kwargs)


def quick_query(query: str, **kwargs) -> str:
    """
    Quick one-shot query without maintaining state.
    
    Args:
        query: Question or task
        **kwargs: Arguments for ChemToolsAgent
    
    Returns:
        Agent response
    
    Example:
        >>> from lang_chain.chemtools_agent import quick_query
        >>> result = quick_query("Recommend conditions for Suzuki coupling")
    """
    agent = ChemToolsAgent(**kwargs)
    return agent.run(query)


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
