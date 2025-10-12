"""
LLM Tools - Advanced Language Model Operations
==============================================

This package provides utilities and tools for integrating large language models
into chemistry workflows, including:

- Multi-provider LLM client management (OpenAI, Aliyun/DeepSeek, etc.)
- Chemistry-specific prompt templates and agents
- Structured output parsing and validation
- Reasoning chains for complex chemistry tasks
- Integration with chemtools for enhanced capabilities

Modules:
    clients: LLM client management and configuration
    prompts: Chemistry-specific prompt templates
    agents: Intelligent agents for chemistry tasks
    parsers: Output parsing and structured extraction
    reasoning: Multi-step reasoning chains
    cache: Response caching and optimization

Usage:
    from llmtools import LLMClient, ChemistryAgent
    
    client = LLMClient(provider="openai", model="gpt-4o")
    agent = ChemistryAgent(client)
    
    result = agent.suggest_conditions(
        reaction="Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
        reaction_type="Buchwald_CN"
    )

Version: 0.1.0
Author: Synthesis-Automation Team
"""

__version__ = "0.1.0"

from llmtools.clients import LLMClient, build_client
from llmtools.agents import ChemistryAgent
from llmtools.reagent_review import (
    LLMReviewOptions,
    build_review_prompt,
    review_taxonomy_proposal,
)

__all__ = [
    "LLMClient",
    "build_client",
    "ChemistryAgent",
    "LLMReviewOptions",
    "build_review_prompt",
    "review_taxonomy_proposal",
    "__version__",
]
