"""
LLM Tools - Advanced Language Model Operations
==============================================

This package provides optional language-model client and review helpers:

- Multi-provider LLM client management (OpenAI, Aliyun/DeepSeek, etc.)
- Chemistry-specific prompt templates
- Structured output parsing and validation
- Reasoning chains for complex chemistry tasks

Modules:
    clients: LLM client management and configuration
    prompts: Chemistry-specific prompt templates
    parsers: Output parsing and structured extraction
    reasoning: Multi-step reasoning chains
    cache: Response caching and optimization

It does not own chemistry analysis or recommendation routing.

Version: 0.1.0
Author: Synthesis-Automation Team
"""

__version__ = "0.1.0"

from llmtools.clients import LLMClient, build_client
from llmtools.reagent_review import (
    LLMReviewOptions,
    build_review_prompt,
    review_taxonomy_proposal,
)

__all__ = [
    "LLMClient",
    "build_client",
    "LLMReviewOptions",
    "build_review_prompt",
    "review_taxonomy_proposal",
    "__version__",
]
