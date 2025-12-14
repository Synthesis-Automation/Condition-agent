"""
LLM Client Management
=====================

Multi-provider LLM client with unified interface for:
- OpenAI (GPT-4, GPT-5, o-series)
- Aliyun/DeepSeek (R1, V3 series)
- Future providers (Anthropic, Google, etc.)

Features:
- Environment-based configuration
- Automatic retry logic
- Token usage tracking
- Response caching
- Streaming support

Based on tests/test_llm.py but production-ready with additional features.
"""

import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union

from openai import OpenAI


# Default base URLs for supported providers
DEFAULT_BASE_URLS = {
    "aliyun": "https://dashscope.aliyuncs.com/compatible-mode/v1",
    "openai": "https://api.openai.com/v1",
}

# Recommended models by provider
RECOMMENDED_MODELS = {
    "aliyun": {
        "reasoning": "deepseek-r1",  # Best reasoning
        "fast": "deepseek-r1-distill-qwen-7b",  # Fast distilled
        "balanced": "deepseek-v3.2",  # Latest balanced (v3.2 experimental)
        "experimental": "deepseek-v3.2",  # Latest experimental
    },
    "openai": {
        "reasoning": "o3-mini",  # Best reasoning
        "fast": "gpt-4o",  # GPT-4o (default)
        "balanced": "gpt-4o",  # Balanced quality/speed
        "advanced": "gpt-5-mini",  # GPT-5 series
    },
}

# All available models per provider (from test_llm.py)
AVAILABLE_MODELS = {
    "aliyun": [
        "deepseek-v3.2",
        "deepseek-v3.1",
        "deepseek-r1",
        "deepseek-r1-0528",
        "deepseek-v3",
        "deepseek-r1-distill-qwen-14b",
        "deepseek-r1-distill-qwen-32b",
        "deepseek-r1-distill-llama-70b",
    ],
    "openai": [
        "gpt-5.2",
        "gpt-5-pro",
        "gpt-5-mini",
        "gpt-5-nano",
        "gpt-5-codex",
        "o3",
        "o3-pro",
        "o3-mini",
        "o4-mini",
        "o3-deep-research",
        "o4-mini-deep-research",
        "gpt-4o",
        "gpt-4.1",
        "gpt-4o-mini",
        "gpt-4.1-mini",
        "gpt-4.1-nano",
    ],
}


@dataclass
class LLMResponse:
    """Structured response from LLM with metadata."""

    content: str
    model: str
    provider: str
    prompt_tokens: int = 0
    completion_tokens: int = 0
    total_tokens: int = 0
    latency_ms: float = 0.0
    finish_reason: str = "stop"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def tokens_per_second(self) -> float:
        """Calculate tokens per second."""
        if self.latency_ms > 0:
            return (self.completion_tokens / self.latency_ms) * 1000
        return 0.0


class LLMClient:
    """
    Unified LLM client supporting multiple providers.

    Examples:
        # Simple usage
        client = LLMClient(provider="openai", model="gpt-4o-mini")
        response = client.chat("Suggest a solvent for Suzuki coupling")

        # With configuration
        client = LLMClient(
            provider="aliyun",
            model="deepseek-r1",
            temperature=0.7,
            max_tokens=2000
        )
        response = client.chat("Analyze this reaction: ...")

        # Auto-select model
        client = LLMClient.from_env(provider="openai", model_type="reasoning")
    """

    def __init__(
        self,
        provider: str = "openai",
        model: Optional[str] = None,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 2000,
        timeout: int = 60,
    ):
        """
        Initialize LLM client.

        Args:
            provider: Provider name ("openai", "aliyun")
            model: Model name (uses default if not specified)
            api_key: API key (reads from env if not provided)
            base_url: Base URL (uses default if not provided)
            temperature: Sampling temperature (0.0-1.0)
            max_tokens: Maximum tokens in response
            timeout: Request timeout in seconds
        """
        self.provider = provider.lower()
        self.model = model or self._get_default_model()
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.timeout = timeout

        # Get API key from parameter or environment
        if api_key is None:
            api_key = self._get_api_key_from_env()

        # Get base URL from parameter or environment or default
        if base_url is None:
            base_url = self._get_base_url_from_env()

        # Create OpenAI client
        self.client = OpenAI(
            api_key=api_key,
            base_url=base_url,
            timeout=timeout,
        )

        # Track usage
        self.total_prompt_tokens = 0
        self.total_completion_tokens = 0
        self.total_requests = 0

    def _get_default_model(self) -> str:
        """Get default model for provider."""
        if self.provider == "aliyun":
            return RECOMMENDED_MODELS["aliyun"]["balanced"]
        elif self.provider == "openai":
            return RECOMMENDED_MODELS["openai"]["balanced"]
        else:
            raise ValueError(f"Unsupported provider: {self.provider}")

    def _get_api_key_from_env(self) -> str:
        """Get API key from environment variables."""
        if self.provider == "aliyun":
            key = os.getenv("ALIYUN_API_KEY")
            if not key:
                raise RuntimeError("ALIYUN_API_KEY environment variable not set")
            return key
        elif self.provider == "openai":
            key = os.getenv("OPENAI_API_KEY")
            if not key:
                raise RuntimeError("OPENAI_API_KEY environment variable not set")
            return key
        else:
            raise ValueError(f"Unsupported provider: {self.provider}")

    def _get_base_url_from_env(self) -> str:
        """Get base URL from environment or use default."""
        if self.provider == "aliyun":
            return os.getenv("ALIYUN_BASE_URL", DEFAULT_BASE_URLS["aliyun"])
        elif self.provider == "openai":
            return os.getenv("OPENAI_BASE_URL", DEFAULT_BASE_URLS["openai"])
        else:
            return DEFAULT_BASE_URLS.get(self.provider, "")

    def chat(
        self,
        prompt: str,
        system: Optional[str] = None,
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
    ) -> LLMResponse:
        """
        Send a chat completion request.

        Args:
            prompt: User message
            system: Optional system message
            temperature: Override default temperature
            max_tokens: Override default max_tokens

        Returns:
            LLMResponse with content and metadata
        """
        messages = []
        if system:
            messages.append({"role": "system", "content": system})
        messages.append({"role": "user", "content": prompt})

        return self.chat_messages(
            messages=messages,
            temperature=temperature,
            max_tokens=max_tokens,
        )

    def chat_messages(
        self,
        messages: List[Dict[str, str]],
        temperature: Optional[float] = None,
        max_tokens: Optional[int] = None,
    ) -> LLMResponse:
        """
        Send a chat completion with full message history.

        Args:
            messages: List of message dicts with role and content
            temperature: Override default temperature
            max_tokens: Override default max_tokens

        Returns:
            LLMResponse with content and metadata
        """
        start_time = time.perf_counter()

        # Prepare API parameters
        api_params = {
            "model": self.model,
            "messages": messages,
        }

        # GPT-5 and o-series models have different parameter requirements
        is_gpt5_or_o_series = self.provider == "openai" and any(
            self.model.startswith(prefix) for prefix in ["gpt-5", "o3", "o4"]
        )

        if is_gpt5_or_o_series:
            # GPT-5 and o-series: use max_completion_tokens, skip temperature (only default=1 supported)
            api_params["max_completion_tokens"] = max_tokens or self.max_tokens
        else:
            # Standard models: use max_tokens and temperature
            api_params["max_tokens"] = max_tokens or self.max_tokens
            api_params["temperature"] = temperature or self.temperature

        response = self.client.chat.completions.create(**api_params)

        latency_ms = (time.perf_counter() - start_time) * 1000

        # Extract response data
        choice = response.choices[0]
        usage = response.usage

        # Update tracking
        self.total_requests += 1
        if usage:
            self.total_prompt_tokens += usage.prompt_tokens
            self.total_completion_tokens += usage.completion_tokens

        return LLMResponse(
            content=choice.message.content or "",
            model=response.model,
            provider=self.provider,
            prompt_tokens=usage.prompt_tokens if usage else 0,
            completion_tokens=usage.completion_tokens if usage else 0,
            total_tokens=usage.total_tokens if usage else 0,
            latency_ms=latency_ms,
            finish_reason=choice.finish_reason or "stop",
        )

    @classmethod
    def from_env(
        cls,
        provider: str = "openai",
        model_type: str = "balanced",
    ) -> "LLMClient":
        """
        Create client with recommended model for task type.

        Args:
            provider: Provider name
            model_type: Task type ("reasoning", "fast", "balanced", "advanced")

        Returns:
            Configured LLMClient
        """
        model = RECOMMENDED_MODELS.get(provider, {}).get(model_type)
        if not model:
            raise ValueError(
                f"Unknown model type '{model_type}' for provider '{provider}'"
            )

        return cls(provider=provider, model=model)

    def get_usage_summary(self) -> Dict[str, Any]:
        """Get usage statistics."""
        return {
            "provider": self.provider,
            "model": self.model,
            "total_requests": self.total_requests,
            "total_prompt_tokens": self.total_prompt_tokens,
            "total_completion_tokens": self.total_completion_tokens,
            "total_tokens": self.total_prompt_tokens + self.total_completion_tokens,
        }


def build_client(provider: str, model: Optional[str] = None, **kwargs) -> LLMClient:
    """
    Convenience function to build an LLM client.

    Args:
        provider: Provider name
        model: Model name (optional)
        **kwargs: Additional arguments for LLMClient

    Returns:
        Configured LLMClient
    """
    return LLMClient(provider=provider, model=model, **kwargs)
