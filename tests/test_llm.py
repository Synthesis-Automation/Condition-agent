"""Quick utilities to sanity check LLM connectivity.

This script intentionally keeps the API surface small so it can double
as a troubleshooting aid when API calls fail from the main application.
Set the appropriate environment variables before running:

Aliyun (DeepSeek-compatible):
    setx ALIYUN_API_KEY "sk-..."
    # optional override
    setx ALIYUN_BASE_URL "https://dashscope.aliyuncs.com/compatible-mode/v1"

OpenAI:
    setx OPENAI_API_KEY "sk-..."
    # optional override
    setx OPENAI_BASE_URL "https://api.openai.com/v1"

Usage examples (PowerShell):
    python test_llm.py aliyun --model deepseek-r1-distill-qwen-7b "9.9和9.11谁大"
    python test_llm.py openai --model gpt-4o-mini "Suggest a solvent"
    python test_llm.py  # choose provider and model interactively
"""

import argparse
import os
import sys

from openai import OpenAI


DEFAULT_ALIYUN_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"


def build_client(provider: str) -> OpenAI:
    """Return a configured OpenAI-compatible client for the requested provider."""

    if provider == "aliyun":
        api_key = os.getenv("ALIYUN_API_KEY")
        base_url = os.getenv("ALIYUN_BASE_URL", DEFAULT_ALIYUN_BASE_URL)
    elif provider == "openai":
        api_key = os.getenv("OPENAI_API_KEY")
        base_url = os.getenv("OPENAI_BASE_URL", DEFAULT_OPENAI_BASE_URL)
    else:
        raise ValueError(f"Unsupported provider '{provider}'")

    if not api_key:
        raise RuntimeError(
            f"Missing API key. Set {'ALIYUN_API_KEY' if provider == 'aliyun' else 'OPENAI_API_KEY'}."
        )

    return OpenAI(api_key=api_key, base_url=base_url)


def run_chat_completion(client: OpenAI, model: str, prompt: str) -> str:
    """Send a single-turn chat completion request and return the text reply."""

    response = client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
    )
    return response.choices[0].message.content


PROVIDERS = ("aliyun", "openai")
AVAILABLE_MODELS = {
    "aliyun": [
        "deepseek-r1-distill-qwen-7b",
        "deepseek-v3.2",
        "deepseek-v3.1",
        "deepseek-r1",
        "deepseek-r1-0528",
        "deepseek-v3",
        "deepseek-r1-distill-qwen-1.5b",
        "deepseek-r1-distill-qwen-14b",
        "deepseek-r1-distill-qwen-32b",
        "deepseek-r1-distill-llama-8b",
        "deepseek-r1-distill-llama-70b",
    ],
    "openai": [
        # GPT-5 series (general/agentic)
        "gpt-5",
        "gpt-5-pro",
        "gpt-5-mini",
        "gpt-5-nano",
        "gpt-5-codex",
        # o-series (reasoning models)
        "o3",
        "o3-pro",
        "o3-mini",
        "o4-mini",
        "o3-deep-research",
        "o4-mini-deep-research",
        # GPT-4.x series (versatile, widely supported)
        "gpt-4o",
        "gpt-4o-mini",
        "gpt-4.1-mini",
        "gpt-4.1-nano"
    ],
}


def prompt_for_provider() -> str:
    """Interactively ask the user to select a provider when none was supplied."""

    print("Select a provider:")
    for idx, option in enumerate(PROVIDERS, start=1):
        print(f"  {idx}. {option}")

    while True:
        choice = input("Enter provider number or name: ").strip().lower()
        if not choice:
            continue

        if choice in PROVIDERS:
            return choice

        if choice.isdigit():
            index = int(choice) - 1
            if 0 <= index < len(PROVIDERS):
                return PROVIDERS[index]

        print("Invalid selection. Please try again.")


def prompt_for_user_message() -> str:
    """Ask the user for a prompt when none was supplied via the CLI."""

    while True:
        text = input("Enter prompt: ").strip()
        if text:
            return text
        print("Prompt cannot be empty.")


def prompt_for_model(provider: str) -> str:
    """Help the user pick a model for the chosen provider."""

    models = AVAILABLE_MODELS.get(provider, [])
    if not models:
        return input("Enter model name: ").strip()

    print("Select a model:")
    for idx, name in enumerate(models, start=1):
        default_suffix = " (default)" if idx == 1 else ""
        print(f"  {idx}. {name}{default_suffix}")
    print("  0. Enter custom model name")

    while True:
        choice = input("Enter model number or name [default 1]: ").strip()
        if not choice:
            return models[0]

        if choice == "0":
            custom = input("Enter custom model name: ").strip()
            if custom:
                return custom
            print("Model name cannot be empty.")
            continue

        if choice.isdigit():
            index = int(choice) - 1
            if 0 <= index < len(models):
                return models[index]

        if choice in models:
            return choice

        print("Invalid selection. Please try again.")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "provider",
        nargs="?",
        choices=PROVIDERS,
        help="Which OpenAI-compatible provider to query.",
    )
    parser.add_argument(
        "prompt",
        nargs="?",
        help="User message to send to the model.",
    )
    parser.add_argument(
        "--model",
        help="Model name to invoke. Skips the interactive model picker when provided.",
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    provider = args.provider or prompt_for_provider()
    prompt = args.prompt or prompt_for_user_message()
    model = args.model or prompt_for_model(provider)

    if not model:
        print("Model name cannot be empty.", file=sys.stderr)
        return 1

    print(f"Requesting {provider} model '{model}'...")

    try:
        client = build_client(provider)
    except Exception as exc:  # narrow failure path to ease CLI debugging
        print(f"Failed to create client: {exc}", file=sys.stderr)
        return 1

    try:
        answer = run_chat_completion(client, model, prompt)
    except Exception as exc:  # catch API errors and bubble them up clearly
        print(f"Model request failed: {exc}", file=sys.stderr)
        return 2

    print("Model response:\n")
    print(answer.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))


