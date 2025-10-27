from typing import List
import json
import random
import string
import os
from datetime import datetime, timedelta
from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, AIMessage, BaseMessage
from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent

from dotenv import load_dotenv

load_dotenv()


# -------- API Configuration --------
DEFAULT_ALIYUN_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
DEFAULT_OPENAI_BASE_URL = "https://api.openai.com/v1"


def get_llm_client():
    """Get LLM client using configured provider and model."""
    # Determine provider from environment
    provider = os.getenv("LLM_PROVIDER", "openai")
    model = os.getenv("LLM_MODEL", "gpt-4o-mini")
    
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

    return ChatOpenAI(
        model=model,
        api_key=api_key,
        base_url=base_url,
        temperature=0
    )


# -------- Tools --------
@tool
def write_json(filepath: str, data: dict) -> str:
    """Write a Python dictionary as JSON to a file with pretty formatting."""
    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        return f"Successfully wrote JSON data to '{filepath}' ({len(json.dumps(data))} characters)."
    except Exception as e:
        return f"Error writing JSON: {str(e)}"


@tool
def read_json(filepath: str) -> str:
    """Read and return the contents of a JSON file."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            data = json.load(f)
        return json.dumps(data, indent=2)
    except FileNotFoundError:
        return f"Error: File '{filepath}' not found."
    except json.JSONDecodeError as e:
        return f"Error: Invalid JSON in file - {str(e)}"
    except Exception as e:
        return f"Error reading JSON: {str(e)}"


@tool
def generate_sample_users(
        first_names: List[str],
        last_names: List[str],
        domains: List[str],
        min_age: int,
        max_age: int
) -> dict:
    """
    Generate sample user data. Count is determined by the length of first_names.

    Args:
        first_names: List of first names (one per user)
        last_names: List of last names (will cycle if fewer than first_names)
        domains: List of email domains (will cycle through)
        min_age: Minimum age for users
        max_age: Maximum age for users

    Returns:
        Dictionary with 'users' array or 'error' message
    """
    # Validation
    if not first_names:
        return {"error": "first_names list cannot be empty"}
    if not last_names:
        return {"error": "last_names list cannot be empty"}
    if not domains:
        return {"error": "domains list cannot be empty"}
    if min_age > max_age:
        return {"error": f"min_age ({min_age}) cannot be greater than max_age ({max_age})"}
    if min_age < 0 or max_age < 0:
        return {"error": "ages must be non-negative"}

    users = []
    count = len(first_names)

    for i in range(count):
        first = first_names[i]
        last = last_names[i % len(last_names)]
        domain = domains[i % len(domains)]
        email = f"{first.lower()}.{last.lower()}@{domain}"

        user = {
            "id": i + 1,
            "firstName": first,
            "lastName": last,
            "email": email,
            "username": f"{first.lower()}{random.randint(100, 999)}",
            "age": random.randint(min_age, max_age),
            "registeredAt": (datetime.now() - timedelta(days=random.randint(1, 365))).isoformat()
        }
        users.append(user)

    return {"users": users, "count": len(users)}


TOOLS = [write_json, read_json, generate_sample_users]

llm = get_llm_client()

SYSTEM_MESSAGE = (
    "You are DataGen, a helpful assistant that generates sample data for applications. "
    "You have access to three tools:\n"
    "1. generate_sample_users: Takes first_names, last_names, domains, min_age, max_age and returns a dictionary with 'users' and 'count' keys\n"
    "2. write_json: Takes a filepath and a data dictionary, writes it to a JSON file\n"
    "3. read_json: Takes a filepath and returns the JSON contents\n\n"
    "To generate and save users:\n"
    "- FIRST: Call generate_sample_users with the required parameters\n"
    "- SECOND: Use the returned data in a write_json call with a filepath\n"
    "Always follow this two-step process. Fill in the user details yourself without asking for them."
)

agent = create_react_agent(llm, TOOLS, prompt=SYSTEM_MESSAGE)


def run_agent(user_input: str, history: List[BaseMessage]) -> AIMessage:
    """Single-turn agent runner with automatic tool execution via LangGraph."""
    try:
        result = agent.invoke(
            {"messages": history + [HumanMessage(content=user_input)]},
            config={"recursion_limit": 15}
        )
        # Return the last AI message
        return result["messages"][-1]
    except Exception as e:
        # Return error as an AI message so the conversation can continue
        return AIMessage(content=f"Error: {str(e)}\n\nPlease try rephrasing your request or provide more specific details.")


if __name__ == "__main__":
    print("=" * 60)
    print("DataGen Agent - Sample Data Generator")
    print("=" * 60)
    print("Generate sample user data and save to JSON files.")
    print()
    print("Examples:")
    print("  - Generate users named John, Jane, Mike and save to users.json")
    print("  - Create users with last names Smith, Jones")
    print("  - Make users aged 25-35 with company.com emails")
    print()
    print("Commands: 'quit' or 'exit' to end")
    print("=" * 60)

    history: List[BaseMessage] = []

    while True:
        user_input = input("You: ").strip()

        # Check for exit commands
        if user_input.lower() in ['quit', 'exit', 'q', ""]:
            print("Goodbye!")
            break

        print("Agent: ", end="", flush=True)
        response = run_agent(user_input, history)
        print(response.content)
        print()

        # Update conversation history
        history += [HumanMessage(content=user_input), response]