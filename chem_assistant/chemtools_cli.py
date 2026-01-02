"""
Interactive CLI for the ChemTools featurization/analysis agent.
"""

import sys
import time
from pathlib import Path
from typing import List, Optional
import threading

from langchain_core.messages import BaseMessage

# Add parent directory to path for imports
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from chem_assistant.chemtools_agent import ChemToolsAgent
from chem_assistant.chemtools_wrapper import get_tool_descriptions


class Colors:
    """ANSI color codes for terminal output."""

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"


class Spinner:
    """Animated spinner to show the agent is working."""

    def __init__(self, message: str = "Processing") -> None:
        self.message = message
        self.spinning = False
        self.thread: Optional[threading.Thread] = None
        self.spinner_chars = ["|", "/", "-", "\\"]
        self.current_idx = 0

    def _spin(self) -> None:
        while self.spinning:
            char = self.spinner_chars[self.current_idx]
            print(f"\r{Colors.OKCYAN}{char} {self.message}...{Colors.ENDC}", end="", flush=True)
            self.current_idx = (self.current_idx + 1) % len(self.spinner_chars)
            time.sleep(0.1)

    def start(self) -> None:
        self.spinning = True
        self.thread = threading.Thread(target=self._spin, daemon=True)
        self.thread.start()

    def stop(self) -> None:
        self.spinning = False
        if self.thread:
            self.thread.join()
        print("\r" + " " * (len(self.message) + 20), end="\r", flush=True)


EXAMPLE_QUERIES = [
    "Normalize this SMILES: c1ccccc1",
    "Featurize molecule: c1ccccc1O",
    "Analyze reaction: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "Detect reaction types: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "Featurize pair: electrophile=Brc1ccccc1 nucleophile=Nc1ccccc1",
]


class ChemToolsCLI:
    """Interactive command-line interface for the ChemTools agent."""

    def __init__(self, verbose: bool = False) -> None:
        self.verbose = verbose
        self.history: List[BaseMessage] = []
        self.agent: Optional[ChemToolsAgent] = None

    def initialize_agent(self) -> None:
        if self.agent is None:
            print(f"{Colors.OKCYAN}Initializing ChemTools Agent...{Colors.ENDC}")
            try:
                self.agent = ChemToolsAgent(verbose=self.verbose)
                print(f"{Colors.OKGREEN}Agent ready!{Colors.ENDC}\n")
            except Exception as exc:
                print(f"{Colors.FAIL}Failed to initialize agent: {exc}{Colors.ENDC}")
                print(f"\n{Colors.WARNING}Please ensure your API keys are configured:{Colors.ENDC}")
                print("  - Set OPENAI_API_KEY or ALIYUN_API_KEY")
                print("  - Optionally set LLM_PROVIDER (openai/aliyun) and LLM_MODEL")
                sys.exit(1)

    def print_header(self) -> None:
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.HEADER}{Colors.BOLD}ChemTools Featurization Agent{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print()
        print("Ask analysis questions and get AI-powered answers backed by ChemTools.")
        print()
        print(f"{Colors.OKBLUE}Commands:{Colors.ENDC}")
        print("  help           - Show example queries")
        print("  tools          - List available tools")
        print("  clear/cls      - Clear conversation history")
        print("  verbose on/off - Toggle detailed output")
        print("  quit/exit/q    - Exit")
        print()
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")

    def print_help(self) -> None:
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Example Queries:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
        for i, example in enumerate(EXAMPLE_QUERIES, 1):
            print(f"{Colors.OKCYAN}{i}.{Colors.ENDC} {example}")
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")

    def print_tools(self) -> None:
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Available Tools:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
        tools = get_tool_descriptions()
        for i, tool in enumerate(tools, 1):
            print(f"{Colors.OKCYAN}{i}. {tool['name']}{Colors.ENDC}")
            desc_line = (tool["description"] or "").strip()
            if desc_line:
                first_line = desc_line.split("\n")[0]
                print(f"   {first_line}")
            parameters = tool.get("parameters") or []
            if not parameters:
                print("   Parameters: none")
            else:
                print("   Parameters:")
                for param in parameters:
                    requirement = "required" if param.get("required") else "optional"
                    param_type = param.get("type") or "any"
                    print(f"     - {param['name']} ({param_type}, {requirement})")
                    if param.get("description"):
                        print(f"       {param['description']}")
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")

    def handle_command(self, user_input: str) -> bool:
        cmd = user_input.lower().strip()
        if cmd in ["quit", "exit", "q"]:
            print(f"\n{Colors.OKGREEN}Goodbye!{Colors.ENDC}\n")
            return True
        if cmd == "help":
            self.print_help()
            return False
        if cmd == "tools":
            self.print_tools()
            return False
        if cmd in ["clear", "cls"]:
            self.history = []
            print(f"{Colors.OKGREEN}Conversation history cleared.{Colors.ENDC}\n")
            return False
        if cmd.startswith("verbose"):
            parts = cmd.split()
            if len(parts) > 1:
                if parts[1] == "on":
                    self.verbose = True
                    if self.agent:
                        self.agent.verbose = True
                    print(f"{Colors.OKGREEN}Verbose mode enabled.{Colors.ENDC}\n")
                elif parts[1] == "off":
                    self.verbose = False
                    if self.agent:
                        self.agent.verbose = False
                    print(f"{Colors.OKGREEN}Verbose mode disabled.{Colors.ENDC}\n")
            else:
                status = "ON" if self.verbose else "OFF"
                print(f"{Colors.OKCYAN}Verbose mode is {status}{Colors.ENDC}\n")
            return False
        return False

    def run(self) -> None:
        self.print_header()
        self.initialize_agent()

        while True:
            try:
                user_input = input(f"{Colors.OKGREEN}You > {Colors.ENDC}").strip()
                if not user_input:
                    continue

                should_exit = self.handle_command(user_input)
                if should_exit:
                    break

                if user_input.lower() in ["help", "tools", "clear", "cls"] or user_input.lower().startswith("verbose"):
                    continue

                spinner = Spinner("Agent is thinking")
                spinner.start()
                try:
                    response, self.history = self.agent.chat(
                        user_input,
                        self.history,
                        recursion_limit=15,
                    )
                finally:
                    spinner.stop()

                print(f"{Colors.OKBLUE}{Colors.BOLD}Agent:{Colors.ENDC} {response}\n")
            except KeyboardInterrupt:
                print(f"\n\n{Colors.WARNING}Interrupted. Type 'quit' to exit.{Colors.ENDC}\n")
            except EOFError:
                print(f"\n{Colors.OKGREEN}Goodbye!{Colors.ENDC}\n")
                break
            except Exception as exc:
                print(f"\n{Colors.FAIL}Error: {exc}{Colors.ENDC}\n")
                if self.verbose:
                    import traceback

                    traceback.print_exc()


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Interactive CLI for ChemTools featurization agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose output",
    )

    args = parser.parse_args()
    cli = ChemToolsCLI(verbose=args.verbose)
    cli.run()


if __name__ == "__main__":
    main()
