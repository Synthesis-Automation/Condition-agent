"""
Interactive CLI for ChemTools LangGraph Agent.

This provides a command-line interface for interacting with the ChemTools
ReAct agent, allowing users to ask chemistry questions and get AI-powered
answers backed by deterministic ChemTools functions.

Usage:
    python -m lang_chain.chemtools_cli
    
    # Or if in lang_chain directory:
    python chemtools_cli.py

Features:
    - Conversational interface with history
    - Chemistry-specific commands
    - Tool execution visibility
    - Example queries for quick start

Commands:
    quit, exit, q  - Exit the CLI
    clear, cls     - Clear conversation history
    help           - Show help and examples
    tools          - List available tools
    builder        - Launch rule-builder wizard (create or edit rule DBs)
    verbose on/off - Toggle verbose mode
"""

import sys
import os
import re
from pathlib import Path
from typing import List, Optional
import threading
import time

# Add parent directory to path for imports
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from langchain_core.messages import BaseMessage
from chem_assistant.chemtools_agent import ChemToolsAgent
from chem_assistant.chemtools_wrapper import (
    get_tool_descriptions,
    clear_recommendation_cache,
    recommendation_cache_stats,
)
from chem_assistant.constraint_parser import (
    ConstraintSpec,
    build_constraint_spec,
    format_constraints_for_prompt,
    merge_specs,
)
from chem_assistant.rule_builder_session import RuleBuilderSession
from chemtools.rule import RuleBuilder

# ============================================================================
# Colors for Terminal Output
# ============================================================================

class Colors:
    """ANSI color codes for terminal output."""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ============================================================================
# Progress Spinner
# ============================================================================

class Spinner:
    """Animated spinner to show agent is working."""
    
    def __init__(self, message="Processing"):
        """Initialize spinner.
        
        Args:
            message: Message to display while spinning
        """
        self.message = message
        self.spinning = False
        self.thread = None
        self.spinner_chars = ['|', '/', '-', '\\']
        self.current_idx = 0
    
    def _spin(self):
        """Internal spinning loop."""
        while self.spinning:
            char = self.spinner_chars[self.current_idx]
            print(f'\r{Colors.OKCYAN}{char} {self.message}...{Colors.ENDC}', end='', flush=True)
            self.current_idx = (self.current_idx + 1) % len(self.spinner_chars)
            time.sleep(0.1)
    
    def start(self):
        """Start the spinner."""
        self.spinning = True
        self.thread = threading.Thread(target=self._spin, daemon=True)
        self.thread.start()
    
    def stop(self):
        """Stop the spinner and clear the line."""
        self.spinning = False
        if self.thread:
            self.thread.join()
        # Clear the spinner line
        print('\r' + ' ' * (len(self.message) + 20), end='\r', flush=True)


# ============================================================================
# Input Prompt Indicator
# ============================================================================

class InputPrompt:
    """Visual indicator for user input with pulsing effect."""
    
    def __init__(self):
        """Initialize input prompt indicator."""
        self.pulsing = False
        self.thread = None
        self.prompt_chars = ['⚛', '⚛️', '⚛ ', '⚛  ']
        self.current_idx = 0
    
    def _pulse(self):
        """Internal pulsing loop."""
        while self.pulsing:
            char = self.prompt_chars[self.current_idx]
            print(f'\r{Colors.OKGREEN}{Colors.BOLD}You {char}▸{Colors.ENDC} ', end='', flush=True)
            self.current_idx = (self.current_idx + 1) % len(self.prompt_chars)
            time.sleep(0.3)
    
    def start(self):
        """Start the pulsing prompt."""
        self.pulsing = True
        self.thread = threading.Thread(target=self._pulse, daemon=True)
        self.thread.start()
    
    def stop(self):
        """Stop pulsing and show final prompt."""
        self.pulsing = False
        if self.thread:
            self.thread.join()
        # Show final static prompt
        print(f'\r{Colors.OKGREEN}{Colors.BOLD}You ▸{Colors.ENDC} ', end='', flush=True)


# ============================================================================
# Example Queries
# ============================================================================

EXAMPLE_QUERIES = [
    "Normalize this SMILES: c1ccccc1",
    "What functional groups are in CCO?",
    "Classify this reactant: Brc1ccccc1",
    "What reaction family is: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "Recommend conditions for: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "Recommend Pd-free conditions for: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
    "Search for precedents of Suzuki coupling",
    "What is the CAS number for Cs2CO3?",
    "Look up information about XPhos ligand",
]


# ============================================================================
# CLI Class
# ============================================================================

class ChemToolsCLI:
    """Interactive command-line interface for ChemTools agent."""
    
    def __init__(self, verbose: bool = False):
        """
        Initialize CLI.
        
        Args:
            verbose: Show detailed tool execution
        """
        self.verbose = verbose
        self.history: List[BaseMessage] = []
        self.agent = None
        self.constraint_spec = ConstraintSpec()
        
    def initialize_agent(self):
        """Initialize the agent (lazy loading)."""
        if self.agent is None:
            print(f"{Colors.OKCYAN}Initializing ChemTools Agent...{Colors.ENDC}")
            try:
                self.agent = ChemToolsAgent(
                    verbose=self.verbose,
                    session_constraints=self.constraint_spec
                )
                print(f"{Colors.OKGREEN}Agent ready!{Colors.ENDC}\n")
            except Exception as e:
                print(f"{Colors.FAIL}Failed to initialize agent: {e}{Colors.ENDC}")
                print(f"\n{Colors.WARNING}Please ensure your API keys are configured:{Colors.ENDC}")
                print("  - Set OPENAI_API_KEY or ALIYUN_API_KEY")
                print("  - Optionally set LLM_PROVIDER (openai/aliyun) and LLM_MODEL")
                sys.exit(1)
    
    def print_header(self):
        """Print CLI header."""
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.HEADER}{Colors.BOLD}ChemTools AI Agent - Interactive Chemistry Assistant{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print()
        print("Ask chemistry questions and get AI-powered answers backed by ChemTools.")
        print()
        print(f"{Colors.OKBLUE}Commands:{Colors.ENDC}")
        print("  help           - Show example queries")
        print("  tools          - List available tools")
        print("  clear/cls      - Clear conversation history")
        print("  verbose on/off - Toggle detailed output")
        print("  constraints    - Manage catalyst/solvent constraints")
        print("  cache          - Show or clear recommendation cache")
        print("  quit/exit/q    - Exit")
        print()
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
    
    def print_help(self):
        """Print help and examples."""
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Example Queries:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
        
        for i, example in enumerate(EXAMPLE_QUERIES, 1):
            print(f"{Colors.OKCYAN}{i}.{Colors.ENDC} {example}")
        
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.WARNING}Tip:{Colors.ENDC} The agent uses multiple tools automatically.")
        print("You can ask follow-up questions - conversation history is maintained.")
        print(f"\n{Colors.OKBLUE}Constraint commands:{Colors.ENDC}")
        print("  constraints show          - Display active constraints")
        print("  constraints set Pd-free   - Replace with new preferences")
          print("  constraints allow Cu Ni   - Override allowed metals")
          print("  constraints cross on      - Enable cross-family search")
          print(f"\n{Colors.OKBLUE}Cache commands:{Colors.ENDC}")
          print("  cache show                - Display cached recommendation entries")
          print("  cache clear               - Flush cached entries")
          print(f"\n{Colors.OKBLUE}Rule builder:{Colors.ENDC}")
          print("  builder                   - Launch guided rule DB wizard")
          print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
    
    def print_tools(self):
        """Print available tools."""
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Available ChemTools:{Colors.ENDC}")
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
    
    def get_user_input_with_indicator(self):
        """Get user input with animated prompt indicator.
        
        Returns:
            str: User input string
        """
        prompt = InputPrompt()
        prompt.start()
        
        # Wait a moment for visual effect, then stop before input
        time.sleep(0.5)
        prompt.stop()
        
        # Now get input with cyan color for user text
        sys.stdout.write(f"{Colors.OKCYAN}")
        sys.stdout.flush()
        user_input = input()
        sys.stdout.write(f"{Colors.ENDC}")
        sys.stdout.flush()
        
        return user_input.strip()
    
    def print_constraints(self):
        """Display the currently active constraint summary."""
        summary = self.constraint_spec.formatted_summary()
        prompt_hint = format_constraints_for_prompt(self.constraint_spec) or ""
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Active Constraints:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{summary if summary != 'none' else 'None'}")
        if self.constraint_spec.constraint_rules:
            rules = ", ".join(sorted(self.constraint_spec.constraint_rules.keys()))
            print(f"Rules enabled: {rules}")
        if prompt_hint:
            print(prompt_hint)
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")

    def print_cache_stats(self) -> None:
        """Display recommendation cache statistics."""
        stats = recommendation_cache_stats()
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Recommendation Cache:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"Entries: {stats['entries']}")
        if stats["entries"]:
            for idx, item in enumerate(stats["items"], 1):
                print(f"\n[{idx}] Reaction key: {item['reaction']}")
                print(f"     k={item['k']} max_variants={item['max_variants']}")
                print(f"     search_all_families={item['search_all_families']}")
                if item["allow_metals"]:
                    print(f"     allow_metals={', '.join(item['allow_metals'])}")
                if item["exclude_metals"]:
                    print(f"     exclude_metals={', '.join(item['exclude_metals'])}")
                if item["prefer_metals"]:
                    print(f"     prefer_metals={', '.join(item['prefer_metals'])}")
                if item["constraint_rules"]:
                    print(f"     constraint_rules={dict(item['constraint_rules'])}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
    
    def sync_agent_constraints(self) -> None:
        """Push the current constraint spec into the agent instance."""
        if self.agent:
            self.agent.session_constraints = self.constraint_spec
    
    @staticmethod
    def _tokenize_metals(payload: str) -> List[str]:
        return [
            token for token in re.split(r"[\\s,]+", payload.strip())
            if token
        ]
    
    def handle_constraints_command(self, raw_input: str) -> None:
        """Parse and execute constraint-prefixed commands."""
        parts = raw_input.strip().split(maxsplit=2)
        if len(parts) == 1:
            self.print_constraints()
            return
        
        action = parts[1].lower()
        payload = parts[2] if len(parts) > 2 else ""
        
        truthy = {"on", "true", "yes", "1", "enable", "enabled", "all"}
        falsy = {"off", "false", "no", "0", "disable", "disabled", "auto", "detect"}
        
        if action in {"show", "status"}:
            self.print_constraints()
            return
        
        if action in {"clear", "reset"}:
            self.constraint_spec = ConstraintSpec()
            self.sync_agent_constraints()
            print(f"{Colors.OKGREEN}Constraints cleared.{Colors.ENDC}\n")
            return
        
        if action == "set":
            spec = build_constraint_spec(text=payload)
            self.constraint_spec = spec
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action == "add":
            spec = build_constraint_spec(text=payload)
            self.constraint_spec = merge_specs(self.constraint_spec, spec)
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action in {"allow", "include"}:
            metals = self._tokenize_metals(payload)
            if not metals:
                print(f"{Colors.WARNING}Provide one or more metal symbols to allow.{Colors.ENDC}\n")
                return
            spec = build_constraint_spec(allow_metals=metals)
            self.constraint_spec.allow_metals = spec.allow_metals
            self.constraint_spec.exclude_metals -= self.constraint_spec.allow_metals
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action in {"exclude", "avoid"}:
            metals = self._tokenize_metals(payload)
            if not metals:
                print(f"{Colors.WARNING}Provide one or more metal symbols to exclude.{Colors.ENDC}\n")
                return
            spec = build_constraint_spec(exclude_metals=metals)
            self.constraint_spec.exclude_metals = spec.exclude_metals
            self.constraint_spec.allow_metals -= self.constraint_spec.exclude_metals
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action == "prefer":
            metals = self._tokenize_metals(payload)
            if not metals:
                print(f"{Colors.WARNING}Provide one or more metal symbols to prefer.{Colors.ENDC}\n")
                return
            spec = build_constraint_spec(prefer_metals=metals)
            self.constraint_spec.prefer_metals = spec.prefer_metals
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action in {"cross", "crossfamily", "families"}:
            state = payload.strip().lower()
            if not state:
                self.constraint_spec.search_all_families = not self.constraint_spec.search_all_families
            elif state in truthy:
                self.constraint_spec.search_all_families = True
            elif state in falsy:
                self.constraint_spec.search_all_families = False
            else:
                print(f"{Colors.WARNING}Specify 'on' or 'off' to control cross-family search.{Colors.ENDC}\n")
                return
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        if action in {"rule", "rules"}:
            if not payload:
                if self.constraint_spec.constraint_rules:
                    rules = ", ".join(f"{k}=on" for k in sorted(self.constraint_spec.constraint_rules.keys()))
                    print(f"{Colors.OKBLUE}Rules active: {rules}{Colors.ENDC}\n")
                else:
                    print(f"{Colors.OKBLUE}No optional rules active.{Colors.ENDC}\n")
                return
            parts_rule = payload.split(maxsplit=1)
            alias = parts_rule[0].lower()
            state = parts_rule[1].lower() if len(parts_rule) > 1 else ""
            mapping = {
                "no_chlorinated": "no_chlorinated",
                "chlorinated": "no_chlorinated",
                "aqueous": "aqueous_only",
                "aqueous_only": "aqueous_only",
                "no_hmpa": "no_HMPA",
                "hmpa": "no_HMPA",
            }
            if alias not in mapping:
                print(f"{Colors.WARNING}Unknown rule '{alias}'. Supported: no_chlorinated, aqueous, no_hmpa.{Colors.ENDC}\n")
                return
            rule_key = mapping[alias]
            if not state or state in truthy:
                self.constraint_spec.constraint_rules[rule_key] = True
            elif state in falsy:
                self.constraint_spec.constraint_rules.pop(rule_key, None)
            else:
                print(f"{Colors.WARNING}Specify 'on' or 'off' for rule state.{Colors.ENDC}\n")
                return
            self.sync_agent_constraints()
            self.print_constraints()
            return
        
        print(
            f"{Colors.WARNING}Unrecognized constraints command. "
            "Usage examples: 'constraints show', 'constraints set Pd-free', "
            "'constraints allow Pd', 'constraints cross on'.{Colors.ENDC}\n"
        )

    def launch_rule_builder(self) -> None:
        """Launch the rule-builder wizard from the CLI."""
        print(f"\n{Colors.HEADER}{Colors.BOLD}Rule Builder Wizard{Colors.ENDC}")
        print(
            "Guided workflow for rule databases.\n"
            "Hint: baseline JSON lives in data/rule_db_v2/.\n"
        )
        choice = input("Load existing file (L) or create new (N)? [L]: ").strip().lower()
        target_path: Optional[Path]

        try:
            if choice in {"", "l", "load"}:
                load_path = input("Path to existing rule DB JSON: ").strip()
                if not load_path:
                    print(f"{Colors.WARNING}No path provided. Aborting builder.{Colors.ENDC}\n")
                    return
                source_path = Path(load_path).expanduser().resolve()
                if not source_path.exists():
                    print(f"{Colors.FAIL}File not found: {source_path}{Colors.ENDC}\n")
                    return
                builder = RuleBuilder.from_file(source_path)
                new_path = input(f"Save changes to path [{source_path}]: ").strip()
                target_path = Path(new_path).expanduser() if new_path else source_path
                print(f"{Colors.OKCYAN}Loaded {source_path}{Colors.ENDC}\n")
            else:
                family = input("Reaction family (e.g., Suzuki_Miyaura): ").strip()
                if not family:
                    print(f"{Colors.WARNING}Family is required to start a new rule DB.{Colors.ENDC}\n")
                    return
                builder = RuleBuilder.new(family)
                suggested = Path("data") / "rule_db_v2" / f"{family}_draft.json"
                out_path = input(f"Output path [{suggested}]: ").strip()
                target_path = Path(out_path).expanduser() if out_path else suggested

            session = RuleBuilderSession(builder)
            session.run_wizard()

            print(f"{Colors.OKCYAN}Running validation...{Colors.ENDC}")
            issues = builder.validate(strict=False)
            if issues:
                print(f"{Colors.WARNING}Validation findings:{Colors.ENDC}")
                for issue in issues:
                    mark = "!" if issue.severity == "error" else "-"
                    print(f"  {mark} [{issue.severity}] {issue.field}: {issue.message}")
                if any(i.severity == "error" for i in issues):
                    print(f"{Colors.FAIL}Errors detected. Resolve them before saving.{Colors.ENDC}\n")
                    return
            else:
                print(f"{Colors.OKGREEN}No validation issues detected.{Colors.ENDC}")

            diff_text = builder.diff()
            if diff_text.strip():
                print(f"\n{Colors.OKCYAN}Change preview:{Colors.ENDC}\n{diff_text}\n")
            else:
                print(f"{Colors.WARNING}No changes detected relative to source.{Colors.ENDC}\n")

            confirm = input(f"Save to {target_path}? (y/n) [y]: ").strip().lower()
            if confirm in {"", "y", "yes"}:
                builder.save(target_path)
                print(f"{Colors.OKGREEN}Saved rule database to {target_path}{Colors.ENDC}\n")
            else:
                print(f"{Colors.WARNING}Save skipped. No files were written.{Colors.ENDC}\n")
        except KeyboardInterrupt:
            print(f"\n{Colors.WARNING}Rule builder cancelled by user.{Colors.ENDC}\n")
        except Exception as exc:
            print(f"\n{Colors.FAIL}Rule builder error: {exc}{Colors.ENDC}\n")
    
    def handle_command(self, user_input: str) -> bool:
        """
        Handle special commands.
        
        Args:
            user_input: User input string
        
        Returns:
            True if command was handled, False otherwise
        """
        cmd = user_input.lower().strip()

        if cmd.startswith('constraint'):
            self.handle_constraints_command(user_input)
            return False

        if cmd.startswith('cache'):
            self.handle_cache_command(user_input)
            return False

        if cmd.startswith('builder'):
            self.launch_rule_builder()
            return False
        
        # Exit commands
        if cmd in ['quit', 'exit', 'q']:
            print(f"\n{Colors.OKGREEN}Goodbye!{Colors.ENDC}\n")
            return True
        
        # Help
        if cmd == 'help':
            self.print_help()
            return False
        
        # Tools
        if cmd == 'tools':
            self.print_tools()
            return False
        
        # Clear history
        if cmd in ['clear', 'cls']:
            self.history = []
            print(f"{Colors.OKGREEN}Conversation history cleared.{Colors.ENDC}\n")
            return False
        
        # Verbose toggle
        if cmd.startswith('verbose'):
            parts = cmd.split()
            if len(parts) > 1:
                if parts[1] == 'on':
                    self.verbose = True
                    if self.agent:
                        self.agent.verbose = True
                    print(f"{Colors.OKGREEN}Verbose mode enabled.{Colors.ENDC}\n")
                elif parts[1] == 'off':
                    self.verbose = False
                    if self.agent:
                        self.agent.verbose = False
                    print(f"{Colors.OKGREEN}Verbose mode disabled.{Colors.ENDC}\n")
            else:
                status = "ON" if self.verbose else "OFF"
                print(f"{Colors.OKCYAN}Verbose mode is {status}{Colors.ENDC}\n")
            return False
        
        return False  # Not a command

    def handle_cache_command(self, raw_input: str) -> None:
        """Parse cache-related commands (show/clear)."""
        parts = raw_input.strip().split()
        if len(parts) == 1 or parts[1].lower() in {"show", "status"}:
            self.print_cache_stats()
            return
        if parts[1].lower() == "clear":
            clear_recommendation_cache()
            print(f"{Colors.OKGREEN}Recommendation cache cleared.{Colors.ENDC}\n")
            return
        print(
            f"{Colors.WARNING}Unrecognized cache command. "
            "Usage: 'cache', 'cache show', or 'cache clear'.{Colors.ENDC}\n"
        )
    
    def run(self):
        """Run the interactive CLI."""
        self.print_header()
        self.initialize_agent()
        
        while True:
            try:
                # Get user input with animated indicator
                user_input = self.get_user_input_with_indicator()
                
                # Empty input
                if not user_input:
                    continue
                # Handle commands
                should_exit = self.handle_command(user_input)
                if should_exit:
                    break
                
                # Skip if it was a command
                if user_input.lower() in ['help', 'tools', 'clear', 'cls'] or \
                   user_input.lower().startswith('verbose'):
                    continue
                
                # Process query with agent
                spinner = Spinner("Agent is thinking")
                spinner.start()
                
                try:
                    response, self.history = self.agent.chat(
                        user_input,
                        self.history,
                        recursion_limit=15,
                        constraints=self.constraint_spec
                    )
                finally:
                    spinner.stop()
                
                print(f"{Colors.OKBLUE}{Colors.BOLD}Agent:{Colors.ENDC} {response}")
                print()
                
            except KeyboardInterrupt:
                print(f"\n\n{Colors.WARNING}Interrupted. Type 'quit' to exit.{Colors.ENDC}\n")
                continue
            except EOFError:
                print(f"\n{Colors.OKGREEN}Goodbye!{Colors.ENDC}\n")
                break
            except Exception as e:
                print(f"\n{Colors.FAIL}Error: {e}{Colors.ENDC}\n")
                if self.verbose:
                    import traceback
                    traceback.print_exc()


# ============================================================================
# Main Entry Point
# ============================================================================

def main():
    """Main entry point for CLI."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Interactive CLI for ChemTools AI Agent",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python chemtools_cli.py
  python chemtools_cli.py --verbose
  python -m lang_chain.chemtools_cli

Environment Variables:
  LLM_PROVIDER    - LLM provider: "openai" or "aliyun" (default: openai)
  LLM_MODEL       - Model name (default: gpt-4o)
  OPENAI_API_KEY  - OpenAI API key (required if provider=openai)
  ALIYUN_API_KEY  - Aliyun API key (required if provider=aliyun)
        """
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output (show tool execution details)'
    )
    
    args = parser.parse_args()
    
    # Create and run CLI
    cli = ChemToolsCLI(verbose=args.verbose)
    cli.run()


if __name__ == "__main__":
    main()
