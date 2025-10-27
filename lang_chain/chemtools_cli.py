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
    verbose on/off - Toggle verbose mode
"""

import sys
import os
from pathlib import Path
from typing import List
import threading
import time

# Add parent directory to path for imports
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from langchain_core.messages import BaseMessage
from lang_chain.chemtools_agent import ChemToolsAgent
from lang_chain.chemtools_wrapper import get_tool_descriptions


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
        self.spinner_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']
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
# Example Queries
# ============================================================================

EXAMPLE_QUERIES = [
    "Normalize this SMILES: c1ccccc1",
    "What functional groups are in CCO?",
    "Classify this reactant: Brc1ccccc1",
    "What reaction family is: Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1",
    "Recommend conditions for: Brc1ccccc1.Nc1ccccc1>>c1ccccc1Nc1ccccc1",
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
        
    def initialize_agent(self):
        """Initialize the agent (lazy loading)."""
        if self.agent is None:
            print(f"{Colors.OKCYAN}⚙️  Initializing ChemTools Agent...{Colors.ENDC}")
            try:
                self.agent = ChemToolsAgent(verbose=self.verbose)
                print(f"{Colors.OKGREEN}✓ Agent ready!{Colors.ENDC}\n")
            except Exception as e:
                print(f"{Colors.FAIL}✗ Failed to initialize agent: {e}{Colors.ENDC}")
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
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
    
    def print_tools(self):
        """Print available tools."""
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{Colors.BOLD}Available ChemTools:{Colors.ENDC}")
        print(f"{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
        
        tools = get_tool_descriptions()
        for i, tool in enumerate(tools, 1):
            print(f"{Colors.OKCYAN}{i}. {tool['name']}{Colors.ENDC}")
            desc = tool['description'].split('\n')[0]  # First line only
            print(f"   {desc[:100]}...")
        
        print(f"\n{Colors.BOLD}{'=' * 70}{Colors.ENDC}\n")
    
    def handle_command(self, user_input: str) -> bool:
        """
        Handle special commands.
        
        Args:
            user_input: User input string
        
        Returns:
            True if command was handled, False otherwise
        """
        cmd = user_input.lower().strip()
        
        # Exit commands
        if cmd in ['quit', 'exit', 'q']:
            print(f"\n{Colors.OKGREEN}Goodbye! 👋{Colors.ENDC}\n")
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
            print(f"{Colors.OKGREEN}✓ Conversation history cleared.{Colors.ENDC}\n")
            return False
        
        # Verbose toggle
        if cmd.startswith('verbose'):
            parts = cmd.split()
            if len(parts) > 1:
                if parts[1] == 'on':
                    self.verbose = True
                    if self.agent:
                        self.agent.verbose = True
                    print(f"{Colors.OKGREEN}✓ Verbose mode enabled.{Colors.ENDC}\n")
                elif parts[1] == 'off':
                    self.verbose = False
                    if self.agent:
                        self.agent.verbose = False
                    print(f"{Colors.OKGREEN}✓ Verbose mode disabled.{Colors.ENDC}\n")
            else:
                status = "ON" if self.verbose else "OFF"
                print(f"{Colors.OKCYAN}Verbose mode is {status}{Colors.ENDC}\n")
            return False
        
        return False  # Not a command
    
    def run(self):
        """Run the interactive CLI."""
        self.print_header()
        self.initialize_agent()
        
        while True:
            try:
                # Get user input
                user_input = input(f"{Colors.BOLD}You:{Colors.ENDC} ").strip()
                
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
                        recursion_limit=15
                    )
                finally:
                    spinner.stop()
                
                print(f"{Colors.OKBLUE}{Colors.BOLD}Agent:{Colors.ENDC} {response}")
                print()
                
            except KeyboardInterrupt:
                print(f"\n\n{Colors.WARNING}Interrupted. Type 'quit' to exit.{Colors.ENDC}\n")
                continue
            except EOFError:
                print(f"\n{Colors.OKGREEN}Goodbye! 👋{Colors.ENDC}\n")
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
