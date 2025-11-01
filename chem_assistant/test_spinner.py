"""
Simple test of the CLI with spinner.

Run this to see the spinner in action without needing an API key.
"""

import time

class Colors:
    OKCYAN = '\033[96m'
    OKBLUE = '\033[94m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'

class Spinner:
    """Animated spinner to show agent is working."""
    
    def __init__(self, message="Processing"):
        self.message = message
        self.spinning = False
        import threading
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
        import threading
        self.thread = threading.Thread(target=self._spin, daemon=True)
        self.thread.start()
    
    def stop(self):
        """Stop the spinner and clear the line."""
        self.spinning = False
        if self.thread:
            self.thread.join()
        # Clear the spinner line
        print('\r' + ' ' * (len(self.message) + 20), end='\r', flush=True)


if __name__ == "__main__":
    print("Testing spinner...")
    
    # Test 1: Short task
    spinner = Spinner("Loading tools")
    spinner.start()
    time.sleep(2)
    spinner.stop()
    print(f"{Colors.OKBLUE}✓ Tools loaded!{Colors.ENDC}\n")
    
    # Test 2: Longer task
    spinner = Spinner("Agent is thinking")
    spinner.start()
    time.sleep(3)
    spinner.stop()
    print(f"{Colors.OKBLUE}{Colors.BOLD}Agent:{Colors.ENDC} Here's my response!")
    
    print("\n✅ Spinner test complete!")
