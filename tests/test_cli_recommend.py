"""Test script for CLI with automated input."""

import sys
from pathlib import Path
from io import StringIO

# Add parent to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Mock input for testing
test_inputs = [
    "Clc1c(C)cccc1C.COc1ccc(B(O)O)cc1>>Cc1cccc(C)c1-c1ccc(OC)cc1",  # Reaction SMILES
    "",  # No requirements (press Enter)
    "yes",  # Confirm submission
]

input_iter = iter(test_inputs)

def mock_input(prompt=""):
    print(prompt, end="")
    value = next(input_iter)
    print(value)  # Echo the input
    return value

# Replace built-in input with mock
import builtins
builtins.input = mock_input

# Now run the CLI
from app.cli_recommend import main

if __name__ == "__main__":
    try:
        main()
    except StopIteration:
        print("\n✅ Test completed (ran out of mock inputs)")
    except SystemExit as e:
        print(f"\n✅ Test completed with exit code: {e.code}")
