# ChemAssistant GUI - Fixed Import Issue

## Problem Fixed

The GUI application (`chem_assistant/gui/app.py`) was using relative imports that didn't work when running the script directly. This has been fixed to support both:
1. Direct execution: `python chem_assistant/gui/app.py`
2. Module execution: `python -m chem_assistant.gui.app`
3. Launcher script: `python launch_gui.py`

## How to Run

### Method 1: Using the Launcher Script (Recommended)
```bash
python launch_gui.py
```

### Method 2: Direct Execution
```bash
python chem_assistant/gui/app.py
```

### Method 3: Module Execution
```bash
python -m chem_assistant.gui.app
```

### Method 4: From Python
```python
from chem_assistant.gui.app import main
main()
```

## What Was Changed

### File: `chem_assistant/gui/app.py`

**Before:**
```python
from .main_window import ChemAssistantWindow
```

**After:**
```python
import sys
from pathlib import Path

# Add project root to path for direct execution
if __name__ == "__main__":
    project_root = Path(__file__).parent.parent.parent
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

# Handle both relative and absolute imports
try:
    from .main_window import ChemAssistantWindow
except ImportError:
    from chem_assistant.gui.main_window import ChemAssistantWindow
```

## Benefits

1. ✅ Works when run directly from any location
2. ✅ Works when imported as a module
3. ✅ Works with the launcher script
4. ✅ Maintains compatibility with existing code
5. ✅ No changes needed to other files

## Requirements

Make sure you have the required dependencies installed:

```bash
pip install PyQt6 langchain-core
```

And ensure the project is in your Python path or run from the project root directory.

## Project Structure

```
Condition-agent/
├── chem_assistant/
│   └── gui/
│       ├── __init__.py
│       ├── app.py          # Main entry point (FIXED)
│       ├── main_window.py  # Main window implementation
│       └── dialogs.py      # Dialog components
├── launch_gui.py           # Convenient launcher (NEW)
└── ...
```

## Troubleshooting

### Issue: "No module named 'PyQt6'"
**Solution:**
```bash
pip install PyQt6
```

### Issue: "No module named 'chem_assistant'"
**Solution:** Make sure you're running from the project root:
```bash
cd C:\Git-softwares\Condition-agent
python launch_gui.py
```

### Issue: Import errors in other modules
**Solution:** The fix in `app.py` ensures the project root is added to `sys.path`, so all imports should work correctly.

## Testing

The GUI application should now:
- ✅ Launch without import errors
- ✅ Display the ChemAssistant main window
- ✅ Load all tools and components correctly
- ✅ Be ready for user interaction

## Related Files

- `chem_assistant/gui/app.py` - Main entry point (fixed imports)
- `chem_assistant/gui/main_window.py` - Main window UI
- `chem_assistant/gui/dialogs.py` - Dialog components
- `launch_gui.py` - Convenient launcher script
- `chem_assistant/chemtools_agent.py` - Agent implementation
- `chem_assistant/chemtools_wrapper.py` - Tool wrappers

## Next Steps

Now that the import issue is fixed, you can:
1. Run the GUI using any of the methods above
2. Test the ChemAssistant features
3. Use the constraint builder
4. Use the rule builder
5. Query chemical reactions and recommendations

**Status:** ✅ **FIXED** - The GUI application now runs successfully!
