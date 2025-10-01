from pathlib import Path
line = Path('app/ui_gradio.py').read_text(encoding='utf-8').splitlines()[1331]
print(repr(line))
