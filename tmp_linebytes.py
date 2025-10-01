from pathlib import Path
line = Path('app/ui_gradio.py').read_bytes().splitlines()[61]
print(line)
