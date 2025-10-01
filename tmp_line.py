from pathlib import Path
line_no = 519
line = Path('app/ui_gradio.py').read_text(encoding='latin-1').splitlines()[line_no-1]
print(line.encode('unicode_escape'))
