from pathlib import Path
text = Path('app/ui_gradio.py').read_text(encoding='latin-1')
for idx, line in enumerate(text.splitlines(), 1):
    if '\ufffd' in line:
        print(idx, line)
