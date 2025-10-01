from pathlib import Path
text = Path('ui_gradio_origin.py').read_text(encoding='utf-8')
print('origin ok length', len(text))
