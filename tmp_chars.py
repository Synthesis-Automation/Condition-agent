from pathlib import Path
text = Path('app/ui_gradio.py').read_text(encoding='latin-1')
chars = sorted({ch for ch in text if ord(ch) >= 128})
print([ch.encode('unicode_escape').decode() for ch in chars])
