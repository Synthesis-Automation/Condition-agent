from pathlib import Path
text = Path('app/ui_gradio.py').read_text(encoding='latin-1')
for idx, line in enumerate(text.splitlines(), 1):
    non_ascii = [ch for ch in line if ord(ch) >= 128]
    if non_ascii:
        encoded = ''.join(non_ascii).encode('unicode_escape').decode()
        print(idx, encoded)
