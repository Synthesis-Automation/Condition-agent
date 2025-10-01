from pathlib import Path
path = Path('app/ui_gradio.py')
data = path.read_bytes()
patterns = [
    (b'\xe9\x88\xa5?', b'- '),
    (b'\xe2\x80?', b'- '),
]
for old, new in patterns:
    if old in data:
        data = data.replace(old, new)
path.write_bytes(data)
