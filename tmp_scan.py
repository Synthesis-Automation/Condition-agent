from pathlib import Path
import re
path = Path('app/ui_gradio.py')
data = path.read_bytes()
for match in re.finditer(b'\xe2\x80.', data):
    seq = match.group(0)
    if seq.endswith(b'?'):
        start = match.start()
        print('bad sequence at', start, seq)
