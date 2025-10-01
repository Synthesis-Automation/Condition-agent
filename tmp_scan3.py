from pathlib import Path
with Path('app/ui_gradio.py').open('rb') as f:
    data = f.read()
for idx, line in enumerate(data.splitlines(), 1):
    if b'\xe2' in line or b'?' in line:
        if b'\xe2' in line and b'?' in line:
            pos = line.find(b'\xe2')
            if pos >=0 and pos+2 < len(line) and line[pos:pos+3] in (b'\xe2\x80?', b'\xe2\x80?'):
                print(idx, line)
