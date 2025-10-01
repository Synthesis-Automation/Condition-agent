from pathlib import Path
path = Path('app/ui_gradio.py')
data = path.read_bytes()
old = '            "?" if rec.get("ok", True) else "??",'.encode('utf-8')
new = '            "OK" if rec.get("ok", True) else "WARN",'.encode('utf-8')
if old in data:
    data = data.replace(old, new)
else:
    raise SystemExit('pattern not found')
path.write_bytes(data)
