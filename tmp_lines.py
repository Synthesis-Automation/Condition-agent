from pathlib import Path
for idx in [61,67,124,744,1285,1286,1829,1830,1856,1857,1860,2177,2178]:
    line = Path('app/ui_gradio.py').read_bytes().splitlines()[idx]
    print(idx+1, line)
