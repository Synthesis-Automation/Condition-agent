import sys
path = sys.argv[1]
sym = sys.argv[2]
with open(path, 'r', encoding='utf-8') as f:
    for i, line in enumerate(f, 1):
        if line.startswith(sym):
            print(i)
            break
