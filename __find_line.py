import re
import sys
path = sys.argv[1]
pat = re.compile(r'^def\s+_enrich_props\(')
with open(path, 'r', encoding='utf-8') as f:
    for i, line in enumerate(f, 1):
        if pat.search(line):
            print(i)
