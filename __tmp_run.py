import json
from chemtools.registry import resolve
print(json.dumps(resolve("108-88-3"), ensure_ascii=False, indent=2))
