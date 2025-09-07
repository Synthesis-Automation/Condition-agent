from typing import Dict, Any, List
def apply_filter(candidates: List[str], rules: Dict[str, Any]) -> Dict[str, Any]:
    rules = rules or {}; inventory=set([x.lower() for x in rules.get("inventory", [])]); blacklist=set([x.lower() for x in rules.get("blacklist", [])])
    allowed=[]; blocked=[]
    for c in candidates:
        lc=c.lower()
        if blacklist and lc in blacklist: blocked.append({"id":c,"reason":"blacklisted"}); continue
        if inventory and lc not in inventory: blocked.append({"id":c,"reason":"not_in_inventory"}); continue
        allowed.append(c)
    return {"allowed": allowed, "blocked": blocked}
