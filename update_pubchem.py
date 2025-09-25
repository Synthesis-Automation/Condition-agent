from pathlib import Path
import re
path = Path('data-processor/reagent_taxonomy_generator.py')
text = path.read_text(encoding='utf-8')
pattern = r"def _resolve_via_pubchem\(session: Any, cas: str, timeout: float\) -> Optional\[Dict\[str, Any\]\]:\r?\n(?:    .*\r?\n)*?(?=def _resolve_via_cactus)"
new_func = '''def _resolve_via_pubchem(session: Any, cas: str, timeout: float) -> Optional[Dict[str, Any]]:
    base = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/{quote(cas)}"
    props = _http_get_json(
        session,
        base + "/property/Title,IUPACName,IsomericSMILES,CanonicalSMILES/JSON",
        timeout,
    )
    entries: List[Dict[str, Any]] = []
    if props:
        entries = props.get("PropertyTable", {}).get("Properties", []) or []

    syn_url = base + "/synonyms/JSON"
    syn_data = _http_get_json(session, syn_url, timeout)
    synonyms_by_cid: Dict[int, List[str]] = {}
    extra_synonyms: List[str] = []
    if syn_data:
        info = syn_data.get("InformationList", {}).get("Information", []) or []
        for block in info:
            names = block.get("Synonym", []) or []
            if not names:
                continue
            cid = block.get("CID")
            if cid is None:
                extra_synonyms.extend(names)
            else:
                synonyms_by_cid.setdefault(int(cid), []).extend(names)

    normalized_target: Optional[str] = None
    try:
        normalized_target = normalize_cas(cas)
    except ValueError:
        normalized_target = None

    def contains_target(synonyms: Sequence[str]) -> bool:
        if not normalized_target:
            return False
        for candidate in synonyms:
            try:
                if normalize_cas(candidate) == normalized_target:
                    return True
            except ValueError:
                continue
        return False

    record: Optional[Dict[str, Any]] = None
    selected_synonyms: List[str] = []
    for candidate in entries:
        cid = candidate.get("CID")
        if cid is None:
            continue
        syns = synonyms_by_cid.get(int(cid), [])
        if contains_target(syns):
            record = candidate
            selected_synonyms = list(syns)
            break

    if record is None and entries:
        record = entries[0]
        cid = record.get("CID")
        if cid is not None:
            selected_synonyms = list(synonyms_by_cid.get(int(cid), []))

    if record is None:
        return None

    smiles = record.get("IsomericSMILES") or record.get("CanonicalSMILES")
    name = record.get("Title") or record.get("IUPACName")

    raw: List[str] = []
    if name:
        raw.append(name)
    raw.extend(selected_synonyms)
    raw.extend(extra_synonyms)
    deduped = dedupe_synonyms(raw)
    if deduped:
        name = name or deduped[0]
    if not name:
        return None
    if deduped and len(deduped) > 16:
        deduped = deduped[:16]
    return {"name": name, "synonyms": deduped or [name], "smiles": smiles}


'''
new_text, count = re.subn(pattern, new_func, text, flags=re.DOTALL)
if count != 1:
    raise SystemExit(f"Replacement count {count}, expected 1")
path.write_text(new_text, encoding='utf-8')
