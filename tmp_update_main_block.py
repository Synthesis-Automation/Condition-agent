from pathlib import Path
from textwrap import dedent

path = Path('data-processor/reagent_taxonomy_generator.py')
text = path.read_text(encoding='utf-8')
old = dedent('''
    if not args.cas or not args.name:
        raise SystemExit("--cas and --name are required unless --list-families is used.")

    cas = normalize_cas(args.cas)
    abbr = args.abbr or args.name
    synonyms = dedupe_synonyms([args.name, abbr, *args.synonym])
    role = args.role
    family_id = args.family
    used_default = False
    family_reason: Optional[List[str]] = None
    role_reason: Optional[str] = None

    if family_id:
        family_role = store.role_for_family(family_id)
''')
new = dedent('''
    if not args.cas:
        raise SystemExit("--cas is required unless --list-families is used.")

    cas = normalize_cas(args.cas)
    resolved_identity: Optional[Dict[str, Any]] = None
    if not args.no_auto_resolve and not args.name:
        resolved_identity = resolve_identity_from_cas(cas, timeout=args.resolver_timeout)
        if args.verbose:
            if resolved_identity:
                print(f"# auto-resolved via {resolved_identity.get('source')}: {resolved_identity.get('name')}")
            else:
                print("# auto-resolve failed to supply a name")

    name = args.name or (resolved_identity.get("name") if resolved_identity else None)
    if not name:
        raise SystemExit("Unable to determine reagent name. Provide --name or use --no-auto-resolve to skip lookup.")

    abbr = args.abbr or name
    resolved_synonyms = resolved_identity.get("synonyms", []) if resolved_identity else []
    synonyms = dedupe_synonyms([name, abbr, *args.synonym, *resolved_synonyms])
    role = args.role
    family_id = args.family
    used_default = False
    family_reason: Optional[List[str]] = None
    role_reason: Optional[str] = None
    auto_resolve_source = resolved_identity.get("source") if resolved_identity else None
    resolved_smiles = resolved_identity.get("smiles") if resolved_identity else None

    if family_id:
        family_role = store.role_for_family(family_id)
''')
if old not in text:
    raise SystemExit('block to replace not found')
text = text.replace(old, new, 1)
path.write_text(text, encoding='utf-8')
