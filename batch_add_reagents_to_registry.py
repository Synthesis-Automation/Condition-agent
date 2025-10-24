#!/usr/bin/env python3
"""
Batch process reagents from the mapped markdown file and add them to the reagent registry.
Searches for CAS numbers using reagent names, skips existing entries.
"""

import json
import re
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

# Add root directory to path
ROOT_DIR = Path(__file__).resolve().parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chemtools.reagent import (
    normalize_cas,
    resolve_identity_from_cas,
    dedupe_synonyms,
)

# Import from reagent_taxonomy_ui
from app.reagent_taxonomy_ui import (
    ReagentRegistryStore,
    generate_taxonomy_entry_llm,
    ROLE_CONFIG,
)

# Check if LLM support is available
try:
    from llmtools.clients import LLMClient
    LLM_AVAILABLE = True
except ImportError:
    LLMClient = None
    LLM_AVAILABLE = False


def search_cas_by_name(name: str, timeout: float = 6.0) -> Optional[Tuple[str, Dict[str, Any]]]:
    """
    Search for CAS number using reagent name via PubChem.
    
    Returns:
        Tuple of (cas, identity_dict) if found, None otherwise
    """
    try:
        # Try direct PubChem search by name
        import requests
        
        # Search PubChem by name
        search_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/IUPACName,MolecularFormula,CanonicalSMILES,InChIKey/JSON"
        
        # Clean the name for URL
        clean_name = name.strip().replace(" ", "%20")
        url = search_url.format(clean_name)
        
        response = requests.get(url, timeout=timeout)
        
        if response.status_code == 200:
            data = response.json()
            if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
                props = data["PropertyTable"]["Properties"][0]
                
                # Now get CAS number using CID
                cid = props.get("CID")
                if cid:
                    cas_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
                    cas_response = requests.get(cas_url, timeout=timeout)
                    
                    if cas_response.status_code == 200:
                        cas_data = cas_response.json()
                        if "InformationList" in cas_data and "Information" in cas_data["InformationList"]:
                            synonyms = cas_data["InformationList"]["Information"][0].get("Synonym", [])
                            
                            # Find CAS number in synonyms
                            cas_pattern = re.compile(r'^\d{2,7}-\d{2}-\d$')
                            for syn in synonyms:
                                if cas_pattern.match(str(syn)):
                                    cas = normalize_cas(str(syn))
                                    
                                    # Build identity dict
                                    identity = {
                                        "cas": cas,
                                        "name": name,
                                        "smiles": props.get("CanonicalSMILES"),
                                        "inchi_key": props.get("InChIKey"),
                                        "synonyms": synonyms[:50],  # Limit to 50
                                        "source": "pubchem_name_search"
                                    }
                                    
                                    return cas, identity
        
        return None
        
    except Exception as e:
        print(f"  Error searching for '{name}': {e}")
        return None


def parse_reagents_from_md(md_file: Path) -> List[Dict[str, Any]]:
    """
    Parse reagents from the mapped markdown file.
    
    Returns:
        List of reagent dicts with name, registry_role, occurrences, reactions
    """
    content = md_file.read_text(encoding='utf-8')
    lines = content.split('\n')
    
    reagents = []
    current_role = None
    current_reagent = None
    
    for line in lines:
        # Detect role section
        if line.startswith("## ") and " reagents)" in line and "**Registry Role:**" not in line:
            # Extract registry role from next few lines
            continue
        
        # Get registry role from description line
        if line.startswith("**Registry Role:** `"):
            match = re.search(r'`([^`]+)`', line)
            if match:
                current_role = match.group(1)
        
        # New reagent
        if line.startswith("### ") and current_role:
            # Save previous reagent
            if current_reagent:
                reagents.append(current_reagent)
            
            reagent_name = line[4:].strip()
            current_reagent = {
                'name': reagent_name,
                'registry_role': current_role,
                'original_role': None,
                'occurrences': 0,
                'reactions': []
            }
        
        # Parse reagent details
        elif line.startswith("- **Original Role:**") and current_reagent:
            current_reagent['original_role'] = line.split("**Original Role:**")[1].strip()
        
        elif line.startswith("- **Occurrences:**") and current_reagent:
            match = re.search(r'\d+', line)
            if match:
                current_reagent['occurrences'] = int(match.group())
        
        elif "**Reaction types:**" in line and current_reagent:
            # Extract reaction types
            reactions_str = line.split("**Reaction types:**")[1].strip()
            # Remove "... and X more" suffix
            reactions_str = re.sub(r',?\s*\.\.\.\s*and\s+\d+\s+more', '', reactions_str)
            current_reagent['reactions'] = [r.strip() for r in reactions_str.split(',')]
    
    # Add last reagent
    if current_reagent:
        reagents.append(current_reagent)
    
    return reagents


def reagent_exists_in_registry(
    store: ReagentRegistryStore,
    name: str,
    role: str,
    cas: Optional[str] = None
) -> bool:
    """
    Check if a reagent already exists in the registry.
    
    Args:
        store: Registry store instance
        name: Reagent name
        role: Registry role
        cas: Optional CAS number to check
    
    Returns:
        True if reagent exists, False otherwise
    """
    # Check by CAS if provided
    if cas:
        normalized_cas = normalize_cas(cas)
        existing = store.find_by_cas(normalized_cas, role=role)
        if existing:
            return True
    
    # Check by name in the role's entries
    entries = store.role_entries.get(role, [])
    name_lower = name.lower()
    
    for entry in entries:
        # Check primary name
        if entry.get("name", "").lower() == name_lower:
            return True
        
        # Check abbreviations
        abbreviations = entry.get("abbreviation", [])
        if isinstance(abbreviations, list):
            if any(abbr.lower() == name_lower for abbr in abbreviations):
                return True
        
        # Check aliases
        aliases = entry.get("aliases", [])
        if isinstance(aliases, list):
            if any(alias.lower() == name_lower for alias in aliases):
                return True
    
    return False


def batch_add_reagents(
    md_file: Path,
    registry_dir: Path,
    llm_client: Optional[Any] = None,
    use_llm: bool = False,
    dry_run: bool = False,
    max_reagents: Optional[int] = None,
    min_occurrences: int = 1,
    skip_roles: Optional[Set[str]] = None,
    delay_seconds: float = 1.0
) -> Dict[str, Any]:
    """
    Batch process reagents from markdown file and add to registry.
    
    Args:
        md_file: Path to reagents_mapped_to_registry_roles.md
        registry_dir: Path to reagent registry directory
        llm_client: Optional LLM client for enhanced processing
        use_llm: Whether to use LLM workflow
        dry_run: If True, don't save to registry
        max_reagents: Maximum number of reagents to process
        min_occurrences: Minimum occurrence count to include
        skip_roles: Set of roles to skip
        delay_seconds: Delay between API calls to avoid rate limiting
    
    Returns:
        Summary dict with statistics
    """
    if not md_file.exists():
        raise FileNotFoundError(f"Markdown file not found: {md_file}")
    
    if not registry_dir.exists():
        raise FileNotFoundError(f"Registry directory not found: {registry_dir}")
    
    print(f"Loading reagents from {md_file}...")
    reagents = parse_reagents_from_md(md_file)
    print(f"Found {len(reagents)} reagents in markdown file")
    
    # Filter by occurrences
    reagents = [r for r in reagents if r['occurrences'] >= min_occurrences]
    print(f"After filtering by min occurrences ({min_occurrences}): {len(reagents)} reagents")
    
    # Filter by roles
    if skip_roles:
        reagents = [r for r in reagents if r['registry_role'] not in skip_roles]
        print(f"After filtering roles: {len(reagents)} reagents")
    
    # Limit max reagents
    if max_reagents:
        reagents = reagents[:max_reagents]
        print(f"Limited to {max_reagents} reagents")
    
    # Load registry
    print(f"\nLoading registry from {registry_dir}...")
    store = ReagentRegistryStore(registry_dir)
    
    # Track results
    results = {
        "total_processed": 0,
        "already_exists": 0,
        "cas_not_found": 0,
        "added_successfully": 0,
        "failed": 0,
        "skipped": 0,
        "details": []
    }
    
    print(f"\n{'='*80}")
    print("Starting batch processing...")
    print(f"{'='*80}\n")
    
    for i, reagent in enumerate(reagents, 1):
        name = reagent['name']
        role = reagent['registry_role']
        occurrences = reagent['occurrences']
        
        print(f"[{i}/{len(reagents)}] Processing: {name} (role: {role}, occurrences: {occurrences})")
        
        results["total_processed"] += 1
        
        # Check if reagent already exists (by name first, without CAS)
        if reagent_exists_in_registry(store, name, role):
            print(f"  ✓ Already exists in registry (by name), skipping")
            results["already_exists"] += 1
            results["details"].append({
                "name": name,
                "role": role,
                "status": "already_exists",
                "reason": "Found by name in registry"
            })
            continue
        
        # Search for CAS number
        print(f"  Searching for CAS number...")
        search_result = search_cas_by_name(name)
        
        if not search_result:
            print(f"  ✗ CAS number not found, skipping")
            results["cas_not_found"] += 1
            results["details"].append({
                "name": name,
                "role": role,
                "status": "cas_not_found",
                "reason": "PubChem search returned no results"
            })
            time.sleep(delay_seconds)  # Rate limiting
            continue
        
        cas, identity = search_result
        print(f"  ✓ Found CAS: {cas}")
        
        # Check again by CAS
        if reagent_exists_in_registry(store, name, role, cas):
            print(f"  ✓ Already exists in registry (by CAS), skipping")
            results["already_exists"] += 1
            results["details"].append({
                "name": name,
                "role": role,
                "cas": cas,
                "status": "already_exists",
                "reason": "Found by CAS in registry"
            })
            time.sleep(delay_seconds)
            continue
        
        # Add to registry
        try:
            if use_llm and llm_client:
                print(f"  Using LLM workflow to generate entry...")
                # Use LLM workflow
                result = generate_taxonomy_entry_llm(
                    cas=cas,
                    registry_dir=registry_dir,
                    llm_client=llm_client,
                    name_override=name,
                    role_override=role,
                    resolver_timeout=6.0
                )
                
                result_status = result.get("status")
                
                if result_status == "ready_to_save" or result_status == "needs_review":
                    entry = result.get("entry")
                    if entry and not dry_run:
                        store.add_entry(role, entry)
                        store.save_role(role)
                        status_msg = "Added successfully (LLM workflow)"
                        if result_status == "needs_review":
                            status_msg += f" [with warnings: {result.get('message', '')}]"
                        print(f"  ✓ {status_msg}")
                        results["added_successfully"] += 1
                        results["details"].append({
                            "name": name,
                            "role": role,
                            "cas": cas,
                            "status": "added",
                            "method": "llm_workflow",
                            "warnings": result.get("message") if result_status == "needs_review" else None
                        })
                    elif dry_run:
                        status_msg = "Would add (dry run, LLM workflow)"
                        if result_status == "needs_review":
                            status_msg += f" [with warnings: {result.get('message', '')}]"
                        print(f"  ✓ {status_msg}")
                        results["added_successfully"] += 1
                        results["details"].append({
                            "name": name,
                            "role": role,
                            "cas": cas,
                            "status": "would_add",
                            "method": "llm_workflow",
                            "warnings": result.get("message") if result_status == "needs_review" else None
                        })
                else:
                    # Get detailed error message
                    result_status = result.get("status", "unknown")
                    error_msg = result.get('error')
                    
                    if not error_msg:
                        # No explicit error - show status and message
                        message = result.get('message', '')
                        error_msg = f"Status '{result_status}'"
                        if message:
                            error_msg += f": {message}"
                        
                        # Try to extract from workflow steps
                        workflow = result.get('workflow', {})
                        for step_name, step_data in workflow.items():
                            if isinstance(step_data, dict) and step_data.get('status') != 'success':
                                step_error = step_data.get('error', 'Failed without details')
                                error_msg += f" | {step_name}: {step_error}"
                    
                    print(f"  ✗ LLM workflow failed: {error_msg}")
                    results["failed"] += 1
                    results["details"].append({
                        "name": name,
                        "role": role,
                        "cas": cas,
                        "status": "failed",
                        "error": error_msg,
                        "llm_result_status": result_status,
                        "llm_message": result.get("message")
                    })
            
            else:
                # Simple deterministic workflow - create minimal entry
                print(f"  Adding with deterministic workflow...")
                
                # Get or infer family
                family_id = store.default_family(role)
                if not family_id:
                    print(f"  ✗ No default family for role '{role}', skipping")
                    results["skipped"] += 1
                    results["details"].append({
                        "name": name,
                        "role": role,
                        "cas": cas,
                        "status": "skipped",
                        "reason": "No default family for role"
                    })
                    continue
                
                # Build role payload
                role_payload = store.build_role_payload(role, family_id)
                
                # Create entry
                entry_id = identity.get("inchi_key") or f"cas-{cas}"
                
                # Infer abbreviations
                synonyms = identity.get("synonyms", [])
                abbreviations = [name] if len(name) <= 10 else []
                
                # Build entry
                entry = {
                    "id": entry_id,
                    "name": name,
                    "abbreviation": abbreviations if abbreviations else None,
                    "aliases": synonyms[:20] if synonyms else [],  # Limit aliases
                    "cas": cas,
                    "inchi_key": identity.get("inchi_key"),
                    "smiles": identity.get("smiles"),
                    "roles": {
                        role: role_payload
                    }
                }
                
                # Build embedding text
                family_entry = store.family_entry(role, family_id)
                if family_entry:
                    from app.reagent_taxonomy_ui import build_embedding_text
                    entry["embedding_text"] = build_embedding_text(
                        role, family_entry, entry, synonyms
                    )
                
                # Remove None values
                entry = {k: v for k, v in entry.items() if v is not None}
                
                if not dry_run:
                    store.add_entry(role, entry)
                    store.save_role(role)
                    print(f"  ✓ Added successfully (deterministic)")
                    results["added_successfully"] += 1
                    results["details"].append({
                        "name": name,
                        "role": role,
                        "cas": cas,
                        "family": family_id,
                        "status": "added",
                        "method": "deterministic"
                    })
                else:
                    print(f"  ✓ Would add (dry run)")
                    results["added_successfully"] += 1
                    results["details"].append({
                        "name": name,
                        "role": role,
                        "cas": cas,
                        "family": family_id,
                        "status": "would_add",
                        "method": "deterministic"
                    })
        
        except Exception as e:
            print(f"  ✗ Error: {e}")
            results["failed"] += 1
            results["details"].append({
                "name": name,
                "role": role,
                "cas": cas if 'cas' in locals() else None,
                "status": "failed",
                "error": str(e)
            })
        
        # Rate limiting delay
        time.sleep(delay_seconds)
        print()
    
    return results


def print_summary(results: Dict[str, Any]) -> None:
    """Print summary statistics."""
    print(f"\n{'='*80}")
    print("BATCH PROCESSING SUMMARY")
    print(f"{'='*80}\n")
    
    print(f"Total processed:      {results['total_processed']}")
    print(f"Already exists:       {results['already_exists']}")
    print(f"CAS not found:        {results['cas_not_found']}")
    print(f"Added successfully:   {results['added_successfully']}")
    print(f"Failed:               {results['failed']}")
    print(f"Skipped:              {results['skipped']}")
    
    print(f"\n{'='*80}\n")


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Batch add reagents from mapped markdown to registry"
    )
    parser.add_argument(
        "--md-file",
        type=Path,
        default=ROOT_DIR / "reagents_mapped_to_registry_roles.md",
        help="Path to mapped reagents markdown file"
    )
    parser.add_argument(
        "--registry-dir",
        type=Path,
        default=ROOT_DIR / "data" / "reagents",
        help="Path to reagent registry directory"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Don't actually save to registry (preview only)"
    )
    parser.add_argument(
        "--max-reagents",
        type=int,
        help="Maximum number of reagents to process"
    )
    parser.add_argument(
        "--min-occurrences",
        type=int,
        default=10,
        help="Minimum occurrence count to include (default: 10)"
    )
    parser.add_argument(
        "--skip-roles",
        nargs="+",
        help="Roles to skip (e.g., solvent ligand)"
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=1.0,
        help="Delay between API calls in seconds (default: 1.0)"
    )
    parser.add_argument(
        "--use-llm",
        action="store_true",
        help="Use LLM workflow for enhanced processing"
    )
    parser.add_argument(
        "--llm-provider",
        default="aliyun",
        help="LLM provider (default: aliyun)"
    )
    parser.add_argument(
        "--llm-model",
        default="deepseek-r1-distill-qwen-7b",
        help="LLM model (default: deepseek-r1-distill-qwen-7b)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Save results to JSON file"
    )
    
    args = parser.parse_args()
    
    # Setup LLM if requested
    llm_client = None
    if args.use_llm:
        if not LLM_AVAILABLE:
            print("Error: LLM support not available. Install llmtools dependencies.")
            sys.exit(1)
        
        print(f"Initializing LLM client (provider: {args.llm_provider}, model: {args.llm_model})...")
        try:
            llm_client = LLMClient(provider=args.llm_provider, model=args.llm_model)
        except Exception as e:
            print(f"Error initializing LLM client: {e}")
            sys.exit(1)
    
    # Convert skip_roles to set
    skip_roles = set(args.skip_roles) if args.skip_roles else None
    
    # Run batch processing
    try:
        results = batch_add_reagents(
            md_file=args.md_file,
            registry_dir=args.registry_dir,
            llm_client=llm_client,
            use_llm=args.use_llm,
            dry_run=args.dry_run,
            max_reagents=args.max_reagents,
            min_occurrences=args.min_occurrences,
            skip_roles=skip_roles,
            delay_seconds=args.delay
        )
        
        # Print summary
        print_summary(results)
        
        # Save results to JSON if requested
        if args.output:
            args.output.write_text(json.dumps(results, indent=2, ensure_ascii=False))
            print(f"Results saved to: {args.output}")
    
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
