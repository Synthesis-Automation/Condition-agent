"""
Migrate protocol database files from v1 to v2.0 schema format.

This script:
1. Reads all protocol JSON files from data/protocol_db/
2. Converts them to v2.0 schema format with:
   - schema_version: "2.0"
   - source_type: "protocol"
   - metadata section (id, name, version, tags)
   - Ensures reaction_smiles exists for DRFP computation
3. Writes converted files to data/protocol_db_v2/

Usage:
    python scripts/migrate_protocols_to_v2.py
"""

import json
from pathlib import Path
from typing import Any, Dict, List
import re


def extract_metadata_from_protocol(protocol: Dict[str, Any]) -> Dict[str, Any]:
    """Extract metadata from protocol v1 format."""
    from datetime import datetime
    
    source = protocol.get("source", {})
    reaction = protocol.get("reaction", {})
    
    # Generate ID from title or family
    title = source.get("title", "")
    family = reaction.get("family", "Unknown")
    
    # Create a clean ID
    if title:
        # Use first 50 chars of title, make lowercase, replace spaces/special chars
        clean_title = re.sub(r'[^a-z0-9]+', '_', title[:50].lower()).strip('_')
        protocol_id = clean_title
    else:
        protocol_id = re.sub(r'[^a-z0-9]+', '_', family.lower()).strip('_')
    
    # Generate name
    name = title if title else family
    
    # Generate version from source info (must match semver: \d+.\d+.\d+)
    version = "1.0.0"
    if source.get("year"):
        # Use year as major version
        year = str(source['year'])
        if source.get("volume"):
            # Use volume as minor version
            vol = str(source['volume']).replace('.', '')  # Remove any dots
            version = f"{year}.{vol}.0"
        else:
            version = f"{year}.1.0"
    
    # Extract tags
    tags = []
    if reaction.get("tags"):
        # Parse tags string
        tags_str = reaction["tags"]
        if isinstance(tags_str, str):
            tags = [t.strip() for t in tags_str.split(";") if t.strip()]
        elif isinstance(tags_str, list):
            tags = tags_str
    
    # Add family as tag if not present
    if family and family not in tags:
        tags.append(family)
    
    # Get current timestamp for dates
    now = datetime.now().isoformat()
    
    return {
        "id": protocol_id,
        "name": name,
        "version": version,
        "description": source.get("title", ""),
        "tags": tags,
        "created_date": now,
        "last_modified": now,
        "source_reference": {
            k: v for k, v in source.items() 
            if v is not None and k not in ["reference_note"]
        }
    }


def ensure_reaction_smiles(protocol: Dict[str, Any]) -> str:
    """Ensure protocol has a valid reaction_smiles field."""
    reaction = protocol.get("reaction", {})
    
    # Check if reaction_smiles already exists
    if reaction.get("reaction_smiles"):
        return reaction["reaction_smiles"]
    
    # If no reaction_smiles, this is a problem for v2.0 protocols
    # Try to construct from reaction_setup if possible
    print(f"  ⚠️  WARNING: No reaction_smiles found in protocol: {protocol.get('source', {}).get('title', 'Unknown')}")
    
    return None


def convert_protocol_to_v2(protocol: Dict[str, Any]) -> Dict[str, Any]:
    """Convert a single protocol from v1 to v2.0 format."""
    
    # Check if already v2
    if protocol.get("schema_version") == "2.0":
        return protocol
    
    metadata = extract_metadata_from_protocol(protocol)
    reaction_smiles = ensure_reaction_smiles(protocol)
    
    # Build v2 structure
    v2_protocol = {
        "schema_version": "2.0",
        "source_type": "protocol",
        "metadata": metadata,
        "reaction": {
            "reaction_smiles": reaction_smiles,
            "family": protocol.get("reaction", {}).get("family", "Unknown"),
            "notes": protocol.get("reaction", {}).get("notes", ""),
        },
        "reaction_setup": protocol.get("reaction_setup", []),
    }
    
    # Add optional fields if present
    reaction = protocol.get("reaction", {})
    if reaction.get("reaction_SMARTS"):
        v2_protocol["reaction"]["reaction_SMARTS"] = reaction["reaction_SMARTS"]
    
    if reaction.get("compatible_functional_groups"):
        v2_protocol["reaction"]["compatible_functional_groups"] = reaction["compatible_functional_groups"]
    
    if reaction.get("incompatible_functional_groups"):
        v2_protocol["reaction"]["incompatible_functional_groups"] = reaction["incompatible_functional_groups"]
    
    # Add workup if present
    if protocol.get("workup_and_purification"):
        v2_protocol["workup_and_purification"] = protocol["workup_and_purification"]
    
    if protocol.get("original_procedure"):
        v2_protocol["original_procedure"] = protocol["original_procedure"]
    
    return v2_protocol


def migrate_protocols():
    """Migrate all protocol files to v2.0 format."""
    
    protocol_db = Path("data/protocol_db")
    output_dir = Path("data/protocol_db_v2")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all JSON files (exclude hidden files and already migrated)
    protocol_files = [
        f for f in protocol_db.glob("*.json")
        if not f.name.startswith(".") and "_v2.json" not in f.name
    ]
    
    print(f"Found {len(protocol_files)} protocol files to migrate\n")
    
    migrated_count = 0
    skipped_count = 0
    error_count = 0
    no_smiles_count = 0
    
    for protocol_file in sorted(protocol_files):
        try:
            print(f"Processing: {protocol_file.name}")
            
            with open(protocol_file, "r", encoding="utf-8") as f:
                data = json.load(f)
            
            # Handle both single protocol and array of protocols
            if isinstance(data, list):
                # Multiple protocols in one file
                v2_protocols = []
                has_no_smiles = False
                
                for idx, protocol in enumerate(data):
                    v2_protocol = convert_protocol_to_v2(protocol)
                    
                    # Check if reaction_smiles is missing
                    if not v2_protocol["reaction"]["reaction_smiles"]:
                        has_no_smiles = True
                        print(f"    Protocol {idx+1}: Missing reaction_smiles")
                    
                    v2_protocols.append(v2_protocol)
                
                # Only write if all protocols have reaction_smiles
                if not has_no_smiles:
                    output_path = output_dir / protocol_file.name
                    with open(output_path, "w", encoding="utf-8") as f:
                        json.dump(v2_protocols, f, indent=2, ensure_ascii=False)
                    
                    migrated_count += 1
                    print(f"  ✅ Migrated {len(v2_protocols)} protocols to: {output_path.name}")
                else:
                    no_smiles_count += 1
                    print(f"  ⏭️  Skipped (missing reaction_smiles): {protocol_file.name}")
            
            else:
                # Single protocol
                v2_protocol = convert_protocol_to_v2(data)
                
                # Check if reaction_smiles is missing
                if not v2_protocol["reaction"]["reaction_smiles"]:
                    no_smiles_count += 1
                    print(f"  ⏭️  Skipped (missing reaction_smiles): {protocol_file.name}")
                else:
                    output_path = output_dir / protocol_file.name
                    with open(output_path, "w", encoding="utf-8") as f:
                        json.dump([v2_protocol], f, indent=2, ensure_ascii=False)
                    
                    migrated_count += 1
                    print(f"  ✅ Migrated to: {output_path.name}")
        
        except Exception as e:
            error_count += 1
            print(f"  ❌ Error: {str(e)}")
    
    print(f"\n{'='*60}")
    print(f"Migration Summary:")
    print(f"  ✅ Migrated: {migrated_count} files")
    print(f"  ⏭️  Skipped (no reaction_smiles): {no_smiles_count} files")
    print(f"  ❌ Errors: {error_count} files")
    print(f"  📁 Total processed: {len(protocol_files)} files")
    print(f"{'='*60}\n")
    
    if no_smiles_count > 0:
        print("⚠️  Note: Files without reaction_smiles need manual review.")
        print("   These protocols cannot be used for DRFP-based matching.")
        print("   Consider adding reaction_smiles to these files manually.\n")


if __name__ == "__main__":
    migrate_protocols()
