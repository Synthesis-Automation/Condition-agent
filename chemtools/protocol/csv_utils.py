from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from chemtools.reagent.constants import ROLE_ALIASES

# Protocol CSV columns matching the standard protocol format (29 columns)
# Based on protocols_01.csv format with integrated reference information
PROTOCOL_CSV_COLUMNS = (
    "reaction_id",              # Unique reaction identifier
    "reaction_type",            # Reaction family/type
    "yield_pct",                # Overall yield percentage
    "temperature_c",            # Reaction temperature in Celsius
    "time_h",                   # Reaction time in hours
    "reaction_smiles",          # Reaction SMILES string
    "reference",                # Full reference with journal, volume, pages, year, DOI, URL
    "reactant_cas",             # CAS numbers of reactants (comma-separated)
    "product_cas",              # CAS numbers of products (comma-separated)
    "reagent_cas",              # CAS numbers of reagents (comma-separated)
    "catalyst_cas",             # CAS numbers of catalysts (comma-separated)
    "solvent_cas",              # CAS numbers of solvents (comma-separated)
    "reactant_amd",             # AMD IDs of reactants (comma-separated)
    "product_amd",              # AMD IDs of products (comma-separated)
    "reagent_amd",              # AMD IDs of reagents (comma-separated)
    "catalyst_amd",             # AMD IDs of catalysts (comma-separated)
    "solvent_amd",              # AMD IDs of solvents (comma-separated)
    "experimental_procedure",   # Full experimental procedure text
    "stages",                   # Number of reaction stages
    "steps",                    # Number of steps
    "product_yield_1",          # Yield for product 1
    "product_yield_2",          # Yield for product 2
    "product_yield_3",          # Yield for product 3
    "product_yield_4",          # Yield for product 4
    "product_yield_5",          # Yield for product 5
    "product_yield_6",          # Yield for product 6
    "product_yield_7",          # Yield for product 7
    "notes",                    # Additional notes
    "reaction_setup_json",      # JSON-encoded full reaction setup
)


def _clean_text(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text


def _split_delimited(value: Any, sep: str) -> List[str]:
    text = _clean_text(value)
    if not text:
        return []
    return [part.strip() for part in text.split(sep) if part.strip()]


def _split_slash(value: Any) -> List[str]:
    text = _clean_text(value)
    if not text:
        return []
    return [part.strip() for part in text.split("/") if part.strip()]


def parse_semicolon_list(value: Any) -> List[str]:
    items = _split_delimited(value, ";")
    return items


def serialize_semicolon_list(items: List[str]) -> str:
    clean_items = [item.strip() for item in items if item and str(item).strip()]
    return ";".join(clean_items)


def parse_source_field(value: Any) -> Dict[str, Any]:
    if isinstance(value, dict):
        return value
    text = _clean_text(value)
    if not text:
        return {}
    try:
        parsed = json.loads(text)
        if isinstance(parsed, dict):
            return parsed
        return {"title": text}
    except Exception:
        return {"title": text}


def extract_source(protocol_data: Dict[str, Any]) -> Dict[str, Any]:
    if not protocol_data:
        return {}
    source = protocol_data.get("source")
    if not source:
        metadata = protocol_data.get("metadata", {})
        source = metadata.get("source_reference") or metadata.get("source")
    return parse_source_field(source)


def extract_tags(protocol_data: Dict[str, Any]) -> List[str]:
    if not protocol_data:
        return []
    reaction = protocol_data.get("reaction", {}) or {}
    tags_raw = reaction.get("tags")
    if not tags_raw:
        tags_raw = protocol_data.get("metadata", {}).get("tags", [])
    if isinstance(tags_raw, list):
        return [str(tag).strip() for tag in tags_raw if str(tag).strip()]
    if isinstance(tags_raw, str):
        return parse_semicolon_list(tags_raw)
    return []


def extract_notes(protocol_data: Dict[str, Any]) -> str:
    if not protocol_data:
        return ""
    reaction = protocol_data.get("reaction", {}) or {}
    notes = reaction.get("notes")
    if notes:
        return _clean_text(notes)
    metadata = protocol_data.get("metadata", {}) or {}
    return _clean_text(metadata.get("description"))


def _parse_json_field(value: Any) -> Optional[Any]:
    if isinstance(value, (dict, list)):
        return value
    text = _clean_text(value)
    if not text:
        return None
    try:
        return json.loads(text)
    except Exception:
        return None


def _format_source_field(source: Dict[str, Any]) -> str:
    if not source:
        return ""
    try:
        return json.dumps(source, ensure_ascii=True, separators=(",", ":"))
    except Exception:
        return ""


def _extract_cas_numbers(chemicals: List[Dict[str, Any]], role_filter: Optional[List[str]] = None) -> str:
    """Extract CAS numbers from chemicals list, optionally filtered by role."""
    cas_numbers = []
    for chem in chemicals:
        role = _clean_text(chem.get("role", ""))
        cas = _clean_text(chem.get("cas", ""))
        
        if cas and (role_filter is None or role in role_filter):
            cas_numbers.append(cas)
    
    return ", ".join(cas_numbers) if cas_numbers else ""


def _extract_temperature(protocol_data: Dict[str, Any]) -> str:
    """Extract temperature from protocol conditions."""
    reaction_setup = protocol_data.get("reaction_setup", [])
    if not reaction_setup or not isinstance(reaction_setup, list):
        return ""
    
    # Get conditions from first setup step
    first_setup = reaction_setup[0] if reaction_setup else {}
    conditions = first_setup.get("conditions", []) if isinstance(first_setup, dict) else []
    
    if conditions and len(conditions) > 0:
        # Use the main reaction conditions (usually the last or longest step)
        main_cond = conditions[-1] if len(conditions) > 1 else conditions[0]
        temp = main_cond.get("temperature_C")
        if temp is not None:
            return str(temp)
    
    return ""


def _extract_time(protocol_data: Dict[str, Any]) -> str:
    """Extract time from protocol conditions."""
    reaction_setup = protocol_data.get("reaction_setup", [])
    if not reaction_setup or not isinstance(reaction_setup, list):
        return ""
    
    # Get conditions from first setup step
    first_setup = reaction_setup[0] if reaction_setup else {}
    conditions = first_setup.get("conditions", []) if isinstance(first_setup, dict) else []
    
    if conditions and len(conditions) > 0:
        # Use the main reaction conditions (usually the last or longest step)
        main_cond = conditions[-1] if len(conditions) > 1 else conditions[0]
        time_h = main_cond.get("time_h")
        if time_h is not None:
            return str(time_h)
    
    return ""


def _count_stages(protocol_data: Dict[str, Any]) -> str:
    """Count number of reaction stages."""
    reaction_setup = protocol_data.get("reaction_setup", [])
    if reaction_setup and isinstance(reaction_setup, list):
        return str(len(reaction_setup))
    return ""


def _count_steps(protocol_data: Dict[str, Any]) -> str:
    """Count total number of steps across all stages."""
    reaction_setup = protocol_data.get("reaction_setup", [])
    if not reaction_setup or not isinstance(reaction_setup, list):
        return ""
    
    total_steps = 0
    for setup in reaction_setup:
        if isinstance(setup, dict):
            conditions = setup.get("conditions", [])
            if conditions and isinstance(conditions, list):
                total_steps += len(conditions)
    
    return str(total_steps) if total_steps > 0 else ""


def _format_integrated_reference(source: Dict[str, Any]) -> str:
    """Format all reference details into a single integrated string."""
    if not source:
        return ""
    
    parts = []
    
    # Title
    title = _clean_text(source.get("title"))
    if title:
        parts.append(title)
    
    # Journal, Volume, Pages (Year)
    journal_parts = []
    journal = _clean_text(source.get("journal"))
    if journal:
        journal_parts.append(journal)
    
    volume = source.get("volume")
    if volume:
        journal_parts.append(f"Vol. {volume}")
    
    pages = _clean_text(source.get("pages"))
    if pages:
        journal_parts.append(f"pp. {pages}")
    
    year = source.get("year")
    if journal_parts and year:
        parts.append(f"{', '.join(journal_parts)} ({year})")
    elif journal_parts:
        parts.append(", ".join(journal_parts))
    elif year:
        parts.append(f"({year})")
    
    # DOI
    doi = _clean_text(source.get("doi"))
    if doi:
        parts.append(f"DOI: {doi}")
    
    # URL
    url = _clean_text(source.get("url"))
    if url:
        parts.append(f"URL: {url}")
    
    return ". ".join(parts) if parts else ""


def row_to_protocol(row: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert CSV row to protocol JSON format.
    
    Args:
        row: CSV row as dictionary
        
    Returns:
        Protocol dictionary in JSON format
    """
    reaction_smiles = _clean_text(row.get("reaction_smiles"))
    family = _clean_text(row.get("reaction_type"))
    notes = _clean_text(row.get("notes"))
    protocol_id = _clean_text(row.get("reaction_id"))
    reference = _clean_text(row.get("reference"))

    # Parse reaction_setup_json if present
    reaction_setup = _parse_json_field(row.get("reaction_setup_json"))
    if not isinstance(reaction_setup, list):
        reaction_setup = []

    # Build source metadata from integrated reference string
    # The reference is now a formatted string, so we store it as title
    # For now, keep it simple - full parsing would be complex
    source: Dict[str, Any] = {}
    if reference:
        source["title"] = reference

    # Build metadata
    metadata: Dict[str, Any] = {
        "id": protocol_id,
        "source_reference": source,
    }

    # Build reaction object
    reaction: Dict[str, Any] = {
        "reaction_smiles": reaction_smiles,
        "family": family,
        "notes": notes,
    }

    # Build protocol data
    protocol_data: Dict[str, Any] = {
        "schema_version": "2.0",
        "source_type": "protocol",
        "metadata": metadata,
        "reaction": reaction,
        "reaction_setup": reaction_setup,
        "original_procedure": _clean_text(row.get("experimental_procedure")),
        "source": source,
    }

    return protocol_data


def protocol_to_row(protocol_data: Dict[str, Any]) -> Dict[str, str]:
    """
    Convert protocol JSON to CSV row matching protocols_01.csv format.
    
    Args:
        protocol_data: Protocol dictionary in JSON format
        
    Returns:
        CSV row as dictionary
    """
    reaction = protocol_data.get("reaction", {}) or {}
    metadata = protocol_data.get("metadata", {}) or {}
    source = extract_source(protocol_data)

    reaction_smiles = _clean_text(reaction.get("reaction_smiles"))
    family = _clean_text(reaction.get("family"))
    notes = extract_notes(protocol_data)
    protocol_id = _clean_text(metadata.get("id")) or _clean_text(metadata.get("protocol_id"))

    # Get reaction setup
    reaction_setup = protocol_data.get("reaction_setup", [])
    reaction_setup_json = ""
    if reaction_setup:
        reaction_setup_json = json.dumps(reaction_setup, ensure_ascii=True, separators=(",", ":"))

    # Extract chemicals from reaction setup
    chemicals = []
    if reaction_setup and isinstance(reaction_setup, list):
        first_setup = reaction_setup[0] if reaction_setup else {}
        chemicals = first_setup.get("chemicals", []) if isinstance(first_setup, dict) else []

    # Extract CAS numbers by role
    reactant_roles = ["starting_material", "reactant"]
    product_roles = ["product"]
    reagent_roles = ["reagent", "base", "acid", "oxidant", "reductant", "additive", "condensation_agent", "other_reagent"]
    catalyst_roles = ["metal_catalyst", "co_catalyst", "catalyst", "ligand"]
    solvent_roles = ["solvent"]

    reactant_cas = _extract_cas_numbers(chemicals, reactant_roles)
    product_cas = _extract_cas_numbers(chemicals, product_roles)
    reagent_cas = _extract_cas_numbers(chemicals, reagent_roles)
    catalyst_cas = _extract_cas_numbers(chemicals, catalyst_roles)
    solvent_cas = _extract_cas_numbers(chemicals, solvent_roles)

    # Extract experimental parameters
    temperature_c = _extract_temperature(protocol_data)
    time_h = _extract_time(protocol_data)
    stages = _count_stages(protocol_data)
    steps = _count_steps(protocol_data)

    # Get experimental procedure
    experimental_procedure = _clean_text(protocol_data.get("original_procedure"))

    # Format integrated reference string from all source details
    integrated_reference = _format_integrated_reference(source)

    # Initialize row with all columns
    row: Dict[str, str] = {col: "" for col in PROTOCOL_CSV_COLUMNS}
    
    # Populate row
    row.update({
        "reaction_id": protocol_id,
        "reaction_type": family,
        "yield_pct": "",  # Not typically in protocol data
        "temperature_c": temperature_c,
        "time_h": time_h,
        "reaction_smiles": reaction_smiles,
        "reference": integrated_reference,
        "reactant_cas": reactant_cas,
        "product_cas": product_cas,
        "reagent_cas": reagent_cas,
        "catalyst_cas": catalyst_cas,
        "solvent_cas": solvent_cas,
        "reactant_amd": "",  # AMD IDs not in protocol data
        "product_amd": "",
        "reagent_amd": "",
        "catalyst_amd": "",
        "solvent_amd": "",
        "experimental_procedure": experimental_procedure,
        "stages": stages,
        "steps": steps,
        "product_yield_1": "",  # Multi-product yields not in protocol data
        "product_yield_2": "",
        "product_yield_3": "",
        "product_yield_4": "",
        "product_yield_5": "",
        "product_yield_6": "",
        "product_yield_7": "",
        "notes": notes,
        "reaction_setup_json": reaction_setup_json,
    })

    return row


def read_protocol_csv(csv_path: str) -> List[Dict[str, Any]]:
    """
    Read protocol CSV file and convert to protocol JSON format.
    
    Args:
        csv_path: Path to CSV file
        
    Returns:
        List of protocol dictionaries in JSON format
    """
    import csv
    
    protocols = []
    with open(csv_path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            protocol = row_to_protocol(row)
            protocols.append(protocol)
    
    return protocols


def write_protocol_csv(
    protocols: List[Dict[str, Any]], 
    csv_path: str, 
    mode: str = "w"
) -> None:
    """
    Write protocols to CSV file.
    
    Args:
        protocols: List of protocol dictionaries in JSON format
        csv_path: Path to output CSV file
        mode: File mode ('w' for write, 'a' for append)
    """
    import csv
    
    if not protocols:
        return
    
    # Convert protocols to rows
    rows = [protocol_to_row(protocol) for protocol in protocols]
    
    # Write CSV
    write_header = mode == "w"
    with open(csv_path, mode, encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PROTOCOL_CSV_COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


def load_multiple_csvs(csv_paths: List[str]) -> List[Dict[str, Any]]:
    """
    Load and combine multiple protocol CSV files.
    
    Args:
        csv_paths: List of paths to CSV files
        
    Returns:
        Combined list of protocol dictionaries
    """
    all_protocols = []
    for csv_path in csv_paths:
        protocols = read_protocol_csv(csv_path)
        all_protocols.extend(protocols)
    
    return all_protocols


def find_protocol_csvs(directory: str) -> List[str]:
    """
    Find all protocol CSV files in a directory (non-recursive).
    
    Args:
        directory: Directory path to search
        
    Returns:
        List of CSV file paths
    """
    from pathlib import Path
    
    dir_path = Path(directory)
    if not dir_path.exists():
        return []
    
    csv_files = []
    for csv_path in dir_path.glob("*.csv"):
        # Skip hidden files and temp files
        if csv_path.name.startswith(".") or csv_path.name.startswith("~"):
            continue
        csv_files.append(str(csv_path))
    
    return sorted(csv_files)

