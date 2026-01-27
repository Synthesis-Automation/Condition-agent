from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from chemtools.reagent.constants import ROLE_ALIASES

STANDARD_LITERATURE_COLUMNS = (
    "reaction_id",
    "detected_reaction_type",
    "reaction_smiles",
    "yield",
    "z_score",
    "reactant_1",
    "reactant_2",
    "reactant_3",
    "formed_motifs",
    "catalyst",
    "ligand",
    "base",
    "acid",
    "oxidant",
    "reductant",
    "additive",
    "condensation_agent",
    "other_reagent",
    "solvent",
    "spectator_groups",
    "reference",
    "Reaction_Key",
)

EXTRA_PROTOCOL_COLUMNS = (
    "protocol_id",
    "reaction_smarts",
    "tags",
    "notes",
    "source",
    "reaction_setup_json",
    "original_procedure",
)

PROTOCOL_CSV_COLUMNS = STANDARD_LITERATURE_COLUMNS + EXTRA_PROTOCOL_COLUMNS

ROLE_COLUMNS = (
    ("catalyst", "metal_catalyst"),
    ("ligand", "ligand"),
    ("base", "base"),
    ("acid", "acid"),
    ("oxidant", "oxidant"),
    ("reductant", "reductant"),
    ("additive", "additive"),
    ("condensation_agent", "condensation_agent"),
    ("other_reagent", "other_reagent"),
    ("solvent", "solvent"),
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


def _append_role_value(bucket: Dict[str, List[str]], role: str, value: str) -> None:
    if not value:
        return
    if role not in bucket:
        bucket[role] = []
    bucket[role].append(value)


def _chemicals_to_role_map(chemicals: List[Dict[str, Any]]) -> Dict[str, List[str]]:
    bucket: Dict[str, List[str]] = {}
    for chem in chemicals:
        role = _clean_text(chem.get("role"))
        name = _clean_text(chem.get("abbreviation")) or _clean_text(chem.get("name"))
        if not role or not name:
            continue
        normalized = ROLE_ALIASES.get(role, role)
        if normalized in ("metal_catalyst", "co_catalyst", "catalyst"):
            _append_role_value(bucket, "catalyst", name)
        elif normalized == "ligand":
            _append_role_value(bucket, "ligand", name)
        elif normalized == "base":
            _append_role_value(bucket, "base", name)
        elif normalized == "acid":
            _append_role_value(bucket, "acid", name)
        elif normalized == "oxidant":
            _append_role_value(bucket, "oxidant", name)
        elif normalized == "reductant":
            _append_role_value(bucket, "reductant", name)
        elif normalized == "additive":
            _append_role_value(bucket, "additive", name)
        elif normalized == "condensation_agent":
            _append_role_value(bucket, "condensation_agent", name)
        elif normalized == "other_reagent":
            _append_role_value(bucket, "other_reagent", name)
        elif normalized == "solvent":
            _append_role_value(bucket, "solvent", name)
    return bucket


def _build_reaction_setup_from_row(row: Dict[str, Any]) -> List[Dict[str, Any]]:
    chemicals: List[Dict[str, Any]] = []
    for column, role in ROLE_COLUMNS:
        value = row.get(column)
        if not value:
            continue
        for item in _split_slash(value):
            chemicals.append({"name": item, "role": role})
    if not chemicals:
        return []
    return [{"chemicals": chemicals, "conditions": []}]


def row_to_protocol(row: Dict[str, Any]) -> Dict[str, Any]:
    reaction_smiles = _clean_text(row.get("reaction_smiles"))
    family = _clean_text(row.get("detected_reaction_type")) or _clean_text(row.get("reaction_family"))
    reaction_smarts = parse_semicolon_list(row.get("reaction_smarts"))
    tags = parse_semicolon_list(row.get("tags"))
    notes = _clean_text(row.get("notes"))
    protocol_id = _clean_text(row.get("protocol_id")) or _clean_text(row.get("reaction_id"))
    source = parse_source_field(row.get("source"))
    reference = _clean_text(row.get("reference"))
    if reference and not source.get("title"):
        source["title"] = reference

    reaction_setup = _parse_json_field(row.get("reaction_setup_json"))
    if not isinstance(reaction_setup, list):
        reaction_setup = _build_reaction_setup_from_row(row)

    original_procedure = _clean_text(row.get("original_procedure"))

    metadata: Dict[str, Any] = {
        "id": protocol_id,
        "source_reference": source,
    }
    if tags:
        metadata["tags"] = tags

    reaction: Dict[str, Any] = {
        "reaction_smiles": reaction_smiles,
        "family": family,
        "reaction_SMARTS": reaction_smarts,
        "notes": notes,
        "tags": tags,
    }

    protocol_data: Dict[str, Any] = {
        "schema_version": "2.0",
        "source_type": "protocol",
        "metadata": metadata,
        "reaction": reaction,
        "reaction_setup": reaction_setup,
        "original_procedure": original_procedure,
        "source": source,
    }

    return protocol_data


def protocol_to_row(protocol_data: Dict[str, Any]) -> Dict[str, str]:
    reaction = protocol_data.get("reaction", {}) or {}
    metadata = protocol_data.get("metadata", {}) or {}
    source = extract_source(protocol_data)

    reaction_smiles = _clean_text(reaction.get("reaction_smiles"))
    family = _clean_text(reaction.get("family"))
    tags = extract_tags(protocol_data)
    notes = extract_notes(protocol_data)
    protocol_id = _clean_text(metadata.get("id")) or _clean_text(metadata.get("protocol_id"))

    reaction_smarts_raw = reaction.get("reaction_SMARTS", [])
    if isinstance(reaction_smarts_raw, str):
        reaction_smarts = parse_semicolon_list(reaction_smarts_raw)
    else:
        reaction_smarts = [str(item).strip() for item in reaction_smarts_raw if str(item).strip()]

    reaction_setup = protocol_data.get("reaction_setup", [])
    reaction_setup_json = ""
    if reaction_setup:
        reaction_setup_json = json.dumps(reaction_setup, ensure_ascii=True, separators=(",", ":"))

    original_procedure = _clean_text(protocol_data.get("original_procedure"))

    chemicals = []
    if reaction_setup and isinstance(reaction_setup, list):
        first_setup = reaction_setup[0] if reaction_setup else {}
        chemicals = first_setup.get("chemicals", []) if isinstance(first_setup, dict) else []

    role_values = _chemicals_to_role_map(chemicals)

    row: Dict[str, str] = {col: "" for col in PROTOCOL_CSV_COLUMNS}
    row.update(
        {
            "reaction_id": protocol_id,
            "detected_reaction_type": family,
            "reaction_smiles": reaction_smiles,
            "reference": _clean_text(source.get("title")),
            "protocol_id": protocol_id,
            "reaction_smarts": serialize_semicolon_list(reaction_smarts),
            "tags": serialize_semicolon_list(tags),
            "notes": notes,
            "source": _format_source_field(source),
            "reaction_setup_json": reaction_setup_json,
            "original_procedure": original_procedure,
        }
    )

    for column, _role in ROLE_COLUMNS:
        values = role_values.get(column, [])
        if values:
            row[column] = "/".join(values)

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

