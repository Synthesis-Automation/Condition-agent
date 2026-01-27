"""
Convert protocol JSON files to CSV format

This CLI tool converts protocol JSON files (in protocol_db_v2 format) to CSV format.
Supports converting individual JSON files or entire directories.

Usage:
    # Convert single JSON file
    python -m chemtools.protocol.convert_json_to_csv input.json output.csv
    
    # Convert entire directory
    python -m chemtools.protocol.convert_json_to_csv data/protocol_db_v2/ data/protocols.csv
    
    # Convert directory to multiple family-specific CSVs
    python -m chemtools.protocol.convert_json_to_csv data/protocol_db_v2/ data/protocol_db_v2_csv/ --split-by-family
    
    # Combine multiple CSVs into one
    python -m chemtools.protocol.convert_json_to_csv --combine data/protocol_db_v2_csv/*.csv data/protocols_combined.csv
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Any, List

from .csv_utils import (
    protocol_to_row,
    write_protocol_csv,
    PROTOCOL_CSV_COLUMNS,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_json_protocol(json_path: Path) -> List[Dict[str, Any]]:
    """
    Load protocol(s) from JSON file.
    
    Handles both single protocol objects and arrays of protocols.
    
    Args:
        json_path: Path to JSON file
        
    Returns:
        List of protocol dictionaries
    """
    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        
        # Handle both single protocol and array of protocols
        if isinstance(data, list):
            return data
        elif isinstance(data, dict):
            return [data]
        else:
            logger.error(f"Invalid JSON structure in {json_path}")
            return []
    except Exception as e:
        logger.error(f"Error loading {json_path}: {e}")
        return []


def convert_json_file_to_csv(json_path: Path, csv_path: Path) -> int:
    """
    Convert a single JSON file to CSV.
    
    Args:
        json_path: Input JSON file path
        csv_path: Output CSV file path
        
    Returns:
        Number of protocols converted
    """
    protocols = load_json_protocol(json_path)
    if not protocols:
        return 0
    
    write_protocol_csv(protocols, str(csv_path), mode="w")
    logger.info(f"Converted {len(protocols)} protocol(s) from {json_path.name} to {csv_path}")
    return len(protocols)


def convert_directory_to_csv(
    input_dir: Path,
    output_path: Path,
    split_by_family: bool = False
) -> int:
    """
    Convert all JSON files in a directory to CSV.
    
    Args:
        input_dir: Directory containing JSON files
        output_path: Output CSV file path (or directory if split_by_family)
        split_by_family: If True, create separate CSV for each reaction family
        
    Returns:
        Total number of protocols converted
    """
    if not input_dir.exists():
        logger.error(f"Input directory not found: {input_dir}")
        return 0
    
    # Find all JSON files
    json_files = []
    for json_path in input_dir.glob("*.json"):
        # Skip hidden files and schema files
        if json_path.name.startswith(".") or "schema" in json_path.name.lower():
            continue
        json_files.append(json_path)
    
    if not json_files:
        logger.error(f"No JSON files found in {input_dir}")
        return 0
    
    logger.info(f"Found {len(json_files)} JSON files")
    
    # Load all protocols
    all_protocols = []
    for json_path in json_files:
        protocols = load_json_protocol(json_path)
        all_protocols.extend(protocols)
    
    if not all_protocols:
        logger.error("No valid protocols found")
        return 0
    
    logger.info(f"Loaded {len(all_protocols)} protocols")
    
    if split_by_family:
        # Group protocols by reaction family
        family_groups: Dict[str, List[Dict[str, Any]]] = {}
        for protocol in all_protocols:
            reaction = protocol.get("reaction", {}) or {}
            family = reaction.get("family", "Unknown").strip()
            if not family or family.lower() == "unknown":
                family = "Unknown"
            
            if family not in family_groups:
                family_groups[family] = []
            family_groups[family].append(protocol)
        
        # Create output directory
        output_dir = output_path
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Write each family to separate CSV
        total_count = 0
        for family, protocols in family_groups.items():
            # Sanitize family name for filename
            safe_family = family.replace("/", "_").replace("\\", "_").replace(" ", "_")
            csv_path = output_dir / f"{safe_family}.csv"
            write_protocol_csv(protocols, str(csv_path), mode="w")
            logger.info(f"  {family}: {len(protocols)} protocols -> {csv_path.name}")
            total_count += len(protocols)
        
        logger.info(f"Created {len(family_groups)} family-specific CSV files")
        return total_count
    else:
        # Write all protocols to single CSV
        write_protocol_csv(all_protocols, str(output_path), mode="w")
        logger.info(f"Converted {len(all_protocols)} protocols to {output_path}")
        return len(all_protocols)


def combine_csv_files(csv_paths: List[Path], output_path: Path) -> int:
    """
    Combine multiple CSV files into one.
    
    Args:
        csv_paths: List of input CSV file paths
        output_path: Output CSV file path
        
    Returns:
        Total number of protocols combined
    """
    from .csv_utils import read_protocol_csv
    
    all_protocols = []
    for csv_path in csv_paths:
        if not csv_path.exists():
            logger.warning(f"CSV file not found: {csv_path}")
            continue
        
        protocols = read_protocol_csv(str(csv_path))
        all_protocols.extend(protocols)
        logger.info(f"Loaded {len(protocols)} protocols from {csv_path.name}")
    
    if not all_protocols:
        logger.error("No protocols found in input CSV files")
        return 0
    
    write_protocol_csv(all_protocols, str(output_path), mode="w")
    logger.info(f"Combined {len(all_protocols)} protocols to {output_path}")
    return len(all_protocols)


def main():
    parser = argparse.ArgumentParser(
        description="Convert protocol JSON files to CSV format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "input",
        help="Input JSON file or directory"
    )
    
    parser.add_argument(
        "output",
        help="Output CSV file or directory (if --split-by-family)"
    )
    
    parser.add_argument(
        "--split-by-family",
        action="store_true",
        help="Create separate CSV file for each reaction family (output must be directory)"
    )
    
    parser.add_argument(
        "--combine",
        action="store_true",
        help="Combine multiple CSV files (input should be glob pattern)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    input_path = Path(args.input)
    output_path = Path(args.output)
    
    try:
        if args.combine:
            # Combine mode: input is glob pattern
            import glob
            csv_paths = [Path(p) for p in glob.glob(args.input)]
            if not csv_paths:
                logger.error(f"No CSV files found matching: {args.input}")
                return 1
            
            count = combine_csv_files(csv_paths, output_path)
            logger.info(f"✓ Successfully combined {count} protocols")
            return 0
        
        elif input_path.is_file():
            # Single file conversion
            if not input_path.suffix == ".json":
                logger.error("Input file must be .json")
                return 1
            
            count = convert_json_file_to_csv(input_path, output_path)
            if count > 0:
                logger.info(f"✓ Successfully converted {count} protocol(s)")
                return 0
            else:
                logger.error("✗ Conversion failed")
                return 1
        
        elif input_path.is_dir():
            # Directory conversion
            count = convert_directory_to_csv(
                input_path,
                output_path,
                split_by_family=args.split_by_family
            )
            if count > 0:
                logger.info(f"✓ Successfully converted {count} protocol(s)")
                return 0
            else:
                logger.error("✗ Conversion failed")
                return 1
        
        else:
            logger.error(f"Input not found: {input_path}")
            return 1
    
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
