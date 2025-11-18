"""
RDF File Reorganization Script

This script reorganizes RDF files from a complex nested structure into a flat structure:
- Target structure: /reaction_type/*.rdf
- Handles year subfolders (2020, 2024, etc.)
- Renames duplicate files to avoid conflicts (adds year/folder prefix)
- Creates a backup mapping file for tracking

Usage:
    python reorganize_rdf_files.py [--dry-run] [--backup]
"""

import os
import shutil
import json
import argparse
from pathlib import Path
from collections import defaultdict
from datetime import datetime


def sanitize_folder_name(name: str) -> str:
    """Convert folder name to safe filename prefix."""
    # Remove special characters and spaces
    safe_name = name.replace(' ', '_').replace('–', '-').replace('—', '-')
    safe_name = ''.join(c for c in safe_name if c.isalnum() or c in ['_', '-'])
    return safe_name


def find_all_rdf_files(base_path: Path) -> dict:
    """
    Scan all RDF files in the directory structure.
    
    Returns:
        dict: {reaction_type: [(source_path, relative_path, filename), ...]}
    """
    rdf_files = defaultdict(list)
    
    for reaction_folder in base_path.iterdir():
        if not reaction_folder.is_dir():
            continue
        
        if reaction_folder.name == 'reorganized':  # Skip if already processed
            continue
            
        reaction_type = reaction_folder.name
        
        # Find all RDF files recursively
        for rdf_file in reaction_folder.rglob('*.rdf'):
            relative_path = rdf_file.relative_to(reaction_folder)
            rdf_files[reaction_type].append((rdf_file, relative_path, rdf_file.name))
    
    return rdf_files


def generate_unique_filename(filename: str, parent_folder: str, year_folder: str, existing_files: set) -> str:
    """
    Generate a unique filename by adding prefixes if needed.
    
    Args:
        filename: Original filename (e.g., "1000.rdf")
        parent_folder: Parent subfolder name (e.g., "Chan-Lam")
        year_folder: Year folder if applicable (e.g., "2020")
        existing_files: Set of already used filenames
    
    Returns:
        Unique filename
    """
    # Remove .rdf extension
    base_name = filename.rsplit('.rdf', 1)[0]
    
    # Strategy 1: Try original filename
    if filename not in existing_files:
        return filename
    
    # Strategy 2: Add year prefix if available
    if year_folder:
        new_name = f"{year_folder}_{base_name}.rdf"
        if new_name not in existing_files:
            return new_name
    
    # Strategy 3: Add parent folder prefix
    if parent_folder:
        safe_parent = sanitize_folder_name(parent_folder)
        new_name = f"{safe_parent}_{base_name}.rdf"
        if new_name not in existing_files:
            return new_name
    
    # Strategy 4: Add both year and parent folder
    if year_folder and parent_folder:
        safe_parent = sanitize_folder_name(parent_folder)
        new_name = f"{year_folder}_{safe_parent}_{base_name}.rdf"
        if new_name not in existing_files:
            return new_name
    
    # Strategy 5: Add sequential number
    counter = 1
    while True:
        new_name = f"{base_name}_{counter}.rdf"
        if new_name not in existing_files:
            return new_name
        counter += 1


def extract_folder_info(relative_path: Path) -> tuple:
    """
    Extract year folder and parent subfolder from relative path.
    
    Returns:
        (year_folder, parent_subfolder)
    """
    parts = relative_path.parts[:-1]  # Exclude filename
    
    year_folder = None
    parent_subfolder = None
    
    for part in parts:
        # Check if it's a year folder (e.g., "2020", "2020-2024")
        if part[0].isdigit() and any(c.isdigit() for c in part[:4]):
            year_folder = part
        else:
            parent_subfolder = part
    
    return year_folder, parent_subfolder


def reorganize_files(base_path: Path, output_base: Path, dry_run: bool = False, create_backup: bool = True) -> dict:
    """
    Reorganize RDF files into flat structure by reaction type.
    
    Args:
        base_path: Source directory (original_dataset)
        output_base: Target directory (reorganized)
        dry_run: If True, only print actions without executing
        create_backup: If True, create mapping JSON file
    
    Returns:
        dict: Mapping of original to new file paths
    """
    print(f"📂 Scanning {base_path}...")
    rdf_files = find_all_rdf_files(base_path)
    
    total_files = sum(len(files) for files in rdf_files.values())
    print(f"✓ Found {total_files} RDF files across {len(rdf_files)} reaction types\n")
    
    mapping = {}
    stats = defaultdict(int)
    
    for reaction_type, files in sorted(rdf_files.items()):
        print(f"Processing: {reaction_type} ({len(files)} files)")
        
        # Create output directory
        output_dir = output_base / reaction_type
        if not dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)
        
        existing_files = set()
        renamed_count = 0
        
        for source_path, relative_path, original_filename in files:
            # Extract folder information
            year_folder, parent_subfolder = extract_folder_info(relative_path)
            
            # Generate unique filename
            new_filename = generate_unique_filename(
                original_filename,
                parent_subfolder,
                year_folder,
                existing_files
            )
            
            existing_files.add(new_filename)
            
            # Track if renamed
            if new_filename != original_filename:
                renamed_count += 1
            
            # Prepare destination path
            dest_path = output_dir / new_filename
            
            # Store mapping
            mapping[str(source_path)] = {
                'new_path': str(dest_path),
                'original_name': original_filename,
                'new_name': new_filename,
                'year_folder': year_folder,
                'parent_folder': parent_subfolder,
                'relative_path': str(relative_path)
            }
            
            # Copy or show action
            if dry_run:
                action = "RENAME" if new_filename != original_filename else "COPY"
                print(f"  [{action}] {relative_path} → {new_filename}")
            else:
                shutil.copy2(source_path, dest_path)
        
        stats[reaction_type] = {
            'total': len(files),
            'renamed': renamed_count
        }
        
        print(f"  ✓ Processed {len(files)} files ({renamed_count} renamed)\n")
    
    # Save mapping file
    if create_backup and not dry_run:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        mapping_file = output_base.parent / f"reorganization_mapping_{timestamp}.json"
        
        with open(mapping_file, 'w', encoding='utf-8') as f:
            json.dump({
                'timestamp': timestamp,
                'source': str(base_path),
                'destination': str(output_base),
                'statistics': stats,
                'mappings': mapping
            }, f, indent=2, ensure_ascii=False)
        
        print(f"💾 Mapping saved to: {mapping_file}")
    
    return mapping


def print_summary(stats: dict):
    """Print reorganization summary."""
    print("\n" + "="*70)
    print("REORGANIZATION SUMMARY")
    print("="*70)
    
    total_files = sum(s['total'] for s in stats.values())
    total_renamed = sum(s['renamed'] for s in stats.values())
    
    print(f"Total reaction types: {len(stats)}")
    print(f"Total files processed: {total_files}")
    print(f"Files renamed: {total_renamed}")
    print(f"Files kept original name: {total_files - total_renamed}")
    
    print("\nTop reaction types by file count:")
    sorted_stats = sorted(stats.items(), key=lambda x: x[1]['total'], reverse=True)
    for reaction_type, counts in sorted_stats[:10]:
        print(f"  {reaction_type}: {counts['total']} files ({counts['renamed']} renamed)")


def main():
    parser = argparse.ArgumentParser(description='Reorganize RDF files by reaction type')
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually copying files'
    )
    parser.add_argument(
        '--no-backup',
        action='store_true',
        help='Do not create mapping backup file'
    )
    parser.add_argument(
        '--source',
        type=str,
        default='original_dataset',
        help='Source directory (default: original_dataset)'
    )
    parser.add_argument(
        '--dest',
        type=str,
        default='reorganized',
        help='Destination directory (default: reorganized)'
    )
    
    args = parser.parse_args()
    
    # Get paths
    script_dir = Path(__file__).parent
    source_path = script_dir / args.source
    dest_path = script_dir / args.dest
    
    if not source_path.exists():
        print(f"❌ Error: Source directory not found: {source_path}")
        return
    
    print("="*70)
    print("RDF FILE REORGANIZATION")
    print("="*70)
    print(f"Source: {source_path}")
    print(f"Destination: {dest_path}")
    print(f"Mode: {'DRY RUN (no changes)' if args.dry_run else 'LIVE (files will be copied)'}")
    print(f"Backup mapping: {'No' if args.no_backup else 'Yes'}")
    print("="*70 + "\n")
    
    if not args.dry_run:
        response = input("⚠️  This will copy files. Continue? (yes/no): ")
        if response.lower() not in ['yes', 'y']:
            print("Aborted.")
            return
    
    # Perform reorganization
    mapping = reorganize_files(
        source_path,
        dest_path,
        dry_run=args.dry_run,
        create_backup=not args.no_backup
    )
    
    # Extract stats from mapping
    stats = defaultdict(lambda: {'total': 0, 'renamed': 0})
    for source, info in mapping.items():
        reaction_type = Path(info['new_path']).parent.name
        stats[reaction_type]['total'] += 1
        if info['original_name'] != info['new_name']:
            stats[reaction_type]['renamed'] += 1
    
    print_summary(dict(stats))
    
    if args.dry_run:
        print("\n💡 This was a dry run. Use without --dry-run to actually copy files.")
    else:
        print(f"\n✅ Reorganization complete! Files copied to: {dest_path}")


if __name__ == "__main__":
    main()
