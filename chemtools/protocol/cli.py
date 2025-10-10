#!/usr/bin/env python3
"""
Protocol Indexer CLI

Command-line tool for building and managing the protocol database index.

Usage:
    # Build index with DRFP fingerprints
    python -m chemtools.protocol.cli build
    
    # Build without DRFP (faster)
    python -m chemtools.protocol.cli build --no-drfp
    
    # Force rebuild (ignore existing index)
    python -m chemtools.protocol.cli build --force
    
    # Show index statistics
    python -m chemtools.protocol.cli stats
    
    # List all families
    python -m chemtools.protocol.cli list-families
    
    # List all tags
    python -m chemtools.protocol.cli list-tags
    
    # Show protocols for a family
    python -m chemtools.protocol.cli show-family Suzuki
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Optional

from .indexer import ProtocolIndexer, build_index

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def cmd_build(args):
    """Build or rebuild the protocol index"""
    print("=" * 60)
    print("Protocol Database Indexer")
    print("=" * 60)
    print()
    
    protocol_dir = Path(args.protocol_dir) if args.protocol_dir else None
    output_path = Path(args.output) if args.output else None
    
    print(f"Building index...")
    if protocol_dir:
        print(f"  Protocol dir: {protocol_dir}")
    if output_path:
        print(f"  Output: {output_path}")
    print(f"  Compute DRFP: {args.drfp}")
    print(f"  Force rebuild: {args.force}")
    print()
    
    try:
        indexer = build_index(
            protocol_dir=protocol_dir,
            output_path=output_path,
            compute_drfp=args.drfp,
            force_rebuild=args.force
        )
        
        print()
        print("=" * 60)
        print("✅ Index built successfully!")
        print("=" * 60)
        print()
        print(f"Indexed {len(indexer.records)} protocols")
        print(f"Index saved to: {indexer.index_path}")
        print()
        
        # Show summary stats
        stats = indexer.get_stats()
        print("Summary:")
        print(f"  Families: {stats['num_families']}")
        print(f"  Tags: {stats['num_tags']}")
        print(f"  DRFP fingerprints: {'Yes' if stats['has_drfp'] else 'No'}")
        print()
        
        if stats['families']:
            print("Top families:")
            for family, count in list(stats['families'].items())[:10]:
                print(f"  {family}: {count}")
        
    except Exception as e:
        print()
        print(f"❌ Error building index: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_stats(args):
    """Show index statistics"""
    try:
        index_path = Path(args.index) if args.index else None
        indexer = ProtocolIndexer.load(index_path)
        
        stats = indexer.get_stats()
        
        print("=" * 60)
        print("Protocol Index Statistics")
        print("=" * 60)
        print()
        print(f"Total protocols: {stats['num_protocols']}")
        print(f"Families: {stats['num_families']}")
        print(f"Tags: {stats['num_tags']}")
        print(f"DRFP fingerprints: {'Yes' if stats['has_drfp'] else 'No'}")
        print()
        print(f"Created: {stats['created_at']}")
        print(f"Updated: {stats['updated_at']}")
        print()
        
        print("Protocols by family:")
        for family, count in stats['families'].items():
            print(f"  {family:50s} {count:3d}")
        
        print()
        print("Top tags:")
        for tag, count in stats['top_tags'].items():
            print(f"  {tag:30s} {count:3d}")
        
    except FileNotFoundError as e:
        print(f"❌ Index not found: {e}")
        print("Run 'build' command first to create the index")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


def cmd_list_families(args):
    """List all reaction families"""
    try:
        index_path = Path(args.index) if args.index else None
        indexer = ProtocolIndexer.load(index_path)
        
        print("=" * 60)
        print("Reaction Families")
        print("=" * 60)
        print()
        
        families = sorted(
            indexer.family_index.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )
        
        for family, filenames in families:
            print(f"{family:50s} ({len(filenames)} protocols)")
            if args.verbose:
                for fname in filenames:
                    record = indexer.records.get(fname)
                    if record:
                        print(f"  - {fname}")
                        print(f"    {record.source_title[:70]}")
                print()
        
        print()
        print(f"Total: {len(families)} families")
        
    except FileNotFoundError as e:
        print(f"❌ Index not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


def cmd_list_tags(args):
    """List all tags"""
    try:
        index_path = Path(args.index) if args.index else None
        indexer = ProtocolIndexer.load(index_path)
        
        print("=" * 60)
        print("Protocol Tags")
        print("=" * 60)
        print()
        
        tags = sorted(
            indexer.tag_index.items(),
            key=lambda x: len(x[1]),
            reverse=True
        )
        
        for tag, filenames in tags:
            print(f"{tag:40s} ({len(filenames)} protocols)")
        
        print()
        print(f"Total: {len(tags)} unique tags")
        
    except FileNotFoundError as e:
        print(f"❌ Index not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


def cmd_show_family(args):
    """Show protocols for a specific family"""
    try:
        index_path = Path(args.index) if args.index else None
        indexer = ProtocolIndexer.load(index_path)
        
        family = args.family
        protocols = indexer.get_by_family(family)
        
        if not protocols:
            print(f"❌ No protocols found for family: {family}")
            print()
            print("Available families:")
            for fam in sorted(indexer.family_index.keys()):
                print(f"  - {fam}")
            sys.exit(1)
        
        print("=" * 60)
        print(f"Protocols for: {family}")
        print("=" * 60)
        print()
        
        for i, protocol in enumerate(protocols, 1):
            print(f"{i}. {protocol.filename}")
            print(f"   Title: {protocol.source_title}")
            print(f"   Journal: {protocol.source_journal} ({protocol.source_year})")
            print(f"   DOI: {protocol.source_doi}")
            print(f"   Reaction: {protocol.reaction_smiles}")
            if protocol.tags:
                print(f"   Tags: {', '.join(protocol.tags)}")
            if protocol.notes:
                print(f"   Notes: {protocol.notes[:100]}...")
            print()
        
        print(f"Total: {len(protocols)} protocols")
        
    except FileNotFoundError as e:
        print(f"❌ Index not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


def cmd_show_tag(args):
    """Show protocols with a specific tag"""
    try:
        index_path = Path(args.index) if args.index else None
        indexer = ProtocolIndexer.load(index_path)
        
        tag = args.tag
        protocols = indexer.get_by_tag(tag)
        
        if not protocols:
            print(f"❌ No protocols found with tag: {tag}")
            print()
            print("Available tags (top 20):")
            tags = sorted(
                indexer.tag_index.items(),
                key=lambda x: len(x[1]),
                reverse=True
            )[:20]
            for t, _ in tags:
                print(f"  - {t}")
            sys.exit(1)
        
        print("=" * 60)
        print(f"Protocols with tag: {tag}")
        print("=" * 60)
        print()
        
        for i, protocol in enumerate(protocols, 1):
            print(f"{i}. {protocol.filename}")
            print(f"   {protocol.source_title}")
            print(f"   Family: {protocol.reaction_family}")
            print()
        
        print(f"Total: {len(protocols)} protocols")
        
    except FileNotFoundError as e:
        print(f"❌ Index not found: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(
        description='Protocol Database Indexer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Build index
  python -m chemtools.protocol.cli build
  
  # Build without DRFP (faster)
  python -m chemtools.protocol.cli build --no-drfp
  
  # Force rebuild
  python -m chemtools.protocol.cli build --force
  
  # Show statistics
  python -m chemtools.protocol.cli stats
  
  # List families
  python -m chemtools.protocol.cli list-families
  
  # Show Suzuki protocols
  python -m chemtools.protocol.cli show-family Suzuki
        """
    )
    
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Build command
    build_parser = subparsers.add_parser('build', help='Build or rebuild the index')
    build_parser.add_argument(
        '--protocol-dir',
        help='Protocol directory (default: data/protocol_db)'
    )
    build_parser.add_argument(
        '--output', '-o',
        help='Output index file path'
    )
    build_parser.add_argument(
        '--no-drfp',
        dest='drfp',
        action='store_false',
        default=True,
        help='Skip DRFP fingerprint computation'
    )
    build_parser.add_argument(
        '--force', '-f',
        action='store_true',
        help='Force rebuild (ignore existing index)'
    )
    build_parser.set_defaults(func=cmd_build)
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Show index statistics')
    stats_parser.add_argument('--index', help='Index file path')
    stats_parser.set_defaults(func=cmd_stats)
    
    # List families command
    list_fam_parser = subparsers.add_parser('list-families', help='List all families')
    list_fam_parser.add_argument('--index', help='Index file path')
    list_fam_parser.set_defaults(func=cmd_list_families)
    
    # List tags command
    list_tags_parser = subparsers.add_parser('list-tags', help='List all tags')
    list_tags_parser.add_argument('--index', help='Index file path')
    list_tags_parser.set_defaults(func=cmd_list_tags)
    
    # Show family command
    show_fam_parser = subparsers.add_parser('show-family', help='Show protocols for a family')
    show_fam_parser.add_argument('family', help='Reaction family name')
    show_fam_parser.add_argument('--index', help='Index file path')
    show_fam_parser.set_defaults(func=cmd_show_family)
    
    # Show tag command
    show_tag_parser = subparsers.add_parser('show-tag', help='Show protocols with a tag')
    show_tag_parser.add_argument('tag', help='Tag name')
    show_tag_parser.add_argument('--index', help='Index file path')
    show_tag_parser.set_defaults(func=cmd_show_tag)
    
    # Parse and execute
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    args.func(args)


if __name__ == '__main__':
    main()
