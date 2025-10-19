"""
Protocol Database Indexer

Builds and maintains an index of protocol JSON files for fast lookup.
The index includes reaction fingerprints for DRFP-based similarity search.

This module focuses on:
1. Scanning protocol_db directory
2. Extracting metadata (reaction, family, tags, source)
3. Computing DRFP fingerprints
4. Building searchable index
5. Incremental updates

Usage:
    from chemtools.protocol.indexer import ProtocolIndexer
    
    # Build index
    indexer = ProtocolIndexer()
    indexer.build_index()
    indexer.save()
    
    # Load existing index
    indexer = ProtocolIndexer.load()
"""

import json
import hashlib
from pathlib import Path
from typing import Dict, Any, List, Optional, Set
from dataclasses import dataclass, asdict, field
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


@dataclass
class ProtocolRecord:
    """
    Single protocol record in the index
    
    Stores essential metadata for fast lookup and similarity comparison
    """
    # File info
    filename: str
    filepath: str  # Relative to protocol_db
    file_hash: str  # MD5 hash to detect changes
    
    # Reaction info
    reaction_smiles: str
    reaction_family: str
    tags: List[str] = field(default_factory=list)
    reaction_smarts: List[str] = field(default_factory=list)  # SMARTS patterns for matching
    notes: str = ""
    
    # Source info
    source_title: str = ""
    source_journal: str = ""
    source_year: int = 0
    source_doi: str = ""
    source_url: str = ""
    
    # Normalized/computed fields
    canonical_reaction: str = ""
    reactants: List[str] = field(default_factory=list)
    products: List[str] = field(default_factory=list)
    
    # DRFP fingerprint (stored as list for JSON serialization)
    drfp_fingerprint: Optional[List[float]] = None
    
    # Metadata
    indexed_at: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ProtocolRecord':
        """Create from dictionary"""
        return cls(**data)


class ProtocolIndexer:
    """
    Build and maintain protocol database index
    
    The indexer:
    - Scans protocol_db/*.json files
    - Extracts reaction info and metadata
    - Computes DRFP fingerprints for similarity search
    - Builds lookup indexes (by family, tags)
    - Supports incremental updates (only reindex changed files)
    """
    
    def __init__(
        self,
        protocol_dir: Optional[Path] = None,
        index_path: Optional[Path] = None
    ):
        """
        Initialize indexer
        
        Args:
            protocol_dir: Directory containing protocol JSON files
            index_path: Path to save index file
        """
        # Set default paths
        if protocol_dir is None:
            project_root = Path(__file__).parent.parent.parent
            protocol_dir = project_root / 'data' / 'protocol_db'
        
        if index_path is None:
            index_path = protocol_dir / '.protocol_index.json'
        
        self.protocol_dir = Path(protocol_dir)
        self.index_path = Path(index_path)
        
        # Index data structures
        self.records: Dict[str, ProtocolRecord] = {}  # filename -> record
        self.family_index: Dict[str, List[str]] = {}  # family -> [filenames]
        self.tag_index: Dict[str, List[str]] = {}  # tag -> [filenames]
        
        # Index metadata
        self.metadata = {
            'version': '1.0',
            'created_at': None,
            'updated_at': None,
            'num_protocols': 0,
            'protocol_dir': str(self.protocol_dir),
            'has_drfp': False  # Whether DRFP fingerprints are computed
        }
    
    def build_index(self, compute_drfp: bool = True, force_rebuild: bool = False):
        """
        Build or update the protocol index
        
        Args:
            compute_drfp: Compute DRFP fingerprints (requires drfp package)
            force_rebuild: Rebuild all protocols (ignore file hashes)
        """
        logger.info(f"Building protocol index from {self.protocol_dir}")
        
        if not self.protocol_dir.exists():
            raise FileNotFoundError(f"Protocol directory not found: {self.protocol_dir}")
        
        # Load existing index if available (for incremental updates)
        existing_records = {}
        if not force_rebuild and self.index_path.exists():
            try:
                existing_indexer = self.load(self.index_path)
                existing_records = existing_indexer.records
                logger.info(f"Loaded existing index with {len(existing_records)} protocols")
            except Exception as e:
                logger.warning(f"Could not load existing index: {e}")
        
        # Find all protocol JSON files
        protocol_files = self._find_protocol_files()
        logger.info(f"Found {len(protocol_files)} protocol files")
        
        # Process each file
        processed = 0
        skipped = 0
        errors = 0
        
        for json_path in protocol_files:
            try:
                filename = json_path.name
                file_hash = self._compute_file_hash(json_path)
                
                # Check if file has changed
                if filename in existing_records:
                    existing_record = existing_records[filename]
                    if existing_record.file_hash == file_hash and not force_rebuild:
                        # File unchanged, reuse existing record
                        self.records[filename] = existing_record
                        skipped += 1
                        continue
                
                # Process the file
                record = self._process_protocol_file(json_path, file_hash, compute_drfp)
                self.records[filename] = record
                processed += 1
                
            except Exception as e:
                logger.error(f"Error processing {json_path.name}: {e}")
                errors += 1
        
        # Build lookup indexes
        self._build_lookup_indexes()
        
        # Update metadata
        if self.metadata['created_at'] is None:
            self.metadata['created_at'] = datetime.utcnow().isoformat()
        self.metadata['updated_at'] = datetime.utcnow().isoformat()
        self.metadata['num_protocols'] = len(self.records)
        self.metadata['has_drfp'] = compute_drfp
        
        logger.info(f"Index build complete:")
        logger.info(f"  Processed: {processed}")
        logger.info(f"  Skipped (unchanged): {skipped}")
        logger.info(f"  Errors: {errors}")
        logger.info(f"  Total protocols: {len(self.records)}")
    
    def _find_protocol_files(self) -> List[Path]:
        """Find all protocol JSON files (excluding examples and hidden files)"""
        json_files = []
        
        for json_path in self.protocol_dir.glob('*.json'):
            # Skip hidden files and schema files
            if json_path.name.startswith('.'):
                continue
            if 'schema' in json_path.name.lower():
                continue
            if json_path.is_file():
                json_files.append(json_path)
        
        return sorted(json_files)
    
    def _compute_file_hash(self, file_path: Path) -> str:
        """Compute MD5 hash of file for change detection"""
        md5 = hashlib.md5()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                md5.update(chunk)
        return md5.hexdigest()
    
    def _process_protocol_file(
        self,
        json_path: Path,
        file_hash: str,
        compute_drfp: bool
    ) -> ProtocolRecord:
        """
        Process a single protocol file and create index record
        
        Args:
            json_path: Path to protocol JSON file
            file_hash: MD5 hash of the file
            compute_drfp: Whether to compute DRFP fingerprint
        
        Returns:
            ProtocolRecord with extracted metadata
        """
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Extract reaction info
        reaction = data.get('reaction', {})
        reaction_smiles = reaction.get('reaction_smiles', '')
        reaction_family = reaction.get('family', '')
        notes = reaction.get('notes', '')
        
        # Extract reaction SMARTS patterns
        reaction_smarts_raw = reaction.get('reaction_SMARTS', [])
        if isinstance(reaction_smarts_raw, list):
            reaction_smarts = [str(s).strip() for s in reaction_smarts_raw if s]
        else:
            reaction_smarts = []
        
        # Extract tags
        tags_raw = reaction.get('tags', '')
        if isinstance(tags_raw, str):
            tags = [t.strip() for t in tags_raw.split(';') if t.strip()]
        elif isinstance(tags_raw, list):
            tags = [str(t).strip() for t in tags_raw if t]
        else:
            tags = []
        
        # Extract source info
        source = data.get('source', {})
        
        # Parse reaction SMILES
        canonical_rxn, reactants, products = self._parse_reaction(reaction_smiles)
        
        # Compute DRFP fingerprint if requested
        drfp_fp = None
        if compute_drfp and reaction_smiles:
            try:
                drfp_fp = self._compute_drfp(reaction_smiles)
            except Exception as e:
                logger.debug(f"Could not compute DRFP for {json_path.name}: {e}")
        
        # Create record
        record = ProtocolRecord(
            filename=json_path.name,
            filepath=str(json_path.relative_to(self.protocol_dir.parent)),
            file_hash=file_hash,
            reaction_smiles=reaction_smiles,
            reaction_smarts=reaction_smarts,
            reaction_family=reaction_family,
            tags=tags,
            notes=notes,
            source_title=source.get('title', ''),
            source_journal=source.get('journal', ''),
            source_year=source.get('year', 0),
            source_doi=source.get('doi', ''),
            source_url=source.get('url', ''),
            canonical_reaction=canonical_rxn,
            reactants=reactants,
            products=products,
            drfp_fingerprint=drfp_fp,
            indexed_at=datetime.utcnow().isoformat()
        )
        
        return record
    
    def _parse_reaction(self, reaction_smiles: str) -> tuple[str, List[str], List[str]]:
        """
        Parse reaction SMILES into canonical form and components
        
        Returns:
            (canonical_reaction, reactants, products)
        """
        if not reaction_smiles or '>>' not in reaction_smiles:
            return reaction_smiles, [], []
        
        try:
            parts = reaction_smiles.split('>>')
            if len(parts) != 2:
                return reaction_smiles, [], []
            
            reactants_str, products_str = parts
            
            # Split components
            reactants = [r.strip() for r in reactants_str.split('.') if r.strip()]
            products = [p.strip() for p in products_str.split('.') if p.strip()]
            
            # Try to canonicalize (optional, depends on RDKit availability)
            try:
                from ..smiles import normalize
                
                canonical_reactants = []
                for r in reactants:
                    norm = normalize(r)
                    canonical_reactants.append(norm.get('canonical', r))
                
                canonical_products = []
                for p in products:
                    norm = normalize(p)
                    canonical_products.append(norm.get('canonical', p))
                
                # Sort for consistency
                canonical_reactants.sort()
                canonical_products.sort()
                
                canonical = '.'.join(canonical_reactants) + '>>' + '.'.join(canonical_products)
                
                return canonical, canonical_reactants, canonical_products
                
            except Exception:
                # If canonicalization fails, return original
                return reaction_smiles, reactants, products
        
        except Exception as e:
            logger.debug(f"Error parsing reaction: {e}")
            return reaction_smiles, [], []
    
    def _compute_drfp(self, reaction_smiles: str) -> Optional[List[float]]:
        """
        Compute DRFP fingerprint for a reaction
        
        Returns:
            List of float values (DRFP vector) or None if computation fails
        """
        try:
            from drfp import DrfpEncoder
            import numpy as np
            
            # Use the correct DrfpEncoder initialization
            encoder = DrfpEncoder()
            fp = encoder.encode([reaction_smiles])[0]
            
            # Convert to list for JSON serialization
            if isinstance(fp, np.ndarray):
                return fp.tolist()
            return list(fp)
            
        except ImportError:
            logger.warning("DRFP package not available. Install with: pip install drfp")
            return None
        except Exception as e:
            logger.debug(f"DRFP computation failed: {e}")
            return None
    
    def _build_lookup_indexes(self):
        """Build family and tag lookup indexes"""
        self.family_index.clear()
        self.tag_index.clear()
        
        for filename, record in self.records.items():
            # Family index
            family = record.reaction_family
            if family:
                if family not in self.family_index:
                    self.family_index[family] = []
                self.family_index[family].append(filename)
            
            # Tag index
            for tag in record.tags:
                tag_lower = tag.lower().strip()
                if tag_lower:
                    if tag_lower not in self.tag_index:
                        self.tag_index[tag_lower] = []
                    self.tag_index[tag_lower].append(filename)
    
    def save(self, output_path: Optional[Path] = None):
        """
        Save index to JSON file
        
        Args:
            output_path: Path to save index (default: self.index_path)
        """
        save_path = output_path or self.index_path
        
        # Prepare data for serialization
        index_data = {
            'metadata': self.metadata,
            'records': {
                filename: record.to_dict()
                for filename, record in self.records.items()
            },
            'family_index': self.family_index,
            'tag_index': self.tag_index
        }
        
        # Save to file
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            json.dump(index_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Saved protocol index to {save_path}")
        logger.info(f"  {len(self.records)} protocols")
        logger.info(f"  {len(self.family_index)} families")
        logger.info(f"  {len(self.tag_index)} tags")
    
    @classmethod
    def load(cls, index_path: Optional[Path] = None) -> 'ProtocolIndexer':
        """
        Load index from JSON file
        
        Args:
            index_path: Path to index file
        
        Returns:
            ProtocolIndexer instance with loaded data
        """
        if index_path is None:
            project_root = Path(__file__).parent.parent.parent
            index_path = project_root / 'data' / 'protocol_db' / '.protocol_index.json'
        
        index_path = Path(index_path)
        
        if not index_path.exists():
            raise FileNotFoundError(f"Index file not found: {index_path}")
        
        with open(index_path, 'r', encoding='utf-8') as f:
            index_data = json.load(f)
        
        # Determine protocol_dir from metadata or index path
        protocol_dir_str = index_data.get('metadata', {}).get('protocol_dir')
        if protocol_dir_str:
            protocol_dir = Path(protocol_dir_str)
        else:
            protocol_dir = index_path.parent
        
        # Create indexer instance
        indexer = cls(protocol_dir=protocol_dir, index_path=index_path)
        
        # Load metadata
        indexer.metadata = index_data.get('metadata', {})
        
        # Load records
        records_data = index_data.get('records', {})
        for filename, record_dict in records_data.items():
            indexer.records[filename] = ProtocolRecord.from_dict(record_dict)
        
        # Load lookup indexes
        indexer.family_index = index_data.get('family_index', {})
        indexer.tag_index = index_data.get('tag_index', {})
        
        logger.info(f"Loaded protocol index from {index_path}")
        logger.info(f"  {len(indexer.records)} protocols")
        
        return indexer
    
    def get_stats(self) -> Dict[str, Any]:
        """Get index statistics"""
        return {
            'num_protocols': len(self.records),
            'num_families': len(self.family_index),
            'num_tags': len(self.tag_index),
            'has_drfp': self.metadata.get('has_drfp', False),
            'created_at': self.metadata.get('created_at'),
            'updated_at': self.metadata.get('updated_at'),
            'families': {
                family: len(files)
                for family, files in sorted(
                    self.family_index.items(),
                    key=lambda x: len(x[1]),
                    reverse=True
                )
            },
            'top_tags': {
                tag: len(files)
                for tag, files in sorted(
                    self.tag_index.items(),
                    key=lambda x: len(x[1]),
                    reverse=True
                )[:20]  # Top 20 tags
            }
        }
    
    def get_by_family(self, family: str) -> List[ProtocolRecord]:
        """Get all protocols for a reaction family"""
        filenames = self.family_index.get(family, [])
        return [self.records[f] for f in filenames if f in self.records]
    
    def get_by_tag(self, tag: str) -> List[ProtocolRecord]:
        """Get all protocols with a specific tag"""
        tag_lower = tag.lower().strip()
        filenames = self.tag_index.get(tag_lower, [])
        return [self.records[f] for f in filenames if f in self.records]
    
    def search_tags(self, tags: List[str]) -> List[ProtocolRecord]:
        """Get protocols matching any of the provided tags"""
        matching_files: Set[str] = set()
        for tag in tags:
            tag_lower = tag.lower().strip()
            matching_files.update(self.tag_index.get(tag_lower, []))
        return [self.records[f] for f in matching_files if f in self.records]


def build_index(
    protocol_dir: Optional[Path] = None,
    output_path: Optional[Path] = None,
    compute_drfp: bool = True,
    force_rebuild: bool = False
) -> ProtocolIndexer:
    """
    Standalone function to build protocol index
    
    Args:
        protocol_dir: Directory containing protocol JSON files
        output_path: Path to save index file
        compute_drfp: Compute DRFP fingerprints
        force_rebuild: Force rebuild all protocols
    
    Returns:
        ProtocolIndexer instance
    
    Example:
        from chemtools.protocol.indexer import build_index
        
        indexer = build_index(compute_drfp=True)
        print(f"Indexed {len(indexer.records)} protocols")
    """
    indexer = ProtocolIndexer(protocol_dir=protocol_dir, index_path=output_path)
    indexer.build_index(compute_drfp=compute_drfp, force_rebuild=force_rebuild)
    indexer.save()
    return indexer
