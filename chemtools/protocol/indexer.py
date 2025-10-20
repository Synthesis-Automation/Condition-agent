"""
Protocol Database Indexer

Builds and maintains an index of protocol JSON files for fast lookup.
The index includes reaction fingerprints for DRFP-based similarity search.

This module focuses on:
1. Scanning protocol_db directory
2. Extracting metadata (reaction, family, tags, source)
3. Computing DRFP fingerprints (stored separately in NPZ format)
4. Building searchable index
5. Incremental updates

DRFP Storage Strategy:
- DRFP fingerprints are stored in a separate .npz file (compressed binary)
- This reduces index JSON size by ~90% (from 4096 floats per protocol)
- Uses uint8 arrays (0-255) for space efficiency
- Protocol index only stores metadata; fingerprints loaded on-demand

Usage:
    from chemtools.protocol.indexer import ProtocolIndexer
    
    # Build index with DRFP fingerprints
    indexer = ProtocolIndexer()
    indexer.build_index(compute_drfp=True)
    indexer.save()
    
    # Load existing index
    indexer = ProtocolIndexer.load()
    
    # Get DRFP fingerprint for a protocol
    fp = indexer.get_drfp_fingerprint("protocol_001.json")
"""

import json
import hashlib
from pathlib import Path
from typing import Dict, Any, List, Optional, Set
from dataclasses import dataclass, asdict, field
from datetime import datetime
import logging

logger = logging.getLogger(__name__)

# Try to import DRFP storage utilities
try:
    from ..util.drfp_storage import DRFPLoader, save_drfp_index
    _DRFP_STORAGE_AVAILABLE = True
except ImportError:
    DRFPLoader = None  # type: ignore
    save_drfp_index = None  # type: ignore
    _DRFP_STORAGE_AVAILABLE = False


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
    
    # DRFP fingerprint reference (NOT stored in JSON - use separate NPZ file)
    # This flag indicates whether DRFP was computed for this protocol
    has_drfp: bool = False
    
    # Metadata
    indexed_at: str = ""
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization (excludes DRFP data)"""
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
    - Computes DRFP fingerprints for similarity search (stored in separate NPZ)
    - Builds lookup indexes (by family, tags)
    - Supports incremental updates (only reindex changed files)
    
    DRFP Storage:
    - Fingerprints stored in separate .protocol_drfp.npz file (compressed binary)
    - Reduces index size by ~90% compared to embedding in JSON
    - Loaded on-demand via DRFPLoader for similarity searches
    """
    
    def __init__(
        self,
        protocol_dir: Optional[Path] = None,
        index_path: Optional[Path] = None,
        drfp_path: Optional[Path] = None
    ):
        """
        Initialize indexer
        
        Args:
            protocol_dir: Directory containing protocol JSON files
            index_path: Path to save index file
            drfp_path: Path to save DRFP fingerprints NPZ file
        """
        # Set default paths
        if protocol_dir is None:
            project_root = Path(__file__).parent.parent.parent
            protocol_dir = project_root / 'data' / 'protocol_db'
        
        if index_path is None:
            index_path = protocol_dir / '.protocol_index.json'
        
        if drfp_path is None:
            drfp_path = protocol_dir / '.protocol_drfp.npz'
        
        self.protocol_dir = Path(protocol_dir)
        self.index_path = Path(index_path)
        self.drfp_path = Path(drfp_path)
        
        # Index data structures
        self.records: Dict[str, ProtocolRecord] = {}  # filename -> record
        self.family_index: Dict[str, List[str]] = {}  # family -> [filenames]
        self.tag_index: Dict[str, List[str]] = {}  # tag -> [filenames]
        
        # DRFP loader (lazy-loaded)
        self._drfp_loader: Optional[Any] = None
        
        # Index metadata
        self.metadata = {
            'version': '1.0',
            'created_at': None,
            'updated_at': None,
            'num_protocols': 0,
            'protocol_dir': str(self.protocol_dir),
            'has_drfp': False,  # Whether DRFP fingerprints are computed
            'drfp_storage': 'npz'  # Storage format for DRFP fingerprints
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
        
        # Process each file and collect DRFP data
        processed = 0
        skipped = 0
        errors = 0
        
        # For DRFP storage
        drfp_fingerprints = []
        drfp_filenames = []
        
        for json_path in protocol_files:
            try:
                filename = json_path.name
                file_hash = self._compute_file_hash(json_path)
                
                # Check if file has changed (check against any record from this file)
                file_unchanged = False
                if not force_rebuild:
                    # Check if any existing record from this file exists
                    for existing_key, existing_record in existing_records.items():
                        # Match by filepath rather than filename
                        if existing_record.filepath.endswith(filename):
                            if existing_record.file_hash == file_hash:
                                file_unchanged = True
                                break
                
                if file_unchanged:
                    # File unchanged, reuse all existing records from this file
                    for existing_key, existing_record in existing_records.items():
                        if existing_record.filepath.endswith(filename):
                            self.records[existing_key] = existing_record
                    skipped += 1
                    continue
                
                # Process the file (returns lists now)
                records_list, drfp_fps_list = self._process_protocol_file(json_path, file_hash, compute_drfp)
                
                # Add all records from this file
                for record, drfp_fp in zip(records_list, drfp_fps_list):
                    self.records[record.filename] = record
                    processed += 1
                    
                    # Collect DRFP fingerprint if computed
                    if drfp_fp is not None:
                        drfp_fingerprints.append(drfp_fp)
                        drfp_filenames.append(record.filename)
                
            except Exception as e:
                logger.error(f"Error processing {json_path.name}: {e}")
                errors += 1
        
        # Save DRFP fingerprints to separate NPZ file
        if compute_drfp and drfp_fingerprints and _DRFP_STORAGE_AVAILABLE:
            try:
                logger.info(f"Saving {len(drfp_fingerprints)} DRFP fingerprints to {self.drfp_path}")
                save_drfp_index(
                    fingerprints=drfp_fingerprints,
                    reaction_ids=drfp_filenames,
                    output_path=str(self.drfp_path),
                    n_bits=4096,
                    radius=3
                )
                self.metadata['has_drfp'] = True
            except Exception as e:
                logger.error(f"Failed to save DRFP fingerprints: {e}")
                self.metadata['has_drfp'] = False
        elif compute_drfp and not _DRFP_STORAGE_AVAILABLE:
            logger.warning("DRFP storage utilities not available. Fingerprints not saved.")
            self.metadata['has_drfp'] = False
        
        # Build lookup indexes
        self._build_lookup_indexes()
        
        # Update metadata
        if self.metadata['created_at'] is None:
            self.metadata['created_at'] = datetime.utcnow().isoformat()
        self.metadata['updated_at'] = datetime.utcnow().isoformat()
        self.metadata['num_protocols'] = len(self.records)
        
        logger.info(f"Index build complete:")
        logger.info(f"  Processed: {processed}")
        logger.info(f"  Skipped (unchanged): {skipped}")
        logger.info(f"  Errors: {errors}")
        logger.info(f"  Total protocols: {len(self.records)}")
        logger.info(f"  DRFP fingerprints: {len(drfp_fingerprints) if drfp_fingerprints else 0}")
    
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
    ) -> tuple[List[ProtocolRecord], List[Optional[Any]]]:
        """
        Process a single protocol file and create index records
        
        Each JSON file can now contain either:
        - A single protocol object (legacy format)
        - An array of protocol objects (new format)
        
        Args:
            json_path: Path to protocol JSON file
            file_hash: MD5 hash of the file
            compute_drfp: Whether to compute DRFP fingerprint
        
        Returns:
            Tuple of (List[ProtocolRecord], List[drfp_fingerprint])
            drfp_fingerprint is None if not computed or computation failed
        """
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        # Handle both array and single protocol formats
        if isinstance(data, list):
            protocols = data
        else:
            protocols = [data]
        
        records = []
        drfp_fps = []
        
        for idx, protocol_data in enumerate(protocols):
            # Generate unique filename for each protocol in array
            if len(protocols) > 1:
                # Multiple protocols: use index suffix
                protocol_id = f"{json_path.stem}[{idx}]"
            else:
                # Single protocol: use original filename
                protocol_id = json_path.stem
            
            record, drfp_fp = self._process_single_protocol(
                protocol_data=protocol_data,
                json_path=json_path,
                protocol_id=protocol_id,
                file_hash=file_hash,
                compute_drfp=compute_drfp
            )
            records.append(record)
            drfp_fps.append(drfp_fp)
        
        return records, drfp_fps
    
    def _process_single_protocol(
        self,
        protocol_data: Dict[str, Any],
        json_path: Path,
        protocol_id: str,
        file_hash: str,
        compute_drfp: bool
    ) -> tuple[ProtocolRecord, Optional[Any]]:
        """
        Process a single protocol object and create index record
        
        Args:
            protocol_data: Protocol JSON object
            json_path: Path to source JSON file
            protocol_id: Unique identifier for this protocol
            file_hash: MD5 hash of the file
            compute_drfp: Whether to compute DRFP fingerprint
        
        Returns:
            Tuple of (ProtocolRecord, drfp_fingerprint)
            drfp_fingerprint is None if not computed or computation failed
        """
        # Extract reaction info
        reaction = protocol_data.get('reaction', {})
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
        source = protocol_data.get('source', {})
        
        # Parse reaction SMILES
        canonical_rxn, reactants, products = self._parse_reaction(reaction_smiles)
        
        # Compute DRFP fingerprint if requested
        drfp_fp = None
        has_drfp = False
        if compute_drfp and reaction_smiles:
            try:
                drfp_fp = self._compute_drfp(reaction_smiles)
                if drfp_fp is not None:
                    has_drfp = True
            except Exception as e:
                logger.debug(f"Could not compute DRFP for {json_path.name}: {e}")
        
        # Create record (without embedding DRFP)
        record = ProtocolRecord(
            filename=protocol_id,  # Use protocol_id instead of json_path.name
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
            has_drfp=has_drfp,
            indexed_at=datetime.utcnow().isoformat()
        )
        
        return record, drfp_fp
    
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
    
    def _compute_drfp(self, reaction_smiles: str) -> Optional[Any]:
        """
        Compute DRFP fingerprint for a reaction
        
        Returns:
            NumPy array (uint8) or None if computation fails
        """
        try:
            from drfp import DrfpEncoder
            import numpy as np
            
            # Use the correct DrfpEncoder initialization
            encoder = DrfpEncoder()
            fp = encoder.encode([reaction_smiles])[0]
            
            # Convert to uint8 numpy array for efficient storage
            if isinstance(fp, np.ndarray):
                return fp.astype('uint8')
            return np.array(fp, dtype='uint8')
            
        except ImportError:
            logger.warning("DRFP package not available. Install with: pip install drfp")
            return None
        except Exception as e:
            logger.debug(f"DRFP computation failed: {e}")
            return None
    
    def get_drfp_fingerprint(self, filename: str) -> Optional[Any]:
        """
        Get DRFP fingerprint for a protocol from the NPZ file
        
        Args:
            filename: Protocol filename (e.g., 'Suzuki_Cu_C(sp3)-C(sp2).json')
        
        Returns:
            NumPy array of DRFP fingerprint or None if not available
        """
        if not self.metadata.get('has_drfp', False):
            return None
        
        if not _DRFP_STORAGE_AVAILABLE or DRFPLoader is None:
            logger.warning("DRFP storage utilities not available")
            return None
        
        # Lazy-load DRFP loader
        if self._drfp_loader is None:
            if not self.drfp_path.exists():
                logger.warning(f"DRFP file not found: {self.drfp_path}")
                return None
            
            try:
                self._drfp_loader = DRFPLoader(str(self.drfp_path))
            except Exception as e:
                logger.error(f"Failed to load DRFP file: {e}")
                return None
        
        # Get fingerprint for this protocol
        try:
            return self._drfp_loader.get_fingerprint(filename)
        except Exception as e:
            logger.debug(f"Could not get DRFP for {filename}: {e}")
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
        
        # Determine DRFP path
        drfp_path = protocol_dir / '.protocol_drfp.npz'
        
        # Create indexer instance
        indexer = cls(protocol_dir=protocol_dir, index_path=index_path, drfp_path=drfp_path)
        
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
        if indexer.metadata.get('has_drfp', False):
            logger.info(f"  DRFP fingerprints available in {indexer.drfp_path}")
        
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
