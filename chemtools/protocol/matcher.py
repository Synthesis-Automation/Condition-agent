"""
Protocol Database Matcher (Legacy)

DEPRECATED: This module is maintained for backward compatibility only.
New code should use chemtools.protocol.recommend.ProtocolRecommender instead,
which provides better SMARTS-based filtering and DRFP similarity search.

This module provides functionality to match user-supplied reactions to standard
procedure protocols stored in the protocol_db directory.

Key features:
- DRFP-based similarity search for protocol matching
- Index-based fast lookup as protocol collection grows
- Tag-based filtering and search
- Reaction family matching
- Extraction of detailed experimental procedures

Usage:
    from chemtools.protocol import ProtocolMatcher
    
    matcher = ProtocolMatcher()
    
    # Find best matching protocol
    result = matcher.find_best_match(
        reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        reaction_family='Suzuki',  # optional
        k=5  # return top 5 matches
    )
    
    # Extract full protocol details
    protocol = matcher.get_protocol(result['matches'][0]['filename'])
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, asdict, field
from datetime import datetime
import logging

from ..smiles import normalize
from ..reaction_similarity import compute_drfp_similarity

logger = logging.getLogger(__name__)


@dataclass
class ProtocolMetadata:
    """Metadata for a single protocol file"""
    filename: str
    reaction_smiles: str
    reaction_family: str
    tags: List[str]
    notes: str
    source_title: str
    source_journal: str
    source_year: int
    source_doi: str
    
    # Computed fields
    reaction_smarts: List[str] = field(default_factory=list)  # SMARTS patterns for matching
    canonical_reaction: str = ""
    drfp_fingerprint: Optional[Any] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization"""
        data = asdict(self)
        # Don't serialize DRFP fingerprint in JSON (too large)
        data.pop('drfp_fingerprint', None)
        return data


@dataclass
class ProtocolMatch:
    """A matched protocol with similarity score"""
    filename: str
    similarity: float
    reaction_smiles: str
    reaction_family: str
    tags: List[str]
    notes: str
    source_title: str
    source_doi: str
    reaction_smarts: List[str] = field(default_factory=list)  # SMARTS patterns for matching
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class ProtocolIndex:
    """
    Index for fast protocol lookup
    
    The index stores:
    - Protocol metadata (reaction SMILES, family, tags, source info)
    - Reaction families -> list of protocol filenames
    - Tags -> list of protocol filenames
    - Precomputed DRFP fingerprints for similarity search
    """
    
    def __init__(self, index_path: Optional[Path] = None):
        self.index_path = index_path
        self.protocols: Dict[str, ProtocolMetadata] = {}
        self.family_index: Dict[str, List[str]] = {}
        self.tag_index: Dict[str, List[str]] = {}
        self.index_metadata: Dict[str, Any] = {
            'created': None,
            'updated': None,
            'num_protocols': 0,
            'version': '1.0'
        }
    
    def add_protocol(self, metadata: ProtocolMetadata):
        """Add a protocol to the index"""
        filename = metadata.filename
        self.protocols[filename] = metadata
        
        # Update family index
        family = metadata.reaction_family
        if family not in self.family_index:
            self.family_index[family] = []
        if filename not in self.family_index[family]:
            self.family_index[family].append(filename)
        
        # Update tag index
        for tag in metadata.tags:
            tag_lower = tag.lower().strip()
            if tag_lower not in self.tag_index:
                self.tag_index[tag_lower] = []
            if filename not in self.tag_index[tag_lower]:
                self.tag_index[tag_lower].append(filename)
    
    def get_by_family(self, family: str) -> List[ProtocolMetadata]:
        """Get all protocols for a reaction family"""
        filenames = self.family_index.get(family, [])
        return [self.protocols[f] for f in filenames if f in self.protocols]
    
    def get_by_tag(self, tag: str) -> List[ProtocolMetadata]:
        """Get all protocols with a specific tag"""
        tag_lower = tag.lower().strip()
        filenames = self.tag_index.get(tag_lower, [])
        return [self.protocols[f] for f in filenames if f in self.protocols]
    
    def search_tags(self, tags: List[str]) -> List[ProtocolMetadata]:
        """Get protocols matching any of the provided tags"""
        matching_files = set()
        for tag in tags:
            tag_lower = tag.lower().strip()
            matching_files.update(self.tag_index.get(tag_lower, []))
        return [self.protocols[f] for f in matching_files if f in self.protocols]
    
    def save(self, path: Optional[Path] = None):
        """Save index to JSON file"""
        save_path = path or self.index_path
        if save_path is None:
            raise ValueError("No index path specified")
        
        self.index_metadata['updated'] = datetime.utcnow().isoformat()
        self.index_metadata['num_protocols'] = len(self.protocols)
        
        index_data = {
            'metadata': self.index_metadata,
            'protocols': {k: v.to_dict() for k, v in self.protocols.items()},
            'family_index': self.family_index,
            'tag_index': self.tag_index
        }
        
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, 'w', encoding='utf-8') as f:
            json.dump(index_data, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Saved protocol index to {save_path} ({len(self.protocols)} protocols)")
    
    @classmethod
    def load(cls, path: Path) -> 'ProtocolIndex':
        """Load index from JSON file"""
        if not path.exists():
            raise FileNotFoundError(f"Index file not found: {path}")
        
        with open(path, 'r', encoding='utf-8') as f:
            index_data = json.load(f)
        
        index = cls(index_path=path)
        index.index_metadata = index_data.get('metadata', {})
        index.family_index = index_data.get('family_index', {})
        index.tag_index = index_data.get('tag_index', {})
        
        # Reconstruct protocol metadata objects
        protocols_data = index_data.get('protocols', {})
        for filename, proto_dict in protocols_data.items():
            metadata = ProtocolMetadata(
                filename=proto_dict['filename'],
                reaction_smiles=proto_dict['reaction_smiles'],
                reaction_family=proto_dict['reaction_family'],
                tags=proto_dict['tags'],
                notes=proto_dict['notes'],
                source_title=proto_dict['source_title'],
                source_journal=proto_dict['source_journal'],
                source_year=proto_dict['source_year'],
                source_doi=proto_dict['source_doi'],
                canonical_reaction=proto_dict.get('canonical_reaction', '')
            )
            index.protocols[filename] = metadata
        
        logger.info(f"Loaded protocol index from {path} ({len(index.protocols)} protocols)")
        return index


class ProtocolMatcher:
    """
    Match user reactions to standard protocols
    
    This class provides the main interface for protocol matching:
    1. Loads the protocol index (or builds it if missing)
    2. Finds best matching protocols using DRFP similarity
    3. Filters by reaction family and tags
    4. Extracts full protocol details
    """
    
    def __init__(
        self,
        protocol_dir: Optional[Path] = None,
        index_path: Optional[Path] = None,
        auto_rebuild_index: bool = True
    ):
        """
        Initialize protocol matcher
        
        Args:
            protocol_dir: Directory containing protocol JSON files
                         (default: data/protocol_db)
            index_path: Path to index file
                       (default: data/protocol_db/.index.json)
            auto_rebuild_index: Rebuild index if outdated or missing
        """
        # Set default paths
        if protocol_dir is None:
            # Assume we're in chemtools/, go up to project root
            project_root = Path(__file__).parent.parent.parent
            protocol_dir = project_root / 'data' / 'protocol_db'
        
        self.protocol_dir = Path(protocol_dir)
        
        if index_path is None:
            index_path = self.protocol_dir / '.protocol_index.json'
        
        self.index_path = Path(index_path)
        
        # Load or build index
        if self.index_path.exists() and not auto_rebuild_index:
            self.index = ProtocolIndex.load(self.index_path)
        else:
            logger.info("Building protocol index...")
            self.index = self._build_index()
            self.index.save(self.index_path)
    
    def _build_index(self) -> ProtocolIndex:
        """Build index from all JSON files in protocol directory"""
        index = ProtocolIndex(index_path=self.index_path)
        index.index_metadata['created'] = datetime.utcnow().isoformat()
        
        if not self.protocol_dir.exists():
            logger.warning(f"Protocol directory not found: {self.protocol_dir}")
            return index
        
        # Scan all JSON files (exclude examples subdirectory)
        json_files = [
            f for f in self.protocol_dir.glob('*.json')
            if f.is_file() and not f.name.startswith('.')
        ]
        
        logger.info(f"Found {len(json_files)} protocol files")
        
        for json_file in json_files:
            try:
                metadata = self._extract_metadata(json_file)
                index.add_protocol(metadata)
            except Exception as e:
                logger.error(f"Error processing {json_file.name}: {e}")
        
        return index
    
    def _extract_metadata(self, json_path: Path) -> ProtocolMetadata:
        """Extract metadata from a protocol JSON file"""
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        reaction = data.get('reaction', {})
        source = data.get('source', {})
        
        # Parse tags
        tags_str = reaction.get('tags', '')
        if isinstance(tags_str, str):
            tags = [t.strip() for t in tags_str.split(';') if t.strip()]
        elif isinstance(tags_str, list):
            tags = tags_str
        else:
            tags = []
        
        # Extract reaction SMARTS patterns
        reaction_smarts_raw = reaction.get('reaction_SMARTS', [])
        if isinstance(reaction_smarts_raw, list):
            reaction_smarts = [str(s).strip() for s in reaction_smarts_raw if s]
        else:
            reaction_smarts = []
        
        # Normalize reaction SMILES
        reaction_smiles = reaction.get('reaction_smiles', '')
        canonical_rxn = ''
        try:
            # Try to canonicalize the reaction
            canonical_rxn = self._canonicalize_reaction(reaction_smiles)
        except Exception as e:
            logger.debug(f"Could not canonicalize {json_path.name}: {e}")
            canonical_rxn = reaction_smiles
        
        metadata = ProtocolMetadata(
            filename=json_path.name,
            reaction_smiles=reaction_smiles,
            reaction_smarts=reaction_smarts,
            reaction_family=reaction.get('family', ''),
            tags=tags,
            notes=reaction.get('notes', ''),
            source_title=source.get('title', ''),
            source_journal=source.get('journal', ''),
            source_year=source.get('year', 0),
            source_doi=source.get('doi', ''),
            canonical_reaction=canonical_rxn
        )
        
        return metadata
    
    def _canonicalize_reaction(self, reaction_smiles: str) -> str:
        """Canonicalize a reaction SMILES"""
        if not reaction_smiles or '>>' not in reaction_smiles:
            return reaction_smiles
        
        parts = reaction_smiles.split('>>')
        if len(parts) != 2:
            return reaction_smiles
        
        reactants, products = parts
        
        # Canonicalize each side
        reactant_list = [normalize(r)['canonical'] for r in reactants.split('.') if r.strip()]
        product_list = [normalize(p)['canonical'] for p in products.split('.') if p.strip()]
        
        # Sort for consistency
        reactant_list.sort()
        product_list.sort()
        
        canonical = '.'.join(reactant_list) + '>>' + '.'.join(product_list)
        return canonical
    
    def find_best_match(
        self,
        reaction_smiles: str,
        reaction_family: Optional[str] = None,
        tags: Optional[List[str]] = None,
        k: int = 5,
        min_similarity: float = 0.0
    ) -> Dict[str, Any]:
        """
        Find best matching protocols for a reaction
        
        Args:
            reaction_smiles: Query reaction SMILES
            reaction_family: Optional reaction family filter
            tags: Optional list of tags to filter by
            k: Number of top matches to return
            min_similarity: Minimum similarity threshold (0.0-1.0)
        
        Returns:
            Dictionary with:
                - matches: List of ProtocolMatch objects
                - query: Query information
                - metadata: Search metadata
        """
        # Canonicalize query reaction
        try:
            canonical_query = self._canonicalize_reaction(reaction_smiles)
        except Exception as e:
            logger.warning(f"Could not canonicalize query: {e}")
            canonical_query = reaction_smiles
        
        # Filter protocols by family and tags
        candidate_protocols = list(self.index.protocols.values())
        
        if reaction_family:
            candidate_protocols = [
                p for p in candidate_protocols
                if p.reaction_family.lower() == reaction_family.lower()
            ]
        
        if tags:
            tag_matches = self.index.search_tags(tags)
            tag_filenames = {p.filename for p in tag_matches}
            candidate_protocols = [
                p for p in candidate_protocols
                if p.filename in tag_filenames
            ]
        
        if not candidate_protocols:
            return {
                'matches': [],
                'query': {
                    'reaction_smiles': reaction_smiles,
                    'canonical': canonical_query,
                    'family': reaction_family,
                    'tags': tags
                },
                'metadata': {
                    'num_candidates': 0,
                    'num_total_protocols': len(self.index.protocols)
                }
            }
        
        # Compute similarities
        matches = []
        for protocol in candidate_protocols:
            try:
                similarity = compute_drfp_similarity(
                    canonical_query,
                    protocol.canonical_reaction or protocol.reaction_smiles
                )
                
                if similarity >= min_similarity:
                    match = ProtocolMatch(
                        filename=protocol.filename,
                        similarity=similarity,
                        reaction_smiles=protocol.reaction_smiles,
                        reaction_smarts=protocol.reaction_smarts,
                        reaction_family=protocol.reaction_family,
                        tags=protocol.tags,
                        notes=protocol.notes,
                        source_title=protocol.source_title,
                        source_doi=protocol.source_doi
                    )
                    matches.append(match)
            except Exception as e:
                logger.debug(f"Error computing similarity for {protocol.filename}: {e}")
        
        # Sort by similarity (descending)
        matches.sort(key=lambda m: m.similarity, reverse=True)
        
        # Take top k
        top_matches = matches[:k]
        
        return {
            'matches': [m.to_dict() for m in top_matches],
            'query': {
                'reaction_smiles': reaction_smiles,
                'canonical': canonical_query,
                'family': reaction_family,
                'tags': tags
            },
            'metadata': {
                'num_candidates': len(candidate_protocols),
                'num_matches': len(matches),
                'num_total_protocols': len(self.index.protocols),
                'min_similarity': min_similarity
            }
        }
    
    def get_protocol(self, filename: str) -> Optional[Dict[str, Any]]:
        """
        Load and return full protocol details
        
        Args:
            filename: Protocol filename (e.g., 'Suzuki_Cu_C(sp3)-C(sp2).json')
        
        Returns:
            Complete protocol dictionary or None if not found
        """
        protocol_path = self.protocol_dir / filename
        
        if not protocol_path.exists():
            logger.error(f"Protocol file not found: {filename}")
            return None
        
        try:
            with open(protocol_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            logger.error(f"Error loading protocol {filename}: {e}")
            return None
    
    def extract_conditions(self, protocol_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract key reaction conditions from protocol
        
        Args:
            protocol_data: Full protocol dictionary
        
        Returns:
            Simplified conditions dictionary with:
                - catalyst
                - ligand
                - base
                - solvent
                - temperature
                - time
                - atmosphere
        """
        conditions = {
            'catalyst': None,
            'ligand': None,
            'base': None,
            'solvent': None,
            'additives': [],
            'temperature_C': None,
            'time_h': None,
            'atmosphere': None
        }
        
        reaction_setup = protocol_data.get('reaction_setup', [])
        if not reaction_setup:
            return conditions
        
        # Get chemicals from first setup step
        chemicals = reaction_setup[0].get('chemicals', [])
        
        for chem in chemicals:
            role = chem.get('role', '')
            name = chem.get('name', '')
            abbrev = chem.get('abbreviation')
            display_name = abbrev if abbrev else name
            
            if role in ('metal_precursor', 'preformed_catalyst', 'co_catalyst'):
                if conditions['catalyst'] is None:
                    conditions['catalyst'] = display_name
                else:
                    conditions['additives'].append(display_name)
            elif role == 'ligand':
                conditions['ligand'] = display_name
            elif role == 'base':
                conditions['base'] = display_name
            elif role == 'solvent':
                conditions['solvent'] = display_name
            elif role == 'additive':
                conditions['additives'].append(display_name)
        
        # Get conditions from first setup step
        cond_list = reaction_setup[0].get('conditions', [])
        if cond_list:
            # Use the main reaction conditions (usually the last or longest step)
            main_cond = cond_list[-1] if len(cond_list) > 1 else cond_list[0]
            conditions['temperature_C'] = main_cond.get('temperature_C')
            conditions['time_h'] = main_cond.get('time_h')
            conditions['atmosphere'] = main_cond.get('atmosphere')
        
        return conditions
    
    def rebuild_index(self):
        """Rebuild the protocol index from scratch"""
        logger.info("Rebuilding protocol index...")
        self.index = self._build_index()
        self.index.save(self.index_path)
        logger.info(f"Index rebuilt with {len(self.index.protocols)} protocols")
    
    def list_families(self) -> List[Tuple[str, int]]:
        """List all reaction families with protocol counts"""
        families = [
            (family, len(files))
            for family, files in self.index.family_index.items()
        ]
        families.sort(key=lambda x: x[1], reverse=True)
        return families
    
    def list_tags(self) -> List[Tuple[str, int]]:
        """List all tags with protocol counts"""
        tags = [
            (tag, len(files))
            for tag, files in self.index.tag_index.items()
        ]
        tags.sort(key=lambda x: x[1], reverse=True)
        return tags


def build_protocol_index(
    protocol_dir: Optional[Path] = None,
    output_path: Optional[Path] = None
) -> ProtocolIndex:
    """
    Standalone function to build protocol index
    
    Usage:
        from chemtools.protocol import build_protocol_index
        index = build_protocol_index()
    """
    matcher = ProtocolMatcher(
        protocol_dir=protocol_dir,
        index_path=output_path,
        auto_rebuild_index=True
    )
    return matcher.index
