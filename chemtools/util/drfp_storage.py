"""
Binary storage utilities for DRFP fingerprints.

This module provides efficient storage and retrieval of DRFP fingerprints
using NumPy's compressed NPZ format instead of embedding 4096-element arrays
in JSONL files.

Storage Benefits:
- Reduces JSONL file size by ~90% (e.g., 670 MB → ~70 MB)
- Faster loading with NumPy binary format
- Uses uint8 arrays (0-255) for compact storage
- Compressed NPZ format for additional space savings

Usage:
    # Save fingerprints to binary file
    save_drfp_index(fingerprints, reaction_ids, "data/C_N_Coupling_Cu_drfp.npz")
    
    # Load fingerprints by reaction_id
    loader = DRFPLoader("data/C_N_Coupling_Cu_drfp.npz")
    fp = loader.get_fingerprint("31-172-CAS-23306204")
"""

from __future__ import annotations

import os
import numpy as np
from typing import Dict, List, Optional, Union
from pathlib import Path


class DRFPLoader:
    """
    Efficient loader for DRFP fingerprints stored in NPZ format.
    
    Loads the entire NPZ file once and provides fast indexed access
    by reaction_id using a hash map.
    """
    
    def __init__(self, npz_path: str | Path):
        """
        Initialize loader with NPZ file path.
        
        Args:
            npz_path: Path to .npz file containing fingerprints
        """
        self.npz_path = str(npz_path)
        self._data: Optional[np.lib.npzfile.NpzFile] = None
        self._fps: Optional[np.ndarray] = None
        self._index: Optional[Dict[str, int]] = None
        self._n_bits: int = 4096
        self._radius: int = 3
        self._loaded = False
    
    def _load(self) -> None:
        """Lazy load NPZ file and build index."""
        if self._loaded:
            return
        
        if not os.path.exists(self.npz_path):
            raise FileNotFoundError(f"DRFP file not found: {self.npz_path}")
        
        # Load compressed NPZ
        data = np.load(self.npz_path, allow_pickle=True)
        
        # Extract arrays
        self._fps = data['fps']  # Shape: (N, 4096), dtype=uint8
        reaction_ids = data['reaction_ids']  # Shape: (N,), dtype=object (strings)
        
        # Store metadata
        if 'n_bits' in data:
            self._n_bits = int(data['n_bits'])
        if 'radius' in data:
            self._radius = int(data['radius'])
        
        # Build index: reaction_id -> row number
        self._index = {str(rid): i for i, rid in enumerate(reaction_ids)}
        
        self._loaded = True
        print(f"Loaded {len(self._index)} DRFP fingerprints from {self.npz_path}")
    
    def get_fingerprint(self, reaction_id: str) -> Optional[np.ndarray]:
        """
        Get fingerprint for a specific reaction_id.
        
        Args:
            reaction_id: Reaction identifier (e.g., "31-172-CAS-23306204")
        
        Returns:
            NumPy array of shape (4096,) with uint8 values, or None if not found
        """
        self._load()
        
        idx = self._index.get(reaction_id)
        if idx is None:
            return None
        
        return self._fps[idx]
    
    def get_all_fingerprints(self) -> np.ndarray:
        """
        Get all fingerprints as a matrix.
        
        Returns:
            NumPy array of shape (N, 4096) with uint8 values
        """
        self._load()
        return self._fps
    
    def get_all_reaction_ids(self) -> List[str]:
        """
        Get list of all reaction IDs in order.
        
        Returns:
            List of reaction ID strings
        """
        self._load()
        return list(self._index.keys())
    
    @property
    def n_bits(self) -> int:
        """Get number of bits in fingerprints (default: 4096)."""
        self._load()
        return self._n_bits
    
    @property
    def radius(self) -> int:
        """Get radius parameter used for fingerprints (default: 3)."""
        self._load()
        return self._radius
    
    def __len__(self) -> int:
        """Return number of fingerprints stored."""
        self._load()
        return len(self._index)


def save_drfp_index(
    fingerprints: List[Union[np.ndarray, List[int]]],
    reaction_ids: List[str],
    output_path: str | Path,
    n_bits: int = 4096,
    radius: int = 3
) -> None:
    """
    Save DRFP fingerprints to compressed NPZ format.
    
    Args:
        fingerprints: List of fingerprint arrays (each 4096 elements)
        reaction_ids: List of reaction IDs (must match fingerprints length)
        output_path: Path to save .npz file
        n_bits: Number of bits in fingerprints (default: 4096)
        radius: Radius parameter used (default: 3)
    
    Example:
        >>> fps = [[0, 1, 0, ...], [0, 0, 2, ...]]  # 4096 elements each
        >>> ids = ["rxn_001", "rxn_002"]
        >>> save_drfp_index(fps, ids, "data/drfp.npz")
    """
    if len(fingerprints) != len(reaction_ids):
        raise ValueError(
            f"Mismatch: {len(fingerprints)} fingerprints but {len(reaction_ids)} reaction_ids"
        )
    
    if not fingerprints:
        raise ValueError("No fingerprints to save")
    
    # Convert to NumPy arrays
    fps_list = []
    for fp in fingerprints:
        if isinstance(fp, np.ndarray):
            fps_list.append(fp.astype('uint8'))
        else:
            fps_list.append(np.array(fp, dtype='uint8'))
    
    # Stack into matrix: (N, 4096)
    fps_matrix = np.vstack(fps_list)
    
    # Convert reaction_ids to object array
    reaction_ids_array = np.array(reaction_ids, dtype=object)
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
    
    # Save as compressed NPZ
    np.savez_compressed(
        output_path,
        fps=fps_matrix,
        reaction_ids=reaction_ids_array,
        n_bits=np.array(n_bits, dtype='int32'),
        radius=np.array(radius, dtype='int32')
    )
    
    # Report file size
    file_size_mb = os.path.getsize(output_path) / (1024 * 1024)
    print(f"Saved {len(reaction_ids)} fingerprints to {output_path}")
    print(f"File size: {file_size_mb:.2f} MB ({file_size_mb/len(reaction_ids)*1000:.1f} KB per reaction)")


def extract_drfp_from_jsonl(
    jsonl_path: str | Path,
    output_npz_path: Optional[str | Path] = None
) -> Dict[str, np.ndarray]:
    """
    Extract DRFP fingerprints from JSONL file's 'precomputed.drfp_fp' field.
    
    Args:
        jsonl_path: Path to JSONL file with embedded drfp_fp arrays
        output_npz_path: Optional path to save NPZ file. If None, only returns dict
    
    Returns:
        Dictionary mapping reaction_id -> fingerprint array
    
    Example:
        >>> fps = extract_drfp_from_jsonl("data/C_N_Coupling_Cu.jsonl", 
        ...                                "data/C_N_Coupling_Cu_drfp.npz")
    """
    import json
    
    fingerprints = []
    reaction_ids = []
    n_bits = 4096
    radius = 3
    
    with open(jsonl_path, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                continue
            
            # Get reaction_id
            reaction_id = record.get('reaction_id')
            if not reaction_id:
                continue
            
            # Extract DRFP fingerprint from precomputed field
            precomputed = record.get('precomputed', {})
            if not isinstance(precomputed, dict):
                continue
            
            drfp_fp = precomputed.get('drfp_fp')
            if drfp_fp is None:
                continue
            
            # Get metadata if available
            if 'drfp_n_bits' in precomputed:
                n_bits = precomputed['drfp_n_bits']
            if 'drfp_radius' in precomputed:
                radius = precomputed['drfp_radius']
            
            fingerprints.append(drfp_fp)
            reaction_ids.append(reaction_id)
    
    print(f"Extracted {len(fingerprints)} DRFP fingerprints from {jsonl_path}")
    
    # Save if output path provided
    if output_npz_path and fingerprints:
        save_drfp_index(fingerprints, reaction_ids, output_npz_path, n_bits, radius)
    
    # Return as dict for convenience
    return {rid: np.array(fp, dtype='uint8') for rid, fp in zip(reaction_ids, fingerprints)}


def get_drfp_path_for_family(family: str, data_dir: str = "data/reaction_dataset") -> str:
    """
    Get the standard DRFP NPZ path for a reaction family.
    
    Args:
        family: Reaction family name (e.g., "C_N_Coupling_Cu")
        data_dir: Directory containing reaction datasets
    
    Returns:
        Path to DRFP NPZ file (e.g., "data/reaction_dataset/C_N_Coupling_Cu_drfp.npz")
    """
    return os.path.join(data_dir, f"{family}_drfp.npz")
