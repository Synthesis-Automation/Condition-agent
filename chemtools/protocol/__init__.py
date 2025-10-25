"""
Protocol Database Package

Provides functionality for matching user reactions to standard experimental protocols
and extracting detailed procedure information.

Main components:
- ProtocolIndexer: Build and maintain protocol database index
- ProtocolRecommender: Find best matching protocols using DRFP similarity
- CLI: Command-line tools for index management

Example usage:
    # Build index from CLI
    python -m chemtools.protocol.cli build
    
    # Get protocol recommendations
    from chemtools.protocol import ProtocolRecommender
    
    recommender = ProtocolRecommender()
    results = recommender.recommend(
        reaction_smiles='CCBr.c1ccccc1B(O)O>>CCc1ccccc1',
        k=5
    )
    
    for match in results['matches']:
        print(f"{match['similarity']:.3f}: {match['source_title']}")
"""

from .indexer import (
    ProtocolIndexer,
    ProtocolRecord,
    build_index,
)

from .recommend import (
    ProtocolRecommender,
    recommend_protocol,
)

from .smarts_generator_cli import (
    SmartsGenerator,
    ReactionSmartsApplicability,
    visualize_smarts_pattern,
    visualize_reaction_smarts,
    visualize_pattern_with_examples,
    main as smarts_generator_main,
    RDKIT_AVAILABLE as PROTOCOL_SMARTS_RDKIT_AVAILABLE,
)

from .batch_update_protocol_smarts import (
    ProtocolSmartsUpdater,
    ProcessingResult,
)

__all__ = [
    'ProtocolIndexer',
    'ProtocolRecord',
    'build_index',
    'ProtocolRecommender',
    'recommend_protocol',
    'SmartsGenerator',
    'ReactionSmartsApplicability',
    'visualize_smarts_pattern',
    'visualize_reaction_smarts',
    'visualize_pattern_with_examples',
    'smarts_generator_main',
    'PROTOCOL_SMARTS_RDKIT_AVAILABLE',
    'ProtocolSmartsUpdater',
    'ProcessingResult',
]
