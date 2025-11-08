"""
Proof of Concept: Similarity-Based Rule Database Selection

Demonstrates using DRFP similarity (like protocol recommendation) to select
the correct rule database instead of relying on brittle family detection.
"""

import json
from pathlib import Path
from typing import List, Dict, Any, Optional
import numpy as np

try:
    from drfp import DrfpEncoder
    DRFP_AVAILABLE = True
except ImportError:
    DRFP_AVAILABLE = False
    print("⚠️  DRFP not available - install with: pip install drfp")


class SimilarityBasedRuleSelector:
    """
    Select the best rule database using DRFP similarity matching
    
    Similar to protocol recommendation but for rule databases.
    """
    
    def __init__(self, rule_db_dir: Path = None):
        if rule_db_dir is None:
            rule_db_dir = Path("data/rule_db")
        
        self.rule_db_dir = Path(rule_db_dir)
        
        # Initialize encoder first
        if DRFP_AVAILABLE:
            self.encoder = DrfpEncoder()  # Correct initialization
        else:
            self.encoder = None
        
        # Then load databases and compute DRFPs
        self.databases = self._load_databases()
        self.reference_drfps = self._precompute_reference_drfps()
    
    def _load_databases(self) -> Dict[str, Dict[str, Any]]:
        """Load all rule database JSON files"""
        databases = {}
        
        for json_file in self.rule_db_dir.glob("*.json"):
            try:
                with open(json_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    databases[json_file.stem] = {
                        'path': json_file,
                        'data': data,
                        'name': data.get('name', json_file.stem)
                    }
            except Exception as e:
                print(f"⚠️  Failed to load {json_file}: {e}")
        
        print(f"✅ Loaded {len(databases)} rule databases")
        return databases
    
    def _get_reference_reactions(self, db_data: dict) -> List[str]:
        """
        Extract or generate reference reactions for a database
        
        For now, manually define representative reactions for each family.
        In full implementation, these would be in the JSON schema.
        """
        db_name = db_data.get('name', '').lower()
        
        # Hardcoded reference reactions for proof-of-concept
        # In production, these would be in the JSON files
        references = {
            'sonogashira': [
                'Ic1ccccc1.C#Cc1ccccc1>>C(#Cc1ccccc1)c1ccccc1',
                'Brc1ccccc1.C#CCCCCC>>C(#CCCCCC)c1ccccc1',
                'Clc1ccc(F)cc1.C#Cc1ccccc1>>Fc1ccc(C#Cc2ccccc2)cc1',
            ],
            'suzuki': [
                'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
                'Ic1cccnc1.c1ccc(B(O)O)cc1>>c1ccc(-c2cccnc2)cc1',
            ],
            'c_o_coupling': [
                'Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1',
                'Ic1ccccc1.CCO>>CCOc1ccccc1',
            ],
            'c_n_coupling_pd': [
                'Brc1ccccc1.c1ccccc1N>>c1ccc(Nc2ccccc2)cc1',
                'Clc1ccccc1.CCN>>CCNc1ccccc1',
            ],
            'c_n_coupling_cu': [
                'Ic1ccccc1.c1ccccc1N>>c1ccc(Nc2ccccc2)cc1',
                'Brc1cccnc1.CCN>>CCNc1cccnc1',
            ],
            'rcm': [
                'C=CCNC(=O)C=C>>C1=CCNC(=O)C1',
                'C=CCCCCCC=C>>C1=CCCCCC1',
            ],
            'amide_formation': [
                'CC(=O)O.CCN>>CC(=O)NCC',
                'c1ccccc1C(=O)O.c1ccccc1N>>O=C(Nc1ccccc1)c1ccccc1',
            ],
            'reductive_amination': [
                'CC(=O)C.CCN>>CC(C)NCC',
                'O=Cc1ccccc1.c1ccccc1N>>c1ccc(CNc2ccccc2)cc1',
            ],
            'snar': [
                'Fc1ccccc1[N+](=O)[O-].CCO>>CCOc1ccccc1[N+](=O)[O-]',
                'Clc1cnccc1.c1ccccc1N>>c1ccc(Nc2cccnc2)cc1',
            ],
        }
        
        # Try to find matching references
        for key in references:
            if key in db_name.lower():
                return references[key]
        
        # Fallback: return empty list (database won't match well)
        print(f"⚠️  No reference reactions defined for {db_name}")
        return []
    
    def _precompute_reference_drfps(self) -> Dict[str, List[np.ndarray]]:
        """Precompute DRFP fingerprints for all reference reactions"""
        if not DRFP_AVAILABLE or self.encoder is None:
            print("⚠️  DRFP not available - cannot precompute fingerprints")
            return {}
        
        reference_drfps = {}
        
        for db_name, db_info in self.databases.items():
            refs = self._get_reference_reactions(db_info['data'])
            
            if not refs:
                continue
            
            drfps = []
            for rxn_smiles in refs:
                try:
                    fp = self.encoder.encode([rxn_smiles])[0]
                    drfps.append(fp)
                except Exception as e:
                    print(f"⚠️  Failed to encode {rxn_smiles}: {e}")
            
            if drfps:
                reference_drfps[db_name] = drfps
                print(f"  {db_name}: {len(drfps)} reference DRFPs")
        
        return reference_drfps
    
    def compute_drfp(self, reaction_smiles: str) -> Optional[np.ndarray]:
        """Compute DRFP fingerprint for a reaction"""
        if not DRFP_AVAILABLE or self.encoder is None:
            return None
        
        try:
            fp = self.encoder.encode([reaction_smiles])[0]
            return fp
        except Exception as e:
            print(f"❌ Failed to compute DRFP: {e}")
            return None
    
    def cosine_similarity(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Compute cosine similarity between two fingerprints"""
        dot_product = np.dot(fp1, fp2)
        norm1 = np.linalg.norm(fp1)
        norm2 = np.linalg.norm(fp2)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        return float(dot_product / (norm1 * norm2))
    
    def select_database(
        self,
        reaction_smiles: str,
        top_k: int = 3
    ) -> List[Dict[str, Any]]:
        """
        Select the best rule database(s) for a reaction using similarity
        
        Returns:
            List of database matches with scores, sorted by similarity
        """
        if not DRFP_AVAILABLE:
            print("❌ DRFP not available - cannot use similarity matching")
            return []
        
        # Compute query DRFP
        query_drfp = self.compute_drfp(reaction_smiles)
        if query_drfp is None:
            return []
        
        # Score each database
        matches = []
        
        for db_name, ref_drfps in self.reference_drfps.items():
            if not ref_drfps:
                continue
            
            # Compute similarity with each reference reaction
            similarities = [
                self.cosine_similarity(query_drfp, ref_drfp)
                for ref_drfp in ref_drfps
            ]
            
            # Aggregate score (weighted: max 70%, avg 30%)
            max_sim = max(similarities)
            avg_sim = np.mean(similarities)
            score = 0.7 * max_sim + 0.3 * avg_sim
            
            matches.append({
                'database': db_name,
                'score': score,
                'max_similarity': max_sim,
                'avg_similarity': avg_sim,
                'num_references': len(similarities),
                'path': str(self.databases[db_name]['path'])
            })
        
        # Sort by score (descending)
        matches.sort(key=lambda x: x['score'], reverse=True)
        
        return matches[:top_k]


def test_similarity_selection():
    """Test the similarity-based selection with example reactions"""
    
    print("=" * 80)
    print("PROOF OF CONCEPT: Similarity-Based Rule Database Selection")
    print("=" * 80)
    print()
    
    if not DRFP_AVAILABLE:
        print("❌ DRFP not installed. Install with: pip install drfp")
        return
    
    # Initialize selector
    selector = SimilarityBasedRuleSelector()
    print()
    
    # Test cases
    test_reactions = [
        {
            'name': 'Sonogashira (User\'s reaction)',
            'smiles': 'Ic1ccccc1C(C)(C)C.C#Cc1ccccc1>>CC(C)(C)c1ccccc1C#Cc1ccccc1',
            'expected': 'sonogashira_db'
        },
        {
            'name': 'RCM (Previously misdetected as SNAr)',
            'smiles': 'C=CCNC(=O)C=C>>C1=CCNC(=O)C1',
            'expected': 'RCM_db'
        },
        {
            'name': 'C-O Coupling',
            'smiles': 'Brc1ccccc1.Oc1ccccc1>>c1ccc(Oc2ccccc2)cc1',
            'expected': 'C_O_coupling_db'
        },
        {
            'name': 'Suzuki',
            'smiles': 'Brc1ccccc1.c1ccc(B(O)O)cc1>>c1ccc(-c2ccccc2)cc1',
            'expected': 'Suzuki_db'
        },
    ]
    
    # Test each reaction
    for i, test in enumerate(test_reactions, 1):
        print(f"\n{'=' * 80}")
        print(f"Test {i}: {test['name']}")
        print(f"{'=' * 80}")
        print(f"Reaction: {test['smiles'][:60]}...")
        print(f"Expected: {test['expected']}")
        print()
        
        matches = selector.select_database(test['smiles'], top_k=5)
        
        if not matches:
            print("❌ No matches found")
            continue
        
        print("Top 5 Database Matches:")
        print("-" * 80)
        
        for rank, match in enumerate(matches, 1):
            marker = "✅" if rank == 1 and match['database'] == test['expected'] else "  "
            print(f"{marker} {rank}. {match['database']:<25} Score: {match['score']:.4f} "
                  f"(max: {match['max_similarity']:.4f}, avg: {match['avg_similarity']:.4f})")
        
        # Check if top match is correct
        top_match = matches[0]['database']
        if top_match == test['expected']:
            print(f"\n✅ CORRECT: Top match is {test['expected']}")
        else:
            print(f"\n❌ INCORRECT: Expected {test['expected']}, got {top_match}")


if __name__ == '__main__':
    test_similarity_selection()
