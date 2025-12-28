"""
Feature-based KNN recommender for reaction conditions.
Uses interpretable chemical features (sterics, electronics, motifs) for similarity.
"""

import json
import numpy as np
from typing import List, Dict, Any, Optional
from pathlib import Path

class FeatureKNNRecommender:
    def __init__(self, dataset_path: str):
        self.dataset_path = Path(dataset_path)
        self.data = []
        self.feature_matrix = []
        self.feature_names = []
        self._load_dataset()

    def _load_dataset(self):
        if not self.dataset_path.exists():
            return

        all_features = []
        with open(self.dataset_path, 'r', encoding='utf-8') as f:
            for line in f:
                entry = json.loads(line)
                payload = entry.get("recommendation_payload", {})
                features = payload.get("features", {})
                if not features:
                    continue
                
                self.data.append(entry)
                all_features.append(features)

        if not all_features:
            return

        # Collect all unique feature keys
        keys = set()
        for f in all_features:
            keys.update(f.keys())
        
        # Filter out non-numeric/boolean features for the matrix
        self.feature_names = sorted([k for k in keys if not k.startswith(("reaction_type", "reaction_category"))])
        
        # Build matrix
        matrix = []
        for f in all_features:
            row = []
            for name in self.feature_names:
                val = f.get(name, 0)
                if isinstance(val, bool):
                    val = 1.0 if val else 0.0
                elif val is None:
                    val = 0.0
                row.append(float(val))
            matrix.append(row)
        
        self.feature_matrix = np.array(matrix)

    def recommend(self, query_features: Dict[str, Any], top_k: int = 5) -> List[Dict[str, Any]]:
        if len(self.feature_matrix) == 0:
            return []

        # Build query vector
        query_vec = []
        for name in self.feature_names:
            val = query_features.get(name, 0)
            if isinstance(val, bool):
                val = 1.0 if val else 0.0
            elif val is None:
                val = 0.0
            query_vec.append(float(val))
        
        query_vec = np.array(query_vec)

        # Compute Euclidean distance (simple for now)
        # In a real scenario, we might want to weight features (e.g. sterics > motifs)
        distances = np.linalg.norm(self.feature_matrix - query_vec, axis=1)
        
        # Get top K indices
        indices = np.argsort(distances)[:top_k]
        
        results = []
        for idx in indices:
            entry = self.data[idx]
            results.append({
                "reaction_id": entry["reaction_id"],
                "reaction_smiles": entry["reaction_smiles"],
                "conditions": entry["conditions"],
                "distance": float(distances[idx]),
                "similarity_score": 1.0 / (1.0 + float(distances[idx])),
                "features": entry["recommendation_payload"]["features"]
            })
            
        return results

def get_recommender(dataset_name: str = "C_N_Coupling_Recommendation_V2.jsonl"):
    base_path = Path(__file__).resolve().parent.parent.parent / "data" / "reaction_dataset"
    return FeatureKNNRecommender(base_path / dataset_name)
