"""
V2 Processing Core for SciFinder RDF data.
Uses taxonomy-based classification and structural featurization.
"""

import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.structural import featurize_reaction
from chemtools.taxonomy.reagent_v2 import ReagentTaxonomyV2, classify_reagent_v2
from chemtools.util.rdkit_helpers import parse_smiles

class V2Processor:
    def __init__(self):
        self.taxonomy = ReagentTaxonomyV2.from_path()
        self.roles = self.taxonomy.roles
        self.families = self.taxonomy.families

    def classify_reagent(self, name: str, cas: str = "", smiles: str = "") -> Dict[str, Any]:
        """Classify a reagent using the V2 taxonomy."""
        record = {"name": name, "cas": cas, "smiles": smiles}
        match = classify_reagent_v2(record, self.roles, self.families)
        
        if match:
            return {
                "name": name,
                "cas": cas,
                "role": match.role_id,
                "family": match.family_id,
                "match_kind": match.match_kind
            }
        return {
            "name": name,
            "cas": cas,
            "role": "UNKNOWN",
            "family": "UNKNOWN",
            "match_kind": None
        }

    def process_reaction(self, reaction_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process a raw reaction from SciFinder into the V2 Recommendation format.
        """
        reaction_smiles = reaction_data.get("reaction_smiles")
        if not reaction_smiles:
            return {}

        # 1. Structural Analysis
        try:
            analysis = featurize_reaction(reaction_smiles)
        except Exception as e:
            print(f"Error featurizing reaction: {e}")
            analysis = {}

        # 2. Extract Recommendation Features
        rec_features = self._extract_recommendation_features(analysis)

        # 3. Standardize Conditions
        raw_reagents = reaction_data.get("reagents", [])
        raw_solvents = reaction_data.get("solvents", [])
        
        conditions = self._standardize_conditions(raw_reagents, raw_solvents, reaction_data)

        # 4. Final V2 Entry (Simplified & Pruned)
        v2_entry = {
            "reaction_id": reaction_data.get("reaction_id"),
            "reaction_smiles": reaction_smiles,
            "conditions": {k: v for k, v in conditions.items() if v not in (None, [], "")},
            "recommendation_payload": {
                "family": rec_features.get("reaction_type"),
                "features": {k: v for k, v in rec_features.items() if v not in (None, [], "", False)}
            }
        }
        
        # Add reference only if it has content
        ref = reaction_data.get("reference")
        if ref:
            v2_entry["reference"] = ref
            
        return v2_entry

    def _extract_recommendation_features(self, analysis: Dict[str, Any]) -> Dict[str, Any]:
        features = {}
        
        # Reaction Type
        detection = analysis.get("detection", {})
        matches = detection.get("matches", [])
        if matches:
            best_match = matches[0]
            features["reaction_type"] = best_match.get("reaction_type")
            features["reaction_category"] = best_match.get("category")
        
        # Sterics (Rounded to 2 decimal places)
        steric = analysis.get("steric", {})
        aryl_sterics = [s.get("result", {}).get("score_0_10", 0) for s in steric.get("aryl", [])]
        alkyl_sterics = [s.get("result", {}).get("score_0_10", 0) for s in steric.get("alkyl", [])]
        features["max_aryl_steric"] = round(max(aryl_sterics), 2) if aryl_sterics else 0
        features["max_alkyl_steric"] = round(max(alkyl_sterics), 2) if alkyl_sterics else 0
        
        # Electronics (Rounded to 2 decimal places)
        electronics = analysis.get("electronics", {})
        aryl_elec = [e.get("result", {}).get("score_0_10", 5.0) for e in electronics.get("aryl", [])]
        features["avg_aryl_electronic"] = round(sum(aryl_elec) / len(aryl_elec), 2) if aryl_elec else 5.0
        
        # Chelators
        nearby = analysis.get("nearby", [])
        has_chelator = False
        for n in nearby:
            for res in n.get("result", []):
                tags = res.get("tags", [])
                if "bidentate_chelator" in tags or "chelator" in tags:
                    has_chelator = True
                    break
        if has_chelator:
            features["has_chelator"] = True
        
        # Motifs (Sparse: only store True)
        motifs = analysis.get("motifs", [])
        for m in motifs:
            features[f"has_{m['compound_id']}"] = True
            
        return features

    def _standardize_conditions(self, raw_reagents: List[Dict[str, str]], raw_solvents: List[Dict[str, str]], reaction_data: Dict[str, Any]) -> Dict[str, Any]:
        conditions = {
            "metal": None,
            "ligand": [],
            "base": [],
            "additive": [],
            "solvent": [],
            "temperature": None,
            "yield": None
        }

        # Process reagents
        for r in raw_reagents:
            name = r.get("name", "")
            cas = r.get("cas", "")
            classification = self.classify_reagent(name, cas)
            role = classification["role"]
            
            if role == "metal_catalyst":
                conditions["metal"] = name
            elif role == "ligand":
                conditions["ligand"].append(name)
            elif role == "base":
                conditions["base"].append(name)
            else:
                conditions["additive"].append(name)

        # Process solvents
        for s in raw_solvents:
            conditions["solvent"].append(s.get("name", ""))

        # Temperature and Yield (Rounded)
        try:
            temp = reaction_data.get("conditions", {}).get("temperature_c")
            if temp is not None:
                conditions["temperature"] = round(float(temp), 1)
        except: pass

        try:
            yld = reaction_data.get("conditions", {}).get("yield_pct")
            if yld is not None:
                conditions["yield"] = round(float(yld), 1)
        except: pass

        return conditions
