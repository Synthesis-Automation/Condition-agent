"""
V2 Processing Core for SciFinder RDF data.
Uses taxonomy-based classification and structural featurization.
"""

import json
import os
import sys
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from chemtools.featurizers.structural import featurize_reaction
from chemtools.reagent.taxonomy_store import RoleHeuristics, TaxonomyStore
from chemtools.util.rdkit_helpers import parse_smiles

class V2Processor:
    def __init__(self):
        self.store = TaxonomyStore()
        self.taxonomy = self.store.taxonomy
        if self.taxonomy:
            self.roles = self.taxonomy.roles
            self.families = self.taxonomy.families
        else:
            self.roles = {}
            self.families = {}
        self.heuristics = RoleHeuristics(self.store)

    def classify_reagent(self, name: str, cas: str = "", smiles: str = "") -> Dict[str, Any]:
        """Classify a reagent using the V2 taxonomy."""
        # If name looks like a CAS and cas is empty, swap them
        if not cas and name and re.match(r'^\d{2,7}-\d{2}-\d$', name):
            cas = name

        match = self.store.classify_reagent(name=name, cas=cas, smiles=smiles)
        
        if match:
            return {
                "name": name,
                "cas": cas,
                "role": match["role_id"],
                "family": match["family_id"],
                "match_kind": match["match_kind"]
            }
            
        # Fallback to heuristics if taxonomy match fails
        # Try both name and cas (if cas was swapped from name)
        inferred = self.heuristics.infer_family(name, [])
        if not inferred:
            inferred_role = self.heuristics.infer_role(name, [])
            if inferred_role:
                role, pattern = inferred_role
                inferred = (role, "INFERRED", [pattern])
                
        if not inferred and cas:
            inferred = self.heuristics.infer_family(cas, [])
            if not inferred:
                inferred_role = self.heuristics.infer_role(cas, [])
                if inferred_role:
                    role, pattern = inferred_role
                    inferred = (role, "INFERRED", [pattern])
            
        if inferred:
            role, family_id, matches = inferred
            return {
                "name": name,
                "cas": cas,
                "role": role,
                "family": family_id,
                "match_kind": "heuristic"
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
        aryl_sterics = []
        for s in steric.get("aryl", []):
            if isinstance(s, dict):
                res = s.get("result", {})
                if isinstance(res, dict):
                    aryl_sterics.append(res.get("score_0_10", 0))
        
        alkyl_sterics = []
        for s in steric.get("alkyl", []):
            if isinstance(s, dict):
                res = s.get("result", {})
                if isinstance(res, dict):
                    alkyl_sterics.append(res.get("score_0_10", 0))
                    
        features["max_aryl_steric"] = round(max(aryl_sterics), 2) if aryl_sterics else 0
        features["max_alkyl_steric"] = round(max(alkyl_sterics), 2) if alkyl_sterics else 0
        
        # Electronics (Rounded to 2 decimal places)
        electronics = analysis.get("electronics", {})
        aryl_elec = []
        for e in electronics.get("aryl", []):
            if isinstance(e, dict):
                res = e.get("result", {})
                if isinstance(res, dict):
                    aryl_elec.append(res.get("score_0_10", 5.0))
                elif isinstance(res, list) and res:
                    # Handle "both" case where result is a list of dicts
                    first_res = res[0]
                    if isinstance(first_res, dict):
                        aryl_elec.append(first_res.get("score_0_10", 5.0))

        features["avg_aryl_electronic"] = round(sum(aryl_elec) / len(aryl_elec), 2) if aryl_elec else 5.0
        
        # Chelators
        nearby = analysis.get("nearby", [])
        has_chelator = False
        for n in nearby:
            if not isinstance(n, dict):
                continue
            result = n.get("result", [])
            if not isinstance(result, list):
                continue
            for res in result:
                if isinstance(res, dict):
                    tags = res.get("tags", [])
                    if "bidentate_chelator" in tags or "chelator" in tags:
                        has_chelator = True
                        break
                elif isinstance(res, str):
                    # Fallback for string labels
                    if "chelator" in res.lower():
                        has_chelator = True
                        break
            if has_chelator:
                break
        if has_chelator:
            features["has_chelator"] = True
        
        # Motifs (Sparse: only store True)
        motifs = analysis.get("motifs", [])
        reactant_types = []
        for m in motifs:
            if isinstance(m, dict) and "compound_id" in m:
                cid = m['compound_id']
                features[f"has_{cid}"] = True
                if cid not in reactant_types:
                    reactant_types.append(cid)
        
        if reactant_types:
            features["reactant_type"] = reactant_types
            
        return features

    def _standardize_conditions(self, raw_reagents: List[Dict[str, str]], raw_solvents: List[Dict[str, str]], reaction_data: Dict[str, Any]) -> Dict[str, Any]:
        conditions = {
            "catalyst": [],
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
                conditions["catalyst"].append(name)
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
