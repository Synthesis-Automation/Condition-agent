import json
import os
from typing import List, Dict, Any, Set, Tuple
from rdkit import Chem

class DynamicClassifier:
    """
    PoC for a score-based chemical motif classifier.
    Uses 'priority' from organic_groups to determine the most specific classification.
    """
    def __init__(self, data_dir: str):
        self.data_dir = data_dir
        self.groups = self._load_json("organic_groups.v1.3.json")["groups"]
        self.compounds = self._load_json("organic_compounds.v1.3.json")["compounds"]
        self.logic = self._load_json("group_logic.json")["group_sets"]
        
        # Map group ID to priority
        self.group_priorities = {g["id"]: g.get("priority", 1) for g in self.groups}
        
        # Add logic set priorities
        for set_id, set_data in self.logic.items():
            if "priority" in set_data:
                self.group_priorities[set_id] = set_data["priority"]
            else:
                # Fallback: min priority of members
                member_priorities = [self.group_priorities.get(m, 1) for m in set_data["members"]]
                if member_priorities:
                    self.group_priorities[set_id] = min(member_priorities)
                else:
                    self.group_priorities[set_id] = 1
        
        self.group_smarts = {g["id"]: g["smarts"] for g in self.groups}
        
    def _load_json(self, filename: str) -> Dict:
        path = os.path.join(self.data_dir, filename)
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)

    def get_compound_score(self, compound: Dict) -> int:
        """Calculate specificity score: Priority(A) + Priority(B)"""
        a_id = compound.get("A")
        b_id = compound.get("B")
        score = self.group_priorities.get(a_id, 1) + self.group_priorities.get(b_id, 1)
        return score

    def classify(self, mol: Chem.Mol) -> List[Dict]:
        """
        Classifies a molecule and returns only the most specific matches per bond.
        """
        all_matches = []
        
        # In a real system, we would use the connection engine.
        # For this PoC, we'll simulate the 'Subsumption' logic.
        
        # 1. Find all matches (Simulated for PoC)
        # In reality, this would run the SMARTS or connection logic.
        # Here we just demonstrate the filtering logic.
        
        # Let's say we found these matches for the same bond (atoms 5 and 6)
        raw_matches = [
            {"id": "Any-X", "atoms": (5, 6), "score": self.get_compound_score({"A": "Any_Scaffold", "B": "X"})},
            {"id": "Csp2-X", "atoms": (5, 6), "score": self.get_compound_score({"A": "sp2_Carbon", "B": "X"})},
            {"id": "Ar-Cl", "atoms": (5, 6), "score": self.get_compound_score({"A": "Ar", "B": "Cl"})},
            {"id": "Ar-X", "atoms": (5, 6), "score": self.get_compound_score({"A": "Ar", "B": "X"})},
        ]
        
        # 2. Group by atom indices (the 'site')
        sites = {}
        for match in raw_matches:
            site = tuple(sorted(match["atoms"]))
            if site not in sites:
                sites[site] = []
            sites[site].append(match)
            
        # 3. For each site, keep only the highest score(s)
        final_results = []
        for site, matches in sites.items():
            max_score = max(m["score"] for m in matches)
            best_matches = [m for m in matches if m["score"] == max_score]
            final_results.extend(best_matches)
            
        return final_results

if __name__ == "__main__":
    # Quick test
    data_path = "chemtools/taxonomy/data"
    classifier = DynamicClassifier(data_path)
    
    print("Specificity Scores Examples:")
    test_cases = [
        {"A": "Any_Scaffold", "B": "X"},      # Any-X
        {"id": "Ar-X", "A": "Ar", "B": "X"},  # Ar-X
        {"id": "Ar-Cl", "A": "Ar", "B": "Cl"}, # Ar-Cl
        {"id": "Bn-Cl", "A": "Bn", "B": "Cl"}, # Bn-Cl
    ]
    
    for tc in test_cases:
        name = tc.get("id", f"{tc['A']}-{tc['B']}")
        score = classifier.get_compound_score(tc)
        print(f"  {name:10} -> Score: {score}")
