
import json
from typing import Dict, Any, List, Optional

class PairingHelper:
    """
    Minimal helper that loads the catalyst/precursor taxonomy and, when a catalyst
    requires a ligand but none is provided, proposes a ligand+conditions pair.
    """

    def __init__(self, taxonomy_path: str):
        with open(taxonomy_path, "r", encoding="utf-8") as f:
            self.tax = json.load(f)
        self._index_families()

    def _index_families(self):
        self.by_family = {f.get("family_id"): f for f in self.tax.get("families", [])}
        self.by_member = {}
        for f in self.tax.get("families", []):
            for em in f.get("example_members", []):
                key = (f.get("family_id"), (em.get("abbr") or em.get("name")))
                self.by_member[key] = (f, em)

    def suggest_for(self, family_id: str, abbr: Optional[str] = None, reaction_hint: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        Return a suggestion dict:
          {
            "ligand": {"family_id": "...", "abbr": "..."},
            "loadings": {"catalyst_mol%": "...", "ligand_mol%": "...", "ratio_L_to_M": 1.2},
            "co_reagents": {"base": [...], "reductant": [...]} or {},
            "solvents": [...],
            "temp_C": [low, high],
            "notes": "...",
            "source": "family" | "member"
          }
        """
        fam = self.by_family.get(family_id)
        if not fam:
            return None

        # Determine if additional ligand is required
        addlig = fam.get("additional_ligand", {})
        usage = (addlig.get("usage") or "optional").lower()

        # Member-level override check
        member = None
        if abbr:
            member = self.by_member.get((family_id, abbr), (None, None))[1]
            if member and member.get("additional_ligand"):
                usage = (member["additional_ligand"].get("usage") or usage).lower()

        if usage == "not_needed":
            return None

        # Candidate pairings
        pairs: List[Dict[str, Any]] = []
        if member and member.get("recommended_pairs"):
            pairs.extend(member["recommended_pairs"])
            source = "member"
        if fam.get("recommended_pairs"):
            pairs.extend(fam["recommended_pairs"])
            source = "family"
        if not pairs:
            return None

        # Simple selection heuristic: filter by reaction_hint when provided
        def match_score(p):
            if reaction_hint and p.get("for_reactions"):
                # score by keyword overlap (case-insensitive substring)
                rh = reaction_hint.lower()
                return max((1 if rh in r.lower() else 0) for r in p["for_reactions"])
            return 0

        pairs = sorted(pairs, key=match_score, reverse=True)
        pick = pairs[0]

        # Pick a ligand abbr example if present
        ligand = {"family_id": None, "abbr": None}
        with_lig = pick.get("with_ligand") or {}
        ligand["family_id"] = with_lig.get("family_id")
        ligand_examples = with_lig.get("examples") or []
        ligand["abbr"] = ligand_examples[0] if ligand_examples else None

        return {
            "ligand": ligand,
            "loadings": pick.get("loadings") or {},
            "co_reagents": pick.get("co_reagents") or {},
            "solvents": pick.get("solvents") or [],
            "temp_C": pick.get("temp_C") or [],
            "notes": pick.get("notes") or "",
            "source": source
        }

# Example:
# helper = PairingHelper("taxonomy_catalysts_precursors_with_pairs.json")
# suggestion = helper.suggest_for("pd_ii_salts", abbr="Pd(OAc)2", reaction_hint="Buchwald–Hartwig")
# print(suggestion)
