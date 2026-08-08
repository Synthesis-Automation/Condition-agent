import csv
import json
from pathlib import Path

from cas_tools.compound_lookup import CompoundLookupResult
from cas_tools.registry_cas_reconciliation import (
    _structure_element_contradiction,
    reconcile_registry_from_cas_csv,
)
from condition_registry import ConditionRegistry
INPUT_FIELDS = (
    "name", "aliases", "roles", "mention_count", "source_columns",
    "reconciliation_status", "normalized_identity", "possible_cas_no",
    "cas_match_status", "cas_confidence", "source_basis", "source_url",
    "match_notes", "cas_checksum_valid",
)


def _write(path: Path, fields, rows=()):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _write_registry(path: Path, records=()):
    path.write_text(
        "".join(json.dumps(record) + "\n" for record in records),
        encoding="utf-8",
    )


def test_reconciliation_adds_existing_alias_and_verified_new_compound(tmp_path):
    substances = tmp_path / "substances.v2.jsonl"
    source = tmp_path / "possible.csv"
    output = tmp_path / "output"
    _write_registry(
        substances,
        ({
            "id": "cas:7732-18-5",
            "name": "Water",
            "cas": "7732-18-5",
            "smiles": "O",
            "aliases": [{"type": "abbreviation", "value": "H2O"}],
            "roles": ["solvent"],
        },),
    )
    common = {
        "roles": "solvent", "mention_count": "2", "source_columns": "Solvent",
        "reconciliation_status": "unresolved", "cas_match_status": "exact_name_match",
        "cas_confidence": "high", "source_basis": "curator check",
        "source_url": "https://example.test", "match_notes": "", "cas_checksum_valid": "yes",
    }
    _write(
        source,
        INPUT_FIELDS,
        (
            {**common, "name": "Aqua reagent", "aliases": "Aqua reagent", "normalized_identity": "Water", "possible_cas_no": "7732-18-5"},
            {**common, "name": "EtOH reagent", "aliases": "EtOH reagent", "normalized_identity": "Ethanol", "possible_cas_no": "64-17-5"},
            {**common, "name": "TCE", "aliases": "TCE", "normalized_identity": "Trichloroethylene or tetrachloroethane", "possible_cas_no": "79-01-6; 79-34-5", "cas_match_status": "ambiguous_multiple_candidates", "cas_confidence": "medium"},
        ),
    )

    def lookup(cas):
        values = {
            "7732-18-5": CompoundLookupResult(cas=cas, status="resolved", canonical_name="Water", smiles="O", formula="H2O", molecular_weight=18.015),
            "64-17-5": CompoundLookupResult(cas=cas, status="resolved", canonical_name="Ethanol", iupac_name="ethanol", smiles="CCO", formula="C2H6O", molecular_weight=46.07, source_ids=("pubchem:cid:702",)),
        }
        return values[cas]

    summary = reconcile_registry_from_cas_csv(
        source,
        output,
        apply_changes=True,
        lookup=lookup,
        lookup_workers=1,
        substances_path=substances,
    )

    assert summary.aliases_added == 1
    assert summary.new_compounds_added == 1
    assert summary.review_rows == 1
    registry = ConditionRegistry(substances_path=substances)
    assert registry.resolve(name="Aqua reagent").substance.substance_id == "cas:7732-18-5"
    assert registry.resolve(name="EtOH reagent").substance.substance_id == "cas:64-17-5"
    with Path(summary.audit_path).open(encoding="utf-8-sig", newline="") as handle:
        decisions = {row["name"]: row["decision"] for row in csv.DictReader(handle)}
    assert decisions == {
        "Aqua reagent": "add_aliases",
        "EtOH reagent": "add_compound",
        "TCE": "review",
    }


def test_shared_cas_with_different_normalized_identities_is_held(tmp_path):
    substances = tmp_path / "substances.v2.jsonl"
    source = tmp_path / "possible.csv"
    _write_registry(substances)
    base = {
        "aliases": "", "roles": "additive", "mention_count": "1",
        "source_columns": "Additive", "reconciliation_status": "unresolved",
        "possible_cas_no": "9002-92-0", "cas_match_status": "commercial_mixture_or_polymer_match",
        "cas_confidence": "medium", "source_basis": "supplier", "source_url": "https://example.test",
        "match_notes": "shared polymer CAS", "cas_checksum_valid": "yes",
    }
    _write(
        source,
        INPUT_FIELDS,
        (
            {**base, "name": "Brij 30", "normalized_identity": "Polyoxyethylene(4) lauryl ether"},
            {**base, "name": "Brij 35", "normalized_identity": "Polyoxyethylene(23) lauryl ether"},
        ),
    )

    result = CompoundLookupResult(
        cas="9002-92-0", status="partial", canonical_name="Laureth", synonyms=("Brij",),
    )
    summary = reconcile_registry_from_cas_csv(
        source,
        tmp_path / "out",
        lookup=lambda _cas: result,
        lookup_workers=1,
        substances_path=substances,
    )

    assert summary.review_rows == 2
    assert summary.new_compounds_added == 0


def test_structure_element_contradiction_detects_wrong_metal_identity():
    assert _structure_element_contradiction(
        ("AlPhos Pd G6 bromide",),
        "C52H67F4OP",
    ) == "Pd"
    assert _structure_element_contradiction(
        ("Ni(COD)(DQ)",),
        "C18H24NiO2",
    ) is None
