import json
from urllib.error import URLError

from cas_tools.compound_lookup import CompoundLookupClient


def _section(heading, values):
    information = []
    for value in values:
        if isinstance(value, tuple):
            number, unit = value
            information.append({"Value": {"Number": [number], "Unit": unit}})
        else:
            information.append(
                {"Value": {"StringWithMarkup": [{"String": value}]}}
            )
    return {
        "Record": {
            "Section": [
                {
                    "TOCHeading": "Chemical and Physical Properties",
                    "Section": [
                        {
                            "TOCHeading": heading,
                            "Information": information,
                        }
                    ],
                }
            ]
        }
    }


def _pubchem_fetch(url: str, _timeout: float) -> bytes:
    if "/property/" in url:
        payload = {
            "PropertyTable": {
                "Properties": [
                    {
                        "CID": 702,
                        "Title": "Ethanol",
                        "IUPACName": "ethanol",
                        "MolecularFormula": "C2H6O",
                        "MolecularWeight": "46.07",
                        "SMILES": "CCO",
                        "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                    }
                ]
            }
        }
    elif "/synonyms/" in url:
        payload = {
            "InformationList": {
                "Information": [
                    {
                        "Synonym": [
                            "ethanol",
                            "ethyl alcohol",
                            "64-17-5",
                            "EtOH",
                            "CHEBI:16236",
                            "grain alcohol",
                        ]
                    }
                ]
            }
        }
    elif "heading=Density" in url:
        payload = _section("Density", ["0.7893 g/cu cm at 20 °C"])
    elif "heading=Boiling%20Point" in url:
        payload = _section("Boiling Point", [(78.2, "°C")])
    elif "heading=Melting%20Point" in url:
        payload = _section("Melting Point", ["-114.1 °C"])
    elif "heading=Physical%20Description" in url:
        payload = _section(
            "Physical Description",
            ["Clear, colorless liquid with a characteristic odor."],
        )
    else:
        raise AssertionError(f"Unexpected URL: {url}")
    return json.dumps(payload).encode("utf-8")


def test_pubchem_lookup_collects_identity_physical_data_and_synonyms() -> None:
    result = CompoundLookupClient(
        fetch_bytes=_pubchem_fetch,
        timeout_seconds=3,
    ).lookup("64-17-5")

    assert result.status == "resolved"
    assert result.canonical_name == "Ethanol"
    assert result.smiles == "CCO"
    assert result.formula == "C2H6O"
    assert result.molecular_weight == 46.07
    assert result.inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert result.abbreviation == "EtOH"
    assert result.synonyms == ("ethyl alcohol", "EtOH", "grain alcohol")
    assert result.density == 0.7893
    assert result.boiling_point_c == 78.2
    assert result.melting_point_c == -114.1
    assert result.substance_kind == "liquid"
    assert result.source_ids == ("pubchem:cid:702",)
    assert result.warnings == ("PUBCHEM_EXPERIMENTAL_FIELDS_REQUIRE_REVIEW",)


def test_invalid_cas_does_not_make_network_request() -> None:
    calls = []

    def fetch(url, _timeout):
        calls.append(url)
        raise AssertionError("network should not be called")

    result = CompoundLookupClient(fetch_bytes=fetch).lookup("64-17-4")

    assert result.status == "invalid_identifier"
    assert result.warnings == ("INVALID_CAS",)
    assert calls == []


def test_pubchem_core_result_survives_optional_endpoint_failures() -> None:
    def fetch(url, timeout):
        if "/property/" in url:
            return _pubchem_fetch(url, timeout)
        raise URLError("optional endpoint unavailable")

    result = CompoundLookupClient(fetch_bytes=fetch).lookup("64-17-5")

    assert result.status == "resolved"
    assert result.canonical_name == "Ethanol"
    assert result.synonyms == ()
    assert result.density is None
    assert any("PUBCHEM_SYNONYM_LOOKUP_FAILED" in item for item in result.warnings)


def test_nci_fallback_supplies_core_fields_when_pubchem_is_unavailable() -> None:
    responses = {
        "/iupac_name": "methanol",
        "/smiles": "CO",
        "/formula": "CH4O",
        "/mw": "32.04",
        "/stdinchiKey": "OKKJLVBELUTLKV-UHFFFAOYSA-N",
        "/names": "Methanol\nMethyl alcohol\nMeOH\n",
    }

    def fetch(url, _timeout):
        if "pubchem" in url:
            raise URLError("PubChem unavailable")
        for suffix, value in responses.items():
            if url.endswith(suffix):
                return value.encode("utf-8")
        raise AssertionError(f"Unexpected URL: {url}")

    result = CompoundLookupClient(fetch_bytes=fetch).lookup("67-56-1")

    assert result.status == "resolved"
    assert result.canonical_name == "Methanol"
    assert result.smiles == "CO"
    assert result.formula == "CH4O"
    assert result.molecular_weight == 32.04
    assert result.abbreviation == "MeOH"
    assert result.source_ids == ("nci:cactus",)
    assert any("PUBCHEM_PROPERTY_LOOKUP_FAILED" in item for item in result.warnings)
