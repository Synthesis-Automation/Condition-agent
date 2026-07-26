from condition_recommender.conversion.references import normalize_reference


def test_doi_identity_is_case_and_url_invariant() -> None:
    first = normalize_reference(
        "Journal Name (2024). https://doi.org/10.1021/acs.joc.4c00001"
    )
    second = normalize_reference("DOI: 10.1021/ACS.JOC.4C00001.")

    assert first.reference_id == second.reference_id
    assert first.doi == "10.1021/acs.joc.4c00001"
    assert first.publication_year == 2024
    assert first.resolution_status == "doi"


def test_patent_identity_is_spacing_invariant() -> None:
    first = normalize_reference(
        "World Intellectual Property Organization, WO2024166133 A1 2024-08-15"
    )
    second = normalize_reference("WO 2024166133-A1")

    assert first.reference_id == second.reference_id
    assert first.patent_number == "WO2024166133A1"
    assert first.resolution_status == "patent_number"


def test_bibliographic_text_normalization_is_deterministic() -> None:
    first = normalize_reference(
        "Journal of Organic Chemistry (2022), 87(15), 10285-10297"
    )
    second = normalize_reference(
        "  JOURNAL OF ORGANIC CHEMISTRY   (2022), 87(15), 10285-10297 "
    )

    assert first.reference_id == second.reference_id
    assert first.resolution_status == "bibliographic_text"
    assert first.warnings == ("REFERENCE_ID_FROM_NORMALIZED_TEXT",)


def test_missing_reference_is_not_assigned_a_shared_identity() -> None:
    identity = normalize_reference("")

    assert identity.reference_id == ""
    assert identity.resolution_status == "missing"
    assert identity.warnings == ("MISSING_REFERENCE",)
