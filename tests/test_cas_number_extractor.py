from __future__ import annotations

import zipfile
from pathlib import Path

from cas_tools.cas_number_extractor import (
    CASMatch,
    discover_candidate_files,
    extract_cas_matches_from_file,
    find_cas_numbers_in_text,
    is_valid_cas_number,
    write_matches_to_csv,
    write_matches_to_markdown,
)


def _write_minimal_docx(path: Path, text: str) -> None:
    document_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<w:document xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">"
        "<w:body><w:p><w:r><w:t>"
        f"{text}"
        "</w:t></w:r></w:p></w:body></w:document>"
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("word/document.xml", document_xml)


def _write_minimal_xlsx(path: Path, text: str) -> None:
    workbook_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
        "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
        "<sheets><sheet name=\"Sheet1\" sheetId=\"1\" r:id=\"rId1\"/></sheets>"
        "</workbook>"
    )
    workbook_rels = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
        "<Relationship Id=\"rId1\" "
        "Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" "
        "Target=\"worksheets/sheet1.xml\"/>"
        "</Relationships>"
    )
    sheet_xml = (
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
        "<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">"
        "<sheetData><row r=\"1\"><c r=\"A1\" t=\"inlineStr\"><is><t>"
        f"{text}"
        "</t></is></c></row></sheetData></worksheet>"
    )
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("xl/workbook.xml", workbook_xml)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        archive.writestr("xl/worksheets/sheet1.xml", sheet_xml)


def test_is_valid_cas_number_uses_checksum() -> None:
    assert is_valid_cas_number("7732-18-5")
    assert is_valid_cas_number("64-17-5")
    assert not is_valid_cas_number("7732-18-4")
    assert not is_valid_cas_number("12-34-5")


def test_find_cas_numbers_in_text_filters_invalid_candidates() -> None:
    text = "Water 7732-18-5, ethanol 64-17-5, and a bad number 7732-18-4."
    assert find_cas_numbers_in_text(text) == ["7732-18-5", "64-17-5"]


def test_discover_candidate_files_accepts_text_like_unknown_extensions(tmp_path: Path) -> None:
    text_file = tmp_path / "notes.custom"
    text_file.write_text("contains 7732-18-5\n", encoding="utf-8")
    binary_file = tmp_path / "image.bin"
    binary_file.write_bytes(b"\x00\x01\x02\x03")

    candidates = discover_candidate_files(tmp_path)

    assert text_file in candidates
    assert binary_file not in candidates


def test_extract_cas_matches_from_csv_reports_cell_location(tmp_path: Path) -> None:
    csv_path = tmp_path / "table.csv"
    csv_path.write_text("name,cas\nwater,7732-18-5\n", encoding="utf-8")

    matches, warnings = extract_cas_matches_from_file(csv_path, base_folder=tmp_path)

    assert not warnings
    assert [match.cas_number for match in matches] == ["7732-18-5"]
    assert matches[0].location == "row 2, col 2"


def test_extract_cas_matches_from_docx_zip_fallback(tmp_path: Path) -> None:
    docx_path = tmp_path / "sample.docx"
    _write_minimal_docx(docx_path, "Sulfuric acid 7664-93-9")

    matches, warnings = extract_cas_matches_from_file(docx_path, base_folder=tmp_path)

    assert not warnings
    assert [match.cas_number for match in matches] == ["7664-93-9"]
    assert matches[0].file_type == "docx"


def test_extract_cas_matches_from_xlsx_zip_fallback(tmp_path: Path) -> None:
    xlsx_path = tmp_path / "sample.xlsx"
    _write_minimal_xlsx(xlsx_path, "Ethanol 64-17-5")

    matches, warnings = extract_cas_matches_from_file(xlsx_path, base_folder=tmp_path)

    assert not warnings
    assert [match.cas_number for match in matches] == ["64-17-5"]
    assert matches[0].location == "sheet Sheet1!A1"


def test_write_matches_to_markdown_groups_results(tmp_path: Path) -> None:
    output_path = tmp_path / "report.md"
    matches = [
        CASMatch(
            source_file="sample.txt",
            relative_path="nested/sample.txt",
            file_type="text",
            location="line 1",
            cas_number="7732-18-5",
            context="Water CAS 7732-18-5",
        ),
        CASMatch(
            source_file="sample.txt",
            relative_path="nested/sample.txt",
            file_type="text",
            location="line 2",
            cas_number="64-17-5",
            context="Ethanol CAS 64-17-5",
        ),
    ]

    write_matches_to_markdown(
        matches,
        output_path,
        source_folder="demo",
        processed_files=["nested/sample.txt"],
        warnings=["Skipped bad.xls"],
    )

    content = output_path.read_text(encoding="utf-8")

    assert "# CAS Number Extraction Report (demo)" in content
    assert "Total files processed: 1" in content
    assert "Unique CAS numbers: 2" in content
    assert "### nested/sample.txt" in content
    assert "- line 1: 7732-18-5 | Water CAS 7732-18-5" in content
    assert "Skipped bad.xls" in content


def test_write_matches_to_csv_appends_only_new_unique_cas_numbers(tmp_path: Path) -> None:
    output_path = tmp_path / "cas_no_all.csv"
    output_path.write_text("cas_number\n7732-18-5\n7732-18-5\n", encoding="utf-8")

    matches = [
        CASMatch(
            source_file="sample_a.txt",
            relative_path="sample_a.txt",
            file_type="text",
            location="line 1",
            cas_number="7732-18-5",
            context="Water CAS 7732-18-5",
        ),
        CASMatch(
            source_file="sample_b.txt",
            relative_path="sample_b.txt",
            file_type="text",
            location="line 2",
            cas_number="64-17-5",
            context="Ethanol CAS 64-17-5",
        ),
        CASMatch(
            source_file="sample_c.txt",
            relative_path="sample_c.txt",
            file_type="text",
            location="line 3",
            cas_number="64-17-5",
            context="Duplicate ethanol CAS 64-17-5",
        ),
    ]

    write_matches_to_csv(matches, output_path)

    content = output_path.read_text(encoding="utf-8").splitlines()

    assert content == ["cas_number", "7732-18-5", "64-17-5"]