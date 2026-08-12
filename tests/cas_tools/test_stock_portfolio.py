"""Supplier stock portfolio regressions."""

from __future__ import annotations

import csv
import json

import pytest

from cas_tools import (
    StockPortfolio,
    StockSourceDefinition,
    build_stock_portfolio,
    load_stock_source_manifest,
    open_stock_lookup,
)


def _write_smi(path, rows) -> None:
    path.write_text(
        "".join(f"{smiles}\t{record_id}\n" for smiles, record_id in rows),
        encoding="utf-8",
    )


def test_portfolio_merges_supplier_offers_and_preserves_provenance(tmp_path) -> None:
    mcule = tmp_path / "mcule.smi"
    enamine = tmp_path / "enamine.csv"
    _write_smi(mcule, (("CCO", "MCULE-1"), ("CCN", "MCULE-2")))
    with enamine.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("SMILES", "ID"))
        writer.writeheader()
        writer.writerow({"SMILES": "OCC", "ID": "EN300-1"})
        writer.writerow({"SMILES": "not-smiles", "ID": "EN-invalid"})
    sources = (
        StockSourceDefinition(
            path=str(mcule),
            supplier="Mcule",
            collection="In Stock",
            snapshot_date="2026-08-02",
            availability_status="in_stock",
            evidence_level="supplier_in_stock",
            terminal_eligible=True,
            source_url="https://mcule.example/source",
            terms_url="https://mcule.example/terms",
        ),
        StockSourceDefinition(
            path=str(enamine),
            supplier="Enamine",
            collection="Authorized Global Stock",
            snapshot_date="2026-07-01",
            availability_status="in_stock",
            evidence_level="supplier_in_stock",
            terminal_eligible=True,
            source_url="https://enamine.example/source",
            terms_url="https://enamine.example/terms",
            format="csv",
            delimiter=",",
            smiles_column="SMILES",
            identifier_column="ID",
        ),
    )
    output = tmp_path / "stock.sqlite"

    report = build_stock_portfolio(sources, output)

    assert report.source_rows == 4
    assert report.accepted_rows == 3
    assert report.invalid_structure_rows == 1
    assert report.unique_molecules == 2
    assert report.offer_count == 3
    with StockPortfolio(output) as portfolio:
        ethanol = portfolio.lookup("C(C)O")
        assert ethanol is not None
        assert ethanol.occurrence_count == 2
        assert {row["supplier"] for row in ethanol.source_records} == {
            "Mcule",
            "Enamine",
        }
        assert all(
            row["source_role"] == "starting_material"
            for row in ethanol.source_records
        )
    with open_stock_lookup(output) as portfolio:
        assert isinstance(portfolio, StockPortfolio)


def test_catalog_only_source_cannot_be_route_terminal(tmp_path) -> None:
    catalog = tmp_path / "catalog.smi"
    _write_smi(catalog, (("CCO", "catalog-1"),))
    source = StockSourceDefinition(
        path=str(catalog),
        supplier="Example",
        collection="Comprehensive Catalog",
        snapshot_date="2026-01-01",
        availability_status="catalog_listed",
        evidence_level="catalog_listed",
        terminal_eligible=False,
        source_url="https://example.test/source",
        terms_url="https://example.test/terms",
    )
    output = tmp_path / "stock.sqlite"

    build_stock_portfolio((source,), output)

    with StockPortfolio(output) as portfolio:
        match = portfolio.lookup("CCO")
        assert match is not None
        assert match.source_records[0]["source_role"] == "catalog_listing"
        assert match.source_records[0]["terminal_eligible"] == "false"


def test_terminal_eligibility_requires_strong_stock_evidence(tmp_path) -> None:
    source = tmp_path / "catalog.smi"
    _write_smi(source, (("CCO", "catalog-1"),))
    definition = StockSourceDefinition(
        path=str(source),
        supplier="Example",
        collection="Catalog",
        snapshot_date="2026-01-01",
        availability_status="catalog_listed",
        evidence_level="catalog_listed",
        terminal_eligible=True,
        source_url="https://example.test/source",
        terms_url="https://example.test/terms",
    )

    with pytest.raises(ValueError, match="only physical or supplier-in-stock"):
        build_stock_portfolio((definition,), tmp_path / "stock.sqlite")


def test_manifest_is_versioned_and_requires_local_sources(tmp_path) -> None:
    source = tmp_path / "stock.smi"
    _write_smi(source, (("CCO", "one"),))
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "sources": [
                    {
                        "path": str(source),
                        "supplier": "Example",
                        "collection": "Stock",
                        "snapshot_date": "2026-01-01",
                        "availability_status": "in_stock",
                        "evidence_level": "supplier_in_stock",
                        "terminal_eligible": True,
                        "source_url": "https://example.test/source",
                        "terms_url": "https://example.test/terms",
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    loaded = load_stock_source_manifest(manifest)

    assert loaded[0].supplier == "Example"
