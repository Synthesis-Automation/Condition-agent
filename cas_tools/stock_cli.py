"""Command-line tools for downloading and compiling stock portfolios."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Any, Sequence

from .stock_portfolio import (
    STOCK_SOURCE_MANIFEST_VERSION,
    StockPortfolio,
    build_stock_portfolio,
    load_stock_source_manifest,
)


MCULE_DATABASE_FILES_URL = "https://mcule.com/api/v1/database-files/"
MCULE_DATABASE_PAGE_URL = "https://mcule.com/database/"
MCULE_TERMS_URL = "https://mcule.com/terms/"


def _download_json(url: str) -> dict[str, Any]:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_mcule_collection(
    output_directory: str | Path,
    *,
    collection: str = "Mcule In Stock",
    file_type: str = "smi.gz",
) -> tuple[Path, Path]:
    """Download and checksum one public Mcule database snapshot and manifest."""

    payload = _download_json(MCULE_DATABASE_FILES_URL)
    matches = [
        item
        for item in payload.get("results") or ()
        if item.get("name") == collection
    ]
    if len(matches) != 1:
        raise ValueError(f"Mcule collection was not uniquely found: {collection}")
    source = matches[0]
    files = [
        item for item in source.get("files") or () if item.get("file_type") == file_type
    ]
    if len(files) != 1:
        raise ValueError(
            f"Mcule collection does not have one {file_type} file: {collection}"
        )
    remote = files[0]
    output_root = Path(output_directory)
    output_root.mkdir(parents=True, exist_ok=True)
    destination = output_root / str(remote["filename"])
    partial = destination.with_suffix(destination.suffix + ".part")
    if destination.is_file() and _sha256(destination) == remote["sha256_checksum"]:
        pass
    else:
        if partial.exists():
            partial.unlink()
        request = urllib.request.Request(
            str(remote["download_url"]),
            headers={"User-Agent": "Condition-agent stock portfolio sync"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            with partial.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
        actual_checksum = _sha256(partial)
        if actual_checksum != remote["sha256_checksum"]:
            partial.unlink()
            raise ValueError("downloaded Mcule checksum does not match API metadata")
        partial.replace(destination)
    manifest = {
        "schema_version": STOCK_SOURCE_MANIFEST_VERSION,
        "sources": [
            {
                "path": str(destination.resolve()),
                "supplier": "Mcule",
                "collection": str(source["name"]),
                "snapshot_date": str(source["last_updated"]),
                "availability_status": "in_stock",
                "evidence_level": "supplier_in_stock",
                "terminal_eligible": True,
                "source_url": MCULE_DATABASE_PAGE_URL,
                "terms_url": MCULE_TERMS_URL,
                "region": "global",
                "format": "smi",
                "delimiter": "\t",
                "smiles_column": "smiles",
                "identifier_column": "mcule_id",
            }
        ],
        "download": {
            "api_url": MCULE_DATABASE_FILES_URL,
            "download_url": remote["download_url"],
            "sha256": remote["sha256_checksum"],
            "entry_count": source["entry_count"],
        },
    }
    manifest_path = output_root / "stock_sources.v1.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return destination, manifest_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    download = commands.add_parser("download-mcule")
    download.add_argument("output_directory")
    download.add_argument("--collection", default="Mcule In Stock")
    download.add_argument("--file-type", default="smi.gz")
    build = commands.add_parser("build")
    build.add_argument("manifest")
    build.add_argument("output")
    build.add_argument("--workers", type=int, default=1)
    build.add_argument("--chunk-size", type=int, default=2_000)
    inspect = commands.add_parser("inspect")
    inspect.add_argument("portfolio")
    inspect.add_argument("smiles", nargs="?")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run stock synchronization, compilation, or inspection."""

    arguments = _parser().parse_args(argv)
    if arguments.command == "download-mcule":
        source, manifest = download_mcule_collection(
            arguments.output_directory,
            collection=arguments.collection,
            file_type=arguments.file_type,
        )
        print(json.dumps({"source": str(source), "manifest": str(manifest)}))
        return 0
    if arguments.command == "build":
        last_reported: dict[str, int] = {}

        def progress(collection: str, rows: int, accepted: int) -> None:
            previous = last_reported.get(collection, 0)
            if rows - previous >= 100_000:
                print(
                    json.dumps(
                        {
                            "collection": collection,
                            "rows": rows,
                            "accepted": accepted,
                        }
                    ),
                    flush=True,
                )
                last_reported[collection] = rows

        report = build_stock_portfolio(
            load_stock_source_manifest(arguments.manifest),
            arguments.output,
            workers=arguments.workers,
            chunk_size=arguments.chunk_size,
            progress=progress,
        )
        print(json.dumps(report.to_dict(), indent=2, sort_keys=True))
        return 0
    with StockPortfolio(arguments.portfolio) as portfolio:
        if arguments.smiles:
            match = portfolio.lookup(arguments.smiles)
            print(
                json.dumps(
                    match.__dict__ if match is not None else None,
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            print(json.dumps(portfolio.metadata(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

