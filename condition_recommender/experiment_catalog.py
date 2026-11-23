"""Import source experiment observations for evidence lookup, including failed controls."""

from __future__ import annotations

import argparse
from collections.abc import Iterator, Mapping
import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
from typing import Any

from .conversion.generic import GenericConversionCache, convert_record
from .conversion.input_schema import RawReactionRecord
from .conversion.intermediate import iter_intermediate_records, raw_record_from_intermediate
from .processed_catalog import build_processed_catalog
from .record_storage import json_bytes


def _inputs(path: Path) -> Iterator[RawReactionRecord]:
    if path.suffix.casefold() != ".json":
        yield from iter_intermediate_records(path)
        return
    values = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(values, list):
        raise ValueError("Expected a JSON array of canonical source observations")
    for index, value in enumerate(values, 1):
        if not isinstance(value, Mapping):
            raise ValueError(f"Expected a source observation object at position {index}")
        yield raw_record_from_intermediate(value)


def build_experiment_catalog(source: str | Path, destination: str | Path) -> dict[str, Any]:
    """Create a separate evidence catalog from canonical source-observation JSON or JSONL.

    Standard conversion preserves source outcomes, target graphs, explicit
    absences, stages, reference IDs and regular recipe admission decisions.
    Each submitted experiment remains addressable by its own observation ID.
    The resulting SQLite artifact is independent of a recommendation index.
    """
    source_path = Path(source).resolve()
    output = Path(destination).resolve()
    if source_path == output:
        raise ValueError("Experiment source and catalog destination must differ")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=output.parent) as temporary:
        canonical = Path(temporary) / "records.jsonl"
        cache = GenericConversionCache(max_entries=512)
        with canonical.open("wb") as handle:
            for record in _inputs(source_path):
                if not record.upstream_observation_id.strip():
                    raise ValueError("Experiment source observation requires a nonempty observation_id")
                converted = convert_record(record, cache=cache)
                handle.write(json_bytes(converted.to_dict()) + b"\n")
        if not canonical.stat().st_size:
            raise ValueError("Experiment source contains no observations")
        return build_processed_catalog(canonical, output, include_local_cores=True)


def main(argv: list[str] | None = None) -> int:
    """Build a searchable experiment catalog from a local source observation file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="source_observation.v1 .json array or .jsonl[.gz] file")
    parser.add_argument("destination", type=Path, help="Output experiment evidence SQLite file")
    args = parser.parse_args(argv)
    try:
        result = build_experiment_catalog(args.source, args.destination)
    except (OSError, ValueError, sqlite3.Error) as exc:
        parser.error(str(exc))
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
