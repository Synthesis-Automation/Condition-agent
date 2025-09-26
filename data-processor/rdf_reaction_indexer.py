#!/usr/bin/env python3
"""Utilities to reorganize SciFinder RDF exports for fast lookup by CAS reaction number.

Given a large concatenated RDF file, this tool can build:
- A JSONL index mapping each CAS reaction number to the byte range
  of its reaction block inside the source RDF file.
- An optional JSONL dump of reaction blocks with metadata, making it
  easier to query with tools like jq or ripgrep.

Examples
--------
    # Build both metadata dump and CAS index
    python rdf_reaction_indexer.py dataset_combined.rdf --jsonl reactions.jsonl --index cas_index.jsonl

    # Create only the CAS index
    python rdf_reaction_indexer.py dataset_combined.rdf --index cas_index.jsonl

    # Dry run for the first 10 reactions
    python rdf_reaction_indexer.py dataset_combined.rdf --jsonl preview.jsonl --limit 10
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from dataclasses import dataclass
from typing import IO, Iterable, Iterator, List, Optional, Tuple

CAS_VARIATION_RE = re.compile(r"VAR\((\d+)\)")


@dataclass
class ReactionBlock:
    """Container for a single `$RFMT` block extracted from the RDF."""

    sequence: int
    offset: int
    end_offset: int
    scheme: Optional[str]
    step: Optional[str]
    lines: List[str]
    cas_entries: List[dict]

    @property
    def length(self) -> int:
        return self.end_offset - self.offset

    @property
    def cas_numbers(self) -> List[str]:
        return [entry["cas"] for entry in self.cas_entries]


def scan_reactions(
    path: str,
    *,
    encoding: str = "latin-1",
    limit: Optional[int] = None,
) -> Iterator[ReactionBlock]:
    """Yield successive reaction blocks from a combined RDF file."""

    sequence = 0
    block_start: Optional[int] = None
    block_lines: List[str] = []
    cas_pairs: List[tuple[Optional[int], str]] = []
    scheme: Optional[str] = None
    step: Optional[str] = None
    pending_variation: Optional[int] = None
    limit_hit = False

    def emit(next_offset: int) -> Optional[ReactionBlock]:
        nonlocal sequence, block_start, block_lines, cas_pairs
        nonlocal scheme, step, pending_variation, limit_hit

        if block_start is None or not block_lines:
            return None

        sequence += 1
        block = ReactionBlock(
            sequence=sequence,
            offset=block_start,
            end_offset=next_offset,
            scheme=scheme,
            step=step,
            lines=list(block_lines),
            cas_entries=[{"variation": var, "cas": cas} for (var, cas) in cas_pairs],
        )

        # Reset holders for the next block
        block_start = None
        block_lines = []
        cas_pairs = []
        scheme = None
        step = None
        pending_variation = None

        if limit is not None and sequence >= limit:
            limit_hit = True

        return block

    with open(path, "rb") as handle:
        while True:
            line_offset = handle.tell()
            raw = handle.readline()
            if not raw:
                block = emit(line_offset)
                if block is not None:
                    yield block
                break

            try:
                text = raw.decode(encoding, errors="replace")
            except Exception:
                text = raw.decode("latin-1", errors="replace")

            if text.startswith("$RFMT"):
                block = emit(line_offset)
                if block is not None:
                    yield block
                    if limit_hit:
                        return

                block_start = line_offset
                block_lines = [text]

                tokens = text.strip().split()
                scheme = next((tok for tok in tokens if tok.startswith("SCHEME")), None)
                step = next((tok for tok in tokens if tok.startswith("STEP")), None)
                continue

            if block_start is None:
                # Skip any preamble before the first $RFMT block
                continue

            block_lines.append(text)

            if text.startswith("$DTYPE") and "CAS_Reaction_Number" in text:
                match = CAS_VARIATION_RE.search(text)
                pending_variation = int(match.group(1)) if match else None
                continue

            if pending_variation is not None and text.startswith("$DATUM"):
                value = text[len("$DATUM"):].strip()
                if value:
                    cas_pairs.append((pending_variation, value))
                pending_variation = None
                continue

        # In case the limit was reached exactly at EOF
        if limit_hit:
            return


def write_outputs(
    blocks: Iterable[ReactionBlock],
    *,
    source_path: str,
    jsonl_path: Optional[str],
    index_path: Optional[str],
    metadata_only: bool,
    quiet: bool,
) -> Tuple[int, int]:
    """Write JSONL dump and/or CAS index from the supplied blocks.

    Returns (reaction_count, cas_entry_count).
    """

    abs_source = os.path.abspath(source_path)
    jsonl_file: Optional[IO[str]] = None
    index_file: Optional[IO[str]] = None

    block_count = 0
    index_count = 0

    try:
        if jsonl_path:
            jsonl_file = open(jsonl_path, "w", encoding="utf-8")
        if index_path:
            index_file = open(index_path, "w", encoding="utf-8")

        for block in blocks:
            block_count += 1
            base_record = {
                "sequence": block.sequence,
                "scheme": block.scheme,
                "step": block.step,
                "offset": block.offset,
                "length": block.length,
                "source_file": abs_source,
            }

            if jsonl_file is not None:
                record = dict(base_record)
                record["cas_numbers"] = block.cas_numbers
                record["variations"] = block.cas_entries
                if not metadata_only:
                    record["raw_rdf"] = "".join(block.lines)
                jsonl_file.write(json.dumps(record, ensure_ascii=False) + "\n")

            if index_file is not None:
                for entry in block.cas_entries:
                    index_record = dict(base_record)
                    index_record["cas"] = entry["cas"]
                    index_record["variation"] = entry["variation"]
                    index_file.write(json.dumps(index_record, ensure_ascii=False) + "\n")
                    index_count += 1

        if not quiet:
            summary = [
                f"Reactions: {block_count}",
                f"Indexed CAS: {index_count}",
            ]
            outputs = []
            if jsonl_path:
                outputs.append(jsonl_path)
            if index_path:
                outputs.append(index_path)
            if outputs:
                summary.append("output -> " + ", ".join(outputs))
            print("; ".join(summary), file=sys.stderr)

    finally:
        if jsonl_file is not None:
            jsonl_file.close()
        if index_file is not None:
            index_file.close()

    return block_count, index_count


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Index SciFinder RDF reactions by CAS number")
    parser.add_argument("input", help="Combined RDF file produced by the GUI combiner")
    parser.add_argument("--jsonl", dest="jsonl_path", help="Write reaction blocks to this JSONL file")
    parser.add_argument("--index", dest="index_path", help="Write CAS->offset index JSONL here")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="When writing --jsonl, omit the raw RDF payload and keep metadata only",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Process at most this many reactions (useful for testing)",
    )
    parser.add_argument(
        "--encoding",
        default="latin-1",
        help="Input character encoding (default: latin-1)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress summary output",
    )

    args = parser.parse_args(argv)
    if not args.jsonl_path and not args.index_path:
        parser.error("At least one of --jsonl or --index must be provided")
    if args.limit is not None and args.limit <= 0:
        parser.error("--limit must be positive")
    return args


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    blocks = scan_reactions(
        args.input,
        encoding=args.encoding,
        limit=args.limit,
    )
    write_outputs(
        blocks,
        source_path=args.input,
        jsonl_path=args.jsonl_path,
        index_path=args.index_path,
        metadata_only=args.metadata_only,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
