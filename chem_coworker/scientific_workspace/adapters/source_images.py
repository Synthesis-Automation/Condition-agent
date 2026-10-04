"""Preserve source scheme/page images without asserting their chemical assignments."""

from __future__ import annotations

import base64
import hashlib
import io
from pathlib import Path
from typing import Any

from PIL import Image

from ..core.store import InvestigationEvent, InvestigationStore
from .literature import load_source_passage

SCHEMA_VERSION = "captured_source_image.v1"
MAX_IMAGE_BYTES = 2 * 1024 * 1024


def _image_metadata(data: bytes) -> dict[str, Any]:
    if not data or len(data) > MAX_IMAGE_BYTES:
        raise ValueError("Source image must contain at most 2 MiB")
    with Image.open(io.BytesIO(data)) as image:
        if image.format not in {"PNG", "JPEG"} or image.width * image.height > 20_000_000:
            raise ValueError("Use a PNG/JPEG image of at most 20 million pixels")
        result = {"mime_type": "image/png" if image.format == "PNG" else "image/jpeg",
                  "width": image.width, "height": image.height}
        image.verify()
    return result


def capture_source_image(
    store: InvestigationStore, file: str, source_ref: str, locator: str,
) -> InvestigationEvent:
    """Hash and preserve an agent-supplied screenshot inside the investigation."""
    _, _, root = load_source_passage(store, source_ref)
    path = (store.root / file).resolve()
    if not path.is_relative_to(store.root.resolve()) or not path.is_file():
        raise ValueError("Source image must be a file inside the investigation")
    if not isinstance(locator, str) or not locator.strip() or len(locator) > 1000:
        raise ValueError("Supply a source page/scheme locator of at most 1000 characters")
    if path.stat().st_size > MAX_IMAGE_BYTES:
        raise ValueError("Source image must contain at most 2 MiB")
    data = path.read_bytes()
    metadata = _image_metadata(data)
    return store.append("derived_file", {
        "schema_version": SCHEMA_VERSION, "source_ref": root, "locator": locator,
        "source_path": str(path.relative_to(store.root.resolve())),
        "source_sha256": hashlib.sha256(data).hexdigest(), "source_bytes": len(data),
        "image_base64": base64.b64encode(data).decode("ascii"), **metadata,
        "origin": "agent_supplied_source_image", "review_status": "unreviewed",
        "assignment_verification": "not_performed",
        "limitations": ["Image origin, page locator and chemical assignments have not been independently verified."],
    }, evidence_refs=(root,))


def load_source_image(store: InvestigationStore, reference: str, source_ref: str) -> dict[str, Any]:
    """Validate recorded bytes and source lineage before exposing an image."""
    if not any(event.kind == "derived_file" and event.artifact_ref == reference for event in store.events()):
        raise ValueError("scheme_refs require recorded captured source images")
    value = store.read_artifact(reference)
    if value.get("schema_version") != SCHEMA_VERSION or value.get("source_ref") != source_ref:
        raise ValueError("Source image must belong to the same captured publication")
    data = base64.b64decode(value["image_base64"], validate=True)
    if hashlib.sha256(data).hexdigest() != value.get("source_sha256") or len(data) != value.get("source_bytes"):
        raise ValueError("Source image bytes disagree with recorded identity")
    metadata = _image_metadata(data)
    if any(value.get(key) != expected for key, expected in metadata.items()):
        raise ValueError("Source image metadata disagrees with its bytes")
    return value
