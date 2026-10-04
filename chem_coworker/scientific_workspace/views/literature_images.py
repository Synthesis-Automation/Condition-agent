"""Bounded presentation of preserved source images, separate from rendered graphs."""

from __future__ import annotations

from typing import Any

from ..adapters.literature import load_source_passage
from ..adapters.literature_reactions import prepared_literature_reaction
from ..adapters.source_images import load_source_image
from ..core.store import InvestigationStore


def answer_literature_images(store: InvestigationStore, answer: dict[str, Any]) -> dict[str, Any]:
    """Read immutable image evidence with a shared 4 MiB/six-image presentation budget."""
    result = {}
    count, total = 0, 0
    for step in answer.get("steps", []):
        for reaction in step.get("literature_reactions", []):
            reference = reaction.get("preparation_ref")
            if not reference or reference in result:
                continue
            entry = {"images": [], "omitted_count": 0, "status": "available"}
            result[reference] = entry
            try:
                prepared_literature_reaction(store, reference)
                preparation = store.read_artifact(reference)["result"]
                _, _, root = load_source_passage(store, preparation["source_ref"])
                for image_ref in preparation.get("scheme_refs", []):
                    record = load_source_image(store, image_ref, root)
                    if count >= 6 or total + record["source_bytes"] > 4 * 1024 * 1024:
                        entry["omitted_count"] += 1
                        continue
                    count += 1
                    total += record["source_bytes"]
                    entry["images"].append({"artifact_ref": image_ref, "locator": record["locator"],
                                            "image_url": "data:" + record["mime_type"] + ";base64," + record["image_base64"],
                                            "assignment_verification": "not_performed"})
            except (OSError, ValueError, KeyError, TypeError):
                entry.update(status="evidence_unavailable", images=[])
    return result
