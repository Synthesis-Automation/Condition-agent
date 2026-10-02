"""Paged inspection of saved JSON with explicit bounds and ancestor cautions."""

from __future__ import annotations

import json
from collections.abc import Mapping
from itertools import islice
from typing import Any

from .projection import _COMMON, _Projection


def _context_fields(view: _Projection, item: Mapping[str, Any], fields: tuple[str, ...], path: str) -> dict[str, Any]:
    """Keep ancestor cautions visible without repeating surrounding result trees."""
    def shallow(value: Any, location: str) -> Any:
        if isinstance(value, (Mapping, list, tuple)):
            view.truncated(location, "ancestor_context_detail", total=len(value))
            return {"summary_omitted": True, "value_type": type(value).__name__, "total": len(value)}
        return view.value(value, location)

    result = {}
    for key in fields:
        if key not in item:
            continue
        value, location = item[key], f"{path}.{key}"
        if isinstance(value, (list, tuple)):
            result[key] = view.preview(value, location, shallow, limit=2, describe_collection=False)
        elif isinstance(value, Mapping) and key == "error":
            result[key] = view.pick(value, ("type", "message"), location)
            if set(value) - {"type", "message"}:
                view.truncated(location, "ancestor_context_detail", total_fields=len(value))
        else:
            result[key] = shallow(value, location)
    return result

def inspect_artifact_payload(
    payload: Any, *, artifact_ref: str, path: tuple[str | int, ...] | list[str | int] = (),
    offset: int = 0, limit: int = 5,
) -> dict[str, Any]:
    """Project one literal JSON path, retaining bounded context and pagination.

    Paths are key/index sequences, never expressions. Collection counts refer to
    this saved artifact; nested previews and context can themselves be truncated.
    """
    if (not isinstance(path, (tuple, list)) or len(path) > 12 or any(
            not (isinstance(part, str) and len(part) <= 80 or type(part) is int and part >= 0)
            for part in path)):
        raise ValueError("path must contain at most 12 literal keys (up to 80 characters) or nonnegative indices")
    if type(offset) is not int or offset < 0 or type(limit) is not int or not 1 <= limit <= 20:
        raise ValueError("offset must be nonnegative and limit must be between 1 and 20")
    value = payload
    selected_path = "$"
    ancestors = [(selected_path, value)]
    for part in path:
        if isinstance(value, Mapping) and isinstance(part, str):
            value = value[part]
        elif isinstance(value, list) and type(part) is int:
            value = value[part]
        else:
            raise ValueError(f"Path component has the wrong container type at {selected_path}")
        selected_path += f"[{json.dumps(part, ensure_ascii=False)}]"
        ancestors.append((selected_path, value))
    if not isinstance(value, (Mapping, list)) and offset:
        raise ValueError("offset is supported only for a selected list or mapping")

    view = _Projection()
    view.remaining_nodes = 120
    view.remaining_text = 4000
    context_view = _Projection()
    context_view.remaining_nodes = 80
    context_view.remaining_text = 1200
    context_view.text_limit = 200
    fields = tuple(dict.fromkeys((*_COMMON, "operation", "execution_status", "errors", "uncertainties",
                                 "evidence_refs", "source_ref", "retrieval_status", "acquisition", "claim_support",
                                 "compatible", "actionable", "admission_eligible", "cautions", "risks",
                                 "unresolved_requirements", "analysis_warnings", "evidence_warnings")))
    context = []
    for index, (location, item) in enumerate(ancestors):
        # A requested warning/error field belongs in the preview, not twice in
        # the context budget as well. Other ancestor cautions stay visible.
        selected_fields = tuple(key for key in fields if index >= len(path) or key != path[index])
        if isinstance(item, Mapping) and any(key in item for key in selected_fields):
            context.append({"path": location, "fields": _context_fields(context_view, item, selected_fields, location)})
    total = len(value) if isinstance(value, (Mapping, list)) else 1
    shown = min(limit, max(0, total - offset)) if isinstance(value, (Mapping, list)) else 1
    page = {"offset": offset, "limit": limit, "total": total, "shown": shown,
            "next_offset": offset + shown if offset + shown < total else None,
            "omitted_before": min(offset, total), "omitted_after": max(0, total - offset - shown),
            "count_scope": "saved_artifact", "value_type": type(value).__name__}
    if isinstance(value, Mapping):
        preview = {key: view.value(value[key], f"{selected_path}[{json.dumps(key, ensure_ascii=False)}]")
                   for key in islice(value, offset, offset + shown)}
    elif isinstance(value, list):
        preview = [view.value(item, f"{selected_path}[{index}]")
                   for index, item in enumerate(value[offset:offset + shown], offset)]
    else:
        preview = view.value(value, selected_path)
    if page["omitted_before"] or page["omitted_after"]:
        view.truncated(selected_path, "collection_page", **page)
    # Selection and ancestor context have separate budgets, so large surrounding
    # warnings cannot consume the requested field's preview allowance.
    view.collections.extend(context_view.collections)
    view.truncations.extend(context_view.truncations[:max(0, 30 - len(view.truncations))])
    view.truncation_count += context_view.truncation_count
    result = {
        "inspection_schema_version": "scientific_artifact_inspection.v1", "artifact_ref": artifact_ref,
        "path": list(path), "json_path": selected_path, "preview": preview, "page": page, "context": context,
        "inspection": {"projection_only": True, "collections": view.collections,
                       "truncations": view.truncations, "truncation_count": view.truncation_count,
                       "text_budget_characters": 4000, "serialized_byte_limit": 24000,
                       "context_text_budget_characters": 1200,
                       "hint": "Inspect a more specific path or the saved artifact for complete evidence and warnings. "
                               "Page counts describe saved collections, not the complete source dataset."},
    }
    # Long JSON keys and deeply nested paths can outweigh the text budget. Report
    # an omitted preview rather than silently returning an oversized console dump.
    def size() -> int:
        return len(json.dumps(result, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))

    if size() > 23500:
        result["preview"] = {"summary_omitted": True}
        page.update({"selected": shown, "shown": 0, "next_offset": offset if offset < total else None,
                     "omitted_after": max(0, total - offset)})
        result["inspection"].update({
            "collections": [], "truncations": [{"path": selected_path, "reason": "serialized_preview_budget"}],
            "truncation_count": view.truncation_count + 1, "metadata_omitted": True,
        })
        while size() > 23500 and len(context) > 1:
            context.pop()
            result["inspection"]["context_omitted"] = True
        if size() > 23500 and context:
            context[0]["fields"] = {key: item for key, item in context[0]["fields"].items()
                                    if item is None or isinstance(item, (bool, int, float, str))}
            result["inspection"]["context_omitted"] = True
    result["inspection"]["truncated"] = bool(result["inspection"]["truncation_count"])
    return result
