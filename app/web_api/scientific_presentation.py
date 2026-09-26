"""Read-only Markdown and molecular drawing views of saved scientific messages."""

from __future__ import annotations

import base64
from functools import lru_cache
import json
import re
from typing import Any, Iterable
from urllib.parse import urlsplit
import xml.etree.ElementTree as ET

from markdown_it import MarkdownIt
from markdown_it.token import Token
from rdkit import rdBase

from visualization import (
    SchemeAnnotation, SchemeMolecule, render_annotated_scheme_svg,
    render_molecule_image_bytes,
)

from chem_coworker.scientific_workspace.answer_contracts import ScientificAnswer


_SMILES_TOKEN = re.compile(r"[A-Za-z0-9@+\[\]()=#$%./\\*:\-]+")
_SMILES_ALPHABET = re.compile(r"(?:Br|Cl|[BCNOPSFIbcnops]|\[[^\]\s]+\]|[0-9@+()=#$%./\\*:\-])+")
_REFERENCE = re.compile(r"sha256:[0-9a-f]{64}")
_MAX_CANDIDATES = 48
_MAX_STRUCTURES = 12


def _walk(tokens: Iterable[Token]) -> Iterable[Token]:
    for token in tokens:
        yield token
        if token.children:
            yield from _walk(token.children)


def _candidates(tokens: list[Token]) -> Iterable[str]:
    """Find possible notation, leaving molecular validity to the existing renderer."""
    seen: set[str] = set()
    for token in _walk(tokens):
        explicit = token.type == "code_inline" or (
            token.type == "fence" and token.info.strip().lower() in {"smiles", "smi"}
        )
        if not explicit and token.type != "text":
            continue
        # URLs, artifact hashes and reaction SMILES are not individual molecules.
        for word in token.content.split():
            if "://" in word or "sha256:" in word or ">" in word:
                continue
            for match in _SMILES_TOKEN.finditer(word):
                candidate = match.group().strip(".,;:")
                if not candidate or len(candidate) > 2000 or candidate in seen:
                    continue
                if not _SMILES_ALPHABET.fullmatch(candidate):
                    continue
                if not explicit and not (
                    re.search(r"[=#[\]()0-9]", candidate)
                    or re.fullmatch(r"(?:Br|Cl|[BCNOPSFI]){2,}", candidate)
                ):
                    continue
                seen.add(candidate)
                yield candidate
                if len(seen) >= _MAX_CANDIDATES:
                    return


@lru_cache(maxsize=64)
def present_message(
    text: str, conversation_id: str,
    references: tuple[tuple[str, str, str], ...] = (),
) -> dict[str, Any]:
    """Render untrusted Markdown safely and depict explicit parseable SMILES.

    This is a disposable view, not scientific evidence. Original messages and
    artifacts are unchanged. Cache work across progress polls and page reloads.
    """
    if not re.fullmatch(r"[0-9a-f]{32}", conversation_id):
        raise ValueError("Invalid conversation identifier")
    markdown = MarkdownIt("commonmark", {"html": False}).enable(["table", "strikethrough"])
    # Generated answers cannot cause external image loads or inject raw HTML.
    markdown.disable("image")
    tokens = markdown.parse(text)
    citations = {ref: (title, url) for ref, title, url in references}

    def citation(reference: str) -> tuple[str, str]:
        if reference not in citations:
            citations[reference] = (
                f"Saved evidence {len(citations) + 1}",
                f"/api/v1/scientific/conversations/{conversation_id}/artifacts/{reference}",
            )
        return citations[reference]

    # Hashes in prose or inline code become readable links. Keep fenced code
    # literal and never nest generated links inside existing Markdown links.
    for token in tokens:
        if not token.children:
            continue
        children = []
        current_reference = None
        in_link = False
        for child in token.children:
            if child.type == "link_open":
                in_link = True
                href = child.attrGet("href") or ""
                current_reference = href if _REFERENCE.fullmatch(href) else None
            elif child.type == "link_close":
                in_link = False
                current_reference = None
            if child.type not in {"text", "code_inline"}:
                children.append(child)
                continue
            if in_link:
                if current_reference and _REFERENCE.fullmatch(child.content):
                    child.type, child.tag = "text", ""
                    child.content = citation(current_reference)[0]
                children.append(child)
                continue
            cursor = 0
            for match in _REFERENCE.finditer(child.content):
                if match.start() > cursor:
                    fragment = Token(child.type, child.tag, 0)
                    fragment.content = child.content[cursor:match.start()]
                    children.append(fragment)
                title, _ = citation(match.group())
                opening = Token("link_open", "a", 1)
                opening.attrSet("href", match.group())
                label = Token("text", "", 0)
                label.content = title
                children.extend((opening, label, Token("link_close", "a", -1)))
                cursor = match.end()
            if cursor == 0:
                children.append(child)
            elif cursor < len(child.content):
                fragment = Token(child.type, child.tag, 0)
                fragment.content = child.content[cursor:]
                children.append(fragment)
        token.children = children
    for token in _walk(tokens):
        if token.type != "link_open":
            continue
        href = token.attrGet("href") or ""
        if _REFERENCE.fullmatch(href):
            _, url = citation(href)
            # Only validated source web URLs or our own scoped artifact path.
            artifact = f"/api/v1/scientific/conversations/{conversation_id}/artifacts/{href}"
            token.attrSet("href", url if urlsplit(url).scheme.lower() in {"http", "https"} else artifact)
        elif urlsplit(href).scheme.lower() not in {"https", "http"}:
            token.attrSet("href", "#")
        token.attrSet("target", "_blank")
        token.attrSet("rel", "noopener noreferrer")
    drawings = []
    for smiles in _candidates(tokens):
        try:
            # Invalid candidates remain in the text, without noisy parser logs.
            with rdBase.BlockLogs():
                svg = render_molecule_image_bytes(smiles, size=(440, 260), image_format="svg")
        except (ValueError, RuntimeError):
            continue
        drawings.append({
            "smiles": smiles,
            "image_url": "data:image/svg+xml;base64," + base64.b64encode(svg).decode("ascii"),
        })
        if len(drawings) == _MAX_STRUCTURES:
            break
    return {
        "html": markdown.renderer.render(tokens, markdown.options, {}),
        "structures": drawings,
    }


def present_conversation(conversation: dict[str, Any]) -> dict[str, Any]:
    """Add presentation fields without mutating saved turns or agent answers."""
    identity = conversation["id"]
    turns = []
    for turn in conversation["turns"]:
        view = {**turn, "question_presentation": present_message(turn["question"], identity)}
        if turn.get("answer"):
            references = ()
            if turn["answer"].get("schema_version") == "scientific_answer.v2":
                view["structured_presentation"] = _structured_view(
                    json.dumps(turn["answer"], sort_keys=True), identity,
                )
                references = tuple(
                    (source["artifact_ref"], source["title"], source["url"] or source["artifact_url"])
                    for source in view["structured_presentation"].get("sources", [])
                )
            view["answer_presentation"] = present_message(turn["answer"]["answer_markdown"], identity, references)
        turns.append(view)
    return {**conversation, "turns": turns}


def _svg_url(svg: bytes) -> str:
    return "data:image/svg+xml;base64," + base64.b64encode(svg).decode("ascii")


def _route_overview(steps: list[dict[str, Any]]) -> str:
    """Draw declared step dependencies, without asserting chemical feasibility."""
    root = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "viewBox": f"0 0 660 {len(steps) * 90 + 20}",
        "role": "img", "aria-label": "Declared route step dependencies",
    })
    ET.SubElement(root, "title").text = "Declared route dependencies; not a feasibility assessment"
    positions = {step["id"]: index for index, step in enumerate(steps)}
    for index, step in enumerate(steps):
        y = 15 + index * 90
        for parent in step["after_step_ids"]:
            parent_y = 15 + positions[parent] * 90 + 48
            lane = 12 + (index % 3) * 8
            ET.SubElement(root, "path", {
                "d": f"M 40 {parent_y} H {lane} V {y + 24} H 37",
                "fill": "none", "stroke": "#668477", "stroke-width": "2",
            })
            ET.SubElement(root, "path", {"d": f"M 32 {y + 19} L 40 {y + 24} L 32 {y + 29}", "fill": "#668477"})
        ET.SubElement(root, "rect", {"x": "40", "y": str(y), "width": "600", "height": "48", "rx": "7", "fill": "#edf3eb", "stroke": "#a8bfae"})
        ET.SubElement(root, "text", {"x": "55", "y": str(y + 29), "font-family": "sans-serif", "font-size": "15", "fill": "#16342f"}).text = f"{step['id']}: {step['title'][:55]} · {step['basis']}"
    return _svg_url(ET.tostring(root, encoding="utf-8"))


@lru_cache(maxsize=32)
def _structured_view(payload: str, identity: str) -> dict[str, Any]:
    """Produce disposable structure/route SVGs from explicit answer objects."""
    raw = json.loads(payload)
    try:
        answer = ScientificAnswer.model_validate({key: raw[key] for key in ScientificAnswer.model_fields if key in raw})
    except ValueError:
        return {"error": "The saved structured answer does not satisfy its schema. Its original explanation remains below."}
    view = answer.model_dump()
    molecules = {item["id"]: item for item in view["molecules"]}
    for molecule in molecules.values():
        try:
            with rdBase.BlockLogs():
                molecule["image_url"] = _svg_url(render_molecule_image_bytes(molecule["smiles"], size=(440, 260), image_format="svg"))
            molecule["drawing_status"] = "drawn"
        except (ValueError, RuntimeError):
            molecule["drawing_status"] = "invalid_or_unsupported_notation"
    steps = {item["id"]: item for item in view["steps"]}
    for step in steps.values():
        step["reaction_smiles"] = (
            ".".join(molecules[key]["smiles"] for key in step["reactant_ids"])
            + ">>" + ".".join(molecules[key]["smiles"] for key in step["product_ids"])
        )
        try:
            with rdBase.BlockLogs():
                svg = render_annotated_scheme_svg(
                    tuple(SchemeMolecule(molecules[key]["name"], molecules[key]["smiles"]) for key in step["reactant_ids"]),
                    tuple(SchemeMolecule(molecules[key]["name"], molecules[key]["smiles"]) for key in step["product_ids"]),
                    title=step["title"], basis=step["basis"],
                    conditions=tuple(SchemeAnnotation(item["text"], item["basis"]) for item in step["conditions"]),
                    yield_info=SchemeAnnotation(step["yield_info"]["text"], step["yield_info"]["basis"]) if step["yield_info"] else None,
                )
                step["image_url"] = _svg_url(svg)
                step["scheme_width"] = float(ET.fromstring(svg).get("width"))
            step["drawing_status"] = "drawn"
        except (ValueError, RuntimeError):
            step["drawing_status"] = "invalid_or_unsupported_notation"
    for route in view["routes"]:
        route_steps = [steps[key] for key in route["step_ids"]]
        route["image_url"] = _route_overview(route_steps)
        produced = {key for step in route_steps for key in step["product_ids"]}
        consumed = {key for step in route_steps for key in step["reactant_ids"]}
        route["unreached_target_ids"] = sorted(set(answer.target_molecule_ids) - (produced - consumed))
    for source in view["sources"]:
        source["artifact_url"] = f"/api/v1/scientific/conversations/{identity}/artifacts/{source['artifact_ref']}"
    return view
