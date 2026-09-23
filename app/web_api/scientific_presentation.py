"""Read-only Markdown and molecular drawing views of saved scientific messages."""

from __future__ import annotations

import base64
from functools import lru_cache
import re
from typing import Any, Iterable
from urllib.parse import urlsplit

from markdown_it import MarkdownIt
from markdown_it.token import Token
from rdkit import rdBase

from visualization import render_molecule_image_bytes


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
def present_message(text: str, conversation_id: str) -> dict[str, Any]:
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
    for token in _walk(tokens):
        if token.type != "link_open":
            continue
        href = token.attrGet("href") or ""
        if _REFERENCE.fullmatch(href):
            token.attrSet("href", f"/api/v1/scientific/conversations/{conversation_id}/artifacts/{href}")
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
            view["answer_presentation"] = present_message(turn["answer"]["answer_markdown"], identity)
        turns.append(view)
    return {**conversation, "turns": turns}
