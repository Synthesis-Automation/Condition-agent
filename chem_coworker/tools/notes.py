"""
Reaction notes tools for ChemCoworker.

Manages a per-reaction-type knowledge base of distilled chemistry wisdom:
warnings, side reactions, substrate scope limitations, and procedural notes.
These notes bridge the gap between structured HTE data (old, no caveats) and
unstructured literature (new, rich in actionable warnings).

Notes folder (in priority order):
  1. CHEM_NOTES_PATH environment variable
  2. {project_root}/notes/

Notes files are plain markdown — human-readable and human-editable.
Broad category files: notes/suzuki_miyaura.md, notes/metal_catalysis.md, etc.
Each note block has a `tags:` line for cross-file keyword retrieval.

Tools:
  - read_reaction_notes : exact lookup by reaction type filename
  - search_notes        : keyword/tag search across ALL notes files
"""
from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional
import os

from ._helpers import _error, _success
from ._base import ToolPlugin

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Folder resolution
# ---------------------------------------------------------------------------

def get_notes_dir() -> Path:
    """Resolve the reaction notes folder."""
    env_path = os.getenv("CHEM_NOTES_PATH")
    if env_path:
        return Path(env_path)
    return Path(__file__).parent.parent.parent / "notes"


# ---------------------------------------------------------------------------
# Reaction type normalization
# ---------------------------------------------------------------------------

# Common aliases → canonical snake_case key used for the notes filename
_ALIASES: Dict[str, str] = {
    # Suzuki
    "suzuki": "suzuki_miyaura",
    "suzuki coupling": "suzuki_miyaura",
    "suzuki-miyaura": "suzuki_miyaura",
    "suzuki miyaura": "suzuki_miyaura",
    # Buchwald-Hartwig
    "buchwald": "buchwald_hartwig",
    "buchwald-hartwig": "buchwald_hartwig",
    "buchwald hartwig": "buchwald_hartwig",
    "bha": "buchwald_hartwig",
    "c-n coupling": "buchwald_hartwig",
    # Negishi
    "negishi": "negishi_coupling",
    # Heck
    "heck": "heck_reaction",
    # Kumada
    "kumada": "kumada_coupling",
    # Sonogashira
    "sonogashira": "sonogashira_coupling",
    # Chan-Lam
    "chan-lam": "chan_lam",
    "chan lam": "chan_lam",
    # SNAr
    "snar": "snar",
    "nucleophilic aromatic substitution": "snar",
    # Mitsunobu
    "mitsunobu": "mitsunobu",
    # Reductive amination
    "reductive amination": "reductive_amination",
    # Amide coupling
    "amide coupling": "amide_coupling",
    "peptide coupling": "amide_coupling",
}


def _normalize_reaction_type(raw: str) -> str:
    """Normalize a reaction type name to a canonical snake_case key."""
    cleaned = raw.strip().lower()
    # Direct alias lookup
    if cleaned in _ALIASES:
        return _ALIASES[cleaned]
    # Already snake_case with no alias needed
    return re.sub(r"[\s\-]+", "_", cleaned)


# ---------------------------------------------------------------------------
# Notes file I/O
# ---------------------------------------------------------------------------

def get_notes_path(reaction_type: str) -> Path:
    """Return the path to the notes file for a reaction type."""
    key = _normalize_reaction_type(reaction_type)
    return get_notes_dir() / f"{key}.md"


def read_notes_text(reaction_type: str) -> Optional[str]:
    """Read the notes file for a reaction type. Returns None if not found."""
    path = get_notes_path(reaction_type)
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8").strip()
    return text if text else None


def append_notes(reaction_type: str, content: str) -> Path:
    """Append extracted notes content to the notes file. Creates it if needed."""
    path = get_notes_path(reaction_type)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists():
        existing = path.read_text(encoding="utf-8")
        path.write_text(existing.rstrip() + "\n\n---\n\n" + content.strip() + "\n",
                        encoding="utf-8")
    else:
        # New file: write a header
        header = f"# {reaction_type.replace('_', ' ').title()} — Reaction Notes\n\n"
        path.write_text(header + content.strip() + "\n", encoding="utf-8")

    return path


# ---------------------------------------------------------------------------
# Tool: read_reaction_notes
# ---------------------------------------------------------------------------

def _read_reaction_notes(reaction_type: str, max_chars: int = 6000) -> Dict[str, Any]:
    """Load accumulated chemistry notes for a reaction type.

    Returns distilled warnings, side reactions, substrate scope limitations,
    and procedural notes extracted from literature and lab experience.
    These notes supplement HTE database conditions with caveats that
    experimental data alone cannot capture.

    Call this in parallel with recommend_conditions to enrich condition
    recommendations with literature-derived warnings.

    Args:
        reaction_type: The reaction type name (e.g. "suzuki_miyaura",
                       "buchwald_hartwig", "snar"). Aliases like "suzuki"
                       or "Suzuki-Miyaura" are also accepted.
        max_chars: Maximum characters to return (default 6000).

    Returns:
        dict with:
          reaction_type: normalized reaction type key
          notes_file: path to the notes file
          content: full notes text (truncated if very long)
          char_count: total characters in notes
          has_notes: True if notes exist for this reaction type
    """
    canonical = _normalize_reaction_type(reaction_type)
    notes_text = read_notes_text(reaction_type)

    if notes_text is None:
        # List what notes do exist, so LLM can suggest alternatives
        notes_dir = get_notes_dir()
        available: List[str] = []
        if notes_dir.exists():
            available = [p.stem for p in sorted(notes_dir.glob("*.md"))
                         if not p.name.startswith("README")]
        return _success({
            "reaction_type": canonical,
            "has_notes": False,
            "content": "",
            "char_count": 0,
            "notes_file": str(get_notes_path(reaction_type)),
            "available_note_files": available,
            "note": (
                f"No notes found for '{canonical}'. "
                + (f"Available: {', '.join(available)}" if available else
                   "Notes folder is empty. Run intake to extract notes from literature.")
            ),
        })

    total_chars = len(notes_text)
    display = notes_text
    if total_chars > max_chars:
        display = notes_text[:max_chars].rsplit("\n\n", 1)[0]
        display += f"\n\n[... truncated at {max_chars} chars; {total_chars} total]"

    return _success({
        "reaction_type": canonical,
        "has_notes": True,
        "content": display,
        "char_count": total_chars,
        "notes_file": str(get_notes_path(reaction_type)),
    })


read_reaction_notes_tool = ToolPlugin(
    name="read_reaction_notes",
    category="notes",
    description=(
        "Load accumulated chemistry notes for a reaction type: solvent warnings, "
        "side reactions, substrate scope limitations, and procedural caveats "
        "extracted from literature. Use when reaction type is confirmed. "
        "Run in parallel with recommend_conditions."
    ),
    prerequisites=[],  # no dependencies — runs in parallel with G1 tools
    fn=_read_reaction_notes,
)


# ---------------------------------------------------------------------------
# Tool: search_notes
# ---------------------------------------------------------------------------

def _search_notes(query: str, top_k: int = 5) -> Dict[str, Any]:
    """Search all reaction notes by keyword or tag across all notes files.

    Use when you don't know which notes file contains the relevant information,
    or when searching for cross-cutting topics (e.g. "copper catalyst",
    "beta-hydride elimination", "sealed vessel", "sp3 coupling").

    Complementary to read_reaction_notes:
      - read_reaction_notes → exact lookup when reaction type is confirmed
      - search_notes        → keyword/tag search when exploring or uncertain

    Tags embedded in note headers (e.g. `tags: copper, sp3_coupling`) receive
    a score bonus, so tag-exact matches rank above body-text matches.

    Args:
        query: Free-text search query. Include catalysts, reagent names,
               bond types, conditions, or problem descriptions.
               e.g. "copper sp3 coupling alkyl halide"
               e.g. "beta hydride elimination side reaction"
               e.g. "pressure vessel sealed toluene"
        top_k: Number of top note sections to return (default 5, max 10).

    Returns:
        dict with:
          found: number of matching sections
          query: the search query used
          results: list of {rank, score, source_file, tags, excerpt}
          notes_searched: list of notes filenames searched
    """
    from .literature import _tokenize, _score_chunk

    notes_dir = get_notes_dir()
    if not notes_dir.exists():
        return _success({
            "found": 0, "query": query, "results": [], "notes_searched": [],
            "note": "Notes folder not found.",
        })

    note_files = [
        f for f in sorted(notes_dir.glob("*.md"))
        if not f.name.startswith("README")
    ]
    if not note_files:
        return _success({
            "found": 0, "query": query, "results": [], "notes_searched": [],
            "note": "No notes files found. Run intake to extract notes from literature.",
        })

    query_tokens = _tokenize(query)
    if not query_tokens:
        return _error("Query is empty or contains only stopwords.")

    scored: List[tuple] = []  # (score, section_text, stem, tags_str)

    for nf in note_files:
        try:
            text = nf.read_text(encoding="utf-8")
        except Exception:
            continue

        # Split into per-source sections (each starts with "## Source:")
        sections = re.split(r"(?=## Source:)", text)

        for section in sections:
            if "## Source:" not in section or len(section.strip()) < 50:
                continue

            # Extract tags line from the section header
            tags_match = re.search(
                r"^tags\s*:\s*([^\n]+)", section, re.IGNORECASE | re.MULTILINE
            )
            tags_str = tags_match.group(1).strip() if tags_match else ""

            # Tag-exact hits get a 2-point bonus each
            tag_tokens = _tokenize(tags_str) if tags_str else []
            tag_bonus = sum(2.0 for qt in query_tokens if qt in tag_tokens)

            body_score = _score_chunk(section, query_tokens)
            total_score = body_score + tag_bonus

            if total_score > 0:
                scored.append((total_score, section, nf.stem, tags_str))

    scored.sort(key=lambda x: x[0], reverse=True)
    top_k = min(top_k, 10)
    top = scored[:top_k]

    if not top:
        return _success({
            "found": 0,
            "query": query,
            "results": [],
            "notes_searched": [f.name for f in note_files],
            "note": (
                f"No sections matched '{query}' across "
                f"{len(note_files)} notes file(s)."
            ),
        })

    results = []
    for rank, (score, section, stem, tags_str) in enumerate(top, 1):
        excerpt = section.strip()
        if len(excerpt) > 800:
            excerpt = excerpt[:800].rsplit("\n", 1)[0] + "\n[...]"
        results.append({
            "rank": rank,
            "score": round(score, 3),
            "source_file": f"{stem}.md",
            "tags": tags_str,
            "excerpt": excerpt,
        })

    return _success({
        "found": len(results),
        "query": query,
        "results": results,
        "notes_searched": [f.name for f in note_files],
    })


search_notes_tool = ToolPlugin(
    name="search_notes",
    category="notes",
    description=(
        "Keyword/tag search across ALL notes files — use when reaction type is "
        "uncertain or when looking for cross-cutting topics (copper catalyst, "
        "beta-hydride elimination, pressure vessel, sp3 coupling, etc.). "
        "Tag-exact matches rank above body-text matches. "
        "Complement to read_reaction_notes (which does exact-type lookup)."
    ),
    prerequisites=[],
    fn=_search_notes,
)


# ---------------------------------------------------------------------------
# All tools in this module
# ---------------------------------------------------------------------------

NOTES_TOOLS = [
    read_reaction_notes_tool,
    search_notes_tool,
]
