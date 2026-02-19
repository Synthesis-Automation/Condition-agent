"""
NotesExtractor — intake pipeline for chemistry documents.

Reads a document (file path, URL, or raw text), calls an LLM to extract
generalizable chemistry notes, and appends them to the appropriate
notes/{reaction_type}.md file.

This is the "write" side of the notes system. The "read" side is the
read_reaction_notes tool in chem_coworker/tools/notes.py.

Usage (Python API):
    from chem_coworker.extractor import NotesExtractor
    ex = NotesExtractor(provider="openai", model="o4-mini")
    result = ex.intake("https://www.orgsyn.org/demo.aspx?prep=v102p0086",
                       reaction_type="suzuki_miyaura")
    print(result["notes_file"])
    print(result["extracted_notes"])

Usage (CLI):
    python chem_coworker/cli.py intake https://www.orgsyn.org/demo.aspx?prep=v102p0086
    python chem_coworker/cli.py intake my_paper.pdf --reaction-type buchwald_hartwig
    python chem_coworker/cli.py intake my_notes.txt
"""
from __future__ import annotations

import logging
import os
import re
from datetime import date
from pathlib import Path
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv

load_dotenv()
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# LLM factory (mirrors agent.py pattern, no circular import)
# ---------------------------------------------------------------------------

def _get_llm(provider: Optional[str] = None, model: Optional[str] = None) -> Any:
    from langchain_openai import ChatOpenAI

    provider = provider or os.getenv("LLM_PROVIDER", "openai")
    model = model or os.getenv("LLM_MODEL", "o4-mini")

    if provider == "aliyun":
        api_key = os.getenv("ALIYUN_API_KEY")
        base_url = os.getenv("ALIYUN_BASE_URL",
                             "https://dashscope.aliyuncs.com/compatible-mode/v1")
    else:
        api_key = os.getenv("OPENAI_API_KEY")
        base_url = os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1")

    _no_temp = ("o1", "o3", "o4", "gpt-5")
    kwargs: Dict[str, Any] = {"model": model, "api_key": api_key, "base_url": base_url}
    if not (provider == "openai" and any(model.startswith(p) for p in _no_temp)):
        kwargs["temperature"] = 0

    return ChatOpenAI(**kwargs)


# ---------------------------------------------------------------------------
# NotesExtractor
# ---------------------------------------------------------------------------

class NotesExtractor:
    """
    Intake pipeline: document → LLM extraction → append to notes/{reaction_type}.md

    Accepts:
      - URL (http/https) — fetches and extracts text
      - File path (.txt, .md, .pdf)
      - Raw text string
    """

    def __init__(
        self,
        provider: Optional[str] = None,
        model: Optional[str] = None,
        verbose: bool = False,
    ):
        self.llm = _get_llm(provider, model)
        self.model_name = model or os.getenv("LLM_MODEL", "o4-mini")
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def intake(
        self,
        source: str,
        reaction_type: str = "",
        save_to_literature: bool = True,
    ) -> Dict[str, Any]:
        """
        Process a source (URL, file path, or raw text) and extract notes.

        Args:
            source: URL, file path, or raw text to process.
            reaction_type: Hint for which notes file to write to.
                           If empty, the LLM will infer it from the document.
            save_to_literature: If True and source is a URL/file, also save
                                the raw text to the literature/ folder.

        Returns:
            dict with source, reaction_types, notes_files, extracted_notes,
            and any warnings.
        """
        # Step 1: Load text
        text, source_name, warnings = self._load_source(source, save_to_literature)
        if not text:
            return {
                "success": False,
                "error": f"Could not load content from: {source!r}",
                "warnings": warnings,
            }

        if self.verbose:
            logger.info(f"[extractor] Loaded {len(text)} chars from {source_name!r}")

        # Step 2: Extract notes via LLM
        # Pass original source so URL is embedded in the citation header
        source_url = source if source.startswith(("http://", "https://")) else ""
        extracted, detected_types = self._extract_with_llm(
            text, source_name, reaction_type, source_url=source_url
        )

        if self.verbose:
            logger.info(f"[extractor] Extracted {len(extracted)} chars, types={detected_types}")

        # Step 3: Append to notes files
        from .tools.notes import append_notes, _normalize_reaction_type

        # Determine which reaction types to file under
        types_to_file: List[str] = []
        if reaction_type:
            types_to_file.append(_normalize_reaction_type(reaction_type))
        types_to_file.extend(
            _normalize_reaction_type(t) for t in detected_types
            if _normalize_reaction_type(t) not in types_to_file
        )
        if not types_to_file:
            types_to_file = ["general"]

        notes_files: List[str] = []
        for rt in types_to_file:
            path = append_notes(rt, extracted)
            notes_files.append(str(path))
            if self.verbose:
                logger.info(f"[extractor] Appended to {path}")

        return {
            "success": True,
            "source": source_name,
            "reaction_types": types_to_file,
            "notes_files": notes_files,
            "extracted_notes": extracted,
            "char_count": len(extracted),
            "warnings": warnings,
        }

    # ------------------------------------------------------------------
    # Internal: load source
    # ------------------------------------------------------------------

    def _load_source(
        self, source: str, save_to_literature: bool
    ) -> tuple[str, str, List[str]]:
        """Returns (text, source_name, warnings)."""
        warnings: List[str] = []

        # URL
        if source.startswith(("http://", "https://")):
            from .tools.literature import _fetch_webpage
            result = _fetch_webpage(
                url=source,
                save_as=self._url_to_filename(source) if save_to_literature else "",
                max_chars=40000,  # intake can handle more than query-time
            )
            if not result.get("success"):
                return "", source, [result.get("error", "fetch failed")]
            text = result.get("text", "")
            title = result.get("title") or source
            if result.get("saved_as"):
                warnings.append(f"Saved to literature/{result['saved_as']}")
            return text, title, warnings

        # File path
        path = Path(source)
        if path.exists():
            suffix = path.suffix.lower()
            if suffix == ".pdf":
                from .tools.literature import _extract_pdf_text
                text = _extract_pdf_text(source=path.name, path=path)
                if save_to_literature:
                    self._maybe_copy_to_literature(path)
            elif suffix in (".txt", ".md"):
                text = path.read_text(encoding="utf-8", errors="replace")
                if save_to_literature:
                    self._maybe_copy_to_literature(path)
            else:
                return "", str(path), [f"Unsupported file type: {suffix}"]
            return text, path.name, warnings

        # Raw text (long strings only — short strings are ambiguous)
        if len(source) > 200:
            return source, "pasted text", warnings

        return "", source, [f"Source not recognized as URL, file path, or text: {source!r}"]

    def _url_to_filename(self, url: str) -> str:
        """Generate a filesystem-safe filename from a URL."""
        # Use the last path component, falling back to domain
        from urllib.parse import urlparse
        parsed = urlparse(url)
        name = parsed.path.rstrip("/").split("/")[-1] or parsed.netloc
        name = re.sub(r"[^\w\-.]", "_", name)
        if not name.endswith((".txt", ".md", ".pdf")):
            name += ".txt"
        return name[:80]  # cap length

    def _maybe_copy_to_literature(self, src: Path) -> None:
        """Copy a local file to the literature/ folder if not already there."""
        from .tools.literature import _get_docs_dir
        docs_dir = _get_docs_dir()
        dest = docs_dir / src.name
        if not dest.exists():
            try:
                docs_dir.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(src.read_bytes())
                logger.info(f"[extractor] Copied {src.name} → literature/")
            except Exception as exc:
                logger.warning(f"[extractor] Could not copy to literature/: {exc}")

    # ------------------------------------------------------------------
    # Internal: LLM extraction
    # ------------------------------------------------------------------

    def _extract_with_llm(
        self,
        text: str,
        source_name: str,
        reaction_type_hint: str,
        source_url: str = "",
    ) -> tuple[str, List[str]]:
        """
        Call LLM to extract notes. Returns (markdown_notes, detected_reaction_types).

        The note block starts with a rich citation header:
          ## Source: {title}
          url: {url}           ← if fetched from URL
          doi: {doi}           ← if LLM found it in the document
          journal: {journal}
          year: {year}
          pages: {vol, pages}
          tags: {tags}
        """
        from .prompts import EXTRACT_NOTES_PROMPT
        from langchain_core.messages import HumanMessage

        # Truncate to keep prompt manageable (LLM context limit)
        MAX_DOC = 16000
        doc_text = text[:MAX_DOC]
        if len(text) > MAX_DOC:
            doc_text += f"\n\n[... document truncated at {MAX_DOC} chars for extraction ...]"

        prompt = EXTRACT_NOTES_PROMPT.format(
            source_name=source_name,
            date=date.today().isoformat(),
            document_text=doc_text,
        )

        try:
            response = self.llm.invoke([HumanMessage(content=prompt)])
            content = getattr(response, "content", str(response))
            if isinstance(content, list):
                content = "\n".join(
                    p.get("text", "") if isinstance(p, dict) else str(p)
                    for p in content
                )
        except Exception as exc:
            logger.error(f"[extractor] LLM call failed: {exc}")
            return f"[Extraction failed: {exc}]\n\nRaw document saved to literature folder.", []

        # ── Parse reaction_types ──────────────────────────────────────────
        detected: List[str] = []
        m = re.search(r"reaction_types?\s*:\s*([^\n]+)", content, re.IGNORECASE)
        if m:
            raw_types = m.group(1).strip()
            raw_list = [re.sub(r"[`']", "", t).strip() for t in re.split(r"[,;]", raw_types)]
            for raw in raw_list:
                raw = re.sub(r"\s*\(.*$", "", raw).strip()
                if raw and re.match(r"^[\w\-]{3,50}$", raw):
                    detected.append(raw)

        if not detected and reaction_type_hint:
            detected = [reaction_type_hint]

        # ── Parse tags ───────────────────────────────────────────────────
        tags_line = ""
        m_tags = re.search(r"^`?tags\s*:\s*([^\n`]+)`?", content, re.IGNORECASE | re.MULTILINE)
        if m_tags:
            raw_tags = m_tags.group(1).strip()
            # Normalize each tag: strip whitespace/backticks, snake_case
            tag_list = [
                re.sub(r"[\s\-]+", "_", re.sub(r"[`']", "", t).strip().lower())
                for t in re.split(r"[,;]", raw_tags)
                if re.sub(r"[`'\s]", "", t)
            ]
            # Keep only clean tokens (≤40 chars, no spaces, no parens)
            tag_list = [t for t in tag_list if re.match(r"^[\w\-]{2,40}$", t)]
            if tag_list:
                tags_line = "tags: " + ", ".join(tag_list) + "\n"

        # ── Parse citation metadata (doi, journal, year, pages) ─────────
        _CITE_FIELDS = ["doi", "journal", "year", "pages"]
        citation: Dict[str, str] = {}
        for field in _CITE_FIELDS:
            m_cite = re.search(
                rf"^`?{field}\s*:\s*([^\n`]+)`?",
                content,
                re.IGNORECASE | re.MULTILINE,
            )
            if m_cite:
                citation[field] = m_cite.group(1).strip()

        # ── Strip ALL metadata lines from body ───────────────────────────
        notes_content = re.sub(r"`?reaction_types?\s*:[^\n]*\n?", "", content)
        notes_content = re.sub(r"`?tags\s*:[^\n]*\n?", "", notes_content, flags=re.IGNORECASE)
        for field in _CITE_FIELDS:
            notes_content = re.sub(
                rf"`?{field}\s*:[^\n]*\n?", "", notes_content, flags=re.IGNORECASE
            )
        notes_content = notes_content.strip()

        # ── Build rich citation header ────────────────────────────────────
        # Line 1: ## Source: {title}  ·  {intake date}
        header_lines = [f"## Source: {source_name}  ·  {date.today().isoformat()}"]
        # Optional URL (always first so it's easiest to click/copy)
        if source_url:
            header_lines.append(f"url: {source_url}")
        # Citation fields in consistent order
        for field in _CITE_FIELDS:
            if field in citation:
                header_lines.append(f"{field}: {citation[field]}")
        # Tags line (last in header block, blank line after)
        if tags_line:
            header_lines.append(tags_line.rstrip())
        header = "\n".join(header_lines) + "\n\n"

        return header + notes_content, detected
