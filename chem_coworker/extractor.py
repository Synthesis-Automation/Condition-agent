"""
NotesExtractor — intake pipeline for chemistry documents.

Reads a document (file path, URL, or raw text), calls an LLM to extract
generalizable chemistry notes, and appends them to the appropriate
knowledge_base/notes/{note_type}/{reaction_type}.md file.

This is the "write" side of the notes system. The "read" side is the
read_reaction_notes tool in chem_coworker/tools/notes.py.

Usage (Python API):
    from chem_coworker.extractor import NotesExtractor
    ex = NotesExtractor(provider="openai", model="gpt-5.6-terra")
    result = ex.intake("https://www.orgsyn.org/demo.aspx?prep=v102p0086",
                       reaction_type="suzuki_miyaura")
    print(result["notes_file"])
    print(result["extracted_notes"])

Usage (CLI):
    python -m chem_coworker._cli.app intake https://www.orgsyn.org/demo.aspx?prep=v102p0086
    python -m chem_coworker._cli.app intake my_paper.pdf --reaction-type buchwald_hartwig
    python -m chem_coworker._cli.app intake my_notes.txt
"""
from __future__ import annotations

import logging
import os
import re
from difflib import get_close_matches
from datetime import date
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from dotenv import load_dotenv

load_dotenv()
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# LLM factory (mirrors agent.py pattern, no circular import)
# ---------------------------------------------------------------------------

def _get_llm(provider: Optional[str] = None, model: Optional[str] = None) -> Any:
    from langchain_openai import ChatOpenAI

    provider = provider or os.getenv("LLM_PROVIDER", "openai")
    model = model or os.getenv("LLM_MODEL", "gpt-5.6-terra")

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
        self.model_name = model or os.getenv("LLM_MODEL", "gpt-5.6-terra")
        self.verbose = verbose

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def intake(
        self,
        source: str,
        reaction_type: str = "",
        note_type: str = "reactions",
        save_to_literature: bool = True,
        mismatch_policy: str = "warn",
        dry_run: bool = False,
        unknown_reaction_policy: str = "general",
        confirm_callback: Optional[Callable[[Dict[str, Any]], bool]] = None,
    ) -> Dict[str, Any]:
        """
        Process a source (URL, file path, or raw text) and extract notes.

        Args:
            source: URL, file path, or raw text to process.
            reaction_type: Hint for which notes file to write to.
                           If empty, the LLM will infer it from the document.
            note_type: Which notes subfolder to write to:
                       "reactions" (default), "mechanisms", "substrates", "protocols".
            save_to_literature: If True and source is a URL/file, also save
                                the raw text to the literature/ folder.
            mismatch_policy: How to handle hint vs detected reaction-type mismatch:
                             warn (default), confirm, reject, force.
            dry_run: If True, do not write notes files (returns planned targets only).
            unknown_reaction_policy: How to handle non-taxonomy labels: general,
                                     quarantine (alias of general), or reject.
            confirm_callback: Optional callback used when mismatch_policy=confirm.
                              Receives a mismatch payload and returns True to proceed.

        Returns:
            dict with source, reaction_types, notes_files, extracted_notes,
            and any warnings.
        """
        # Step 1: Load text
        text, source_name, warnings, saved_filename = self._load_source(source, save_to_literature)
        warning_details = [
            self._warning_detail("info", "source_notice", msg)
            for msg in warnings
        ]
        if not text:
            return {
                "success": False,
                "error": f"Could not load content from: {source!r}",
                "warnings": warnings,
                "warning_details": warning_details,
            }

        if self.verbose:
            logger.info(f"[extractor] Loaded {len(text)} chars from {source_name!r}")

        # Step 2: Extract notes via LLM
        # Pass original source so URL is embedded in the citation header
        source_url = source if source.startswith(("http://", "https://")) else ""
        extracted, detected_types = self._extract_with_llm(
            text, source_name, reaction_type, source_url=source_url,
            note_type=note_type, saved_filename=saved_filename,
        )

        if self.verbose:
            logger.info(f"[extractor] Extracted {len(extracted)} chars, types={detected_types}")

        # Step 3: Taxonomy-canonicalize intake labels + apply mismatch policy
        label_plan = self._plan_reaction_type_filing(
            reaction_type_hint=reaction_type,
            detected_types=detected_types,
            mismatch_policy=mismatch_policy,
            unknown_reaction_policy=unknown_reaction_policy,
            confirm_callback=confirm_callback,
        )
        warnings.extend(label_plan["warnings"])
        warning_details.extend(label_plan.get("warning_details", []))

        if label_plan.get("error"):
            return {
                "success": False,
                "error": label_plan["error"],
                "warnings": warnings,
                "warning_details": warning_details,
                "source": source_name,
                "note_type": note_type,
                "reaction_type_hint_raw": reaction_type or "",
                "reaction_type_hint_canonical": label_plan.get("reaction_type_hint_canonical"),
                "reaction_types_detected_raw": label_plan.get("reaction_types_detected_raw", []),
                "reaction_types_detected_canonical": label_plan.get("reaction_types_detected_canonical", []),
                "reaction_types_unknown": label_plan.get("reaction_types_unknown", []),
                "reaction_type_suggestions": label_plan.get("reaction_type_suggestions", {}),
                "mismatch": bool(label_plan.get("mismatch")),
                "mismatch_policy": mismatch_policy,
                "dry_run": dry_run,
            }

        types_to_file: List[str] = label_plan["reaction_types"]
        quarantine_target: Optional[Dict[str, Any]] = label_plan.get("quarantine_target")

        # Step 4: Append to notes files (or compute planned target files for dry-run)
        from .tools.notes import (
            append_notes,
            append_quarantine_notes,
            get_notes_path,
            get_quarantine_notes_path,
        )
        notes_files: List[str] = []
        for rt in types_to_file:
            if dry_run:
                path = get_notes_path(rt, note_type)
            else:
                path = append_notes(rt, extracted, note_type)
                if self.verbose:
                    logger.info(f"[extractor] Appended to {path}")
            notes_files.append(str(path))

        quarantine_file = ""
        if quarantine_target:
            quarantine_bucket = str(quarantine_target.get("bucket") or "unknown_reaction_labels")
            quarantine_labels = list(quarantine_target.get("labels") or [])
            if dry_run:
                qpath = get_quarantine_notes_path(note_type=note_type, bucket=quarantine_bucket)
            else:
                qpath = append_quarantine_notes(
                    extracted,
                    note_type=note_type,
                    bucket=quarantine_bucket,
                    labels=quarantine_labels,
                )
                if self.verbose:
                    logger.info(f"[extractor] Quarantined intake to {qpath}")
            quarantine_file = str(qpath)
            notes_files.append(quarantine_file)

        if not dry_run:
            # Rebuild index so new notes are immediately discoverable
            try:
                from knowledge_base.notes.build_index import build_index
                build_index()
                if self.verbose:
                    logger.info("[extractor] Rebuilt knowledge_base/notes/_index.json")
            except Exception as exc:
                logger.debug(f"[extractor] Could not rebuild index: {exc}")

        return {
            "success": True,
            "source": source_name,
            "note_type": note_type,
            "reaction_types": types_to_file,
            "reaction_type_hint_raw": reaction_type or "",
            "reaction_type_hint_canonical": label_plan.get("reaction_type_hint_canonical"),
            "reaction_types_detected_raw": label_plan.get("reaction_types_detected_raw", []),
            "reaction_types_detected_canonical": label_plan.get("reaction_types_detected_canonical", []),
            "reaction_types_unknown": label_plan.get("reaction_types_unknown", []),
            "reaction_type_suggestions": label_plan.get("reaction_type_suggestions", {}),
            "mismatch": bool(label_plan.get("mismatch")),
            "mismatch_policy": mismatch_policy,
            "unknown_reaction_policy": unknown_reaction_policy,
            "dry_run": dry_run,
            "write_performed": not dry_run,
            "quarantined": bool(quarantine_target),
            "quarantine_file": quarantine_file,
            "notes_files": notes_files,
            "extracted_notes": extracted,
            "char_count": len(extracted),
            "warnings": warnings,
            "warning_details": warning_details,
        }

    @staticmethod
    def _warning_detail(
        severity: str,
        code: str,
        message: str,
        **extra: Any,
    ) -> Dict[str, Any]:
        detail: Dict[str, Any] = {
            "severity": severity,
            "code": code,
            "message": message,
        }
        if extra:
            detail.update(extra)
        return detail

    @staticmethod
    def _suggest_reaction_type_ids(label: str, limit: int = 5) -> List[str]:
        from chemtools.taxonomy.reaction_catalog import list_reaction_type_ids

        token = (label or "").strip().lower()
        if not token:
            return []
        ids = [rid.strip().lower() for rid in list_reaction_type_ids()]
        return get_close_matches(token, ids, n=limit, cutoff=0.45)

    def _plan_reaction_type_filing(
        self,
        reaction_type_hint: str,
        detected_types: List[str],
        mismatch_policy: str,
        unknown_reaction_policy: str,
        confirm_callback: Optional[Callable[[Dict[str, Any]], bool]] = None,
    ) -> Dict[str, Any]:
        """
        Canonicalize taxonomy labels and compute which reaction-type note IDs to file.

        Returns a dict containing:
          reaction_types, reaction_type_hint_canonical, reaction_types_detected_raw,
          reaction_types_detected_canonical, reaction_types_unknown,
          reaction_type_suggestions, mismatch, warnings, warning_details,
          optional quarantine_target, and optional error.
        """
        from chemtools.taxonomy.reaction_catalog import resolve_reaction_type

        mismatch_policy = (mismatch_policy or "warn").strip().lower()
        unknown_reaction_policy = (unknown_reaction_policy or "general").strip().lower()
        valid_mismatch = {"warn", "confirm", "reject", "force"}
        valid_unknown = {"general", "quarantine", "reject"}
        if mismatch_policy not in valid_mismatch:
            return {
                "reaction_types": [],
                "warnings": [],
                "warning_details": [],
                "error": f"Invalid mismatch policy: {mismatch_policy!r}",
            }
        if unknown_reaction_policy not in valid_unknown:
            return {
                "reaction_types": [],
                "warnings": [],
                "warning_details": [],
                "error": f"Invalid unknown reaction policy: {unknown_reaction_policy!r}",
            }

        warnings: List[str] = []
        warning_details: List[Dict[str, Any]] = []
        reaction_type_suggestions: Dict[str, List[str]] = {}

        def _push_warning(
            message: str,
            *,
            severity: str = "warn",
            code: str = "intake_warning",
            **extra: Any,
        ) -> None:
            warnings.append(message)
            warning_details.append(self._warning_detail(severity, code, message, **extra))

        def _canonicalize_many(labels: List[str]) -> tuple[List[str], List[str]]:
            canonical: List[str] = []
            unknown: List[str] = []
            for raw in labels:
                token = (raw or "").strip()
                if not token:
                    continue
                resolved = resolve_reaction_type(token)
                if resolved:
                    canonical_id = re.sub(r"[\s\-]+", "_", resolved.strip().lower())
                    if canonical_id not in canonical:
                        canonical.append(canonical_id)
                elif token not in unknown:
                    unknown.append(token)
            return canonical, unknown

        hint_raw = (reaction_type_hint or "").strip()
        hint_canonical = None
        if hint_raw:
            resolved_hint = resolve_reaction_type(hint_raw)
            if resolved_hint:
                hint_canonical = re.sub(r"[\s\-]+", "_", resolved_hint.strip().lower())
        hint_unknown = [hint_raw] if hint_raw and not hint_canonical else []
        detected_raw = [t for t in detected_types if (t or "").strip()]
        detected_canonical, detected_unknown = _canonicalize_many(detected_raw)

        unknown_labels = []
        for raw in hint_unknown + detected_unknown:
            if raw not in unknown_labels:
                unknown_labels.append(raw)

        if unknown_labels:
            msg = (
                "Unknown reaction label(s) not found in taxonomy: "
                + ", ".join(unknown_labels)
            )
            for raw in unknown_labels:
                reaction_type_suggestions[raw] = self._suggest_reaction_type_ids(raw)
            if unknown_reaction_policy == "reject":
                return {
                    "reaction_types": [],
                    "reaction_type_hint_canonical": hint_canonical,
                    "reaction_types_detected_raw": detected_raw,
                    "reaction_types_detected_canonical": detected_canonical,
                    "reaction_types_unknown": unknown_labels,
                    "reaction_type_suggestions": reaction_type_suggestions,
                    "mismatch": False,
                    "warnings": warnings,
                    "warning_details": warning_details,
                    "error": msg,
                }
            if unknown_reaction_policy == "quarantine":
                _push_warning(
                    msg + " (quarantine enabled)",
                    code="unknown_reaction_type",
                    labels=unknown_labels,
                    suggestions=reaction_type_suggestions,
                )
            else:
                _push_warning(
                    msg + " (routing unknown labels to general)",
                    code="unknown_reaction_type",
                    labels=unknown_labels,
                    suggestions=reaction_type_suggestions,
                )

        mismatch = False
        if hint_canonical and detected_canonical and hint_canonical not in detected_canonical:
            mismatch = True
            mismatch_msg = (
                "Reaction type hint mismatch: "
                f"hint={hint_canonical} vs detected={', '.join(detected_canonical)}"
            )
            if mismatch_policy == "reject":
                return {
                    "reaction_types": [],
                    "reaction_type_hint_canonical": hint_canonical,
                    "reaction_types_detected_raw": detected_raw,
                    "reaction_types_detected_canonical": detected_canonical,
                    "reaction_types_unknown": unknown_labels,
                    "reaction_type_suggestions": reaction_type_suggestions,
                    "mismatch": True,
                    "warnings": warnings,
                    "warning_details": warning_details,
                    "error": mismatch_msg,
                }
            if mismatch_policy == "confirm":
                if confirm_callback is None:
                    return {
                        "reaction_types": [],
                        "reaction_type_hint_canonical": hint_canonical,
                        "reaction_types_detected_raw": detected_raw,
                        "reaction_types_detected_canonical": detected_canonical,
                        "reaction_types_unknown": unknown_labels,
                        "reaction_type_suggestions": reaction_type_suggestions,
                        "mismatch": True,
                        "warnings": warnings,
                        "warning_details": warning_details,
                        "error": mismatch_msg + " (confirmation required)",
                    }
                proceed = bool(confirm_callback({
                    "hint": hint_canonical,
                    "detected": detected_canonical,
                    "unknown_labels": unknown_labels,
                    "message": mismatch_msg,
                }))
                if not proceed:
                    return {
                        "reaction_types": [],
                        "reaction_type_hint_canonical": hint_canonical,
                        "reaction_types_detected_raw": detected_raw,
                        "reaction_types_detected_canonical": detected_canonical,
                        "reaction_types_unknown": unknown_labels,
                        "reaction_type_suggestions": reaction_type_suggestions,
                        "mismatch": True,
                        "warnings": warnings,
                        "warning_details": warning_details,
                        "error": mismatch_msg + " (user declined confirmation)",
                    }
                _push_warning(
                    mismatch_msg + " (confirmed)",
                    code="reaction_type_mismatch",
                    hint=hint_canonical,
                    detected=detected_canonical,
                )
            elif mismatch_policy == "warn":
                _push_warning(
                    mismatch_msg,
                    code="reaction_type_mismatch",
                    hint=hint_canonical,
                    detected=detected_canonical,
                )

        types_to_file: List[str] = []
        if mismatch_policy == "force" and hint_canonical:
            types_to_file = [hint_canonical]
            if mismatch:
                _push_warning(
                    "Mismatch policy=force: using hinted reaction type and ignoring detected labels",
                    code="reaction_type_mismatch_forced",
                    hint=hint_canonical,
                    detected=detected_canonical,
                )
        else:
            if hint_canonical:
                types_to_file.append(hint_canonical)
            for rt in detected_canonical:
                if rt not in types_to_file:
                    types_to_file.append(rt)

        quarantine_target = None
        if not types_to_file:
            if unknown_labels and unknown_reaction_policy == "quarantine":
                quarantine_target = {
                    "bucket": "unknown_reaction_labels",
                    "labels": unknown_labels,
                }
                _push_warning(
                    "No taxonomy-canonical reaction type resolved; filing to quarantine review store",
                    code="quarantine_routing",
                    labels=unknown_labels,
                )
            else:
                types_to_file = ["general"]
                _push_warning(
                    "No taxonomy-canonical reaction type resolved; filing under general",
                    code="general_fallback",
                )

        return {
            "reaction_types": types_to_file,
            "reaction_type_hint_canonical": hint_canonical,
            "reaction_types_detected_raw": detected_raw,
            "reaction_types_detected_canonical": detected_canonical,
            "reaction_types_unknown": unknown_labels,
            "reaction_type_suggestions": reaction_type_suggestions,
            "mismatch": mismatch,
            "warnings": warnings,
            "warning_details": warning_details,
            "quarantine_target": quarantine_target,
        }

    # ------------------------------------------------------------------
    # Internal: load source
    # ------------------------------------------------------------------

    def _load_source(
        self, source: str, save_to_literature: bool
    ) -> tuple[str, str, List[str], str]:
        """Returns (text, source_name, warnings, saved_filename).

        saved_filename is the filename written to knowledge_base/sources/ (empty string if not saved).
        """
        warnings: List[str] = []

        # URL
        if source.startswith(("http://", "https://")):
            from .tools.literature import _fetch_webpage
            result = _fetch_webpage(
                url=source,
                save_as="__auto__" if save_to_literature else "",
                max_chars=40000,  # intake can handle more than query-time
            )
            if not result.get("success"):
                return "", source, [result.get("error", "fetch failed")], ""
            text = result.get("text", "")
            title = result.get("title") or source
            saved_filename = result.get("saved_as", "")
            if saved_filename:
                warnings.append(f"Saved to knowledge_base/sources/{saved_filename}")
            return text, title, warnings, saved_filename

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
                return "", str(path), [f"Unsupported file type: {suffix}"], ""
            saved_filename = path.name if save_to_literature else ""
            return text, path.name, warnings, saved_filename

        # Raw text (long strings only — short strings are ambiguous)
        if len(source) > 200:
            return source, "pasted text", warnings, ""

        return "", source, [f"Source not recognized as URL, file path, or text: {source!r}"], ""

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
        """Copy a local file to the knowledge_base/sources/ folder if not already there."""
        from .tools.literature import _get_docs_dir
        docs_dir = _get_docs_dir()
        dest = docs_dir / src.name
        if not dest.exists():
            try:
                docs_dir.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(src.read_bytes())
                logger.info(f"[extractor] Copied {src.name} → knowledge_base/sources/")
            except Exception as exc:
                logger.warning(f"[extractor] Could not copy to knowledge_base/sources/: {exc}")

    # ------------------------------------------------------------------
    # Internal: LLM extraction
    # ------------------------------------------------------------------

    def _extract_with_llm(
        self,
        text: str,
        source_name: str,
        reaction_type_hint: str,
        source_url: str = "",
        note_type: str = "reactions",
        saved_filename: str = "",
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
        from .prompts import build_extract_prompt
        from langchain_core.messages import HumanMessage

        # Truncate to keep prompt manageable (LLM context limit)
        MAX_DOC = 16000
        doc_text = text[:MAX_DOC]
        if len(text) > MAX_DOC:
            doc_text += f"\n\n[... document truncated at {MAX_DOC} chars for extraction ...]"

        prompt = build_extract_prompt(note_type).format(
            source_name=source_name,
            source_url=source_url or "(not available)",
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
        # Local source file — enables read_literature_source look-up
        if saved_filename:
            header_lines.append(f"source_file: {saved_filename}")
        # Citation fields in consistent order
        for field in _CITE_FIELDS:
            if field in citation:
                header_lines.append(f"{field}: {citation[field]}")
        # Tags line (last in header block, blank line after)
        if tags_line:
            header_lines.append(tags_line.rstrip())
        header = "\n".join(header_lines) + "\n\n"

        return header + notes_content, detected
