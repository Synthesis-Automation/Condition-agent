"""Conservative display compaction of explicitly supplied scheme annotations.

These lexical helpers do not resolve substances, infer reagents, or assess chemistry.
Callers must retain the complete original annotations in the accompanying details.
"""

from __future__ import annotations

from collections.abc import Iterable
import re


_NUMBER = r"(?:\d+(?:\.\d+)?|\.\d+)"
_AMOUNT = re.compile(
    rf"(?<![\w]){_NUMBER}\s*(?:mol\s*%|wt\s*%|vol\s*%|"
    r"mmol|[µμu]mol|mol|kg|mg|[µμu]g|g|mL|ml|[µμu]L|L|mM|M|N|"
    r"equivalents?|equiv\.?|eq\.?)(?![A-Za-z0-9])"
)
_OPERATING = re.compile(
    rf"(?<![\w])(?:[-−]?{_NUMBER}(?:\s*[–-]\s*{_NUMBER})?\s*"
    r"(?:°\s*[CF]|℃|K|hours?|hrs?|h|minutes?|mins?|min|seconds?|sec|s|atm|bar))(?![A-Za-z0-9])",
    re.IGNORECASE,
)
_OPERATING_SUFFIX = re.compile(
    r"\s+(?:at\s+)?(?:reflux|overnight|room\s+temperature|ambient\s+temperature)\s*\.?$", re.I,
)
_LOCANTS = re.compile(r"\b(?:\d+|[NOPS])['′]*(?:,(?:\d+|[NOPS])['′]*)+(?=-)")
_HEADING = re.compile(r"^(?:initial screen|screening|reagents?|catalysts?|solvents?|conditions?)\s*:\s*", re.I)
_PREFIX = re.compile(r"^(?:(?:add(?:ed)?|in|with|using|dry|anhydrous|aqueous|aq\.)\s+)+", re.I)
_WORKUP = re.compile(
    r"\b(?:work[ -]?up|quench(?:ed)?|brine|ice|extract(?:ion|ed|ing)?|wash(?:ed|ing)?|"
    r"drying|dried|silica|chromatograph\w*|eluent|evaporat\w*|purif\w*|"
    r"recrystalli\w*|triturat\w*|filtrat\w*|isolat\w*)\b", re.I,
)
_PROSE = re.compile(
    r"\b(?:unknown|unreported|unspecified|unavailable|none|no|not|without|avoid|omit|"
    r"see|refer|according|patent|literature|example|reported|proposed|hypothesis|"
    r"yield|conversion|purity|gave|afforded|obtained|heated|stirred|stir|heat|cool|"
    r"reaction|mixture|solution|screen|requires?|suggest\w*|consider|may|might|"
    r"should|could|would|evidence|verified|feasibility|necessary|sufficient|"
    r"if|then|after|before|until|possible|expected|results?|details|of|was|were|is|are)\b", re.I,
)
_GENERIC = re.compile(
    r"(?:core|stannane|substrate|starting material|compound|product|intermediate|"
    r"reagent|catalyst|solvent)(?:\s*[-#]?\s*[A-Za-z]?\d+(?:\.\d+)?)?|"
    r"(?:substrate\s+|reaction\s+)?(?:concentration|temperature|pressure|time)|"
    r"[CP]\d+(?:\.\d+)?|room temperature|ambient temperature|r\.?t\.?|"
    r"reflux|overnight|n/?a|tbd", re.I,
)


def _clauses(text: str) -> Iterable[str]:
    """Split ingredient lists while retaining ligand groups and name locants."""
    protected = {index for match in _LOCANTS.finditer(text)
                 for index in range(match.start(), match.end()) if text[index] == ","}
    depth = 0
    start = 0
    index = 0
    while index < len(text):
        character = text[index]
        if character in "([{":
            depth += 1
        elif character in ")]}":
            depth -= 1
            if depth < 0:
                return
        if depth == 0:
            separator = re.match(r"\s+(?:in|with|and|plus)\s+", text[index:], re.I)
            if character in ";\n" or (character == "," and index not in protected) or separator:
                yield text[start:index]
                index += len(separator.group()) if separator else 1
                start = index
                continue
        index += 1
    if depth == 0:
        yield text[start:]


def compact_condition_labels(texts: Iterable[str], reactant_names: Iterable[str] = ()) -> tuple[str, ...]:
    """Keep short explicit ingredient labels, omitting amounts and procedural prose.

    Exact reactant names, generic substrate identifiers, workup clauses and unknown
    conditions are omitted. Spelling and ligand/formula numerals are preserved;
    retained text is a display label, not a recognized or validated chemical identity.
    This intentionally does not extract names from arbitrary narrative procedures.
    """
    reactants: set[str] = set()
    for name in reactant_names:
        normalized = " ".join(name.split()).casefold()
        reactants.add(normalized)
        # A supplied caption such as "B: ethyl ester" declares its own alias;
        # exclude that exact identifier without guessing chemical synonyms.
        alias = re.fullmatch(r"([a-z]\d*):\s*(.+)", normalized)
        if alias:
            reactants.update(alias.groups())
    labels: list[str] = []
    seen: set[str] = set()
    for text in texts:
        supplied = _HEADING.sub("", text.strip())
        # Do not detach an ingredient fragment from a negated or tentative sentence.
        if _PROSE.search(supplied) or _WORKUP.search(supplied):
            continue
        for raw in _clauses(supplied):
            if _WORKUP.search(raw) or _PROSE.search(raw):
                continue
            label = _PREFIX.sub("", raw.strip())
            label = _OPERATING_SUFFIX.sub("", label)
            label = _AMOUNT.sub("", label)
            label = _OPERATING.sub("", label)
            label = re.sub(r"\(\s*[,;]*\s*\)", "", label)
            label = " ".join(label.strip(" \t\r\n,;.").split())
            # After removing a concentration, an explicit wrapper may be exposed.
            label = _PREFIX.sub("", label).strip()
            identity = label.casefold()
            if (not label or identity in reactants or identity in seen or _GENERIC.fullmatch(label)
                    or len(label) > 72 or len(label.split()) > 4
                    or not any(character.isalpha() for character in label)
                    or any(character in label for character in ":%=<>!?\n")):
                continue
            labels.append(label)
            seen.add(identity)
    return tuple(labels)


def compact_yield_label(text: str | None) -> str | None:
    """Return one explicitly stated percentage or range; omit ambiguous/non-yield text."""
    if not text or re.search(r"\b(?:not|unknown|unreported|unspecified|purity|conversion|ee)\b", text, re.I):
        return None
    matches = re.findall(
        rf"(?<![\w.])(?:[<>≤≥≈~]\s*)?{_NUMBER}(?:\s*[–-]\s*{_NUMBER})?\s*%", text,
    )
    values = tuple(dict.fromkeys(re.sub(r"\s+", "", value) for value in matches))
    return values[0] if len(values) == 1 else None
