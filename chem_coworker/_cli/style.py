"""ANSI style helpers for the terminal CLI."""

from __future__ import annotations


class C:
    """ANSI color codes named for semantic roles."""

    R = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    ITALIC = "\033[3m"

    RED = "\033[91m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    BLUE = "\033[94m"
    MAGENTA = "\033[95m"
    CYAN = "\033[96m"
    WHITE = "\033[97m"

    PROMPT = "\033[1;94m"
    LABEL = "\033[1;96m"
    META = "\033[2m"
    TOOL = "\033[93m"
    ANSWER = "\033[97m"
    HYPO = "\033[96m"
    CONF = "\033[2;92m"
    WARN = "\033[91m"
    OK = "\033[92m"
    ERR = "\033[91m"
    SECTION = "\033[2m"


_W = 60
SEP = f"{C.SECTION}{'─' * _W}{C.R}"
SEP_FAT = f"{C.SECTION}{'━' * _W}{C.R}"


def label(icon: str, text: str) -> str:
    """Render a bold-cyan label line fragment."""
    return f"{C.LABEL}{icon} {text}{C.R}"
