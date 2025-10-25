"""Rule-based matcher package supporting scheme and selector databases."""
from .loader import load_db
from .matcher import match
from .types import (
	MatchResult,
	RuleDB,
	SchemeConditionDB,
	SchemeEntry,
	SelectorRule,
	SelectorRuleDB,
)

__all__ = [
	"load_db",
	"match",
	"MatchResult",
	"RuleDB",
	"SchemeConditionDB",
	"SelectorRuleDB",
	"SchemeEntry",
	"SelectorRule",
]
