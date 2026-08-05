from .cas_number_extractor import (
    CASMatch,
    discover_candidate_files,
    extract_cas_matches_from_file,
    find_cas_numbers_in_text,
    is_valid_cas_number,
    write_matches_to_csv,
    write_matches_to_markdown,
)
from .compound_lookup import (
    CompoundLookupClient,
    CompoundLookupResult,
    lookup_compound_by_cas,
)

__all__ = [
    "CASMatch",
    "CompoundLookupClient",
    "CompoundLookupResult",
    "discover_candidate_files",
    "extract_cas_matches_from_file",
    "find_cas_numbers_in_text",
    "is_valid_cas_number",
    "lookup_compound_by_cas",
    "write_matches_to_csv",
    "write_matches_to_markdown",
]
