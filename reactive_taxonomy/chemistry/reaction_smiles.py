"""Lexical reaction fields shared by readers of RDKit reaction SMILES."""

import re


_FIELD_SEPARATOR = re.compile(r"(?<!-)>")
REACTION_FIELD_PARSER_VERSION = "reaction_smiles_fields.v2"


def split_reaction_fields(value: str) -> list[str]:
    """Split reactants, agents and products while preserving a molecular dative bond (->).

    Return all fields, including empty fields. The caller validates arity and
    the required sides for its contract.
    """
    return _FIELD_SEPARATOR.split(value) if "->" in value else value.split(">")
