"""
chemtools.eval — RDKit-based reaction evaluation and validation.

Provides computable sanity checks for reaction SMILES produced by retrosynthesis
or condition prediction. Designed to complement LLM chemical reasoning, not replace it.

Public API:
    from chemtools.eval.reaction_evaluator import evaluate_reaction, CheckResult, EvalReport
"""
from .reaction_evaluator import evaluate_reaction, CheckResult, EvalReport

__all__ = ["evaluate_reaction", "CheckResult", "EvalReport"]
