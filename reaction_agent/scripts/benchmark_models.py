"""
Benchmark multiple LLM models on the ReactionReasoningAgent.

Tests accuracy (correct reaction_type from taxonomy) and speed across all
models available via openai/aliyun providers.

Usage:
    python benchmark_models.py
    python benchmark_models.py --models deepseek-v3.2 kimi-k2.5
    python benchmark_models.py --reactions 0 1      # run only cases 0 and 1
    python benchmark_models.py --no-color
"""

from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", message=".*pkg_resources.*")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------
class C:
    HEADER = "\033[95m"; BLUE = "\033[94m"; CYAN = "\033[96m"
    GREEN  = "\033[92m"; YELLOW = "\033[93m"; RED = "\033[91m"
    BOLD   = "\033[1m";  END = "\033[0m"

    @classmethod
    def disable(cls):
        for a in ("HEADER","BLUE","CYAN","GREEN","YELLOW","RED","BOLD","END"):
            setattr(cls, a, "")


# ---------------------------------------------------------------------------
# Ground-truth test cases
# ---------------------------------------------------------------------------
TEST_CASES = [
    {
        "name": "Suzuki coupling",
        "smiles": "Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1",
        "expected": "Suzuki_miyaura",
    },
    {
        "name": "C-N coupling (Ar-I + formamide)",
        "smiles": "Ic1ccccc1.O=CNC1=CCCCC1>>O=CN(C1=CCCCC1)c1ccccc1",
        "expected": "C_N_Coupling",
    },
    {
        "name": "Amide formation",
        "smiles": "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1",
        "expected": "Amide_formation",
    },
    {
        "name": "Buchwald-Hartwig (Ar-Br + amine)",
        "smiles": "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1",
        "expected": "C_N_Coupling",
    },
    {
        "name": "C-O coupling (Ar-Br + phenol)",
        "smiles": "Brc1ccccc1.Oc1ccccc1>>c1ccc(-c2ccccc2)cc1",   # simplified; real product is Ar-O-Ar
        "expected": "C_O_Coupling",
    },
    {
        "name": "Miyaura borylation",
        "smiles": "Brc1ccccc1.B1OC(C)(C)C(C)(C)O1>>B1(c2ccccc2)OC(C)(C)C(C)(C)O1",
        "expected": "Miyaura_borylation",
    },
]

# Models to benchmark: (model_id, provider)
DEFAULT_MODELS = [
    ("deepseek-v3.2",  "aliyun"),
    ("kimi-k2.5",      "aliyun"),
    ("glm-4.7",        "aliyun"),
    ("MiniMax-M2.1",   "aliyun"),
    ("gpt-5.2",        "openai"),
]


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
class Result:
    def __init__(self, model, provider, case_name, expected,
                 got, correct, elapsed, tool_calls, error=""):
        self.model = model
        self.provider = provider
        self.case_name = case_name
        self.expected = expected
        self.got = got
        self.correct = correct
        self.elapsed = elapsed
        self.tool_calls = tool_calls
        self.error = error


# ---------------------------------------------------------------------------
# Run one model × one reaction
# ---------------------------------------------------------------------------
def run_one(model: str, provider: str, case: dict) -> Result:
    from reaction_agent.reasoning_agent import ReactionReasoningAgent

    agent = ReactionReasoningAgent(provider=provider, model=model)
    t0 = time.time()
    try:
        profile = agent.analyze(case["smiles"])
        elapsed = time.time() - t0
        got = profile.reaction_type or ""
        correct = (got == case["expected"])
        return Result(
            model=model, provider=provider,
            case_name=case["name"], expected=case["expected"],
            got=got, correct=correct, elapsed=elapsed,
            tool_calls=profile.total_tool_calls,
            error="; ".join(profile.warnings) if profile.warnings else "",
        )
    except Exception as exc:
        elapsed = time.time() - t0
        return Result(
            model=model, provider=provider,
            case_name=case["name"], expected=case["expected"],
            got="ERROR", correct=False, elapsed=elapsed, tool_calls=0,
            error=str(exc)[:120],
        )


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------
def print_result(r: Result):
    tick = f"{C.GREEN}✓{C.END}" if r.correct else f"{C.RED}✗{C.END}"
    got_col = C.GREEN if r.correct else C.RED
    print(f"  {tick} {r.case_name:<35s} "
          f"got={got_col}{r.got or 'None':<25}{C.END} "
          f"{r.elapsed:5.1f}s  tools={r.tool_calls}")
    if r.error and not r.correct:
        print(f"    {C.YELLOW}⚠ {r.error[:100]}{C.END}")


def print_summary_table(all_results: list[Result], models: list):
    """Print an accuracy × speed table across all models."""
    print(f"\n{C.BOLD}{C.HEADER}{'='*72}{C.END}")
    print(f"{C.BOLD}{C.HEADER}  SUMMARY TABLE{C.END}")
    print(f"{C.BOLD}{C.HEADER}{'='*72}{C.END}")

    # Header
    print(f"\n  {'Model':<22} {'Prov':<7} {'Acc':>5} {'Correct':>8} "
          f"{'Avg s':>7} {'Avg tools':>10}")
    print(f"  {'-'*22} {'-'*7} {'-'*5} {'-'*8} {'-'*7} {'-'*10}")

    for (model, provider) in models:
        rows = [r for r in all_results if r.model == model]
        if not rows:
            continue
        n_correct = sum(1 for r in rows if r.correct)
        n_total = len(rows)
        acc = n_correct / n_total * 100
        avg_s = sum(r.elapsed for r in rows) / n_total
        avg_tools = sum(r.tool_calls for r in rows) / n_total
        acc_col = C.GREEN if acc >= 80 else (C.YELLOW if acc >= 50 else C.RED)
        print(f"  {model:<22} {provider:<7} "
              f"{acc_col}{acc:5.0f}%{C.END} "
              f"  {n_correct}/{n_total}     "
              f"{avg_s:6.1f}s  {avg_tools:9.1f}")

    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Benchmark LLM models on the ReactionReasoningAgent"
    )
    parser.add_argument(
        "--models", nargs="+", default=None,
        help="Model IDs to test (default: all). E.g. --models deepseek-v3.2 gpt-5.2"
    )
    parser.add_argument(
        "--reactions", nargs="+", type=int, default=None,
        help="Indices of test cases to run (default: all). E.g. --reactions 0 1"
    )
    parser.add_argument("--no-color", action="store_true")
    args = parser.parse_args()

    if args.no_color:
        C.disable()

    # Select models
    models_to_run = DEFAULT_MODELS
    if args.models:
        model_map = {m: p for m, p in DEFAULT_MODELS}
        models_to_run = [
            (m, model_map.get(m, "aliyun")) for m in args.models
        ]

    # Select test cases
    cases = TEST_CASES
    if args.reactions is not None:
        cases = [TEST_CASES[i] for i in args.reactions if i < len(TEST_CASES)]

    print(f"\n{C.BOLD}{C.HEADER}{'='*72}{C.END}")
    print(f"{C.BOLD}{C.HEADER}  ReasoningAgent Model Benchmark{C.END}")
    print(f"{C.BOLD}{C.HEADER}  {len(models_to_run)} models × {len(cases)} reactions{C.END}")
    print(f"{C.BOLD}{C.HEADER}{'='*72}{C.END}\n")

    all_results: list[Result] = []

    for model, provider in models_to_run:
        print(f"{C.BOLD}{C.CYAN}── {model} ({provider}) ──{C.END}")
        for case in cases:
            print(f"  {C.BLUE}→ {case['name']}{C.END}", flush=True)
            r = run_one(model, provider, case)
            all_results.append(r)
            print_result(r)
        print()

    print_summary_table(all_results, models_to_run)

    # Per-reaction accuracy across models
    print(f"{C.BOLD}Per-reaction accuracy:{C.END}")
    for case in cases:
        rows = [r for r in all_results if r.case_name == case["name"]]
        n_correct = sum(1 for r in rows if r.correct)
        col = C.GREEN if n_correct == len(rows) else (C.YELLOW if n_correct > 0 else C.RED)
        print(f"  {col}{n_correct}/{len(rows)}{C.END}  {case['name']}  "
              f"(expected: {case['expected']})")
    print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{C.CYAN}Interrupted.{C.END}")
        sys.exit(0)
