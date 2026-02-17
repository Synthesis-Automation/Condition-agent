"""
Interactive CLI tester for the ReactionPipeline (5-stage SMILES pipeline).

Shows each stage result in detail so you can see exactly what the pipeline
does and whether the LLM fallback fires.

Usage (no LLM — deterministic stages only):
    python test_pipeline_interactive.py

Usage (with LLM fallback, needs OPENAI_API_KEY):
    python test_pipeline_interactive.py --llm
    python test_pipeline_interactive.py --llm --model gpt-4o-mini
"""

import argparse
import os
import sys
import time
from pathlib import Path

# Ensure project root is on sys.path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))


# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------
class C:
    HEADER  = '\033[95m'
    BLUE    = '\033[94m'
    CYAN    = '\033[96m'
    GREEN   = '\033[92m'
    YELLOW  = '\033[93m'
    RED     = '\033[91m'
    BOLD    = '\033[1m'
    END     = '\033[0m'

    @classmethod
    def disable(cls):
        for attr in ('HEADER','BLUE','CYAN','GREEN','YELLOW','RED','BOLD','END'):
            setattr(cls, attr, '')


def hdr(title: str):
    print(f"\n{C.BOLD}{C.CYAN}{'='*72}{C.END}")
    print(f"{C.BOLD}{C.CYAN}  {title}{C.END}")
    print(f"{C.BOLD}{C.CYAN}{'='*72}{C.END}")

def sub(msg: str):
    print(f"{C.BLUE}  → {msg}{C.END}")

def ok(msg: str):
    print(f"{C.GREEN}  ✓ {msg}{C.END}")

def warn(msg: str):
    print(f"{C.YELLOW}  ⚠ {msg}{C.END}")

def err(msg: str):
    print(f"{C.RED}  ✗ {msg}{C.END}")

def kv(key: str, val):
    print(f"  {C.BOLD}{key}:{C.END} {val}")


# ---------------------------------------------------------------------------
# Preset reactions
# ---------------------------------------------------------------------------
PRESETS = {
    '1': ("Suzuki-Miyaura coupling",
          "Brc1ccccc1.B(O)(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1"),
    '2': ("Amide formation",
          "CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1"),
    '3': ("C-N coupling (Buchwald-Hartwig)",
          "Brc1ccccc1.NCc1ccccc1>>c1ccc(NCc2ccccc2)cc1"),
    '4': ("Complex Suzuki (with protecting group)",
          "CC1(C)OB(c2cnn(CCOC3CCCCO3)c2)OC1(C)C."
          "Cc1nc(-c2cn3c(n2)-c2ccc(Br)cc2OCC3)n(C)n1>>"
          "Cc1nc(-c2cn3c(n2)-c2ccc(-c4cnn(CCO)c4)cc2OCC3)n(C)n1"),
    '5': ("Simple substitution (edge case)",
          "Brc1ccccc1.O>>Oc1ccccc1"),
}


def show_menu(use_llm: bool, model: str):
    print(f"\n{C.BOLD}{C.HEADER}{'='*72}{C.END}")
    print(f"{C.BOLD}{C.HEADER}  ReactionPipeline Interactive Tester{C.END}")
    mode = f"LLM fallback ON ({model})" if use_llm else "deterministic only (no LLM)"
    print(f"{C.BOLD}{C.HEADER}  Mode: {mode}{C.END}")
    print(f"{C.BOLD}{C.HEADER}{'='*72}{C.END}")
    print("\n  Select a reaction:\n")
    for k, (name, smiles) in PRESETS.items():
        print(f"  {C.BOLD}{k}.{C.END} {name}")
        print(f"     {smiles[:78]}{'...' if len(smiles)>78 else ''}")
    print(f"\n  {C.BOLD}6.{C.END} Custom SMILES")
    print(f"  {C.BOLD}0.{C.END} Exit\n")


# ---------------------------------------------------------------------------
# Stage display helpers
# ---------------------------------------------------------------------------

def show_stage1(norm):
    hdr("STAGE 1 — Normalization & Pre-check")
    if norm.success:
        ok("Normalization passed")
        kv("Normalized SMILES", norm.normalized_smiles)
        kv("Reactants", norm.reactants)
        kv("Product", norm.product)
        if norm.warnings:
            for w in norm.warnings:
                warn(w)
    else:
        err(f"Normalization FAILED: {norm.error}")
        for w in norm.warnings:
            warn(w)


def show_stage2(feat):
    hdr("STAGE 2 — Deterministic Featurization")
    if feat.success:
        ok("featurize_reaction() succeeded")
        kv("Reaction type",   feat.reaction_type or "None")
        kv("Confidence",      f"{feat.reaction_type_confidence:.3f}")
        kv("Reacted motifs",  list(feat.reacted_motifs))
        kv("Formed motifs",   list(feat.formed_motifs))
        kv("Spectator motifs",list(feat.spectator_motifs))
        kv("Reaction key",    (feat.reaction_key[:80] + "...") if feat.reaction_key and len(feat.reaction_key) > 80 else feat.reaction_key)
        kv("Unclassified reactant", feat.has_unclassified_reactant)
        kv("Reactant motif count",  feat.reactant_motif_count)
        for w in feat.warnings:
            warn(w)
    else:
        err("featurize_reaction() FAILED")
        for w in feat.warnings:
            err(w)


def show_stage3(quality):
    hdr("STAGE 3 — Quality Gate")
    if quality.passed:
        ok("All quality criteria passed — no LLM fallback needed")
    else:
        warn(f"Quality gate FAILED ({len(quality.issues)} issue(s)) — LLM fallback will be triggered")
        for issue in quality.issues:
            err(issue)
    print(f"\n  {C.BOLD}Raw scores:{C.END}")
    for k, v in quality.scores.items():
        print(f"    {k}: {v}")


def show_stage4(llm_result, skipped_reason: str = ""):
    hdr("STAGE 4 — LLM Taxonomy Fallback")
    if skipped_reason:
        warn(f"Skipped: {skipped_reason}")
        return
    if llm_result is None:
        warn("Not triggered (quality gate passed)")
        return
    if llm_result.success:
        ok("LLM fallback succeeded")
        kv("Reaction type",   llm_result.reaction_type)
        kv("Confidence",      f"{llm_result.confidence:.3f}")
        kv("Reacted motifs",  list(llm_result.reacted_motifs))
        kv("Formed motifs",   list(llm_result.formed_motifs))
        kv("Reasoning",       llm_result.reasoning)
        if llm_result.invalid_ids_found:
            warn(f"Invalid IDs rejected: {llm_result.invalid_ids_found}")
    else:
        err("LLM fallback failed or produced incomplete output")
    for w in llm_result.warnings:
        warn(w)


def show_stage5(result):
    hdr("STAGE 5 — Merged Result")
    source = "LLM fallback" if result.used_llm_fallback else "featurize_reaction()"
    ok(f"Final fields sourced from: {source}")
    kv("reaction_type",           result.reaction_type or "None")
    kv("reaction_type_confidence", f"{result.reaction_type_confidence:.3f}")
    kv("reacted_motifs",          list(result.reacted_motifs))
    kv("formed_motifs",           list(result.formed_motifs))
    kv("spectator_motifs",        list(result.spectator_motifs))
    kv("reaction_key",            (result.reaction_key[:80] + "...") if result.reaction_key and len(result.reaction_key)>80 else result.reaction_key)
    kv("used_llm_fallback",       result.used_llm_fallback)
    if result.pipeline_warnings:
        print(f"\n  {C.BOLD}Pipeline warnings:{C.END}")
        for w in result.pipeline_warnings:
            warn(w)


# ---------------------------------------------------------------------------
# Main test runner
# ---------------------------------------------------------------------------

def run_pipeline(smiles: str, use_llm: bool, model: str):
    from reaction_agent.smiles_pipeline import ReactionPipeline
    from reaction_agent.pipeline_eval import QualityConfig

    print(f"\n{C.BOLD}Input SMILES:{C.END} {smiles}")

    # Build pipeline
    llm_client = None
    if use_llm:
        sub("Creating LLM client...")
        try:
            from llmtools.clients import LLMClient
            llm_client = LLMClient(provider="openai", model=model)
            ok(f"LLM client ready ({model})")
        except Exception as e:
            err(f"Could not create LLM client: {e}")
            err("Running in deterministic-only mode")

    pipeline = ReactionPipeline(llm_client=llm_client)

    t0 = time.time()

    # Run stage by stage (calling individual methods for visibility)
    norm = pipeline.normalize(smiles)
    show_stage1(norm)

    if not norm.success:
        err("Pipeline stopped at Stage 1")
        return

    feat = pipeline.featurize(norm.normalized_smiles)
    show_stage2(feat)

    quality = pipeline.evaluate(feat)
    show_stage3(quality)

    llm_result = None
    llm_skip = ""
    if quality.needs_llm_fallback:
        if llm_client is not None:
            sub("Running LLM fallback...")
            t_llm = time.time()
            llm_result = pipeline.llm_fallback(norm.normalized_smiles, feat, quality)
            ok(f"LLM fallback done ({time.time()-t_llm:.2f}s)")
        else:
            llm_skip = "No LLM client provided (run with --llm to enable)"
    show_stage4(llm_result, skipped_reason=llm_skip)

    result = pipeline.merge(norm, feat, quality, llm_result)
    show_stage5(result)

    elapsed = time.time() - t0
    hdr("SUMMARY")
    kv("Total time",       f"{elapsed:.2f}s")
    kv("Quality passed",   result.quality.passed if result.quality else "N/A")
    kv("LLM fallback used",result.used_llm_fallback)
    kv("Final rxn type",   result.reaction_type or "None")
    kv("Final motifs",     f"{list(result.reacted_motifs)} → {list(result.formed_motifs)}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive tester for the 5-stage ReactionPipeline"
    )
    parser.add_argument("--llm", action="store_true",
                        help="Enable LLM fallback (requires OPENAI_API_KEY)")
    parser.add_argument("--model", default="gpt-4o",
                        help="LLM model to use for fallback (default: gpt-4o)")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable terminal colours")
    args = parser.parse_args()

    if args.no_color:
        C.disable()

    if args.llm and not os.getenv("OPENAI_API_KEY"):
        print(f"{C.RED}ERROR: --llm requires OPENAI_API_KEY to be set{C.END}")
        sys.exit(1)

    while True:
        show_menu(args.llm, args.model)
        choice = input(f"{C.BOLD}Choice (0-6): {C.END}").strip()

        if choice == '0':
            print(f"\n{C.CYAN}Goodbye!{C.END}")
            break
        elif choice in PRESETS:
            name, smiles = PRESETS[choice]
            print(f"\n{C.BOLD}Selected:{C.END} {name}")
            run_pipeline(smiles, args.llm, args.model)
        elif choice == '6':
            smiles = input(f"\n{C.BOLD}Enter reaction SMILES (reactants>>product): {C.END}").strip()
            if not smiles:
                warn("No input — returning to menu")
                continue
            run_pipeline(smiles, args.llm, args.model)
        else:
            err("Invalid choice — enter 0-6")
            continue

        again = input(f"\n{C.BOLD}Test another? (y/n): {C.END}").strip().lower()
        if again not in ('y', 'yes', ''):
            print(f"\n{C.CYAN}Goodbye!{C.END}")
            break


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{C.CYAN}Interrupted. Goodbye!{C.END}")
        sys.exit(0)
