"""
Interactive CLI tester for the ReactionPipeline (5-stage SMILES pipeline).

Shows each stage result in detail so you can see exactly what the pipeline
does and whether the LLM fallback fires.

Usage:
    python test_pipeline_interactive.py
    python test_pipeline_interactive.py --no-llm
    python test_pipeline_interactive.py --model gpt-4o-mini
"""

import argparse
import os
import sys
import time
import warnings
from pathlib import Path

# Suppress noisy third-party deprecation warnings that clutter CLI output
warnings.filterwarnings("ignore", category=UserWarning, module="pkg_resources")
warnings.filterwarnings("ignore", message=".*pkg_resources.*")
warnings.filterwarnings("ignore", category=DeprecationWarning, module="rxnmapper")

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
# Model selector (mirrors reaction_agent/cli.py)
# ---------------------------------------------------------------------------

SELECTABLE_MODELS = [
    {"name": "o4-mini",       "provider": "openai"},   # 1 - default
    {"name": "gpt-5.2",       "provider": "openai"},   # 2
    {"name": "glm-5",         "provider": "aliyun"},   # 3
    {"name": "glm-4.7",       "provider": "aliyun"},   # 4
    {"name": "MiniMax-M2.1",  "provider": "aliyun"},   # 5
    {"name": "deepseek-v3.2", "provider": "aliyun"},   # 6
]

_ALIYUN_MODELS = {m["name"] for m in SELECTABLE_MODELS if m["provider"] == "aliyun"}


def select_model_interactive() -> tuple:
    """Print numbered model menu and return (model_name, provider)."""
    print("\nSelect LLM model:")
    for i, m in enumerate(SELECTABLE_MODELS, 1):
        tag = "  ← default" if i == 1 else ""
        print(f"  [{i}] {m['name']:<18} ({m['provider']}){tag}")

    while True:
        try:
            raw = input("\nEnter number [1]: ").strip()
        except (EOFError, KeyboardInterrupt):
            raw = ""

        if raw == "":
            choice = 1
        else:
            try:
                choice = int(raw)
            except ValueError:
                print(f"  Please enter a number between 1 and {len(SELECTABLE_MODELS)}.")
                continue

        if 1 <= choice <= len(SELECTABLE_MODELS):
            selected = SELECTABLE_MODELS[choice - 1]
            return selected["name"], selected["provider"]
        print(f"  Please enter a number between 1 and {len(SELECTABLE_MODELS)}.")


def show_header(use_llm: bool, model: str, use_reasoning: bool = True):
    print(f"\n{C.BOLD}{C.HEADER}{'='*72}{C.END}")
    print(f"{C.BOLD}{C.HEADER}  ReactionPipeline Interactive Tester{C.END}")
    if use_reasoning:
        mode = f"Reasoning Agent ON ({model})"
    elif use_llm:
        mode = f"Simple LLM fallback ({model})"
    else:
        mode = "deterministic only (no LLM)"
    print(f"{C.BOLD}{C.HEADER}  Mode: {mode}{C.END}")
    print(f"{C.BOLD}{C.HEADER}{'='*72}{C.END}")


# ---------------------------------------------------------------------------
# Stage display helpers
# ---------------------------------------------------------------------------

def _smi(s: str, n: int = 60) -> str:
    """Truncate a SMILES string for display."""
    if not s:
        return "(none)"
    return s if len(s) <= n else s[:n] + "…"


def show_stage1(norm):
    hdr("STAGE 1 — Normalization & Pre-check")
    if norm.success:
        ok("Normalization passed")
        # Reactants on one line, product on next — skip full normalized SMILES
        # (already printed as Input SMILES above)
        rcts = "  |  ".join(_smi(r) for r in norm.reactants)
        kv("Reactants", rcts)
        kv("Product  ", _smi(norm.product))
        # Deduplicated warnings only
        seen: set = set()
        for w in (norm.warnings or []):
            if w not in seen:
                seen.add(w)
                warn(w)
    else:
        err(f"Normalization FAILED: {norm.error}")
        seen: set = set()
        for w in (norm.warnings or []):
            if w not in seen:
                seen.add(w)
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
        kv("Reaction key",    _smi(feat.reaction_key, 80))
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

    # Show rich reasoning profile if available
    profile = getattr(llm_result, 'reactivity_profile', None)
    if profile:
        show_reasoning_profile(profile)


def show_reasoning_profile(profile):
    """Display the rich ReactivityProfile from the reasoning agent."""
    print(f"\n  {C.BOLD}{C.HEADER}--- Reasoning Agent Analysis ---{C.END}")

    # Tools called
    if profile.tools_called:
        unique_tools = sorted(set(profile.tools_called))
        kv("Tools used", f"{len(profile.tools_called)} calls ({', '.join(unique_tools)})")

    # Named reaction (new)
    if profile.named_reaction:
        kv("Named reaction", f"{C.CYAN}{profile.named_reaction}{C.END}")

    # All component roles (new)
    if profile.all_roles:
        print(f"\n  {C.BOLD}Component roles:{C.END}")
        for frag, role in profile.all_roles.items():
            print(f"    {C.CYAN}{_smi(frag, 40)}{C.END}  →  {role}")

    # Electrophile
    e = profile.electrophile
    if e.hybridization or e.leaving_group:
        print(f"\n  {C.BOLD}Electrophile:{C.END}")
        if e.center_atom:
            kv("    Center", e.center_atom)
        if e.hybridization:
            kv("    Hybridization", e.hybridization)
        if e.leaving_group:
            kv("    Leaving group", f"{e.leaving_group} ({e.leaving_group_quality})")
        if e.electronic_class:
            kv("    Electronics", f"score={e.electronic_score:.1f} ({e.electronic_class})")
        if e.steric_class:
            kv("    Sterics", f"score={e.steric_score:.1f} ({e.steric_class})")
        if e.activation_tags:
            kv("    Activation tags", e.activation_tags)

    # Nucleophile
    n = profile.nucleophile
    if n.identity or n.attacking_atom:
        print(f"\n  {C.BOLD}Nucleophile:{C.END}")
        if n.identity:
            kv("    Identity", n.identity)
        if n.attacking_atom:
            kv("    Attacking atom", n.attacking_atom)
        if n.hardsoft:
            kv("    Hard/soft", n.hardsoft)
        if n.is_also_base:
            warn("    Also acts as base (elimination risk)")
        if n.steric_bulk:
            kv("    Steric bulk", n.steric_bulk)

    # Reaction pattern type
    if profile.reaction_pattern_type:
        kv("Pattern type", profile.reaction_pattern_type)

    # Mechanism
    m = profile.mechanism
    if m.primary_class:
        print(f"\n  {C.BOLD}Mechanism:{C.END}")
        kv("    Primary class", m.primary_class)
        kv("    Confidence", f"{m.confidence:.2f}")
        if m.requires_catalyst:
            kv("    Catalyst metals", m.likely_catalyst_metals)
        if m.key_intermediates:
            kv("    Key intermediates", m.key_intermediates)
        if m.stepwise:
            print(f"    {C.BOLD}Stepwise pathway:{C.END}")
            for step in m.stepwise:
                print(f"      {step}")
        if m.alternative_mechanisms:
            kv("    Alternatives", m.alternative_mechanisms)
        if m.evidence:
            print(f"    {C.BOLD}Evidence:{C.END}")
            for ev in m.evidence[:5]:
                print(f"      • {ev}")

    # Transformation
    t = profile.transformation
    if t.bonds_broken or t.bonds_formed:
        print(f"\n  {C.BOLD}Transformation:{C.END}")
        if t.key_bond_type:
            kv("    Key bond type", t.key_bond_type)
        if t.bonds_broken:
            kv("    Bonds broken", t.bonds_broken)
        if t.bonds_formed:
            kv("    Bonds formed", t.bonds_formed)
        if t.fg_removed:
            kv("    FGs removed", t.fg_removed)
        if t.fg_formed:
            kv("    FGs formed", t.fg_formed)
        if t.redox_change:
            kv("    Redox", t.redox_change)

    # Tandem / reactive streams (new)
    if profile.is_tandem:
        warn("Tandem/multi-step reaction detected")
        if profile.reactive_streams:
            print(f"\n  {C.BOLD}Reactive streams:{C.END}")
            for s in profile.reactive_streams:
                sid = s.get("stream_id", "?")
                desc = s.get("description", "")
                prod = s.get("product_fragment", "")
                print(f"    Stream {sid}: {desc}")
                if prod:
                    print(f"      → {prod}")

    # Product verification (new)
    pv = profile.product_verification
    if pv:
        score = pv.get("verification_score", "")
        confirmed = pv.get("confirmed_in_product", [])
        expected = pv.get("expected_motifs", [])
        score_color = C.GREEN if score == "high" else C.YELLOW if score == "medium" else C.RED
        print(f"\n  {C.BOLD}Product verification:{C.END} {score_color}{score}{C.END}")
        if expected != confirmed:
            missing = set(expected) - set(confirmed)
            if missing:
                warn(f"  Expected but not confirmed: {list(missing)}")

    # Missing conditions (new)
    if profile.missing_conditions:
        print(f"\n  {C.BOLD}Missing conditions (inferred):{C.END}")
        for mc in profile.missing_conditions:
            warn(f"  {mc}")

    # Extended analysis
    if profile.selectivity_risks:
        print(f"\n  {C.BOLD}Selectivity risks:{C.END}")
        for risk in profile.selectivity_risks:
            warn(f"  {risk}")

    if profile.condition_implications:
        print(f"\n  {C.BOLD}Condition implications:{C.END}")
        for ck, cv in profile.condition_implications.items():
            kv(f"    {ck}", cv)

    print(f"  {C.BOLD}{C.HEADER}--- End Reasoning Profile ---{C.END}")


def show_stage5(result):
    hdr("STAGE 5 — Merged Result")
    source = "LLM fallback" if result.used_llm_fallback else "featurize_reaction()"
    ok(f"Final fields sourced from: {source}")
    kv("reaction_type",           result.reaction_type or "None")
    kv("reaction_type_confidence", f"{result.reaction_type_confidence:.3f}")
    kv("reacted_motifs",          list(result.reacted_motifs))
    kv("formed_motifs",           list(result.formed_motifs))
    kv("spectator_motifs",        list(result.spectator_motifs))
    kv("reaction_key",            _smi(result.reaction_key, 80))
    kv("used_llm_fallback",       result.used_llm_fallback)
    if result.pipeline_warnings:
        print(f"\n  {C.BOLD}Pipeline warnings:{C.END}")
        for w in result.pipeline_warnings:
            warn(w)


def show_stage6(recs, top_k: int, elapsed: float):
    if recs is None:
        return

    kv("Detected reaction type", recs.predicted_reaction_type or "Unknown")
    kv("Reaction type confidence", f"{recs.reaction_type_confidence:.3f}")
    kv("Reacted motifs (recommender)", list(recs.reacted_motifs or []))
    kv("Formed motifs  (recommender)", list(recs.formed_motifs or []))
    kv("Total matching experiments", recs.total_matching_experiments)
    kv("Fallback match", recs.is_fallback_match)
    kv("Query time", f"{elapsed:.3f}s")

    if not recs.recommendations:
        warn("No recommendations found")
        print("  Possible reasons:")
        print("    - Reaction type / motifs not represented in HTE experiments database")
        print("    - Minimum experiment count threshold not met")
        return

    print(f"\n  {C.BOLD}Top {min(top_k, len(recs.recommendations))} condition sets:{C.END}\n")

    for i, rec in enumerate(recs.recommendations[:top_k], 1):
        if rec.avg_z_score > 1.0:
            score_col, label = C.GREEN, "Excellent"
        elif rec.avg_z_score > 0.0:
            score_col, label = C.YELLOW, "Good"
        else:
            score_col, label = C.RED, "Poor"

        print(f"  {C.BOLD}Rank {i}:{C.END}  "
              f"{score_col}Z={rec.avg_z_score:.2f} ({label}){C.END}  "
              f"conf={rec.confidence_score:.0f}%  n={rec.num_experiments}")

        parts = []
        if rec.catalyst:  parts.append(f"cat={rec.catalyst}")
        if rec.ligand:    parts.append(f"lig={rec.ligand}")
        if rec.base:      parts.append(f"base={rec.base}")
        if rec.solvent:   parts.append(f"solv={rec.solvent}")
        if rec.additive:  parts.append(f"add={rec.additive}")
        print(f"           {' | '.join(parts)}")
        print(f"           success={rec.success_rate:.1f}%  "
              f"avg_yield={rec.avg_yield:.1f}%  "
              f"median_yield={rec.median_yield:.1f}%")


# ---------------------------------------------------------------------------
# Lazy HTERecommender cache — loaded at most once per session, on demand
# ---------------------------------------------------------------------------
_REC_CACHE: dict = {"instance": None}

def _get_recommender(db_path: str):
    """Return a cached HTERecommender, loading it on first call."""
    if _REC_CACHE["instance"] is None:
        print(f"\n  Loading HTE database from: {db_path}")
        print(f"  (first load may take a moment; cached for the rest of the session)")
        t0 = time.time()
        from chemtools.recommend.recommender import HTERecommender
        _REC_CACHE["instance"] = HTERecommender(db_path)
        ok(f"HTE database ready ({time.time()-t0:.1f}s)")
    return _REC_CACHE["instance"]


# ---------------------------------------------------------------------------
# Main test runner
# ---------------------------------------------------------------------------

def run_pipeline(
    smiles: str,
    use_llm: bool,
    model: str,
    top_k: int,
    db_path: str,
    use_reasoning: bool = True,
    reasoning_model: str = None,
    provider: str = "openai",
):
    from reaction_agent.smiles_pipeline import ReactionPipeline

    print(f"\n{C.BOLD}Input SMILES:{C.END} {smiles}")

    # Build pipeline
    llm_client = None
    if use_llm and not use_reasoning:
        sub("Creating LLM client...")
        try:
            from llmtools.clients import LLMClient
            llm_client = LLMClient(provider=provider, model=model)
            ok(f"LLM client ready ({model})")
        except Exception as e:
            err(f"Could not create LLM client: {e}")
            err("Running in deterministic-only mode")

    pipeline = ReactionPipeline(
        llm_client=llm_client,
        use_reasoning_agent=use_reasoning,
        reasoning_provider=provider,
        reasoning_model=reasoning_model or model,
    )

    t0 = time.time()

    # Stage 1 — normalize
    norm = pipeline.normalize(smiles)
    show_stage1(norm)
    if not norm.success:
        err("Pipeline stopped at Stage 1")
        return

    # Stage 2 — featurize
    feat = pipeline.featurize(norm.normalized_smiles)
    show_stage2(feat)

    # Stage 3 — quality gate
    quality = pipeline.evaluate(feat)
    show_stage3(quality)

    # Stage 4 — LLM fallback (conditional)
    llm_result = None
    llm_skip = ""
    if quality.needs_llm_fallback:
        if use_reasoning:
            sub(f"Running reasoning agent ({reasoning_model or model})...")
            t_llm = time.time()
            llm_result = pipeline.reasoning_fallback(norm.normalized_smiles, feat, quality)
            ok(f"Reasoning agent done ({time.time()-t_llm:.2f}s)")
        elif llm_client is not None:
            sub("Running LLM fallback...")
            t_llm = time.time()
            llm_result = pipeline.llm_fallback(norm.normalized_smiles, feat, quality)
            ok(f"LLM fallback done ({time.time()-t_llm:.2f}s)")
        else:
            llm_skip = "No LLM client or reasoning agent"
    show_stage4(llm_result, skipped_reason=llm_skip)

    # Stage 5 — merge
    result = pipeline.merge(norm, feat, quality, llm_result)
    show_stage5(result)

    # Stage 6 — HTE recommender (on demand)
    recs = None
    rec_elapsed = 0.0
    reactants = norm.reactants
    product = norm.product
    hdr(f"STAGE 6 — HTE Condition Recommendations (experiments, top {top_k})")
    if not product or not reactants:
        warn("Stage 6 skipped: no reactants/product available")
    else:
        ans = input(f"\n{C.BOLD}  Run HTE recommendation? (y/n): {C.END}").strip().lower()
        if ans in ('y', 'yes'):
            try:
                rec = _get_recommender(db_path)
                sub(f"Querying HTE experiments database (top {top_k})...")
                t_rec = time.time()
                recs = rec.recommend(
                    reactant_a_smiles=reactants[0],
                    reactant_b_smiles=reactants[1] if len(reactants) > 1 else None,
                    product_smiles=product,
                    top_k=top_k,
                    min_experiments=2,
                    reaction_type_filter=result.reaction_type,
                    source_group="experiments",
                )
                rec_elapsed = time.time() - t_rec
            except Exception as e:
                err(f"HTERecommender raised: {e}")
        else:
            warn("Stage 6 skipped by user")
    show_stage6(recs, top_k, rec_elapsed)

    elapsed = time.time() - t0
    hdr("SUMMARY")
    kv("Total time",           f"{elapsed:.2f}s")
    kv("Quality passed",       result.quality.passed if result.quality else "N/A")
    kv("LLM fallback used",    result.used_llm_fallback)
    if result.reactivity_profile:
        kv("Reasoning agent",  f"YES ({result.reactivity_profile.total_tool_calls} tool calls)")
        kv("Mechanism",        result.reactivity_profile.mechanism.primary_class or "N/A")
    kv("Final rxn type",       result.reaction_type or "None")
    kv("Final motifs",         f"{list(result.reacted_motifs)} → {list(result.formed_motifs)}")
    kv("Recommendations found", len(recs.recommendations) if recs else 0)
    kv("Total experiments",     recs.total_matching_experiments if recs else 0)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive tester for the 5-stage ReactionPipeline"
    )
    parser.add_argument("--no-llm", action="store_true",
                        help="Disable ALL LLM fallback (deterministic stages only)")
    parser.add_argument("--no-reasoning", action="store_true",
                        help="Disable reasoning agent, use simple prompt fallback instead")
    parser.add_argument("--model", default=None,
                        help="LLM model to use (default: interactive selection)")
    parser.add_argument("--provider", default=None,
                        help="LLM provider (default: auto-detected from model)")
    parser.add_argument("--reasoning-model", default=None,
                        help="Override model for reasoning agent (default: same as --model)")
    parser.add_argument("--db-path", default="data/HTE_db",
                        help="Path to HTE database directory (default: data/HTE_db)")
    parser.add_argument("--top-k", type=int, default=5,
                        help="Number of recommendations to show (default: 5)")
    parser.add_argument("--no-color", action="store_true",
                        help="Disable terminal colours")
    args = parser.parse_args()

    if args.no_color:
        C.disable()

    # Model / provider selection
    if args.no_llm:
        args.model = args.model or "none"
        args.provider = args.provider or "none"
    elif args.model is None:
        args.model, args.provider = select_model_interactive()
    else:
        if args.provider is None:
            args.provider = "aliyun" if args.model in _ALIYUN_MODELS else "openai"

    use_llm = not args.no_llm
    use_reasoning = use_llm and not args.no_reasoning

    # Check the right API key for the selected provider
    api_key_env = f"{args.provider.upper()}_API_KEY"
    if use_llm and not os.getenv(api_key_env):
        print(f"{C.YELLOW}WARNING: {api_key_env} not set — LLM fallback disabled{C.END}")
        use_llm = False
        use_reasoning = False

    show_header(use_llm, args.model, use_reasoning=use_reasoning)

    while True:
        try:
            smiles = input(f"\n{C.BOLD}Reaction SMILES (q to exit): {C.END}").strip()
        except EOFError:
            break

        if smiles.lower() in ('q', 'quit', 'exit'):
            print(f"\n{C.CYAN}Goodbye!{C.END}")
            break
        if not smiles:
            continue

        run_pipeline(
            smiles, use_llm, args.model, args.top_k, args.db_path,
            use_reasoning=use_reasoning,
            reasoning_model=args.reasoning_model,
            provider=args.provider,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(f"\n\n{C.CYAN}Interrupted. Goodbye!{C.END}")
        sys.exit(0)
