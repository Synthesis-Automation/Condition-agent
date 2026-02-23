"""
CLI for HTE Analytics Tools

Provides command-line access to HTE database analytics:
- List reactant pairs
- Analyze catalysts
- View reaction type summaries
- Export filtered datasets
- Backtest deterministic HTE recommender on held-out rows

Usage (examples)
----------------
Run a reaction-type summary on the default database (`data/HTE_db`):

    python app/HTE_analytics_cli.py reactions --top 30 --compact

Save a CSV result and automatically generate a Markdown sidecar summary:

    python app/HTE_analytics_cli.py reactions --top 50 -o results/reaction_summary.csv

This writes:
- `results/reaction_summary.csv`
- `results/reaction_summary.md`  (generated automatically)

List reactant pairs with filters and save CSV + Markdown:

    python app/HTE_analytics_cli.py pairs --reaction C_S_Coupling --catalyst Pd --top 20 -o results/cs_pairs.csv

Export a filtered subset (CSV + Markdown export summary):

    python app/HTE_analytics_cli.py export results/pd_cn_subset.csv --reaction C_N_Coupling --catalyst Pd

Backtest the recommender and save JSON + Markdown report:

    python app/HTE_analytics_cli.py backtest --input data/HTE_db/experiments/HTE_canonical.csv --output results/backtest.json

This writes:
- `results/backtest.json`
- `results/backtest.md`  (generated automatically)
"""

import argparse
import json
import random
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import pandas as pd

# Ensure repo root is on sys.path for direct execution.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chemtools.recommend import HTERecommender
from chemtools.recommend.analytics import HTEAnalytics


def _markdown_sidecar_path(output_path: str | Path) -> Path:
    path = Path(output_path)
    return path.with_suffix(".md") if path.suffix else Path(f"{path}.md")


def _md_escape(value: Any) -> str:
    text = "" if value is None else str(value)
    return text.replace("|", r"\|").replace("\n", " ").strip()


def _dataframe_to_markdown(df: pd.DataFrame, *, max_rows: Optional[int] = None) -> str:
    table = df.head(max_rows).copy() if max_rows is not None else df.copy()
    columns = [str(c) for c in table.columns]
    if not columns:
        return "_No columns_"

    rows: List[List[str]] = []
    for _, row in table.iterrows():
        rows.append([_md_escape(row[col]) for col in table.columns])

    header = "| " + " | ".join(_md_escape(c) for c in columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    body = "\n".join("| " + " | ".join(r) + " |" for r in rows) if rows else ""
    return "\n".join([header, divider, body]).rstrip()


def _write_markdown_sidecar_for_dataframe(
    output_path: str | Path,
    *,
    title: str,
    df: pd.DataFrame,
    top: Optional[int] = None,
    notes: Optional[List[str]] = None,
) -> Path:
    md_path = _markdown_sidecar_path(output_path)
    lines = [f"# {title}", ""]
    if notes:
        lines.extend([line for line in notes if line])
        lines.append("")
    lines.append(f"Rows: {len(df):,}")
    if top is not None and len(df) > top:
        lines.append(f"Showing top {top:,} rows.")
    lines.append("")
    lines.append(_dataframe_to_markdown(df, max_rows=top))
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return md_path


def _write_markdown_sidecar_for_backtest_report(output_path: str | Path, report: Dict[str, Any]) -> Path:
    md_path = _markdown_sidecar_path(output_path)
    metrics = report.get("metrics") or {}
    diag = report.get("diagnostics") or {}
    settings = report.get("settings") or {}
    lines = [
        "# HTE Recommender Backtest",
        "",
        "## Dataset",
        "",
        f"- Input CSV: `{report.get('input_csv', '')}`",
        f"- Rows after filters: {int(report.get('rows_after_filters', 0)):,}",
        f"- Train rows: {int(report.get('train_rows', 0)):,}",
        f"- Test rows: {int(report.get('test_rows', 0)):,}",
        f"- Evaluated rows: {int(report.get('evaluated_rows', 0)):,}",
        "",
        "## Metrics",
        "",
        f"- Query coverage: {100.0 * float(metrics.get('query_coverage', 0.0)):.1f}%",
        f"- MRR: {float(metrics.get('mrr', 0.0)):.4f}",
        f"- Avg rank (hits only): {metrics.get('avg_rank') if metrics.get('avg_rank') is not None else 'N/A'}",
    ]
    for key in sorted(metrics):
        if key.startswith("hit@"):
            lines.append(f"- {key}: {100.0 * float(metrics.get(key, 0.0)):.1f}%")
    lines.extend(
        [
            "",
            "## Diagnostics",
            "",
            f"- Avg DB coverage: {float(diag.get('avg_database_coverage_pct', 0.0)):.2f}%",
            f"- Avg matching experiments: {float(diag.get('avg_matching_experiments', 0.0)):.2f}",
            f"- No recommendation count: {int(report.get('no_recommendation_count', 0))}",
            f"- Skipped invalid query: {int(report.get('skipped_invalid_query', 0))}",
            f"- Skipped empty condition: {int(report.get('skipped_empty_condition', 0))}",
            f"- Skipped unseen condition: {int(report.get('skipped_unseen_condition', 0))}",
            "",
            "## Settings",
            "",
        ]
    )
    for key in sorted(settings):
        lines.append(f"- {key}: `{settings[key]}`")
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return md_path


def _prompt(text: str, default: str | None = None) -> str:
    if default:
        prompt = f"{text} [{default}]: "
    else:
        prompt = f"{text}: "
    value = input(prompt).strip()
    return value if value else (default or "")


def _prompt_int(text: str, default: int) -> int:
    value = _prompt(text, str(default))
    try:
        return int(value)
    except ValueError:
        return default


def _prompt_float(text: str, default: float | None = None) -> float | None:
    value = _prompt(text, "" if default is None else str(default))
    if not value:
        return default
    try:
        return float(value)
    except ValueError:
        return default


def _prompt_yes_no(text: str, default: bool = False) -> bool:
    suffix = "Y/n" if default else "y/N"
    value = _prompt(f"{text} ({suffix})", "y" if default else "n").lower()
    if value in {"y", "yes"}:
        return True
    if value in {"n", "no"}:
        return False
    return default


def _run_default_summary(db_path: str) -> int:
    args_reactions = argparse.Namespace(
        db_path=db_path,
        top=None,
        compact=True,
        output=None,
    )
    args_metals = argparse.Namespace(
        db_path=db_path,
        detailed=False,
        output=None,
    )
    cmd_reactions(args_reactions)
    cmd_metals(args_metals)
    return 0


def _run_wizard(db_path: str) -> int:
    print("\nInteractive HTE Analytics Wizard")
    print("=" * 40)
    db_path = _prompt("HTE database path", db_path)

    menu = (
        "\nChoose a command:\n"
        "  1) List reactant pairs\n"
        "  2) Analyze catalysts\n"
        "  3) Reaction type summary\n"
        "  4) Metal usage\n"
        "  5) Export filtered dataset\n"
        "  6) Backtest HTE recommender\n"
        "  q) Quit\n"
    )
    while True:
        choice = _prompt(menu, "1").lower()
        if choice in {"q", "quit", "exit"}:
            return 0
        if choice == "1":
            args = argparse.Namespace(
                db_path=db_path,
                reaction=_prompt("Reaction type filter", ""),
                catalyst=_prompt("Catalyst filter", ""),
                min_experiments=_prompt_int("Minimum experiments", 1),
                sort=_prompt("Sort by (count/success_rate)", "count"),
                top=_prompt_int("Top N results", 20),
                compact=_prompt_yes_no("Compact output", True),
                output=_prompt("Save CSV to (optional)", ""),
            )
            if not args.output:
                args.output = None
            cmd_list_pairs(args)
        elif choice == "2":
            args = argparse.Namespace(
                db_path=db_path,
                reaction=_prompt("Reaction type filter", ""),
                reactant_a=_prompt("Reactant A type", ""),
                reactant_b=_prompt("Reactant B type", ""),
                top=_prompt_int("Top N results", 20),
                compact=_prompt_yes_no("Compact output", True),
                output=_prompt("Save CSV to (optional)", ""),
            )
            if not args.output:
                args.output = None
            cmd_catalysts(args)
        elif choice == "3":
            args = argparse.Namespace(
                db_path=db_path,
                reaction=_prompt("Reaction type filter (optional)", ""),
                detailed_map=_prompt_yes_no("Include detailed map in CSV/report", False),
                detail_top=_prompt_int("Detailed map top N pairs/catalysts", 5),
                top=_prompt_int("Top N results", 20),
                compact=_prompt_yes_no("Compact output", True),
                output=_prompt("Save CSV to (optional)", ""),
            )
            if not args.reaction:
                args.reaction = None
            if not args.output:
                args.output = None
            cmd_reactions(args)
        elif choice == "4":
            args = argparse.Namespace(
                db_path=db_path,
                detailed=_prompt_yes_no("Detailed breakdown", False),
                output=_prompt("Save CSV to (optional)", ""),
            )
            if not args.output:
                args.output = None
            cmd_metals(args)
        elif choice == "5":
            output = _prompt("Output CSV path", "")
            if not output:
                print("Output path is required.")
                continue
            args = argparse.Namespace(
                db_path=db_path,
                output=output,
                reaction=_prompt("Reaction type filter", ""),
                catalyst=_prompt("Catalyst filter", ""),
                reactant_a=_prompt("Reactant A type", ""),
                reactant_b=_prompt("Reactant B type", ""),
                min_yield=_prompt_float("Minimum yield", None),
            )
            cmd_export(args)
        elif choice == "6":
            default_input = db_path
            db_path_obj = Path(db_path)
            if db_path_obj.is_dir():
                default_input = str(db_path_obj / "experiments" / "HTE_canonical.csv")
            args = argparse.Namespace(
                input=_prompt("Input CSV path for backtest", default_input),
                reaction=_prompt("Reaction type filter (optional)", ""),
                catalyst=_prompt("Catalyst filter (optional)", ""),
                source_group=_prompt("Source group filter (optional)", ""),
                test_fraction=_prompt_float("Test fraction", 0.2) or 0.2,
                test_size=_prompt_int("Exact test size (0 to disable)", 0),
                seed=_prompt_int("Random seed", 13),
                top_k=_prompt_int("Recommendation top-k", 10),
                hit_ks=_prompt("Hit@k list (comma-separated)", "1,3,5,10"),
                min_experiments=_prompt_int("Minimum experiments per condition", 1),
                reaction_key_only=_prompt_yes_no("Use reaction-key-only matching", False),
                use_spectator_groups=_prompt_yes_no("Use spectator group weighting", True),
                allow_unseen_conditions=_prompt_yes_no("Include unseen train conditions as misses", False),
                train_output=_prompt("Save train split CSV path (optional)", ""),
                output=_prompt(
                    "Output JSON path",
                    "results/hte_recommender_backtest.json",
                ),
                per_row_output=_prompt("Per-row CSV output path (optional)", ""),
            )
            if not args.reaction:
                args.reaction = None
            if not args.catalyst:
                args.catalyst = None
            if not args.source_group:
                args.source_group = None
            if int(args.test_size) <= 0:
                args.test_size = None
            if not args.train_output:
                args.train_output = None
            if not args.per_row_output:
                args.per_row_output = None
            cmd_backtest(args)
        else:
            print("Invalid choice. Try again.")


def _parse_reaction_smiles(reaction_smiles: str) -> Tuple[str, Optional[str], Optional[str]]:
    text = str(reaction_smiles or "").strip()
    if not text:
        return "", None, None
    if ">>" in text:
        reactants_part, product = text.split(">>", 1)
        reactants = [r for r in reactants_part.split(".") if r]
        reactant_a = reactants[0] if reactants else ""
        reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
        product_smiles = product.strip() or None
        return reactant_a, reactant_b, product_smiles
    reactants = [r for r in text.split(".") if r]
    reactant_a = reactants[0] if reactants else ""
    reactant_b = ".".join(reactants[1:]) if len(reactants) > 1 else None
    return reactant_a, reactant_b, None


def _normalize_text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return ""
    return text


def _first_present_value(row: pd.Series, candidates: Iterable[str]) -> Any:
    for col in candidates:
        if col in row.index:
            return row[col]
    return None


def _extract_query_from_row(row: pd.Series) -> Tuple[str, Optional[str], Optional[str], str]:
    reaction_smiles = _normalize_text(
        _first_present_value(
            row,
            ("reaction_smiles", "Reaction_SMILES", "reactionSmiles"),
        )
    )
    if reaction_smiles:
        reactant_a, reactant_b, product = _parse_reaction_smiles(reaction_smiles)
        return reactant_a, reactant_b, product, reaction_smiles

    reactant_a = _normalize_text(
        _first_present_value(
            row,
            ("reactant_a_smiles", "Reactant_A_SMILES", "reactant_a", "Reactant_A"),
        )
    )
    reactant_b = _normalize_text(
        _first_present_value(
            row,
            ("reactant_b_smiles", "Reactant_B_SMILES", "reactant_b", "Reactant_B"),
        )
    )
    product = _normalize_text(
        _first_present_value(
            row,
            ("product_smiles", "Product_SMILES", "product", "Product"),
        )
    )
    return reactant_a, (reactant_b or None), (product or None), ""


def _condition_key(
    catalyst: Any,
    ligand: Any,
    base: Any,
    solvent: Any,
    secondary_solvent: Any,
    additive: Any,
    coupling_reagent: Any,
) -> Tuple[str, str, str, str, str, str, str]:
    return (
        _normalize_text(catalyst),
        _normalize_text(ligand),
        _normalize_text(base),
        _normalize_text(solvent),
        _normalize_text(secondary_solvent),
        _normalize_text(additive),
        _normalize_text(coupling_reagent),
    )


def _row_condition_key(row: pd.Series) -> Tuple[str, str, str, str, str, str, str]:
    return _condition_key(
        _first_present_value(row, ("Catalyst", "catalyst")),
        _first_present_value(row, ("Ligand", "ligand")),
        _first_present_value(row, ("Base", "base")),
        _first_present_value(row, ("Solvent", "solvent")),
        _first_present_value(row, ("Secondary Solvent", "secondary_solvent", "Secondary_Solvent")),
        _first_present_value(row, ("Additive", "additive")),
        _first_present_value(row, ("Coupling Reagent", "coupling_reagent", "Coupling_Reagent")),
    )


def _rec_condition_key(rec: Any) -> Tuple[str, str, str, str, str, str, str]:
    return _condition_key(
        getattr(rec, "catalyst", ""),
        getattr(rec, "ligand", ""),
        getattr(rec, "base", ""),
        getattr(rec, "solvent", ""),
        getattr(rec, "secondary_solvent", ""),
        getattr(rec, "additive", ""),
        getattr(rec, "coupling_reagent", ""),
    )


def _rank_for_key(recs: Iterable[Any], target_key: Tuple[str, str, str, str, str, str, str]) -> Optional[int]:
    for idx, rec in enumerate(recs, start=1):
        if _rec_condition_key(rec) == target_key:
            return idx
    return None


def _mean_or_zero(values: List[float]) -> float:
    if not values:
        return 0.0
    return float(sum(values) / len(values))


def _compute_backtest_metrics(ranks: List[Optional[int]], ks: Iterable[int]) -> Dict[str, Any]:
    total = len(ranks)
    if total == 0:
        out: Dict[str, Any] = {
            "count": 0,
            "query_coverage": 0.0,
            "mrr": 0.0,
            "avg_rank": None,
        }
        for k in ks:
            out[f"hit@{k}"] = 0.0
        return out

    hits = [r for r in ranks if r is not None]
    out = {
        "count": total,
        "query_coverage": float(len(hits) / total),
        "mrr": float(sum(1.0 / r for r in hits) / total) if total else 0.0,
        "avg_rank": float(sum(hits) / len(hits)) if hits else None,
    }
    for k in ks:
        out[f"hit@{k}"] = float(sum(1 for r in hits if r <= k) / total)
    return out


def _prepare_split(
    df: pd.DataFrame,
    *,
    test_fraction: float,
    test_size: Optional[int],
    seed: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    indices = list(df.index)
    if not indices:
        return df.iloc[0:0].copy(), df.iloc[0:0].copy()
    rng = random.Random(seed)
    rng.shuffle(indices)
    if test_size is None:
        test_count = max(1, int(len(indices) * test_fraction))
    else:
        test_count = min(max(1, int(test_size)), len(indices))
    test_set = set(indices[:test_count])
    test_df = df.loc[list(test_set)].copy()
    train_df = df.loc[[idx for idx in indices if idx not in test_set]].copy()
    return train_df, test_df


def _filter_rows_by_contains(
    df: pd.DataFrame,
    *,
    value: Optional[str],
    candidates: Iterable[str],
) -> pd.DataFrame:
    text = _normalize_text(value)
    if not text:
        return df
    for col in candidates:
        if col in df.columns:
            mask = df[col].fillna("").astype(str).str.contains(text, case=False, regex=False)
            return df[mask].copy()
    return df


def cmd_list_pairs(args):
    """List reactant pairs command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("REACTANT PAIR ANALYSIS")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.catalyst:
        print(f"Catalyst Filter: {args.catalyst}")
    if args.min_experiments > 1:
        print(f"Min Experiments: {args.min_experiments}")
    print()
    
    df = analytics.list_reactant_pairs(
        reaction_type=args.reaction,
        catalyst_filter=args.catalyst,
        min_experiments=args.min_experiments,
        sort_by=args.sort
    )
    
    if len(df) == 0:
        print("No matching reactant pairs found")
        return
    
    print(f"Found {len(df)} reactant pair combinations\n")
    
    # Display results
    if args.compact:
        # Compact format
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Reactant_A_Type']} + {row['Reactant_B_Type']}")
            print(f"   Reaction: {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            print(f"   Top Catalyst: {row['Top_Catalyst']}")
            print()
    else:
        # Full table format
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 50
        }
        
        import pandas as pd
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved results to {args.output}")
        md_path = _write_markdown_sidecar_for_dataframe(
            args.output,
            title="HTE Reactant Pair Analysis",
            df=df,
            top=args.top,
            notes=[
                f"Reaction filter: {args.reaction or 'None'}",
                f"Catalyst filter: {args.catalyst or 'None'}",
                f"Min experiments: {args.min_experiments}",
                f"Sort by: {args.sort}",
            ],
        )
        print(f"Saved markdown summary to {md_path}")


def cmd_catalysts(args):
    """Analyze catalysts command"""
    analytics = HTEAnalytics(args.db_path)
    import pandas as pd
    
    print("\n" + "="*80)
    print("CATALYST ANALYSIS")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.reactant_a:
        print(f"Reactant A Type: {args.reactant_a}")
    if args.reactant_b:
        print(f"Reactant B Type: {args.reactant_b}")
    print()
    
    df = analytics.get_catalyst_stats(
        reaction_type=args.reaction,
        reactant_a_type=args.reactant_a,
        reactant_b_type=args.reactant_b
    )
    
    if len(df) == 0:
        print("No catalysts found matching criteria")
        return
    
    print(f"Found {len(df)} catalysts\n")
    
    # Display results
    if args.compact:
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Catalyst']} [{row['Metal']}]")
            print(f"   Experiments: {row['Num_Experiments']}, "
                  f"Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            if pd.notna(row['Reaction_Types']) and len(row['Reaction_Types']) < 100:
                print(f"   Reactions: {row['Reaction_Types']}")
            print()
    else:
        import pandas as pd
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 60
        }
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved results to {args.output}")
        md_path = _write_markdown_sidecar_for_dataframe(
            args.output,
            title="HTE Catalyst Analysis",
            df=df,
            top=args.top,
            notes=[
                f"Reaction filter: {args.reaction or 'None'}",
                f"Reactant A filter: {args.reactant_a or 'None'}",
                f"Reactant B filter: {args.reactant_b or 'None'}",
            ],
        )
        print(f"Saved markdown summary to {md_path}")


def cmd_reactions(args):
    """Analyze reaction types command"""
    analytics = HTEAnalytics(args.db_path)
    reaction_filter = getattr(args, "reaction", None)
    detailed_map_requested = bool(getattr(args, "detailed_map", False))
    detail_top = max(1, int(getattr(args, "detail_top", 5)))
    include_detailed_map = detailed_map_requested or bool(reaction_filter)
    
    print("\n" + "="*80)
    print("REACTION TYPE SUMMARY")
    print("="*80)
    if reaction_filter:
        print(f"Reaction Type Filter: {reaction_filter}")
    if include_detailed_map:
        print(f"Detailed Map: enabled (top {detail_top})")
    print()
    
    df = analytics.get_reaction_type_summary(
        reaction_type=reaction_filter,
        include_detailed_map=include_detailed_map,
        detail_top_k=detail_top,
    )
    
    if reaction_filter:
        print(f"Found {len(df)} reaction types matching filter\n")
    else:
        print(f"Found {len(df)} reaction types in database\n")
    
    # Display results
    if args.compact:
        for i, row in df.head(args.top).iterrows():
            print(f"{i+1}. {row['Reaction_Type']}")
            print(f"   Experiments: {row['Num_Experiments']:,}, "
                  f"Pairs: {row['Num_Reactant_Pairs']}, "
                  f"Catalysts: {row['Num_Catalysts']}")
            print(f"   Avg Yield: {row['Avg_Yield']:.1f}%, "
                  f"Success Rate: {row['Success_Rate']:.1f}%")
            print(f"   Top Catalyst: {row['Top_Catalyst']}")
            print(f"   Top Pair: {row['Top_Reactant_Pair']}")
            print()
    else:
        import pandas as pd
        pd_options = {
            'display.max_rows': args.top,
            'display.max_columns': None,
            'display.width': None,
            'display.max_colwidth': 40
        }
        with pd.option_context(*[item for pair in pd_options.items() for item in pair]):
            print(df.head(args.top).to_string(index=False))
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved results to {args.output}")
        md_path = _write_markdown_sidecar_for_dataframe(
            args.output,
            title="HTE Reaction Type Summary",
            df=df,
            top=args.top,
            notes=[
                f"Reaction filter: {reaction_filter or 'None'}",
                f"Detailed map included: {'yes' if include_detailed_map else 'no'}",
                f"Detailed map top entries: {detail_top}",
            ],
        )
        print(f"Saved markdown summary to {md_path}")


def cmd_metals(args):
    """Analyze metal usage command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("METAL USAGE ANALYSIS")
    print("="*80)
    print()
    
    result = analytics.analyze_metal_usage()
    
    print(f"Total Experiments: {result['total_experiments']:,}\n")
    
    print("Metal Distribution:")
    print("-" * 50)
    
    import pandas as pd
    df = result['metal_distribution']
    
    for _, row in df.iterrows():
        metal = row['Metal']
        count = row['Num_Experiments']
        pct = row['Percentage']
        bar = "#" * int(pct / 2)
        print(f"{metal:>4}: {bar:<35} {count:>6,} ({pct:>5.1f}%)")
    
    if args.detailed:
        print("\n\nMetal Usage by Reaction Type:")
        print("-" * 50)
        
        for metal, reactions in sorted(result['by_reaction_type'].items(), 
                                      key=lambda x: sum(x[1].values()), reverse=True):
            if metal and metal != 'Other':
                print(f"\n{metal}:")
                for rxn, count in sorted(reactions.items(), key=lambda x: x[1], reverse=True)[:5]:
                    print(f"  {rxn}: {count:,}")
    
    if args.output:
        df.to_csv(args.output, index=False)
        print(f"\nSaved metal distribution to {args.output}")
        notes = [f"Detailed breakdown included in terminal output: {'yes' if args.detailed else 'no'}"]
        md_path = _write_markdown_sidecar_for_dataframe(
            args.output,
            title="HTE Metal Usage Analysis",
            df=df,
            top=None,
            notes=notes,
        )
        print(f"Saved markdown summary to {md_path}")


def cmd_export(args):
    """Export filtered dataset command"""
    analytics = HTEAnalytics(args.db_path)
    
    print("\n" + "="*80)
    print("EXPORT FILTERED DATASET")
    print("="*80)
    
    if args.reaction:
        print(f"Reaction Type: {args.reaction}")
    if args.catalyst:
        print(f"Catalyst Filter: {args.catalyst}")
    if args.reactant_a:
        print(f"Reactant A Type: {args.reactant_a}")
    if args.reactant_b:
        print(f"Reactant B Type: {args.reactant_b}")
    if args.min_yield:
        print(f"Min Yield: {args.min_yield}%")
    print()
    
    count = analytics.export_subset(
        output_path=args.output,
        reaction_type=args.reaction,
        catalyst_filter=args.catalyst,
        reactant_a_type=args.reactant_a,
        reactant_b_type=args.reactant_b,
        min_yield=args.min_yield
    )
    
    print(f"\nExport complete: {count:,} experiments")
    md_path = _markdown_sidecar_path(args.output)
    md_lines = [
        "# HTE Export Summary",
        "",
        f"- Output CSV: `{args.output}`",
        f"- Exported rows: {count:,}",
        f"- Reaction filter: {args.reaction or 'None'}",
        f"- Catalyst filter: {args.catalyst or 'None'}",
        f"- Reactant A filter: {args.reactant_a or 'None'}",
        f"- Reactant B filter: {args.reactant_b or 'None'}",
        f"- Minimum yield: {args.min_yield if args.min_yield is not None else 'None'}",
        "",
    ]
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"Saved markdown summary to {md_path}")


def cmd_backtest(args):
    """Backtest HTERecommender on held-out rows from an input CSV."""
    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    df = pd.read_csv(input_path)
    df = _filter_rows_by_contains(
        df,
        value=getattr(args, "reaction", None),
        candidates=("Reaction_Type_Standardized", "reaction_type"),
    )
    df = _filter_rows_by_contains(
        df,
        value=getattr(args, "catalyst", None),
        candidates=("Catalyst", "catalyst"),
    )
    df = _filter_rows_by_contains(
        df,
        value=getattr(args, "source_group", None),
        candidates=("Source_Group", "source_group"),
    )

    if df.empty:
        raise ValueError("No rows left after filters. Adjust --reaction/--catalyst/--source-group.")

    train_df, test_df = _prepare_split(
        df,
        test_fraction=float(args.test_fraction),
        test_size=args.test_size,
        seed=int(args.seed),
    )
    if train_df.empty or test_df.empty:
        raise ValueError("Train/test split failed. Increase dataset size or adjust split parameters.")

    temp_train_path: Optional[Path] = None
    train_path: Path
    if args.train_output:
        train_path = Path(args.train_output)
        train_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        handle = tempfile.NamedTemporaryFile(delete=False, suffix=".csv")
        handle.close()
        temp_train_path = Path(handle.name)
        train_path = temp_train_path
    train_df.to_csv(train_path, index=False)

    ks = [int(v.strip()) for v in str(args.hit_ks).split(",") if str(v).strip()]
    if not ks:
        ks = [1, 3, 5, 10]
    top_k_eval = max(int(args.top_k), max(ks))

    recommender = HTERecommender(str(train_path))
    train_condition_keys = {_row_condition_key(row) for _, row in train_df.iterrows()}

    ranks: List[Optional[int]] = []
    db_coverages: List[float] = []
    match_counts: List[float] = []
    per_row_records: List[Dict[str, Any]] = []

    skipped_invalid_query = 0
    skipped_empty_condition = 0
    skipped_unseen_condition = 0
    no_recommendation_count = 0

    try:
        for _, row in test_df.iterrows():
            target_key = _row_condition_key(row)
            if not any(target_key):
                skipped_empty_condition += 1
                continue
            if (not args.allow_unseen_conditions) and target_key not in train_condition_keys:
                skipped_unseen_condition += 1
                continue

            reactant_a, reactant_b, product_smiles, query_key = _extract_query_from_row(row)
            if not reactant_a:
                skipped_invalid_query += 1
                continue

            result = recommender.recommend(
                reactant_a_smiles=reactant_a,
                reactant_b_smiles=reactant_b,
                product_smiles=product_smiles,
                top_k=top_k_eval,
                min_experiments=int(args.min_experiments),
                reaction_key_only=bool(args.reaction_key_only),
                use_spectator_groups=bool(args.use_spectator_groups),
            )
            recs = list(getattr(result, "recommendations", []) or [])
            if not recs:
                no_recommendation_count += 1
            rank = _rank_for_key(recs, target_key)
            ranks.append(rank)
            db_coverages.append(float(getattr(result, "database_coverage", 0.0) or 0.0))
            match_counts.append(float(getattr(result, "total_matching_experiments", 0.0) or 0.0))

            if args.per_row_output:
                per_row_records.append(
                    {
                        "query": query_key or reactant_a,
                        "rank": rank,
                        "found": rank is not None,
                        "recommendations_returned": len(recs),
                        "database_coverage_pct": float(getattr(result, "database_coverage", 0.0) or 0.0),
                        "matching_experiments": int(getattr(result, "total_matching_experiments", 0) or 0),
                    }
                )

        metrics = _compute_backtest_metrics(ranks, ks)
        report: Dict[str, Any] = {
            "input_csv": str(input_path),
            "rows_after_filters": int(len(df)),
            "train_rows": int(len(train_df)),
            "test_rows": int(len(test_df)),
            "evaluated_rows": int(len(ranks)),
            "skipped_invalid_query": int(skipped_invalid_query),
            "skipped_empty_condition": int(skipped_empty_condition),
            "skipped_unseen_condition": int(skipped_unseen_condition),
            "no_recommendation_count": int(no_recommendation_count),
            "settings": {
                "reaction_filter": getattr(args, "reaction", None),
                "catalyst_filter": getattr(args, "catalyst", None),
                "source_group_filter": getattr(args, "source_group", None),
                "test_fraction": float(args.test_fraction),
                "test_size": args.test_size,
                "seed": int(args.seed),
                "min_experiments": int(args.min_experiments),
                "top_k_eval": int(top_k_eval),
                "hit_ks": ks,
                "reaction_key_only": bool(args.reaction_key_only),
                "use_spectator_groups": bool(args.use_spectator_groups),
                "allow_unseen_conditions": bool(args.allow_unseen_conditions),
            },
            "metrics": metrics,
            "diagnostics": {
                "avg_database_coverage_pct": _mean_or_zero(db_coverages),
                "avg_matching_experiments": _mean_or_zero(match_counts),
            },
        }

        if args.output:
            output_path = Path(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
            md_path = _write_markdown_sidecar_for_backtest_report(args.output, report)

        if args.per_row_output and per_row_records:
            per_row_path = Path(args.per_row_output)
            per_row_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(per_row_records).to_csv(per_row_path, index=False)

        print("\n" + "=" * 80)
        print("HTE RECOMMENDER BACKTEST")
        print("=" * 80)
        print(f"Input rows (after filters): {report['rows_after_filters']:,}")
        print(f"Train/Test/Evaluated: {report['train_rows']:,} / {report['test_rows']:,} / {report['evaluated_rows']:,}")
        print(f"Query coverage: {metrics['query_coverage'] * 100:.1f}%")
        print(f"MRR: {metrics['mrr']:.4f}")
        if metrics["avg_rank"] is not None:
            print(f"Avg rank (hits only): {metrics['avg_rank']:.2f}")
        else:
            print("Avg rank (hits only): N/A")
        for k in ks:
            print(f"Hit@{k}: {metrics[f'hit@{k}'] * 100:.1f}%")
        print(f"Avg DB coverage: {report['diagnostics']['avg_database_coverage_pct']:.2f}%")
        print(f"Avg matching experiments: {report['diagnostics']['avg_matching_experiments']:.2f}")
        if args.output:
            print(f"\nSaved report to {args.output}")
            print(f"Saved markdown summary to {md_path}")
        if args.per_row_output and per_row_records:
            print(f"Saved per-row details to {args.per_row_output}")
    finally:
        if temp_train_path is not None:
            try:
                temp_train_path.unlink(missing_ok=True)
            except Exception:
                pass


def main():
    if len(sys.argv) == 1:
        return _run_default_summary("data/HTE_db")

    parser = argparse.ArgumentParser(
        description="HTE Database Analytics Tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all Suzuki reactant pairs with Pd catalysts
  python -m chemtools.recommend.analytics pairs --reaction Suzuki --catalyst Pd --top 10
  
  # Analyze Cu catalysts in C-N coupling
  python -m chemtools.recommend.analytics catalysts --reaction "C-N" --catalyst Cu --compact
  
  # View reaction type summary
  python -m chemtools.recommend.analytics reactions --top 20
  
  # Analyze metal usage
  python -m chemtools.recommend.analytics metals --detailed
  
  # Export Pd-catalyzed Suzuki data
  python -m chemtools.recommend.analytics export suzuki_pd.csv --reaction Suzuki --catalyst Pd

  # Backtest deterministic HTE recommender on held-out rows
  python app/HTE_analytics_cli.py backtest --input data/HTE_db/experiments/HTE_canonical.csv --top-k 10 --hit-ks 1,3,5,10
        """
    )
    
    parser.add_argument(
        '--db-path',
        default='data/HTE_db',
        help='Path to HTE database CSV/JSONL or directory (default: data/HTE_db)'
    )
    parser.add_argument(
        '--interactive',
        action='store_true',
        help='Run interactive wizard mode'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Analytics command')
    
    # List pairs command
    pairs_parser = subparsers.add_parser('pairs', help='List reactant pairs')
    pairs_parser.add_argument('--reaction', help='Filter by reaction type')
    pairs_parser.add_argument('--catalyst', help='Filter by catalyst metal (e.g., Pd, Cu)')
    pairs_parser.add_argument('--min-experiments', type=int, default=1,
                             help='Minimum number of experiments (default: 1)')
    pairs_parser.add_argument('--sort', choices=['count', 'success_rate'], default='count',
                             help='Sort by count or success_rate (default: count)')
    pairs_parser.add_argument('--top', type=int, default=20,
                             help='Number of results to show (default: 20)')
    pairs_parser.add_argument('--compact', action='store_true',
                             help='Use compact output format')
    pairs_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Catalysts command
    cat_parser = subparsers.add_parser('catalysts', help='Analyze catalysts')
    cat_parser.add_argument('--reaction', help='Filter by reaction type')
    cat_parser.add_argument('--reactant-a', help='Filter by reactant A type')
    cat_parser.add_argument('--reactant-b', help='Filter by reactant B type')
    cat_parser.add_argument('--top', type=int, default=20,
                           help='Number of results to show (default: 20)')
    cat_parser.add_argument('--compact', action='store_true',
                           help='Use compact output format')
    cat_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Reactions command
    rxn_parser = subparsers.add_parser('reactions', help='Analyze reaction types')
    rxn_parser.add_argument('--reaction', help='Filter by reaction type')
    rxn_parser.add_argument(
        '--detailed-map',
        action='store_true',
        help='Include serialized Detailed_Map JSON column (top pairs/catalysts) in the report',
    )
    rxn_parser.add_argument(
        '--detail-top',
        type=int,
        default=5,
        help='Top N reactant pairs/catalysts to include in Detailed_Map (default: 5)',
    )
    rxn_parser.add_argument('--top', type=int, default=20,
                           help='Number of results to show (default: 20)')
    rxn_parser.add_argument('--compact', action='store_true',
                           help='Use compact output format')
    rxn_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Metals command
    metals_parser = subparsers.add_parser('metals', help='Analyze metal usage')
    metals_parser.add_argument('--detailed', action='store_true',
                              help='Show detailed breakdown by reaction type')
    metals_parser.add_argument('-o', '--output', help='Save results to CSV')
    
    # Export command
    export_parser = subparsers.add_parser('export', help='Export filtered dataset')
    export_parser.add_argument('output', help='Output CSV file path')
    export_parser.add_argument('--reaction', help='Filter by reaction type')
    export_parser.add_argument('--catalyst', help='Filter by catalyst metal')
    export_parser.add_argument('--reactant-a', help='Filter by reactant A type')
    export_parser.add_argument('--reactant-b', help='Filter by reactant B type')
    export_parser.add_argument('--min-yield', type=float, help='Minimum yield threshold')

    # Backtest command
    backtest_parser = subparsers.add_parser('backtest', help='Backtest HTE recommender offline')
    backtest_parser.add_argument(
        '--input',
        default='data/HTE_db/experiments/HTE_canonical.csv',
        help='Input HTE CSV for train/test backtest (default: data/HTE_db/experiments/HTE_canonical.csv)',
    )
    backtest_parser.add_argument('--reaction', help='Optional reaction type filter before split')
    backtest_parser.add_argument('--catalyst', help='Optional catalyst filter before split')
    backtest_parser.add_argument('--source-group', help='Optional source group filter before split')
    backtest_parser.add_argument(
        '--test-fraction',
        type=float,
        default=0.2,
        help='Fraction of rows for test split when --test-size is not set (default: 0.2)',
    )
    backtest_parser.add_argument(
        '--test-size',
        type=int,
        default=None,
        help='Exact number of test rows (overrides --test-fraction)',
    )
    backtest_parser.add_argument('--seed', type=int, default=13, help='Random seed for split reproducibility')
    backtest_parser.add_argument('--top-k', type=int, default=10, help='Top-k recommendations to request')
    backtest_parser.add_argument(
        '--hit-ks',
        default='1,3,5,10',
        help='Comma-separated k values for hit@k reporting (default: 1,3,5,10)',
    )
    backtest_parser.add_argument(
        '--min-experiments',
        type=int,
        default=1,
        help='Minimum experiments required per condition in recommender (default: 1)',
    )
    backtest_parser.add_argument(
        '--reaction-key-only',
        action='store_true',
        help='Only match by reaction key/signatures (disable reactant-type fallback)',
    )
    backtest_parser.add_argument(
        '--without-spectator-groups',
        action='store_true',
        help='Disable spectator-group weighting during backtest queries',
    )
    backtest_parser.add_argument(
        '--allow-unseen-conditions',
        action='store_true',
        help='Include test rows whose exact condition key is unseen in train (counted as misses if not retrieved)',
    )
    backtest_parser.add_argument(
        '--train-output',
        default=None,
        help='Optional path to save train split CSV (default: temporary file)',
    )
    backtest_parser.add_argument(
        '--output',
        default='results/hte_recommender_backtest.json',
        help='Path to write JSON report (default: results/hte_recommender_backtest.json)',
    )
    backtest_parser.add_argument(
        '--per-row-output',
        default=None,
        help='Optional CSV path for per-row ranks and diagnostics',
    )
    
    args = parser.parse_args()
    
    if args.interactive:
        return _run_wizard(args.db_path)

    if not args.command:
        parser.print_help()
        return 1
    
    try:
        if args.command == 'pairs':
            cmd_list_pairs(args)
        elif args.command == 'catalysts':
            cmd_catalysts(args)
        elif args.command == 'reactions':
            cmd_reactions(args)
        elif args.command == 'metals':
            cmd_metals(args)
        elif args.command == 'export':
            cmd_export(args)
        elif args.command == 'backtest':
            args.use_spectator_groups = not bool(args.without_spectator_groups)
            cmd_backtest(args)
        return 0
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
