"""ssbench CLI entrypoint."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from .config import load_config
from .dataset.bgsu_representatives import (
    DatasetFilters,
    download_representatives,
    split_by_rfam,
    stratify_buckets,
    write_manifest,
)
from .dataset.pdb_extract_chain import extract_chain
from .dataset.rcsb_fetch import fetch_pdb
from .metrics.ensemble_metrics import compute_ensemble_metrics
from .metrics.helix_metrics import compute_helix_metrics
from .metrics.pair_metrics import compute_pair_metrics
from .metrics.pk_metrics import compute_pk_metrics
from .predict.parse_dotbracket import pairs_from_dotbracket, parse_dotbracket
from .predict.run_command_predictor import run_predictor
from .report.summarize import summarize_metrics
from .sequence import normalize_rna_sequence, sanitize_rna_sequence
from .truth.barnaba_extract import build_truth
from .truth.pk_detect import split_nested_pk


def cmd_dataset_download(args: argparse.Namespace) -> None:
    cfg = load_config(args.config)
    filters = DatasetFilters(
        release_id=args.release_id or cfg.release_id,
        resolution=args.resolution or cfg.resolution,
        single_chain_only=True,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        limit=args.limit,
    )
    df = download_representatives(filters, full=args.full)
    df = split_by_rfam(
        df,
        seed=cfg.split_seed,
        drop_missing=cfg.drop_missing_rfam,
        fallback_to_ec=cfg.fallback_to_ec,
    )
    if "nts_observed" in df.columns:
        df["bucket"] = stratify_buckets(df["nts_observed"], cfg.length_buckets)
    write_manifest(df, args.out_manifest)
    print(f"Wrote manifest to {args.out_manifest} ({len(df)} rows)")


def cmd_dataset_fetch_pdb(args: argparse.Namespace) -> None:
    df = pd.read_csv(args.manifest)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    kept_rows = []
    fetched = 0
    skipped = 0
    for _, row in df.iterrows():
        pdb_id = str(row["pdb_id"]).lower()
        chain_id = str(row["chain_id"])
        pdb_path = out_dir / f"{pdb_id}.pdb"
        chain_path = out_dir / f"{pdb_id}_{chain_id}.pdb"
        try:
            if not pdb_path.exists() and not pdb_path.with_suffix(".cif").exists():
                pdb_path = fetch_pdb(pdb_id, pdb_path)
            elif pdb_path.with_suffix(".cif").exists():
                pdb_path = pdb_path.with_suffix(".cif")
            if not chain_path.exists():
                extract_chain(pdb_path, chain_id, chain_path)
            fetched += 1
            kept_rows.append(row)
        except Exception as exc:
            print(f"[WARN] Skipping {pdb_id}:{chain_id} due to fetch error: {exc}")
            skipped += 1

    if args.out_manifest is not None:
        pd.DataFrame(kept_rows).to_csv(args.out_manifest, index=False)
        print(f"Wrote fetched manifest to {args.out_manifest} ({len(kept_rows)} rows)")
    print(f"Fetch summary: ok={fetched} skipped={skipped}")


def cmd_truth_build(args: argparse.Namespace) -> None:
    df = pd.read_csv(args.manifest)
    pdb_dir = Path(args.pdb_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    built = 0
    already = 0
    skipped = 0
    for _, row in df.iterrows():
        pdb_id = str(row["pdb_id"]).lower()
        chain_id = str(row["chain_id"])
        target_id = str(row["target_id"])
        out_path = out_dir / f"{target_id}.json"
        if out_path.exists():
            already += 1
            continue
        pdb_path = pdb_dir / f"{pdb_id}_{chain_id}.pdb"
        if not pdb_path.exists():
            print(f"[WARN] Missing PDB chain for {target_id}, skipping")
            continue
        try:
            record = build_truth(pdb_path)
        except SystemExit as exc:
            print(f"[WARN] Truth build failed for {target_id}: {exc}")
            skipped += 1
            continue
        except Exception as exc:
            print(f"[WARN] Truth build failed for {target_id}: {exc}")
            skipped += 1
            continue
        out_path.write_text(json.dumps(record.__dict__, indent=2))
        built += 1
    if built == 0 and already == 0:
        raise SystemExit("No truth records built. Check PDB fetch step.")
    print(f"Truth records built: {built} already={already} skipped={skipped}")


def cmd_predict_run(args: argparse.Namespace) -> None:
    cfg = load_config(args.config)
    df = pd.read_csv(args.manifest)
    truth_dir = Path(args.truth)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    k = args.k if args.k is not None else cfg.k
    produced = 0
    already = 0
    for _, row in df.iterrows():
        target_id = str(row["target_id"])
        safe_id = target_id.replace("|", "_")
        truth_path = truth_dir / f"{target_id}.json"
        if not truth_path.exists():
            continue
        data = json.loads(truth_path.read_text())
        seq_raw = data["sequence"]
        seq = normalize_rna_sequence(seq_raw)
        # Preserve ambiguity/modified bases (e.g. 'X') in the FASTA input.
        # Substituting them (e.g. X->U) can substantially change CaCoFold scaffolds and harm
        # downstream performance (notably for tRNA targets).
        #
        # Some truth sequences can include a Barnaba fragment separator '&'. Predictors that
        # call external tools should sanitize that internally while preserving length.
        if any(ch not in "ACGUXN&" for ch in seq):
            seq = sanitize_rna_sequence(seq, unknown="N")
            print(
                f"[WARN] {target_id}: non-(ACGU/X/N/& ) bases in truth sequence; substituting unknown->N for prediction input"
            )
        elif any(ch not in "ACGU" for ch in seq):
            # Keep as-is; predictors may handle ambiguity better than a forced substitution.
            print(f"[WARN] {target_id}: non-ACGU bases in truth sequence; leaving unchanged for prediction input")
        fasta_path = out_dir / f"{safe_id}.fa"
        fasta_path.write_text(f">{target_id}\n{seq}\n")
        pred_dir = out_dir / safe_id
        pred_dir.mkdir(parents=True, exist_ok=True)
        out_path = pred_dir / "predictions.db"
        if out_path.exists():
            try:
                _seq0, dotbrs0 = parse_dotbracket(out_path.read_text())
            except Exception:
                dotbrs0 = []
            if len(dotbrs0) >= int(k):
                already += 1
                continue
        try:
            dotbrackets = run_predictor(args.predictor_cmd, fasta_path, pred_dir, k)
        except Exception as exc:
            print(f"[WARN] Predictor failed for {target_id}: {exc}")
            dotbrackets = []
        flat_dotbrackets: list[str] = []
        for item in dotbrackets:
            if isinstance(item, str):
                flat_dotbrackets.append(item)
            elif isinstance(item, list):
                flat_dotbrackets.extend([s for s in item if isinstance(s, str)])
        dotbrackets = flat_dotbrackets
        if not dotbrackets:
            dotbrackets = ["." * len(seq)]
        out_path.write_text("\n".join([seq] + dotbrackets) + "\n")
        produced += 1
    if produced == 0 and already == 0:
        raise SystemExit("No predictions produced. Check truth or predictor output.")
    print(f"Predictions produced: {produced} already={already}")


def cmd_score(args: argparse.Namespace) -> None:
    cfg = load_config(args.config)
    df = pd.read_csv(args.manifest)
    truth_dir = Path(args.truth)
    pred_dir = Path(args.predictions)

    helix_overlap = (
        args.helix_overlap if args.helix_overlap is not None else cfg.helix_overlap_threshold
    )
    rows = []
    for _, row in df.iterrows():
        target_id = str(row["target_id"])
        truth_path = truth_dir / f"{target_id}.json"
        safe_id = target_id.replace("|", "_")
        pred_path = pred_dir / safe_id / "predictions.db"
        if not truth_path.exists() or not pred_path.exists():
            continue
        truth_data = json.loads(truth_path.read_text())
        seq = str(truth_data["sequence"])
        L = len(seq)

        def _filter_pairs(pairs: set[tuple[int, int]]) -> set[tuple[int, int]]:
            out: set[tuple[int, int]] = set()
            for i, j in pairs:
                if 0 <= i < j < L:
                    out.add((i, j))
            return out

        ref_pairs = _filter_pairs({tuple(p) for p in truth_data["canonical_pairs"]})
        ref_nested = _filter_pairs({tuple(p) for p in truth_data["nested_pairs"]})
        ref_pk = _filter_pairs({tuple(p) for p in truth_data["pk_pairs"]})

        pred_seq, dotbrs = parse_dotbracket(pred_path.read_text())
        if pred_seq and len(pred_seq) != L:
            print(
                f"[WARN] {target_id}: prediction sequence length {len(pred_seq)} != truth length {L}; "
                "dropping mismatched structures"
            )
        dotbrs = [s for s in dotbrs if len(s) == L]

        pred_pairs_list = [set(pairs_from_dotbracket(s)) for s in dotbrs]
        if not pred_pairs_list:
            pred_pairs_list = [set()]

        best_pred = pred_pairs_list[0]

        # QC metrics that reflect *top-100 yield* and exploration/diversity without requiring
        # any particular ranking to be correct.
        def jaccard_distance(a: set[tuple[int, int]], b: set[tuple[int, int]]) -> float:
            if not a and not b:
                return 0.0
            inter = len(a & b)
            union = len(a | b)
            return 1.0 - (inter / union if union else 0.0)

        topn = min(100, len(pred_pairs_list))
        top_pairs = pred_pairs_list[:topn]
        union_pairs: set[tuple[int, int]] = set()
        for ps in top_pairs:
            union_pairs |= ps
        union_recall_100 = 0.0
        if ref_pairs:
            union_recall_100 = len(union_pairs & ref_pairs) / len(ref_pairs)
        pair_counts = [len(ps) for ps in top_pairs]
        pair_counts_sorted = sorted(pair_counts)
        mean_pred_pairs_100 = float(sum(pair_counts) / len(pair_counts)) if pair_counts else 0.0
        median_pred_pairs_100 = float(pair_counts_sorted[len(pair_counts_sorted) // 2]) if pair_counts else 0.0
        mean_jaccard_to_top1_100 = (
            float(sum(jaccard_distance(best_pred, ps) for ps in top_pairs[1:]) / (topn - 1))
            if topn > 1
            else 0.0
        )
        n_unique_structs_100 = int(len(set(dotbrs[:topn]))) if dotbrs else 0

        pair_metrics = compute_pair_metrics(best_pred, ref_pairs, L)
        helix_metrics = compute_helix_metrics(list(best_pred), list(ref_pairs), helix_overlap)
        nested_pred, pk_pred = split_nested_pk(sorted(best_pred))
        pk_metrics = compute_pk_metrics(
            set(nested_pred),
            set(pk_pred),
            ref_nested,
            ref_pk,
            L,
        )
        ensemble_metrics = compute_ensemble_metrics(pred_pairs_list, ref_pairs, L)

        rows.append(
            {
                "target_id": target_id,
                "split": row.get("split", "train"),
                "bucket": row.get("bucket", "unknown"),
                "precision": pair_metrics.precision,
                "recall": pair_metrics.recall,
                "f1": pair_metrics.f1,
                "mcc": pair_metrics.mcc,
                "helix_precision": helix_metrics.precision,
                "helix_recall": helix_metrics.recall,
                "lonely_pair_rate": helix_metrics.lonely_pair_rate,
                "nested_f1": pk_metrics.nested.f1,
                "pk_f1": pk_metrics.pk.f1,
                "best_of_k_f1": ensemble_metrics.best_f1,
                "top1_f1": ensemble_metrics.top1_f1,
                "union_recall_100": union_recall_100,
                "union_pairs_100": int(len(union_pairs)),
                "mean_pred_pairs_100": mean_pred_pairs_100,
                "median_pred_pairs_100": median_pred_pairs_100,
                "mean_jaccard_to_top1_100": mean_jaccard_to_top1_100,
                "n_unique_structs_100": n_unique_structs_100,
                **{f"coverage_{k}": v for k, v in ensemble_metrics.coverage.items()},
            }
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise SystemExit("No rows scored. Check truth and predictions paths.")
    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_path, index=False)
    print(f"Wrote metrics to {out_path} ({len(rows)} rows)")
    if args.min_f1 is not None:
        mean_f1 = df_out["f1"].mean()
        if mean_f1 < args.min_f1:
            raise SystemExit(
                f"Mean F1 {mean_f1:.3f} below threshold {args.min_f1:.3f}. "
                "Consider increasing top-k, sampling steps, or using covariation-aware predictions."
            )


def cmd_report(args: argparse.Namespace) -> None:
    summarize_metrics(args.metrics, args.out)
    print(f"Wrote report to {args.out}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ssbench", description="RNA structure benchmark runner")
    sub = parser.add_subparsers(dest="command")

    dataset = sub.add_parser("dataset")
    dataset_sub = dataset.add_subparsers(dest="dataset_cmd")

    dl = dataset_sub.add_parser("download")
    dl.add_argument("--release-id", default="current")
    dl.add_argument("--resolution", default="4.0")
    dl.add_argument("--out-manifest", required=True)
    dl.add_argument("--limit", type=int, default=None)
    dl.add_argument("--full", action="store_true")
    dl.add_argument("--config", default=None)
    dl.set_defaults(func=cmd_dataset_download)

    fetch = dataset_sub.add_parser("fetch-pdb")
    fetch.add_argument("--manifest", required=True)
    fetch.add_argument("--out-dir", required=True)
    fetch.add_argument("--out-manifest", default=None)
    fetch.set_defaults(func=cmd_dataset_fetch_pdb)

    truth = sub.add_parser("truth")
    truth_sub = truth.add_subparsers(dest="truth_cmd")
    tb = truth_sub.add_parser("build")
    tb.add_argument("--manifest", required=True)
    tb.add_argument("--pdb-dir", required=True)
    tb.add_argument("--out-dir", required=True)
    tb.set_defaults(func=cmd_truth_build)

    predict = sub.add_parser("predict")
    predict_sub = predict.add_subparsers(dest="predict_cmd")
    pr = predict_sub.add_parser("run")
    pr.add_argument("--manifest", required=True)
    pr.add_argument("--truth", required=True)
    pr.add_argument("--predictor-cmd", required=True)
    pr.add_argument("--out-dir", required=True)
    pr.add_argument("--k", type=int, default=None)
    pr.add_argument("--config", default=None)
    pr.set_defaults(func=cmd_predict_run)

    score = sub.add_parser("score")
    score.add_argument("--manifest", required=True)
    score.add_argument("--truth", required=True)
    score.add_argument("--predictions", required=True)
    score.add_argument("--out", required=True)
    score.add_argument("--helix-overlap", type=float, default=None)
    score.add_argument("--config", default=None)
    score.add_argument("--min-f1", type=float, default=None)
    score.set_defaults(func=cmd_score)

    report = sub.add_parser("report")
    report.add_argument("--metrics", required=True)
    report.add_argument("--out", required=True)
    report.set_defaults(func=cmd_report)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if not hasattr(args, "func"):
        parser.print_help()
        raise SystemExit(1)
    args.func(args)


if __name__ == "__main__":
    main()
