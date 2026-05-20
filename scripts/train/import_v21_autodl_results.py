#!/usr/bin/env python3
"""Import and summarize v21 AutoDL results."""
from __future__ import annotations

import argparse
import json
import shutil
import tarfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path("/home/zdq-as/mouse_methyl_work")
DEFAULT_OUT = ROOT / "results" / "ralph_v21_raid_raw_clock" / "autodl_import"
REPORT = ROOT / "doc" / "20_analysis" / "53_20260521_v21_autodl_import_report.md"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str), encoding="utf-8")


def materialize_source(source: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        target = out_dir / source.name
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(source, target)
        return target
    if source.suffixes[-2:] == [".tar", ".gz"] or source.suffix == ".tgz":
        target = out_dir / source.stem.replace(".tar", "")
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)
        with tarfile.open(source, "r:gz") as tar:
            tar.extractall(target)
        return target
    raise SystemExit(f"Unsupported AutoDL result source: {source}")


def find_first(root: Path, name: str) -> Path | None:
    matches = sorted(root.rglob(name))
    return matches[0] if matches else None


def read_csv_if_exists(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    imported_root = materialize_source(args.source, args.out_dir)
    summary_path = find_first(imported_root, "v21_dl_summary.csv")
    lodo_path = find_first(imported_root, "v21_lodo_summary.csv")
    random_path = find_first(imported_root, "result.json")
    summary = read_csv_if_exists(summary_path)
    lodo = read_csv_if_exists(lodo_path)

    best = {}
    if not summary.empty:
        ranked = summary.sort_values(["mae_weeks", "rmse_weeks"], ascending=True)
        best = ranked.iloc[0].to_dict()
    lodo_metrics = {}
    if not lodo.empty and "mae_weeks" in lodo.columns:
        lodo_metrics = {
            "lodo_mean_mae_weeks": round(float(lodo["mae_weeks"].mean()), 3),
            "lodo_worst_mae_weeks": round(float(lodo["mae_weeks"].max()), 3),
            "lodo_rows": int(len(lodo)),
        }

    payload = {
        "timestamp": utc_now(),
        "status": "completed" if not summary.empty else "no_summary_found",
        "source": str(args.source),
        "imported_root": str(imported_root),
        "summary_path": str(summary_path) if summary_path else "",
        "lodo_path": str(lodo_path) if lodo_path else "",
        "n_summary_rows": int(len(summary)),
        "best_config": best.get("config_id"),
        "best_group_mae_weeks": best.get("mae_weeks"),
        "best_group_pearson_r": best.get("pearson_r"),
        **lodo_metrics,
        "training_authorized_local": False,
        "autoresearch_authorized": False,
    }
    write_json(args.out_dir / "v21_autodl_import_summary.json", payload)

    lines = [
        "# v21 AutoDL Import Report",
        "",
        f"Date: {utc_now()}",
        "",
        f"- status: `{payload['status']}`",
        f"- source: `{args.source}`",
        f"- best config: `{payload.get('best_config')}`",
        f"- best GroupKFold MAE: `{payload.get('best_group_mae_weeks')}`",
        f"- LODO mean/worst MAE: `{payload.get('lodo_mean_mae_weeks')}` / `{payload.get('lodo_worst_mae_weeks')}`",
        "",
        "Local DL training remains unauthorized; these are imported AutoDL artifacts.",
    ]
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
