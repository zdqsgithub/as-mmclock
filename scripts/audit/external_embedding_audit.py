#!/usr/bin/env python3
"""
Audit external methylation foundation models for this mouse RRBS project.

This script is intentionally non-integrating: it checks installation/runtime
feasibility and feature compatibility, then writes a report. It does not use
external model outputs in benchmark metrics.
"""
from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
from pathlib import Path

import pandas as pd

ROOT = Path("/home/zdq-as/mouse_methyl_work")


def package_status(package: str) -> dict:
    spec = importlib.util.find_spec(package)
    installed = spec is not None
    version = None
    if installed:
        try:
            version = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            version = "unknown"
    return {"installed": installed, "version": version}


def torch_status() -> dict:
    spec = importlib.util.find_spec("torch")
    if spec is None:
        return {"installed": False, "cuda_available": False, "device_count": 0}
    import torch

    return {
        "installed": True,
        "version": torch.__version__,
        "cuda_available": bool(torch.cuda.is_available()),
        "device_count": int(torch.cuda.device_count()),
    }


def infer_feature_space(matrix_path: Path) -> dict:
    matrix = pd.read_parquet(matrix_path)
    features = [str(idx) for idx in matrix.index[:1000]]
    region_like = sum((":" in f and "-" in f) for f in features)
    cpg_like = sum(("_" in f and f.rsplit("_", 1)[-1].isdigit()) for f in features)
    if region_like > cpg_like:
        kind = "mouse_region_bins"
    elif cpg_like > 0:
        kind = "mouse_cpg_coordinates"
    else:
        kind = "unknown"
    return {
        "matrix_path": str(matrix_path),
        "shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "feature_space": kind,
        "first_features": features[:5],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_path", default=str(ROOT / "results" / "phase0" / "region_matrix_5kb.parquet"))
    parser.add_argument("--out_dir", default=str(ROOT / "results" / "external_embedding_audit"))
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_info = infer_feature_space(Path(args.matrix_path))
    torch_info = torch_status()
    methylgpt = package_status("methylgpt")
    cpgpt = package_status("cpgpt")

    findings = {
        "methylgpt": {
            **methylgpt,
            "known_capabilities": [
                "embedding_extraction",
                "age_prediction",
                "imputation",
                "cpg_selection",
            ],
            "compatibility": "blocked_for_benchmark"
            if matrix_info["feature_space"] == "mouse_region_bins"
            else "requires_vocabulary_mapping",
            "reason": "Current project matrix is mouse 5kb RRBS regions, while public MethylGPT examples use human methylation CpG/probe vocabularies.",
            "source": "https://github.com/albert-ying/MethylGPT",
        },
        "cpgpt": {
            **cpgpt,
            "known_capabilities": [
                "methylation_foundation_model",
                "imputation",
                "age_prediction",
                "attention_based_feature_interpretability",
            ],
            "compatibility": "blocked_for_benchmark"
            if matrix_info["feature_space"] == "mouse_region_bins"
            else "requires_vocabulary_mapping",
            "reason": "Current project matrix is mouse 5kb RRBS regions; CpGPT compatibility requires a supported CpG vocabulary and feature ordering.",
            "source": "https://spacefrontiers.org/r/10.1101/2024.10.24.619766",
        },
    }

    payload = {
        "project": "mouse_methyl_work",
        "policy": "external embeddings are gated POC only and must not enter benchmark until compatible with mouse RRBS features",
        "runtime": {"torch": torch_info},
        "matrix": matrix_info,
        "models": findings,
        "decision": "do_not_use_external_embeddings_in_v4_benchmark",
    }
    (out_dir / "external_embedding_audit.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# External Embedding Model Audit",
        "",
        "External foundation models are audited as a gated POC only. They are not used in v4 benchmark metrics.",
        "",
        f"- Matrix: `{matrix_info['matrix_path']}`",
        f"- Shape: `{matrix_info['shape'][0]} features x {matrix_info['shape'][1]} samples`",
        f"- Feature space: `{matrix_info['feature_space']}`",
        f"- Torch installed: `{torch_info['installed']}`; CUDA available: `{torch_info['cuda_available']}`",
        "",
        "| Model | Installed | Compatibility | Decision |",
        "|---|---:|---|---|",
    ]
    for name, model in findings.items():
        lines.append(
            f"| {name} | {model['installed']} | {model['compatibility']} | not used in v4 benchmark |"
        )
    lines.extend(
        [
            "",
            "## Rationale",
            "",
            "The current v4 input is mouse RRBS 5kb region bins. MethylGPT/CpGPT public interfaces are built around methylation CpG/probe vocabularies, mainly human-scale pretraining. Direct use would require a validated mouse vocabulary or mapping layer. AS-DS-Ops forbids using human clock CpGs or unvalidated cross-species mappings as mouse clock features.",
        ]
    )
    (out_dir / "external_embedding_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[Audit] Wrote {out_dir / 'external_embedding_audit.json'}")
    print(f"[Audit] Wrote {out_dir / 'external_embedding_audit.md'}")


if __name__ == "__main__":
    main()
