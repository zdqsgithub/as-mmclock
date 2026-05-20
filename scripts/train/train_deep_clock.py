#!/usr/bin/env python3
"""Leakage-safe PyTorch methylation clock training for v20 DL autoresearch."""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.impute import SimpleImputer
from sklearn.model_selection import GroupKFold, KFold, train_test_split
from sklearn.preprocessing import QuantileTransformer, RobustScaler, StandardScaler

import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


ROOT = Path("/home/zdq-as/mouse_methyl_work")
GSM_RE = re.compile(r"(GSM\d+)")
GSE60012_TILE_RE = re.compile(r"^(GSE60012_tile_\d{3})")


def resolve_matrix_sample_id(value: object) -> str:
    text = str(value)
    tile = GSE60012_TILE_RE.match(text)
    if tile:
        return tile.group(1)
    gsm = GSM_RE.search(text)
    if gsm:
        return gsm.group(1)
    return text


def seed_everything(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = True


def compute_metrics(y_true_weeks: np.ndarray, y_pred_weeks: np.ndarray) -> dict:
    if len(y_true_weeks) < 3 or np.std(y_pred_weeks) == 0 or np.std(y_true_weeks) == 0:
        r, pval = 0.0, 1.0
    else:
        r, pval = stats.pearsonr(y_true_weeks, y_pred_weeks)
        if np.isnan(r):
            r, pval = 0.0, 1.0
    mae = float(np.mean(np.abs(y_true_weeks - y_pred_weeks)))
    medae = float(np.median(np.abs(y_true_weeks - y_pred_weeks)))
    rmse = float(np.sqrt(np.mean((y_true_weeks - y_pred_weeks) ** 2)))
    ss_res = float(np.sum((y_true_weeks - y_pred_weeks) ** 2))
    ss_tot = float(np.sum((y_true_weeks - y_true_weeks.mean()) ** 2))
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    return {
        "pearson_r": round(float(r), 4),
        "pearson_pval": float(pval),
        "mae_weeks": round(mae, 3),
        "medae_weeks": round(medae, 3),
        "rmse_weeks": round(rmse, 3),
        "r2": round(float(r2), 4),
    }


def fast_pearson_correlations(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    X_mean = np.nanmean(X, axis=0)
    y_mean = float(np.mean(y))
    X_centered = X - X_mean
    y_centered = y - y_mean
    ss_X = np.nansum(X_centered**2, axis=0)
    ss_y = float(np.sum(y_centered**2))
    cov = np.nansum(X_centered * y_centered[:, np.newaxis], axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / np.sqrt(ss_X * ss_y)
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)


def select_features(X_train: np.ndarray, y_train: np.ndarray, n_features: int | str) -> np.ndarray:
    if str(n_features).lower() == "all":
        n_features = X_train.shape[1]
    n_features = min(int(n_features), X_train.shape[1])
    corrs = fast_pearson_correlations(X_train, y_train)
    return np.argsort(np.abs(corrs))[-n_features:]


def train_dataset_stability_mask(
    X_train: np.ndarray,
    groups_train: np.ndarray,
    min_group_presence: float | None,
    max_group_mean_shift: float | None,
) -> np.ndarray:
    mask = np.ones(X_train.shape[1], dtype=bool)
    unique_groups = [group for group in pd.unique(groups_train) if pd.notna(group)]
    group_means = []
    for group in unique_groups:
        group_X = X_train[groups_train == group]
        finite = np.isfinite(group_X)
        if min_group_presence is not None:
            mask &= finite.mean(axis=0) >= min_group_presence
        if max_group_mean_shift is not None:
            with np.errstate(invalid="ignore"):
                group_means.append(np.nanmean(group_X, axis=0))
    if max_group_mean_shift is not None and group_means:
        means = np.vstack(group_means)
        finite_means = np.isfinite(means).all(axis=0)
        with np.errstate(invalid="ignore"):
            shift = np.nanmax(means, axis=0) - np.nanmin(means, axis=0)
        mask &= finite_means & (shift <= max_group_mean_shift)
    return mask


def train_agebin_presence_mask(
    X_train: np.ndarray,
    age_weeks_train: np.ndarray,
    min_agebin_presence: float | None,
    min_agebin_samples: int,
) -> np.ndarray:
    mask = np.ones(X_train.shape[1], dtype=bool)
    if min_agebin_presence is None:
        return mask
    bins = np.array([0, 4, 13, 26, 52, 104, np.inf], dtype=float)
    bin_ids = np.digitize(age_weeks_train, bins[1:-1], right=False)
    populated = [bin_id for bin_id in np.unique(bin_ids) if np.sum(bin_ids == bin_id) >= min_agebin_samples]
    if len(populated) < 2:
        return mask
    finite = np.isfinite(X_train)
    for bin_id in populated:
        mask &= finite[bin_ids == bin_id].mean(axis=0) >= min_agebin_presence
    return mask


def make_scaler(preprocess: str, seed: int, n_train: int):
    if preprocess == "standard":
        return StandardScaler()
    if preprocess == "robust":
        return RobustScaler()
    if preprocess == "quantile_uniform":
        return QuantileTransformer(
            n_quantiles=max(10, min(1000, n_train)),
            output_distribution="uniform",
            random_state=seed,
        )
    if preprocess == "none":
        return None
    raise ValueError(f"Unsupported preprocess: {preprocess}")


def build_target(age_days: np.ndarray, transform: str) -> np.ndarray:
    if transform == "log1p_days":
        return np.log1p(age_days).astype(np.float32)
    if transform == "linear_weeks":
        return (age_days / 7).astype(np.float32)
    raise ValueError(f"Unsupported target transform: {transform}")


def inverse_target(pred: np.ndarray, transform: str, max_age_days: float) -> tuple[np.ndarray, np.ndarray]:
    if transform == "log1p_days":
        pred = np.clip(pred, 0, np.log1p(max(max_age_days * 1.5, 1.0)))
        pred_days = np.expm1(pred)
        return pred_days, pred_days / 7
    if transform == "linear_weeks":
        pred_weeks = np.clip(pred, 0, max(max_age_days / 7 * 1.5, 1.0))
        return pred_weeks * 7, pred_weeks
    raise ValueError(f"Unsupported target transform: {transform}")


def parse_region_start(region_id: str) -> tuple[str, int]:
    match = re.match(r"([^:]+):(\d+)-", str(region_id))
    if not match:
        return ("zz", 10**18)
    return (match.group(1), int(match.group(2)))


def dataset_filter(meta: pd.DataFrame, selector: str) -> pd.Series:
    selector = selector.strip()
    if selector.startswith("all_except:"):
        excluded = {item.strip() for item in selector.split(":", 1)[1].split(",") if item.strip()}
        return ~meta["dataset_batch"].isin(excluded)
    if selector in {"all", "*"}:
        return pd.Series(True, index=meta.index)
    allowed = {item.strip() for item in selector.split(",") if item.strip()}
    return meta["dataset_batch"].isin(allowed)


def load_data(matrix_path: Path, metadata_path: Path, dataset_selector: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not matrix_path.exists():
        raise SystemExit(f"Matrix missing: {matrix_path}")
    if not metadata_path.exists():
        raise SystemExit(f"Metadata missing: {metadata_path}")
    matrix = pd.read_parquet(matrix_path)
    matrix.columns = [resolve_matrix_sample_id(col) for col in matrix.columns]
    if pd.Index(matrix.columns).duplicated().any():
        duplicates = sorted(pd.Index(matrix.columns)[pd.Index(matrix.columns).duplicated()].unique())
        raise SystemExit(f"Duplicate resolved matrix sample ids: {duplicates[:10]}")
    meta = pd.read_csv(metadata_path)
    meta = meta[dataset_filter(meta, dataset_selector) & meta["age_days"].notna()].copy()
    common = [sample for sample in meta["sample_id"] if sample in matrix.columns]
    meta = meta[meta["sample_id"].isin(common)].reset_index(drop=True)
    if not common:
        raise SystemExit("No common samples between matrix and metadata.")
    return matrix[common], meta


class ResidualBlock(nn.Module):
    def __init__(self, width: int, dropout: float):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(width),
            nn.Linear(width, width),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(width, width),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class ResidualMLP(nn.Module):
    def __init__(self, n_features: int, width: int, depth: int, dropout: float):
        super().__init__()
        layers: list[nn.Module] = [nn.Linear(n_features, width), nn.GELU(), nn.Dropout(dropout)]
        layers += [ResidualBlock(width, dropout) for _ in range(max(depth - 1, 0))]
        layers += [nn.LayerNorm(width), nn.Linear(width, 1)]
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x).squeeze(-1)


class WideDeepMLP(nn.Module):
    def __init__(self, n_features: int, width: int, bottleneck: int, depth: int, dropout: float):
        super().__init__()
        deep: list[nn.Module] = [nn.Linear(n_features, width), nn.GELU(), nn.Dropout(dropout)]
        for _ in range(max(depth - 1, 0)):
            deep += [nn.Linear(width, width), nn.GELU(), nn.Dropout(dropout)]
        deep += [nn.Linear(width, bottleneck), nn.GELU()]
        self.deep = nn.Sequential(*deep)
        self.head = nn.Sequential(nn.LayerNorm(n_features + bottleneck), nn.Linear(n_features + bottleneck, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(torch.cat([x, self.deep(x)], dim=1)).squeeze(-1)


class GenomicCNN(nn.Module):
    def __init__(self, n_features: int, channels: int, kernel: int, dropout: float):
        super().__init__()
        pad = kernel // 2
        self.net = nn.Sequential(
            nn.Conv1d(1, channels, kernel, padding=pad),
            nn.BatchNorm1d(channels),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel, padding=pad),
            nn.BatchNorm1d(channels),
            nn.GELU(),
            nn.AdaptiveAvgPool1d(1),
            nn.Flatten(),
            nn.Linear(channels, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x.unsqueeze(1)).squeeze(-1)


class TabTransformerLite(nn.Module):
    def __init__(self, n_features: int, token_dim: int, heads: int, layers: int, dropout: float):
        super().__init__()
        self.value_proj = nn.Linear(1, token_dim)
        self.pos = nn.Parameter(torch.zeros(1, n_features, token_dim))
        enc_layer = nn.TransformerEncoderLayer(
            d_model=token_dim,
            nhead=heads,
            dim_feedforward=token_dim * 4,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=layers)
        self.head = nn.Sequential(nn.LayerNorm(token_dim), nn.Linear(token_dim, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.value_proj(x.unsqueeze(-1)) + self.pos[:, : x.shape[1], :]
        pooled = self.encoder(tokens).mean(dim=1)
        return self.head(pooled).squeeze(-1)


class DenoisingAutoEncoderRegressor(nn.Module):
    def __init__(self, n_features: int, width: int, latent: int, dropout: float):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_features, width),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(width, latent),
            nn.GELU(),
        )
        self.decoder = nn.Sequential(nn.Linear(latent, width), nn.GELU(), nn.Linear(width, n_features))
        self.head = nn.Sequential(nn.LayerNorm(latent), nn.Linear(latent, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.encoder(x)).squeeze(-1)

    def reconstruct(self, x: torch.Tensor) -> torch.Tensor:
        return self.decoder(self.encoder(x))


@dataclass
class ModelConfig:
    architecture: str
    width: int
    depth: int
    dropout: float
    bottleneck: int
    cnn_channels: int
    cnn_kernel: int
    token_dim: int
    transformer_heads: int
    transformer_layers: int
    dae_latent: int


def make_model(n_features: int, cfg: ModelConfig) -> nn.Module:
    if cfg.architecture == "res_mlp":
        return ResidualMLP(n_features, cfg.width, cfg.depth, cfg.dropout)
    if cfg.architecture == "wide_deep":
        return WideDeepMLP(n_features, cfg.width, cfg.bottleneck, cfg.depth, cfg.dropout)
    if cfg.architecture == "cnn1d":
        return GenomicCNN(n_features, cfg.cnn_channels, cfg.cnn_kernel, cfg.dropout)
    if cfg.architecture == "tab_transformer":
        return TabTransformerLite(
            n_features, cfg.token_dim, cfg.transformer_heads, cfg.transformer_layers, cfg.dropout
        )
    if cfg.architecture == "dae_mlp":
        return DenoisingAutoEncoderRegressor(n_features, cfg.width, cfg.dae_latent, cfg.dropout)
    raise ValueError(f"Unsupported architecture: {cfg.architecture}")


def train_val_split(train_idx: np.ndarray, y: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    bins = pd.qcut(pd.Series(y[train_idx]), q=min(5, len(train_idx) // 20), duplicates="drop", labels=False)
    stratify = bins.to_numpy() if len(np.unique(bins)) > 1 else None
    tr, val = train_test_split(train_idx, test_size=0.18, random_state=seed, stratify=stratify)
    return np.asarray(tr), np.asarray(val)


def tensor_loader(X: np.ndarray, y: np.ndarray | None, batch_size: int, shuffle: bool) -> DataLoader:
    X_t = torch.from_numpy(X.astype(np.float32))
    if y is None:
        return DataLoader(TensorDataset(X_t), batch_size=batch_size, shuffle=shuffle)
    y_t = torch.from_numpy(y.astype(np.float32))
    return DataLoader(TensorDataset(X_t, y_t), batch_size=batch_size, shuffle=shuffle)


def predict(model: nn.Module, X: np.ndarray, batch_size: int, device: torch.device) -> np.ndarray:
    model.eval()
    preds = []
    with torch.no_grad():
        for (xb,) in tensor_loader(X, None, batch_size, False):
            preds.append(model(xb.to(device)).detach().cpu().numpy())
    return np.concatenate(preds).astype(np.float32)


def pretrain_dae(
    model: nn.Module,
    X_train: np.ndarray,
    args: argparse.Namespace,
    device: torch.device,
) -> None:
    if not isinstance(model, DenoisingAutoEncoderRegressor) or args.dae_pretrain_epochs <= 0:
        return
    model.train()
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loss_fn = nn.MSELoss()
    for _ in range(args.dae_pretrain_epochs):
        for (xb,) in tensor_loader(X_train, None, args.batch_size, True):
            xb = xb.to(device)
            noise_mask = torch.rand_like(xb) > args.dae_noise
            recon = model.reconstruct(xb * noise_mask)
            loss = loss_fn(recon, xb)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
            opt.step()


def train_regressor(
    model: nn.Module,
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_val: np.ndarray,
    y_val: np.ndarray,
    args: argparse.Namespace,
    device: torch.device,
) -> dict:
    model.to(device)
    pretrain_dae(model, X_train, args, device)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    loss_fn = nn.SmoothL1Loss(beta=args.huber_beta)
    scaler = torch.amp.GradScaler("cuda", enabled=args.amp and device.type == "cuda")
    best_state = None
    best_val = float("inf")
    stale = 0
    history = []
    for epoch in range(1, args.epochs + 1):
        model.train()
        train_losses = []
        for xb, yb in tensor_loader(X_train, y_train, args.batch_size, True):
            xb = xb.to(device)
            yb = yb.to(device)
            opt.zero_grad(set_to_none=True)
            with torch.amp.autocast("cuda", enabled=args.amp and device.type == "cuda"):
                loss = loss_fn(model(xb), yb)
            scaler.scale(loss).backward()
            scaler.unscale_(opt)
            nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
            scaler.step(opt)
            scaler.update()
            train_losses.append(float(loss.detach().cpu()))
        val_pred = predict(model, X_val, args.batch_size, device)
        val_loss = float(np.mean(np.abs(val_pred - y_val)))
        history.append({"epoch": epoch, "train_loss": float(np.mean(train_losses)), "val_mae_target": val_loss})
        if val_loss < best_val - args.min_delta:
            best_val = val_loss
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            stale = 0
        else:
            stale += 1
        if stale >= args.patience:
            break
    if best_state is not None:
        model.load_state_dict(best_state)
    return {"best_val_mae_target": best_val, "epochs_ran": len(history), "history": history}


def build_cv(groups: np.ndarray, seed: int):
    unique_groups = np.unique(groups)
    if len(unique_groups) >= 2:
        cv = GroupKFold(n_splits=min(5, len(unique_groups)))
        return "GroupKFold_dataset_batch", list(cv.split(np.zeros(len(groups)), groups=groups))
    cv = KFold(n_splits=5, shuffle=True, random_state=seed)
    return "KFold", list(cv.split(np.zeros(len(groups))))


def preprocess_fold(
    X_raw: np.ndarray,
    y_fit: np.ndarray,
    y_days_filter: np.ndarray,
    groups: np.ndarray,
    feature_ids: np.ndarray,
    train_idx: np.ndarray,
    apply_idx: np.ndarray,
    args: argparse.Namespace,
    seed: int,
):
    train_presence = np.isfinite(X_raw[train_idx]).mean(axis=0) >= args.min_train_feature_presence
    stability = train_dataset_stability_mask(
        X_raw[train_idx],
        groups[train_idx],
        args.min_train_group_feature_presence,
        args.max_train_dataset_mean_shift,
    )
    agebin = train_agebin_presence_mask(
        X_raw[train_idx],
        y_days_filter[train_idx] / 7,
        args.min_train_agebin_feature_presence,
        args.min_train_agebin_samples,
    )
    mask = train_presence & stability & agebin
    presence_idx = np.flatnonzero(mask)
    if len(presence_idx) == 0:
        raise RuntimeError("No features pass fold filters.")
    top_local_idx = select_features(X_raw[train_idx][:, presence_idx], y_fit[train_idx], args.n_feature_prefilter)
    top_idx = presence_idx[top_local_idx]
    if args.architecture == "cnn1d":
        order = sorted(range(len(top_idx)), key=lambda pos: parse_region_start(feature_ids[top_idx[pos]]))
        top_idx = top_idx[order]
    imputer = SimpleImputer(strategy=args.imputation)
    X_train = imputer.fit_transform(X_raw[train_idx][:, top_idx]).astype(np.float32)
    X_apply = imputer.transform(X_raw[apply_idx][:, top_idx]).astype(np.float32)
    scaler = make_scaler(args.preprocess, seed, len(train_idx))
    if scaler is not None:
        X_train = scaler.fit_transform(X_train).astype(np.float32)
        X_apply = scaler.transform(X_apply).astype(np.float32)
    artifacts = {
        "selected_feature_indices": top_idx,
        "selected_feature_ids": feature_ids[top_idx].astype(str).tolist(),
        "imputer": imputer,
        "scaler": scaler,
    }
    counts = {
        "n_features_presence": int(train_presence.sum()),
        "n_features_stability": int(stability.sum()),
        "n_features_agebin": int(agebin.sum()),
        "n_features_after_all_filters": int(mask.sum()),
        "n_features_selected": int(len(top_idx)),
    }
    return X_train, X_apply, artifacts, counts


def run_groupkfold(args: argparse.Namespace, out_dir: Path) -> tuple[dict, pd.DataFrame]:
    matrix, meta = load_data(Path(args.matrix_path), Path(args.metadata_path), "all")
    feature_ids = matrix.index.astype(str).to_numpy()
    X_raw = matrix.T.values.astype(np.float32)
    y_days = meta["age_days"].values.astype(np.float32)
    y_fit = build_target(y_days, args.target_transform)
    y_days_filter = y_days.copy()
    if args.randomize_labels:
        rng = np.random.default_rng(args.random_seed)
        perm = rng.permutation(len(y_fit))
        y_fit = y_fit[perm]
        y_days_filter = y_days_filter[perm]
    groups = meta["dataset_batch"].fillna("unknown").astype(str).values
    cv_strategy, cv_iter = build_cv(groups, args.random_seed)
    y_pred_fit = np.full(len(meta), np.nan, dtype=np.float32)
    fold_rows = []
    for fold, (train_idx, test_idx) in enumerate(cv_iter, start=1):
        seed_everything(args.random_seed + fold)
        X_train_all, X_test, _, counts = preprocess_fold(
            X_raw, y_fit, y_days_filter, groups, feature_ids, train_idx, test_idx, args, args.random_seed + fold
        )
        tr_idx_local, val_idx_local = train_val_split(np.arange(len(train_idx)), y_fit[train_idx], args.random_seed + fold)
        model_cfg = model_config_from_args(args)
        model = make_model(X_train_all.shape[1], model_cfg)
        train_info = train_regressor(
            model,
            X_train_all[tr_idx_local],
            y_fit[train_idx][tr_idx_local],
            X_train_all[val_idx_local],
            y_fit[train_idx][val_idx_local],
            args,
            torch.device(args.device),
        )
        y_pred_fit[test_idx] = predict(model, X_test, args.batch_size, torch.device(args.device))
        fold_rows.append({"fold": fold, **counts, **{k: v for k, v in train_info.items() if k != "history"}})
        print(f"Fold {fold}: selected={counts['n_features_selected']} epochs={train_info['epochs_ran']}")
    y_pred_days, y_pred_weeks = inverse_target(y_pred_fit, args.target_transform, float(y_days.max()))
    metrics = compute_metrics(y_days / 7, y_pred_weeks)
    pred_df = prediction_frame(meta, y_days, y_pred_days)
    result = common_result(args, metrics, X_raw.shape, cv_strategy, fold_rows)
    return result, pred_df


def run_heldout(args: argparse.Namespace, out_dir: Path) -> tuple[dict, pd.DataFrame]:
    train_matrix, train_meta = load_data(Path(args.matrix_path), Path(args.metadata_path), f"all_except:{args.heldout_dataset}")
    test_matrix, test_meta = load_data(Path(args.matrix_path), Path(args.metadata_path), args.heldout_dataset)
    train_matrix, test_matrix = train_matrix.align(test_matrix, join="inner", axis=0)
    feature_ids = train_matrix.index.astype(str).to_numpy()
    X_train_raw = train_matrix.T.values.astype(np.float32)
    X_test_raw = test_matrix.T.values.astype(np.float32)
    X_all_raw = np.vstack([X_train_raw, X_test_raw])
    train_idx = np.arange(len(train_meta))
    test_idx = np.arange(len(train_meta), len(train_meta) + len(test_meta))
    meta_all = pd.concat([train_meta, test_meta], ignore_index=True)
    y_days_all = meta_all["age_days"].values.astype(np.float32)
    y_fit_all = build_target(y_days_all, args.target_transform)
    groups = meta_all["dataset_batch"].fillna("unknown").astype(str).values
    seed_everything(args.random_seed)
    X_train_all, X_test, _, counts = preprocess_fold(
        X_all_raw, y_fit_all, y_days_all, groups, feature_ids, train_idx, test_idx, args, args.random_seed
    )
    tr_idx_local, val_idx_local = train_val_split(np.arange(len(train_idx)), y_fit_all[train_idx], args.random_seed)
    model_cfg = model_config_from_args(args)
    model = make_model(X_train_all.shape[1], model_cfg)
    train_info = train_regressor(
        model,
        X_train_all[tr_idx_local],
        y_fit_all[train_idx][tr_idx_local],
        X_train_all[val_idx_local],
        y_fit_all[train_idx][val_idx_local],
        args,
        torch.device(args.device),
    )
    pred_fit = predict(model, X_test, args.batch_size, torch.device(args.device))
    y_pred_days, y_pred_weeks = inverse_target(pred_fit, args.target_transform, float(y_days_all.max()))
    metrics = compute_metrics(test_meta["age_days"].values.astype(np.float32) / 7, y_pred_weeks)
    pred_df = prediction_frame(test_meta, test_meta["age_days"].values.astype(np.float32), y_pred_days)
    result = common_result(args, metrics, X_all_raw.shape, f"LODO_{args.heldout_dataset}", [{**counts, **train_info}])
    return result, pred_df


def run_final(args: argparse.Namespace, out_dir: Path) -> tuple[dict, pd.DataFrame]:
    matrix, meta = load_data(Path(args.matrix_path), Path(args.metadata_path), "all")
    feature_ids = matrix.index.astype(str).to_numpy()
    X_raw = matrix.T.values.astype(np.float32)
    y_days = meta["age_days"].values.astype(np.float32)
    y_fit = build_target(y_days, args.target_transform)
    groups = meta["dataset_batch"].fillna("unknown").astype(str).values
    all_idx = np.arange(len(meta))
    X_train_all, _, artifacts, counts = preprocess_fold(
        X_raw, y_fit, y_days, groups, feature_ids, all_idx, all_idx, args, args.random_seed
    )
    tr_idx_local, val_idx_local = train_val_split(all_idx, y_fit, args.random_seed)
    seed_everything(args.random_seed)
    model_cfg = model_config_from_args(args)
    model = make_model(X_train_all.shape[1], model_cfg)
    train_info = train_regressor(
        model,
        X_train_all[tr_idx_local],
        y_fit[tr_idx_local],
        X_train_all[val_idx_local],
        y_fit[val_idx_local],
        args,
        torch.device(args.device),
    )
    pred_fit = predict(model, X_train_all, args.batch_size, torch.device(args.device))
    y_pred_days, y_pred_weeks = inverse_target(pred_fit, args.target_transform, float(y_days.max()))
    metrics = compute_metrics(y_days / 7, y_pred_weeks)
    torch.save(
        {
            "state_dict": model.state_dict(),
            "model_config": asdict(model_cfg),
            "args": vars(args),
            "selected_feature_ids": artifacts["selected_feature_ids"],
        },
        out_dir / "final_deep_clock_model.pt",
    )
    joblib.dump(artifacts, out_dir / "final_preprocess_artifacts.joblib")
    pd.DataFrame(
        {"feature_id": artifacts["selected_feature_ids"], "feature_index": artifacts["selected_feature_indices"]}
    ).to_csv(out_dir / "selected_features.csv", index=False)
    pred_df = prediction_frame(meta, y_days, y_pred_days)
    result = common_result(args, metrics, X_raw.shape, "final_all_data_apparent", [{**counts, **train_info}])
    result["warning"] = "Apparent all-data fit metrics are optimistic and not held-out evidence."
    return result, pred_df


def prediction_frame(meta: pd.DataFrame, y_days: np.ndarray, y_pred_days: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": meta["sample_id"].values,
            "dataset_batch": meta.get("dataset_batch", pd.Series(["unknown"] * len(meta))).values,
            "tissue": meta.get("tissue", pd.Series(["unknown"] * len(meta))).values,
            "intervention": meta.get("intervention", pd.Series(["control"] * len(meta))).values,
            "age_days_true": y_days,
            "age_days_pred": y_pred_days,
            "age_weeks_true": y_days / 7,
            "age_weeks_pred": y_pred_days / 7,
            "residual_weeks": y_days / 7 - y_pred_days / 7,
        }
    )


def model_config_from_args(args: argparse.Namespace) -> ModelConfig:
    return ModelConfig(
        architecture=args.architecture,
        width=args.width,
        depth=args.depth,
        dropout=args.dropout,
        bottleneck=args.bottleneck,
        cnn_channels=args.cnn_channels,
        cnn_kernel=args.cnn_kernel,
        token_dim=args.token_dim,
        transformer_heads=args.transformer_heads,
        transformer_layers=args.transformer_layers,
        dae_latent=args.dae_latent,
    )


def common_result(
    args: argparse.Namespace,
    metrics: dict,
    matrix_shape: tuple[int, int],
    cv_strategy: str,
    fold_rows: list[dict],
) -> dict:
    return {
        "architecture": args.architecture,
        "n_feature_prefilter": args.n_feature_prefilter,
        "imputation": args.imputation,
        "preprocess": args.preprocess,
        "target_transform": args.target_transform,
        "min_train_feature_presence": args.min_train_feature_presence,
        "min_train_group_feature_presence": args.min_train_group_feature_presence,
        "max_train_dataset_mean_shift": args.max_train_dataset_mean_shift,
        "min_train_agebin_feature_presence": args.min_train_agebin_feature_presence,
        "model_config": asdict(model_config_from_args(args)),
        "training": {
            "epochs": args.epochs,
            "patience": args.patience,
            "batch_size": args.batch_size,
            "lr": args.lr,
            "weight_decay": args.weight_decay,
            "amp": args.amp,
            "device": args.device,
        },
        "cv_strategy": cv_strategy,
        "random_seed": args.random_seed,
        "randomize_labels": args.randomize_labels,
        "eval_mode": args.eval_mode,
        "heldout_dataset": args.heldout_dataset,
        "n_samples": int(matrix_shape[0]),
        "n_features_total": int(matrix_shape[1]),
        "fold_rows": fold_rows,
        "leakage_controls": [
            "cv_split_before_feature_selection",
            "feature_presence_filter_fit_train_fold_only",
            "dataset_stability_filter_fit_train_fold_only",
            "agebin_presence_filter_fit_train_fold_only",
            "feature_selection_fit_train_fold_only",
            "imputer_fit_train_fold_only",
            "scaler_fit_train_fold_only",
            "model_fit_train_fold_only",
        ],
        **metrics,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix_path", required=True)
    parser.add_argument("--metadata_path", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--eval_mode", default="groupkfold", choices=["groupkfold", "heldout", "final"])
    parser.add_argument("--heldout_dataset", default=None)
    parser.add_argument("--architecture", choices=["res_mlp", "wide_deep", "cnn1d", "tab_transformer", "dae_mlp"], default="res_mlp")
    parser.add_argument("--n_feature_prefilter", default="1000")
    parser.add_argument("--imputation", default="median", choices=["median", "mean"])
    parser.add_argument("--preprocess", default="standard", choices=["standard", "robust", "quantile_uniform", "none"])
    parser.add_argument("--min_train_feature_presence", type=float, default=0.8)
    parser.add_argument("--min_train_group_feature_presence", type=float, default=0.8)
    parser.add_argument("--max_train_dataset_mean_shift", type=float, default=0.125)
    parser.add_argument("--min_train_agebin_feature_presence", type=float, default=0.8)
    parser.add_argument("--min_train_agebin_samples", type=int, default=8)
    parser.add_argument("--target_transform", default="log1p_days", choices=["log1p_days", "linear_weeks"])
    parser.add_argument("--width", type=int, default=256)
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--dropout", type=float, default=0.25)
    parser.add_argument("--bottleneck", type=int, default=64)
    parser.add_argument("--cnn_channels", type=int, default=64)
    parser.add_argument("--cnn_kernel", type=int, default=15)
    parser.add_argument("--token_dim", type=int, default=16)
    parser.add_argument("--transformer_heads", type=int, default=4)
    parser.add_argument("--transformer_layers", type=int, default=2)
    parser.add_argument("--dae_latent", type=int, default=64)
    parser.add_argument("--dae_noise", type=float, default=0.1)
    parser.add_argument("--dae_pretrain_epochs", type=int, default=20)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--patience", type=int, default=30)
    parser.add_argument("--min_delta", type=float, default=1e-4)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-3)
    parser.add_argument("--huber_beta", type=float, default=0.5)
    parser.add_argument("--grad_clip", type=float, default=1.0)
    parser.add_argument("--random_seed", type=int, default=42)
    parser.add_argument("--randomize_labels", action="store_true")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    parser.add_argument("--amp", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    if args.eval_mode == "heldout" and not args.heldout_dataset:
        raise SystemExit("--heldout_dataset is required for heldout mode")
    if args.device == "cuda" and not torch.cuda.is_available():
        raise SystemExit("CUDA requested but torch.cuda.is_available() is false")
    if args.architecture == "tab_transformer" and str(args.n_feature_prefilter).lower() == "all":
        raise SystemExit("tab_transformer requires finite n_feature_prefilter")

    if args.eval_mode == "groupkfold":
        result, pred_df = run_groupkfold(args, out_dir)
    elif args.eval_mode == "heldout":
        result, pred_df = run_heldout(args, out_dir)
    else:
        result, pred_df = run_final(args, out_dir)
    result["exec_time_sec"] = round(time.time() - started, 1)
    result["torch_version"] = torch.__version__
    result["cuda_available"] = torch.cuda.is_available()
    result["cuda_device"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    (out_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str), encoding="utf-8")
    pred_df.to_csv(out_dir / "predictions.csv", index=False)
    print(
        f"[{args.eval_mode}:{args.architecture}] R={result['pearson_r']} "
        f"MAE={result['mae_weeks']} R2={result['r2']} out={out_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
