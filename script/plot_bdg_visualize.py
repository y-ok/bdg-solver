import argparse
import glob
import math
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go


KIND_TO_PATTERN = {
    "density": "bdg_density_*.csv",
    "orderParameter": "bdg_orderParameter_*.csv",
    "magnetization": "bdg_magnetization_*.csv",
}

KIND_TO_COLUMN = {
    "density": "density",
    "orderParameter": "orderParameter",
    "magnetization": "magnetization",
}


def pick_latest_file(directory: Path, pattern: str) -> Path:
    paths = [Path(p) for p in glob.glob(str(directory / pattern))]
    if not paths:
        raise FileNotFoundError(f"Not found: {pattern} in {directory}")
    return max(paths, key=lambda p: p.stat().st_mtime)


def to_grid(df: pd.DataFrame, value_col: str):
    if "x" not in df.columns or "y" not in df.columns:
        raise ValueError(f"CSV must have x,y columns. got={list(df.columns)}")

    xs = df["x"].to_numpy(dtype=int)
    ys = df["y"].to_numpy(dtype=int)
    vals = df[value_col].to_numpy(dtype=float)

    lx = int(xs.max()) + 1
    ly = int(ys.max()) + 1

    grid = np.full((ly, lx), np.nan, dtype=float)  # [y, x]
    grid[ys, xs] = vals
    return grid, lx, ly


def robust_range(values: np.ndarray, low_q: float = 0.02, high_q: float = 0.98):
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return -1.0, 1.0
    lo = float(np.quantile(finite, low_q))
    hi = float(np.quantile(finite, high_q))
    if lo == hi:
        eps = 1e-12 if lo == 0.0 else abs(lo) * 1e-3
        return lo - eps, hi + eps
    return lo, hi


def nice_step_125(x: float) -> float:
    """x を 1-2-5 系のステップに切り上げ"""
    if x <= 0.0 or not math.isfinite(x):
        return 1.0
    exp = 10.0 ** math.floor(math.log10(x))
    f = x / exp
    if f <= 1.0:
        return 1.0 * exp
    if f <= 2.0:
        return 2.0 * exp
    if f <= 5.0:
        return 5.0 * exp
    return 10.0 * exp


def choose_zrange(kind: str, grid: np.ndarray, robust: bool = True):
    finite = grid[np.isfinite(grid)]
    if finite.size == 0:
        # 値が無い場合の保険
        return 0.0, 1.0

    max_val = float(finite.max())

    # ★ここだけ変更：zmax を「次の目盛り」まで伸ばす
    # 目盛りがだいたい 8 個くらいになるようにステップを決める（0.7なら step=0.1 になりやすい）
    step = nice_step_125(max_val / 7.0)
    zmax = step * math.ceil(max_val / step)

    # max_val が目盛りぴったりなら、次の目盛り（例: 0.7 -> 0.8）を出すために1段上げる
    if math.isclose(zmax, max_val, rel_tol=1e-12, abs_tol=1e-12):
        zmax += step

    # 下限ルール（既存のまま）
    if kind in ("density", "magnetization"):
        zmin = 0.0
    else:
        zmin = float(np.quantile(finite, 0.02)) if robust else float(finite.min())

    # 範囲が潰れるのを防ぐ（既存のまま）
    if zmin == zmax:
        eps = 1e-12 if zmax == 0.0 else abs(zmax) * 1e-3
        zmin -= eps
        zmax += eps

    return zmin, zmax

def write_surface_html(kind: str, csv_path: Path, out_dir: Path, robust: bool):
    df = pd.read_csv(csv_path)
    col = KIND_TO_COLUMN[kind]
    if col not in df.columns:
        raise ValueError(f"Missing column '{col}' in {csv_path.name}. got={list(df.columns)}")

    grid, lx, ly = to_grid(df, col)
    zmin, zmax = choose_zrange(kind, grid, robust=robust)

    x = np.arange(lx)
    y = np.arange(ly)

    fig = go.Figure()
    fig.add_trace(
        go.Surface(
            x=x,
            y=y,
            z=grid,
            cmin=zmin,
            cmax=zmax,
            colorbar=dict(title=col),
            contours={"z": {"show": True, "usecolormap": True, "project_z": True}},
            hovertemplate="x=%{x}<br>y=%{y}<br>z=%{z}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"{kind} 3D (from {csv_path.name})",
        scene=dict(
            xaxis_title="x",
            yaxis_title="y",
            zaxis_title=col,
            zaxis=dict(range=[zmin, zmax]),
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )

    out_path = make_out_path(out_dir, "3d", csv_path)
    fig.write_html(out_path, include_plotlyjs="cdn")
    return out_path, (zmin, zmax), (lx, ly)


def write_heatmap_html(kind: str, csv_path: Path, out_dir: Path, robust: bool):
    df = pd.read_csv(csv_path)
    col = KIND_TO_COLUMN[kind]
    if col not in df.columns:
        raise ValueError(f"Missing column '{col}' in {csv_path.name}. got={list(df.columns)}")

    grid, lx, ly = to_grid(df, col)
    zmin, zmax = choose_zrange(kind, grid, robust=robust)

    fig = go.Figure(
        data=go.Heatmap(
            z=grid,
            zmin=zmin,
            zmax=zmax,
            hovertemplate="x=%{x}<br>y=%{y}<br>z=%{z}<extra></extra>",
            colorbar=dict(title=col),
        )
    )
    fig.update_layout(
        title=f"{kind} 2D (from {csv_path.name})",
        xaxis_title="x",
        yaxis_title="y",
        yaxis=dict(scaleanchor="x", scaleratio=1),
        margin=dict(l=40, r=10, t=40, b=40),
    )

    out_path = make_out_path(out_dir, "2d", csv_path)
    fig.write_html(out_path, include_plotlyjs="cdn")
    return out_path, (zmin, zmax), (lx, ly)

def make_out_path(out_dir: Path, dim: str, csv_path: Path) -> Path:
    stem = csv_path.stem  # 例: "bdg_density_N=100_h=0.00"
    if stem.startswith("bdg_"):
        stem = stem[len("bdg_"):]  # 例: "density_N=100_h=0.00"
    return out_dir / f"bdg_plot_{dim}_{stem}.html"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="JavaのCSV出力フォルダ")
    ap.add_argument("--out", default=None, help="出力フォルダ（省略時は --dir と同じ）")
    ap.add_argument("--no-robust", action="store_true",
                    help="ロバスト推定を無効化（min/maxを使用）")
    ap.add_argument("--only3d", action="store_true", help="3Dのみ出力")
    ap.add_argument("--only2d", action="store_true", help="2Dのみ出力")
    args = ap.parse_args()

    if args.only3d and args.only2d:
        raise ValueError("--only3d と --only2d は同時に指定できません")

    in_dir = Path(args.dir)
    out_dir = Path(args.out) if args.out else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    robust = not args.no_robust

    chosen = {k: pick_latest_file(in_dir, pat) for k, pat in KIND_TO_PATTERN.items()}

    print("Selected files:")
    for k, p in chosen.items():
        print(f"  {k}: {p.name}")

    print("\nGenerating interactive HTML...")
    for kind, csv_path in chosen.items():
        if not args.only2d:
            out3d, (zmin, zmax), (lx, ly) = write_surface_html(kind, csv_path, out_dir, robust)
            print(f"  saved {out3d.name}  zrange=({zmin:.6g}, {zmax:.6g})  size={lx}x{ly}")

        if not args.only3d:
            out2d, (zmin, zmax), (lx, ly) = write_heatmap_html(kind, csv_path, out_dir, robust)
            print(f"  saved {out2d.name}  zrange=({zmin:.6g}, {zmax:.6g})  size={lx}x{ly}")

    print("\nDone.")


if __name__ == "__main__":
    main()
