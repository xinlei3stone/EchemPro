"""提取 CS Studio 扫速 CV txt 的最后一圈, 存 CSV, 并按 LSV 对比图样式绘图.

与 RuNiMoHER/scripts/extract_scanrate_cv.py 对齐; 数据格式相同
(CS Studio 导出: 首行元数据 + 表头 E(V)/i(A/cm²)/T(s))。
电流已是 A/cm² (面积 1 cm²), 无需面积换算。

数据默认位于 ../RuHeterojunctionHER/DATA/EC/<日期>/<样品>-扫速/,
推荐用命令行参数指定 (不必改源码):
    ../RuHeterojunctionHER/.venv/bin/python extract_scanrate_cv.py \
        --date 20260921 --sample Ru_NiP_annealing
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from echem_process import COLORS, origin_axes, format_ticks_legend, read_cs_scan_txt

# 数据根目录: 默认指向同级的 RuHeterojunctionHER (实验数据所在项目)
DEFAULT_ROOT = Path(__file__).resolve().parents[1] / "RuHeterojunctionHER"
DATE = "YYYYMMDD"        # 实验日期文件夹
SAMPLE = "样品名"        # 样品名 (文件夹为 <SAMPLE>-扫速)
FILES = [20, 40, 60, 80, 100]  # 扫速 mV/s, 按实际文件名 (如 20.txt)
CURVE_COLOR = "#1f77b4"

PH = 14.0       # 1 M KOH
EREF = 0.098    # Hg/HgO vs. SHE / V


def last_loop(df: pd.DataFrame) -> pd.DataFrame:
    """按 E 转向点识别循环, 取最后一个低电位转折点之后的完整一圈"""
    E = df["E(V)"].to_numpy()
    d = np.sign(np.diff(E))
    turns = np.where(np.diff(d) != 0)[0] + 1
    lo_turns = [t for t in turns if E[t] < E.min() + 0.05]
    start = lo_turns[-1] if lo_turns else 0
    return df.iloc[start:].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help=f"数据根目录 (默认: {DEFAULT_ROOT})")
    ap.add_argument("--date", default=DATE,
                    help="实验日期文件夹, 如 20260921")
    ap.add_argument("--sample", default=SAMPLE,
                    help="样品名 (文件夹为 <SAMPLE>-扫速)")
    ap.add_argument("--files", type=int, nargs="+", default=FILES,
                    help="扫速列表 mV/s, 对应 <rate>.txt")
    args = ap.parse_args()

    data_dir = args.root / f"DATA/EC/{args.date}/{args.sample}-扫速"
    if not data_dir.is_dir():
        raise SystemExit(
            f"数据目录不存在: {data_dir}\n"
            "请用 --date/--sample 指定 (或用 --root 指定数据根目录)")

    curves = []
    for rate in args.files:
        txt = data_dir / f"{rate}.txt"
        df = read_cs_scan_txt(txt)
        seg = last_loop(df)
        out_csv = txt.with_name(f"{rate}_last_loop.csv")
        seg.to_csv(out_csv, index=False)
        print(f"Saved: {out_csv}  ({len(seg)} points)")

        E_RHE = seg["E(V)"] + 0.0592 * PH + EREF
        j = seg["i(A/cm²)"] * 1000.0  # A/cm² -> mA/cm²
        curves.append((E_RHE, j, f"{rate} mV s⁻¹",
                       COLORS[len(curves) % len(COLORS)]))

    fig, ax = origin_axes()
    for E_RHE, j, label, color in curves:
        ax.plot(E_RHE, j, color=color, lw=2.5, label=label)
    ax.axhline(0, color="0.5", lw=1.0)
    ax.set_xlabel("Potential (V vs. RHE)", fontsize=36, fontweight="bold")
    ax.set_ylabel("Current density (mA cm⁻²)", fontsize=36, fontweight="bold")
    # y 轴按数据自适应: 最大 |j| 加 15% 边距, 图例区再留顶部空间
    j_max = max(float(j.abs().max()) for _, j, _, _ in curves)
    ax.set_ylim(-1.15 * j_max, 1.55 * j_max)
    ax.minorticks_on()
    # 三列图例置于顶部, 避免与曲线重叠
    ax.tick_params(labelsize=24)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    leg = ax.legend(frameon=False, fontsize=18, loc="upper center",
                    ncol=3, handlelength=1.2, handletextpad=0.5,
                    columnspacing=0.8, borderaxespad=0.3)
    for text in leg.get_texts():
        text.set_fontweight("bold")
    out_png = data_dir.parent / f"CV_scanrate_comparison_{args.sample}.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_png}")

    # ---- H 脱附峰 Ep vs. 扫速 (定性比较, 参照论文图 4-7c 样式) ----
    print("\nH desorption peak (anodic branch):")
    eps = []
    for (E_RHE, j, label, color), rate in zip(curves, args.files):
        # 阳极分支: E 上升段
        anodic = j.index[np.diff(E_RHE.to_numpy(), prepend=E_RHE.iloc[0]) > 0]
        E_a, j_a = E_RHE.loc[anodic], j.loc[anodic]
        win = (E_a > 0.1) & (E_a < 0.4)  # 脱附峰窗口 / V vs. RHE
        j_sm = j_a[win].rolling(15, center=True).mean()
        Ep = E_a[win].iloc[int(np.nanargmax(j_sm.to_numpy()))]
        eps.append(Ep)
        print(f"  {label}: Ep = {Ep:.4f} V vs. RHE")

    # 峰位钉在窗口边缘 = 窗口内单调无峰, 跳过 Ep 分析
    if any(abs(ep - 0.1) < 0.01 or abs(ep - 0.4) < 0.01 for ep in eps):
        print("警告: 峰位钉在窗口边缘, 该样品无可分辨 H 脱附峰, 跳过 Ep 图")
        return

    slope, intercept = np.polyfit(args.files, eps, 1)
    r = np.corrcoef(args.files, eps)[0, 1]
    print(f"\nSlope = {slope:.2e} V per (mV/s), R² = {r**2:.4f}")

    summary = pd.DataFrame({
        "v/(mV s-1)": args.files, "Ep/V vs RHE": eps})
    out_csv = data_dir.parent / f"Hdes_peak_Ep_vs_scanrate_{args.sample}.csv"
    summary.to_csv(out_csv, index=False)
    print(f"Saved: {out_csv}")

    fig, ax = origin_axes()
    ax.plot(args.files, eps, "o", color=CURVE_COLOR, markersize=10,
            markeredgecolor="white", markeredgewidth=1.0, label=args.sample)
    x_fit = np.linspace(min(args.files) - 20, max(args.files) + 20, 50)
    ax.plot(x_fit, intercept + slope * x_fit, "-", color=CURVE_COLOR, lw=2.0)
    exp = int(np.floor(np.log10(abs(slope))))
    mant = slope / 10**exp
    ax.text(0.40, 0.25, f"Slope = {mant:.2f}*10$^{{{exp}}}$",
            transform=ax.transAxes, fontsize=24, fontweight="bold",
            color=CURVE_COLOR)
    ax.set_xlabel("Scan rate (mV s⁻¹)", fontsize=36, fontweight="bold")
    ax.set_ylabel("Hydrogen Desorption\nPeak Position (V)",
                  fontsize=36, fontweight="bold")
    ax.set_xlim(min(args.files) - 20, max(args.files) + 20)
    ax.minorticks_on()
    format_ticks_legend(ax, "upper left")
    out_png = data_dir.parent / f"Hdes_Ep_vs_scanrate_{args.sample}.png"
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
