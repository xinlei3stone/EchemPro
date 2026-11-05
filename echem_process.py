#!/usr/bin/env python3
"""通用电化学数据处理脚本: CS Studio txt 读取 + EIS / LSV / Tafel 批量对比

与 RuNiMoHER/scripts/echem_process.py 对齐, 输入层改为本项目的
CS Studio 文本导出 (20260921 起的数据均为该格式):
  - CV / LSV: <样品>-CV.txt / <样品>-LSV.txt, 首行为 CSStudioFile 元数据,
    其后表头 E(V) / i(A/cm²) / T(s) (tab 分隔, i 已是电流密度, 面积 1 cm²)
  - EIS: <样品>-EIS.txt, 表头 Freq(Hz) / Z'(Ohm.cm²) / Z''(Ohm.cm²) / ...,
    Z'' 电容性为负 (与 BioLogic -Im(Z) 同号约定)
  - 旧 BioLogic .mpr 仍可转换 (需 .venv 装 galvani, 未装时给出提示)

功能 (对每个样品文件夹):
  1. 自动识别 -EIS / -LSV 的 .txt (或 .mpr 转 .csv;
     已内置 galvani 列 ID 182 = "step time/s" <f8 补丁, 见 echemdata/galvani#112)
  2. EIS: DRT (Tikhonov 正则化 + NNLS) 反演, 得到 R∞ / R_pol / γ(τ) 与峰位
     (λ 默认 1e-2, 可留空走 L-curve 或手动指定) -> 单样品 Nyquist/Bode
     (实测点+DRT 拟合线) 与 DRT 图, 样品间 Nyquist 对比图 + DRT 汇总/峰位 CSV
  3. LSV: 取回程段 (扫描顶点后), 电位换算 RHE (Hg/HgO, 默认 pH=14,
     90% iR 补偿, Rs 取该样品 DRT 的高频电阻 R∞, EIS 阻抗已是 Ω·cm² 与
     i(A/cm²) 直接相乘得 iR 压降), 保存 *_RHE.csv, 画 LSV 对比图
  4. Tafel: 在指定 Log|J| 窗口内拟合并画对比图 (默认 0.8-1.6)

LSV 说明: 本组 LSV 常以 CV 方式记录完整循环 (去程到顶点后返回),
脚本自动检测回程段作为 LSV 数据; 单程数据则全取。

用法 (在 EchemPro-src/ 下, 复用 RuHeterojunctionHER 的 venv):
  ../RuHeterojunctionHER/.venv/bin/python echem_process.py \
      ../RuHeterojunctionHER/DATA/EC/20260921 Ru_NiP_annealing
  ../RuHeterojunctionHER/.venv/bin/python echem_process.py \
      ../RuHeterojunctionHER/DATA/EC/20260921 Ru_NiP_annealing \
      --labels "Ru/NiP annealing" --tag NiP-annealing
  ../RuHeterojunctionHER/.venv/bin/python echem_process.py \
      ../RuHeterojunctionHER/DATA/EC/<date> S1 S4 \
      --rs 1.414 1.460          # 手动指定 Rs, 跳过 DRT 的 R∞ 取值
"""

import argparse
import os
import re
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from drt_analysis import run_drt

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

plt.rcParams['savefig.dpi'] = 300

COLORS = ["black", "red", "blue", "green", "orange", "purple", "brown", "cyan"]

TECHNIQUE_TOKENS = {
    "eis": {"EIS", "PEIS", "GEIS"},
    "lsv": {"LSV"},
    "cv": {"CV"},
}


def detect_technique(stem: str):
    """从文件名(不含扩展名)识别电化学技术, 返回类别或 None"""
    tokens = set(re.split(r"[^A-Za-z0-9]+", stem.upper()))
    for category, keys in TECHNIQUE_TOKENS.items():
        if tokens & keys:
            return category
    return None


# ---------------- CS Studio txt 读取 ----------------

def read_cs_scan_txt(path: Path) -> pd.DataFrame:
    """读 CS Studio CV/LSV txt: 首行 CSStudioFile 元数据, 其后表头
    E(V) / i(A/cm²) / T(s); 逐行数值解析, 自动跳过元数据/表头/空行,
    强转数值失败的行为 NaN 后丢弃 (兼容 "N CYCLE" 边界标记)"""
    rows = []
    with open(path, "r", encoding="utf-8-sig") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"{path.name} 无有效数据行")
    return pd.DataFrame(rows, columns=["E(V)", "i(A/cm²)", "T(s)"])


def read_cs_eis_txt(path: Path) -> pd.DataFrame:
    """读 CS Studio EIS txt, 归一化为 BioLogic 风格列名:
    freq/Hz, Re(Z)/Ohm, -Im(Z)/Ohm (CS Z'' 电容性为负, 与 BioLogic
    -Im(Z) 同为正值表示电容, 故直接取负)。感抗点 (Z''>0) 保留,
    DRT 模型含电感项 L 可吸收。列名按表头文本匹配, 不依赖列序"""
    rows = []
    with open(path, "r", encoding="utf-8-sig") as f:
        col_map = None
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if col_map is None:
                stripped = [p.strip() for p in parts]
                if stripped and stripped[0].lower().startswith("freq"):
                    col_map = {}
                    for k, name in enumerate(stripped):
                        n = name.lower()
                        if n.startswith("freq"):
                            col_map["f"] = k
                        elif n.startswith("z''"):   # 必须先判 Z'' (它也是 Z' 前缀)
                            col_map["im"] = k
                        elif n.startswith("z'"):
                            col_map["re"] = k
                    missing = {"f", "re", "im"} - col_map.keys()
                    if missing:
                        raise ValueError(f"{path.name} 表头缺列: {missing}")
                continue
            try:
                rows.append((float(parts[col_map["f"]]),
                             float(parts[col_map["re"]]),
                             float(parts[col_map["im"]])))
            except (ValueError, IndexError):
                continue  # 跳过表头残行/空行
    if not rows:
        raise ValueError(f"{path.name} 无有效数据行")
    df = pd.DataFrame(rows, columns=["freq/Hz", "Re(Z)/Ohm", "Zim_raw"])
    df["-Im(Z)/Ohm"] = -df["Zim_raw"]
    return df[["freq/Hz", "Re(Z)/Ohm", "-Im(Z)/Ohm"]]


def mpr_to_csv(mpr_path: Path) -> Path:
    """mpr -> CSV (CSV 已存在且更新时跳过; galvani 未安装时给出安装提示)"""
    csv_path = mpr_path.with_suffix(".csv")
    if csv_path.exists() and csv_path.stat().st_mtime >= mpr_path.stat().st_mtime:
        return csv_path
    try:
        from galvani import BioLogic
    except ImportError:
        raise SystemExit(
            f"ERROR: 转换 {mpr_path.name} 需要 galvani: .venv/bin/pip install galvani\n"
            "       或直接从 EC-Lab/CS Studio 导出 txt 再运行")
    # galvani <= 0.5.0 不认识新版 EC-Lab 的列 ID 182 (step time/s), 按 master 定义补上
    BioLogic.VMPdata_colID_dtype_map.setdefault(182, ("step time/s", "<f8"))
    try:
        df = pd.DataFrame(BioLogic.MPRfile(str(mpr_path)).data)
    except Exception:
        if csv_path.exists():
            print(f"  WARNING: {mpr_path.name} 转换失败, 使用已有 CSV")
            return csv_path
        raise
    df.to_csv(csv_path, index=False)
    return csv_path


def load_eis(path: Path) -> pd.DataFrame:
    """按扩展名读取 EIS 数据并归一化为 BioLogic 风格三列 DataFrame"""
    if path.suffix.lower() == ".txt":
        return read_cs_eis_txt(path)
    return pd.read_csv(path)


def load_scan(path: Path) -> pd.DataFrame:
    """按扩展名读取 CV/LSV 扫描数据, 归一化为 E(V)/i(A/cm²)/T(s);
    galvani 转出的 CSV 用 Ewe/V 与 <I>/mA 列换算 (面积 1 cm²)"""
    if path.suffix.lower() == ".txt":
        return read_cs_scan_txt(path)
    df = pd.read_csv(path)
    return pd.DataFrame({
        "E(V)": df["Ewe/V"],
        "i(A/cm²)": df["<I>/mA"] / 1000.0,
        "T(s)": df["time/s"] if "time/s" in df.columns else np.arange(len(df)),
    })


def find_one(folder: Path, technique: str):
    """在文件夹内找某技术的 -EIS/-LSV 文件: .txt (CS Studio 导出) 优先,
    其次 .csv (galvani 转出) 与 .mpr (BioLogic 原始, 运行时转换);
    多个时取排序后第一个; 跳过 macOS AppleDouble 垃圾文件 (._ 前缀)"""
    for ext in ("*.txt", "*.csv", "*.mpr"):
        hits = sorted(p for p in folder.glob(ext)
                      if not p.name.startswith("._")
                      and detect_technique(p.stem) == technique)
        if hits:
            return hits[0]
    return None


# ---------------- EIS: DRT (Distribution of Relaxation Times) ----------------

def drt_plot_text(drt):
    return (f"R∞ = {drt['r_inf']:.2f} Ω·cm²\n"
            f"R_pol = {drt['r_pol']:.2f} Ω·cm²\n"
            f"λ = {drt['lam']:.1e},  rel. err = {drt['rel_err']:.1%}")


def plot_nyquist(df, drt, out, label=None):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(df["Re(Z)/Ohm"], df["-Im(Z)/Ohm"], "o", color="0.45",
            markersize=5, label="measured")
    ax.plot(drt["z_fit_re"], -drt["z_fit_im"], "b-", linewidth=2,
            label="DRT fit")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Re(Z) / Ω·cm²")
    ax.set_ylabel("-Im(Z) / Ω·cm²")
    if label:
        ax.set_title(f"Nyquist - {label}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    ax.text(0.03, 0.03, drt_plot_text(drt), transform=ax.transAxes,
            va="bottom", fontsize=12,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_bode(df, drt, out):
    z_meas = (df["Re(Z)/Ohm"].to_numpy(dtype=float)
              - 1j * df["-Im(Z)/Ohm"].to_numpy(dtype=float))
    z_fit = drt["z_fit_re"] + 1j * drt["z_fit_im"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    ax1.semilogx(df["freq/Hz"], np.abs(z_meas), "o", color="0.45", markersize=5)
    phase = np.degrees(np.arctan2(z_meas.imag, z_meas.real))
    ax2.semilogx(df["freq/Hz"], phase, "o", color="0.45", markersize=5)
    ax1.semilogx(drt["f"], np.abs(z_fit), "b-", linewidth=2)
    ax2.semilogx(drt["f"], np.degrees(np.angle(z_fit)), "b-", linewidth=2)
    ax1.set_xlabel("Frequency / Hz")
    ax1.set_ylabel("|Z| / Ω·cm²")
    ax1.grid(True, which="both", alpha=0.3)
    ax2.set_xlabel("Frequency / Hz")
    ax2.set_ylabel("Phase / deg")
    ax2.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_drt(drt, out, label=None):
    tau, gamma = drt["tau"], drt["gamma"]
    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(tau, gamma, color="tab:blue", linewidth=1.8)
    ax.fill_between(tau, 0, gamma, color="tab:blue", alpha=0.25)
    for k, i in enumerate(drt["peaks"]):
        left = (k % 2 == 0)
        ax.annotate(f"τ={tau[i]:.2e} s\n(f*={1/(2*np.pi*tau[i]):.2g} Hz)",
                    (tau[i], gamma[i]), textcoords="offset points",
                    xytext=(-14 if left else 14, 14 + 12 * (k % 3)),
                    ha="right" if left else "left", va="bottom", fontsize=8,
                    arrowprops=dict(arrowstyle="-", color="0.5", lw=0.6,
                                    shrinkA=0, shrinkB=3))
    ax.set_xscale("log")
    ax.set_xlabel("τ / s")
    ax.set_ylabel("γ / Ω·cm²")
    ax.set_title(f"DRT - {label}  (λ = {drt['lam']:.1e})"
                 if label else f"DRT  (λ = {drt['lam']:.1e})")
    ax.set_xlim(tau[0], tau[-1])
    ax.set_ylim(0, max(gamma.max() * 1.45, 1e-12))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


# ---------------- Origin 风格对比图 ----------------
# Page 250.02 x 208.41 mm; layer Left 25.08%, Top 8.349%, W 68.18%, H 71.79%

MM = 1.0 / 25.4


def origin_axes():
    plt.rcParams.update({
        "font.family": "STIX Two Text",
        "font.sans-serif": ["STIX Two Text"],
        "font.size": 22,
        "axes.linewidth": 2.0,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 8, "ytick.major.size": 8,
        "xtick.minor.size": 4.5, "ytick.minor.size": 4.5,
        "xtick.major.width": 2.0, "ytick.major.width": 2.0,
        "xtick.minor.width": 1.5, "ytick.minor.width": 1.5,
    })
    fig = plt.figure(figsize=(250.02 * MM, 208.41 * MM))
    ax = fig.add_axes([0.2508, 1.0 - 0.08349 - 0.7179, 0.6818, 0.7179])
    return fig, ax


def format_ticks_legend(ax, legend_loc, bbox=None):
    """刻度标签 24pt 加粗; 图例 20pt 加粗"""
    ax.tick_params(labelsize=24)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")
    leg = ax.legend(frameon=False, fontsize=20, loc=legend_loc,
                    handlelength=1.8, bbox_to_anchor=bbox)
    for text in leg.get_texts():
        text.set_fontweight("bold")


def cathodic_sweep(df: pd.DataFrame) -> pd.DataFrame:
    """LSV 若为完整循环 (去程到顶点后返回), 取顶点后的回程段; 否则全取"""
    E = df["E(V)"].to_numpy()
    v = int(np.argmin(E))
    if E[-1] > E[v] + 0.05:
        return df.iloc[v + 1:]
    return df


# ---------------- 主流程 ----------------

def process(root, folders, labels=None, tag=None, ph=14.0, eref=0.098,
            ir_comp=0.9, rs_manual=None, tafel_window=(0.8, 1.6),
            drt_lambda=None, eis_model=None, log=print):
    """对一组样品文件夹执行 EIS DRT 分析 + LSV RHE/iR 补偿 + Tafel 对比。

    folders: 样品文件夹路径列表 (相对 root 或绝对路径均可)。
    labels / tag 为 None 时按文件夹名自动取值; rs_manual 为各样品手动
    Rs 列表 (None 则用 DRT 的 R∞)。drt_lambda 为 None 时由 L-curve 自动
    选取正则化参数。log 为日志回调 (GUI 里接 Qt 信号)。
    eis_model 仅为旧调用兼容保留: EIS 已无等效电路拟合, 传值会被忽略。
    返回 {"tag": ..., "outputs": [...]} 供 GUI 预览生成的文件。"""
    if eis_model is not None:
        log("  NOTE: eis_model 已废弃 (EIS 改用 DRT 分析), 忽略该参数")
    root = Path(root)
    if not root.is_dir():
        raise ValueError(f"not a directory: {root}")
    folders = [Path(f) for f in folders]
    for f in folders:
        if not f.is_dir():
            raise ValueError(f"sample folder not found: {f}")

    if labels is None:
        labels = [f.name for f in folders]
    if len(labels) != len(folders):
        raise ValueError("--labels 数量须与样品数一致")
    if rs_manual is not None and len(rs_manual) != len(folders):
        raise ValueError("--rs 数量须与样品数一致")

    if tag is None:
        names = [f.name for f in folders]
        prefix = os.path.commonprefix(names)
        short = [n[len(prefix):] for n in names]
        tag = "_vs_".join(short if all(short) else names)

    lo, hi = tafel_window
    outputs = []

    # ---- 1. 读取 + EIS DRT 分析 ----
    log(f"Folder: {root}   samples: {', '.join(labels)}   tag: {tag}")
    samples = []
    for folder, label in zip(folders, labels):
        s = {"folder": folder, "label": label}

        peis = find_one(folder, "eis")
        if peis:
            try:
                csv = mpr_to_csv(peis) if peis.suffix == ".mpr" else peis
                df = load_eis(csv)
                f = df["freq/Hz"].to_numpy(dtype=float)
                z_re = df["Re(Z)/Ohm"].to_numpy(dtype=float)
                z_im = -df["-Im(Z)/Ohm"].to_numpy(dtype=float)
                drt = run_drt(f, z_re, z_im, lam=drt_lambda)
                s["drt"] = drt
                s["eis_df"] = df
                stem = csv.with_suffix("")
                nyquist_out = Path(f"{stem}_Nyquist.png")
                bode_out = Path(f"{stem}_Bode.png")
                drt_out = Path(f"{stem}_DRT.png")
                plot_nyquist(df, drt, nyquist_out, label=label)
                plot_bode(df, drt, bode_out)
                plot_drt(drt, drt_out, label=label)
                outputs.extend([nyquist_out, bode_out, drt_out])
                log(f"  {label}: DRT {csv.name} "
                    f"[λ={drt['lam']:.1e}, {drt['lam_desc']}] -> "
                    f"R∞={drt['r_inf']:.3f} Ω·cm², R_pol={drt['r_pol']:.3f} Ω·cm², "
                    f"peaks={len(drt['peaks'])}, rel.err={drt['rel_err']:.1%}")
                if drt["rel_err"] > 0.10:
                    log(f"  WARNING: {label} DRT 相对残差偏大 "
                        f"({drt['rel_err']:.1%} > 10%), 检查原始谱/低频离群点, "
                        "可调整 DRT λ")
            except Exception:
                log(f"  WARNING: EIS DRT analysis failed for {label}")
                traceback.print_exc(limit=1)
        else:
            log(f"  {label}: no EIS found")

        lsv = find_one(folder, "lsv")
        if lsv:
            s["lsv_path"] = mpr_to_csv(lsv) if lsv.suffix == ".mpr" else lsv
        else:
            log(f"  {label}: no LSV found")
        samples.append(s)

    # ---- 2. EIS DRT 汇总 + Nyquist 对比 ----
    fits = [s for s in samples if "drt" in s]
    if fits:
        summary = pd.DataFrame(
            [(s["label"], f["lam"], f["r_inf"], f["r_pol"],
              f["rel_err"] * 100, len(f["peaks"]))
             for s in fits for f in [s["drt"]]],
            columns=["sample", "lambda", "R_inf/Ohm.cm2", "R_pol/Ohm.cm2",
                     "rel_err/%", "n_peaks"])
        out_csv = root / f"EIS_summary_{tag}.csv"
        summary.to_csv(out_csv, index=False, float_format="%.4g")
        outputs.append(out_csv)
        log(f"Saved: {out_csv}")

        peak_rows = [(s["label"], f["tau"][i], 1 / (2 * np.pi * f["tau"][i]),
                      f["gamma"][i])
                     for s in fits for f in [s["drt"]] for i in f["peaks"]]
        if peak_rows:
            peaks_csv = root / f"EIS_DRT_peaks_{tag}.csv"
            pd.DataFrame(peak_rows, columns=["sample", "tau/s", "f*/Hz",
                                             "gamma/Ohm.cm2"]).to_csv(
                peaks_csv, index=False, float_format="%.4g")
            outputs.append(peaks_csv)
            log(f"Saved: {peaks_csv}")

        fig, ax = origin_axes()
        for i, s in enumerate(fits):
            c = COLORS[i % len(COLORS)]
            drt = s["drt"]
            ax.plot(s["eis_df"]["Re(Z)/Ohm"], s["eis_df"]["-Im(Z)/Ohm"],
                    "o", markerfacecolor="none", markeredgecolor=c,
                    markeredgewidth=1.5, markersize=6)
            ax.plot(drt["z_fit_re"], -drt["z_fit_im"], "-", color=c, lw=2.5,
                    label=f"{s['label']} ($R_{{pol}}$={drt['r_pol']:.2f})")
        xs = np.concatenate([s["eis_df"]["Re(Z)/Ohm"].to_numpy() for s in fits])
        xrange_ = xs.max() * 1.05  # x 从 0 起, 顶端留 5%
        box_ratio = (0.7179 * 208.41) / (0.6818 * 250.02)  # 图层框高宽比
        ax.set_xlim(0, xrange_)
        ax.set_ylim(0, xrange_ * box_ratio)  # y 从 0 起, 范围与等比例坐标匹配
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_xlabel("Z' (Ω·cm²)", fontsize=36, fontweight="bold")
        ax.set_ylabel("-Z'' (Ω·cm²)", fontsize=36, fontweight="bold")
        ax.minorticks_on()
        format_ticks_legend(ax, "upper center", bbox=(0.55, 0.99))
        out_png = root / f"Nyquist_comparison_{tag}.png"
        fig.savefig(out_png, dpi=300)
        plt.close(fig)
        outputs.append(out_png)
        log(f"Saved: {out_png}")

    # ---- 3. LSV: RHE 换算 + iR 补偿 + 对比图 ----
    curves = []
    for i, s in enumerate(samples):
        if "lsv_path" not in s:
            continue
        if rs_manual is not None:
            rs = rs_manual[i]
        elif "drt" in s:
            rs = s["drt"]["r_inf"]
        else:
            rs = 0.0
            log(f"  WARNING: {s['label']} 无 Rs, 未做 iR 补偿")
        s["rs_used"] = rs

        seg = cathodic_sweep(load_scan(s["lsv_path"]))
        Ewe = seg["E(V)"].to_numpy()
        i_A = seg["i(A/cm²)"].to_numpy()      # 已是电流密度 (面积 1 cm²)
        j_mA = i_A * 1000.0
        # iR 压降 = i(A/cm²) × Rs(Ω·cm²); Rs 取 DRT 的 R∞ (高频截距),
        # EIS 阻抗已是面积归一化值, 直接相乘
        E_RHE = Ewe + 0.0592 * ph + eref - i_A * ir_comp * rs
        out = pd.DataFrame({
            "T/s": seg["T(s)"], "E/V": Ewe, "i/(A cm-2)": i_A,
            "j/(mA cm-2)": j_mA, "E_RHE_iR90/V": E_RHE,
        })
        out_csv = s["lsv_path"].with_name(f"{s['lsv_path'].stem}_RHE.csv")
        out.to_csv(out_csv, index=False)
        outputs.append(out_csv)
        log(f"Saved: {out_csv}  (Rs={rs:.3f} Ω·cm²)")
        curves.append((out, s["label"], COLORS[i % len(COLORS)]))

    if curves:
        fig, ax = origin_axes()
        for out, label, color in curves:
            ax.plot(out["E_RHE_iR90/V"], out["j/(mA cm-2)"],
                    color=color, lw=2.5, label=label)
        ax.axhline(0, color="0.5", lw=1.0)
        ax.set_xlabel("Potential (V vs. RHE)", fontsize=36, fontweight="bold")
        ax.set_ylabel("Current density (mA cm⁻²)", fontsize=36, fontweight="bold")
        ax.minorticks_on()
        format_ticks_legend(ax, "lower right")
        out_png = root / f"LSV_comparison_{tag}.png"
        fig.savefig(out_png, dpi=300)
        plt.close(fig)
        outputs.append(out_png)
        log(f"Saved: {out_png}")

        # ---- 4. Tafel 对比 ----
        fig, ax = origin_axes()
        results = []
        y_values = []
        for out, label, color in curves:
            mask = out["j/(mA cm-2)"] < 0
            eta = -out.loc[mask, "E_RHE_iR90/V"].to_numpy()
            J = -out.loc[mask, "j/(mA cm-2)"].to_numpy()
            logJ = np.log10(J)
            window = (logJ >= lo) & (logJ <= hi) & np.isfinite(logJ)
            if window.sum() < 5:
                log(f"  WARNING: {label} 在 Log|J|={lo}-{hi} 窗口内点数不足, 跳过 Tafel")
                continue
            slope, intercept = np.polyfit(logJ[window], eta[window], 1)
            tafel = slope * 1000.0
            results.append((label, tafel, color))
            ax.plot(logJ[window], eta[window], "o", markerfacecolor="none",
                    markeredgecolor=color, markeredgewidth=1.5, markersize=6,
                    label=label)
            y_values.append(eta[window])
            x_line = np.linspace(lo, hi, 200)
            y_line = slope * x_line + intercept
            ax.plot(x_line, y_line, "-", color=color, lw=2)
            y_values.append(y_line)

        for i, (label, tafel, color) in enumerate(results):
            # 第一条曲线标注在右下, 其余从右上向下排
            y = 0.03 if i == 0 else 0.90 - (i - 1) * 0.09
            ax.text(0.95, y, f"{tafel:.0f} mV dec⁻¹", transform=ax.transAxes,
                    fontsize=24, ha="right", va="bottom", color=color)

        ax.set_xlabel("Log (|J| (mA cm⁻²))", fontsize=36, fontweight="bold")
        ax.set_ylabel("Overpotential (V vs. RHE)", fontsize=36, fontweight="bold")
        ax.set_xlim(lo - 0.05, hi + 0.05)
        if y_values:
            y_min = min(np.min(y) for y in y_values)
            y_max = max(np.max(y) for y in y_values)
            margin = max(0.02, 0.05 * (y_max - y_min))
            ax.set_ylim(y_min - margin, y_max + margin)
        ax.minorticks_on()
        format_ticks_legend(ax, "upper left")
        out_png = root / f"Tafel_comparison_{tag}.png"
        fig.savefig(out_png, dpi=300)
        plt.close(fig)
        outputs.append(out_png)
        log(f"Saved: {out_png}")
        for label, tafel, _ in results:
            log(f"  {label}: Tafel slope = {tafel:.1f} mV/dec (Log|J| = {lo}-{hi})")

    log("Done.")
    return {"tag": tag, "outputs": outputs}


def main():
    ap = argparse.ArgumentParser(
        description="通用电化学数据处理: CS Studio txt + EIS/LSV/Tafel 对比图")
    ap.add_argument("root", help="日期文件夹, 如 DATA/EC/20260921")
    ap.add_argument("samples", nargs="+",
                    help="样品子文件夹名 (相对 root), 也可直接给文件夹路径")
    ap.add_argument("--labels", nargs="+", default=None, help="图例标签, 与样品一一对应")
    ap.add_argument("--tag", default=None,
                    help="输出文件名后缀, 默认取样品名去公共前缀后用 _vs_ 连接")
    ap.add_argument("--ph", type=float, default=14.0, help="电解液 pH (默认 14, 1 M KOH)")
    ap.add_argument("--eref", type=float, default=0.098,
                    help="参比电极电位 vs. SHE / V (默认 0.098, Hg/HgO)")
    ap.add_argument("--ir-comp", type=float, default=0.9, help="iR 补偿比例 (默认 0.9)")
    ap.add_argument("--rs", nargs="+", type=float, default=None,
                    help="手动指定各样品 Rs / Ω·cm², 跳过 DRT 的 R∞ 取值")
    ap.add_argument("--tafel-window", nargs=2, type=float, default=[0.8, 1.6],
                    metavar=("LOGJ_MIN", "LOGJ_MAX"), help="Tafel 拟合窗口 (默认 0.8 1.6)")
    ap.add_argument("--drt-lambda", type=float, default=1e-2,
                    help="DRT 正则化参数 λ (默认 1e-2; 用 --drt-lcurve 改为自动)")
    ap.add_argument("--drt-lcurve", action="store_true",
                    help="用 L-curve 自动选取 DRT λ (覆盖 --drt-lambda)")
    ap.add_argument("--eis-model", default=None, help=argparse.SUPPRESS)
    args = ap.parse_args()

    folders = []
    for s in args.samples:
        p = Path(s)
        folders.append(p if p.is_dir() else Path(args.root) / s)

    try:
        process(root=args.root, folders=folders, labels=args.labels,
                tag=args.tag, ph=args.ph, eref=args.eref, ir_comp=args.ir_comp,
                rs_manual=args.rs, tafel_window=args.tafel_window,
                drt_lambda=None if args.drt_lcurve else args.drt_lambda,
                eis_model=args.eis_model)
    except ValueError as e:
        sys.exit(f"ERROR: {e}")


if __name__ == "__main__":
    main()
