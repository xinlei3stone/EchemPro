# -*- coding: utf-8 -*-
"""DRT (Distribution of Relaxation Times) 分析脚本

与 RuNiMoHER/scripts/drt_analysis.py 对齐, 输入层同时支持:
  - CS Studio 导出的 *-EIS.txt (Freq(Hz) / Z'(Ohm.cm²) / Z''(Ohm.cm²), 见
    echem_process.read_cs_eis_txt)
  - galvani 从 .mpr 转出的 PEIS csv (需含 freq/Hz, Re(Z)/Ohm, -Im(Z)/Ohm 三列)

方法参照 Ciucci 组 DRTtools (https://github.com/Mycroft2333/DRTtools.com 所述方法):
  - 阻抗模型: Z(ω) = R∞ + iωL + ∫ γ(τ) / (1 + iωτ) dlnτ
  - Tikhonov 正则化 (一阶差分罚) + 非负最小二乘 (NNLS, γ ≥ 0)
  - 正则化参数 λ 通过 L-curve 拐点自动选取, 也可用 --lambda 手动指定

用法 (在 EchemPro-src/ 下):
  ../RuHeterojunctionHER/.venv/bin/python drt_analysis.py <eis.txt|peis_csv> \
      [--lambda 1e-3] [--out-dir DIR]

输出: <文件名>_DRT.png (Nyquist / Bode / DRT 三联图) 和 <文件名>_DRT.csv (τ, γ)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import nnls

import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def load_eis(path):
    """读取 EIS 数据, 返回 f, Z', Z'' (Z'' 为真实虚部, 电容性为负)"""
    path = Path(path)
    if path.suffix.lower() == ".txt":
        # 同目录导入, 避免包名冲突
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from echem_process import read_cs_eis_txt
        df = read_cs_eis_txt(path)
    else:
        df = pd.read_csv(path)
    f = df["freq/Hz"].to_numpy(dtype=float)
    z_re = df["Re(Z)/Ohm"].to_numpy(dtype=float)
    z_im = -df["-Im(Z)/Ohm"].to_numpy(dtype=float)
    order = np.argsort(f)  # 按频率升序, 方便绘图
    return f[order], z_re[order], z_im[order]


def build_matrices(f, tau):
    """构造 DRT 线性系统的实部/虚部矩阵 (含 R∞ 和 L 两列)"""
    omega = 2 * np.pi * f
    dln_tau = np.diff(np.log(tau)).mean()
    wt = np.outer(omega, tau)                     # (Nf, Nt)
    k_re = 1.0 / (1.0 + wt**2) * dln_tau          # Re 核: 1/(1+ω²τ²)
    k_im = -wt / (1.0 + wt**2) * dln_tau          # Im 核: -ωτ/(1+ω²τ²)
    n_tau = len(tau)
    a_re = np.hstack([np.ones((len(f), 1)), np.zeros((len(f), 1)), k_re])
    a_im = np.hstack([np.zeros((len(f), 1)), omega[:, None], k_im])
    A = np.vstack([a_re, a_im])
    # 一阶差分罚矩阵, 只作用于 γ 部分 (R∞, L 不正则化)
    D = np.zeros((n_tau - 1, n_tau + 2))
    for i in range(n_tau - 1):
        D[i, 2 + i] = -1.0
        D[i, 2 + i + 1] = 1.0
    return A, D


def solve_drt(A, D, b, lam):
    """Tikhonov + NNLS: min ||Ax-b||² + λ²||Dγ||², s.t. x ≥ 0"""
    A_aug = np.vstack([A, lam * D])
    b_aug = np.concatenate([b, np.zeros(D.shape[0])])
    x, _ = nnls(A_aug, b_aug)
    return x


def lcurve_lambda(A, D, b, lambdas=None):
    """L-curve 拐点自动选 λ (log 残差范数 vs log 解半范数, 最大曲率)"""
    if lambdas is None:
        lambdas = np.logspace(-7, 1, 17)
    res, norm = [], []
    for lam in lambdas:
        x = solve_drt(A, D, b, lam)
        res.append(np.linalg.norm(A @ x - b))
        norm.append(np.linalg.norm(D @ x))
    lr, ln = np.log10(res), np.log10(norm)
    # 参数曲线 (lr(λ), ln(λ)) 的曲率
    d1r, d1n = np.gradient(lr), np.gradient(ln)
    d2r, d2n = np.gradient(d1r), np.gradient(d1n)
    curv = np.abs(d1r * d2n - d1n * d2r) / (d1r**2 + d1n**2) ** 1.5
    i = int(np.argmax(curv[1:-1])) + 1  # 端点曲率不可靠, 只在内部取
    return lambdas[i], lambdas, lr, ln


def find_peaks(tau, gamma):
    """在 γ(τ) 上找显著局部极大值 (>5% 峰高)"""
    peaks = []
    for i in range(1, len(gamma) - 1):
        if gamma[i] > gamma[i - 1] and gamma[i] >= gamma[i + 1] and gamma[i] > 0.05 * gamma.max():
            peaks.append(i)
    return peaks


def main():
    ap = argparse.ArgumentParser(description="DRT 分析 (Tikhonov + NNLS)")
    ap.add_argument("csv", help="EIS 文件路径 (CS Studio *-EIS.txt 或 galvani 转出 csv)")
    ap.add_argument("--lambda", dest="lam", type=float, default=None,
                    help="正则化参数 (默认 L-curve 自动选取)")
    ap.add_argument("--out-dir", default=None, help="输出目录 (默认与输入同目录)")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    out_dir = Path(args.out_dir) if args.out_dir else csv_path.parent
    stem = csv_path.stem

    f, z_re, z_im = load_eis(csv_path)
    # τ 网格: 覆盖数据频率对应的时间常数范围, 两端各延伸 1 个数量级
    tau_min = 1 / (2 * np.pi * f.max()) * 0.1
    tau_max = 1 / (2 * np.pi * f.min()) * 10
    n_dec = np.log10(tau_max / tau_min)
    tau = np.logspace(np.log10(tau_min), np.log10(tau_max), int(n_dec * 10))

    A, D = build_matrices(f, tau)
    b = np.concatenate([z_re, z_im])

    if args.lam is not None:
        lam = args.lam
    else:
        lam, _, _, _ = lcurve_lambda(A, D, b)
    x = solve_drt(A, D, b, lam)

    r_inf, L, gamma = x[0], x[1], x[2:]
    z_fit = A @ x
    n = len(f)
    z_fit_re, z_fit_im = z_fit[:n], z_fit[n:]
    rel_err = np.linalg.norm(z_fit - b) / np.linalg.norm(b)

    # 汇总输出
    r_pol = gamma.sum() * np.diff(np.log(tau)).mean()
    peaks = find_peaks(tau, gamma)
    print(f"λ = {lam:.2e} (L-curve 自动)" if args.lam is None else f"λ = {lam:.2e} (手动)")
    print(f"R∞ = {r_inf:.4f} Ω·cm²,  L = {L:.3e} H,  R_pol(Σγ) ≈ {r_pol:.4f} Ω·cm²")
    print(f"拟合相对残差 = {rel_err:.2%}")
    for i in peaks:
        print(f"DRT 峰: τ = {tau[i]:.3e} s  (f* = {1/(2*np.pi*tau[i]):.3g} Hz),  γ = {gamma[i]:.4f} Ω·cm²")

    # 导出 γ(τ)
    drt_csv = out_dir / f"{stem}_DRT.csv"
    pd.DataFrame({"tau/s": tau, "gamma/Ohm.cm2": gamma}).to_csv(drt_csv, index=False)

    # ---- 三联图: Nyquist / Bode / DRT ----
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    fig.suptitle(f"{stem}  |  λ = {lam:.1e}", fontsize=11)

    ax = axes[0]
    ax.scatter(z_re, -z_im, c="k", s=22, label="Original EIS Data", zorder=3)
    ax.plot(z_fit_re, -z_fit_im, c="tab:red", lw=1.5, label="DRT Fit")
    ax.set_xlabel("Z' (Ω·cm²)")
    ax.set_ylabel("-Z'' (Ω·cm²)")
    ax.set_title("Nyquist Plot")
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend()
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax2 = ax.twinx()
    ax.scatter(f, np.abs(z_re + 1j * z_im), c="k", s=18, label="|Z|")
    ax.plot(f, np.abs(z_fit_re + 1j * z_fit_im), c="k", lw=1.2)
    # 与 DRTtools 网页一致: 电容性相位取正 (即 -Phase 约定)
    ax2.scatter(f, -np.degrees(np.arctan2(z_im, z_re)), c="mediumslateblue", s=18, label="Phase")
    ax2.plot(f, -np.degrees(np.arctan2(z_fit_im, z_fit_re)), c="mediumslateblue", lw=1.2)
    ax.set_xscale("log")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("|Z| (Ω·cm²)")
    ax2.set_ylabel("Phase (°)")
    ax.set_title("Bode Plot")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right")
    ax.grid(alpha=0.3)

    ax = axes[3 - 1]
    ax.plot(tau, gamma, c="tab:blue", lw=1.5)
    ax.fill_between(tau, 0, gamma, color="tab:blue", alpha=0.25)
    for i in peaks:
        ax.annotate(f"τ={tau[i]:.2e} s\n(f*={1/(2*np.pi*tau[i]):.2g} Hz)",
                    (tau[i], gamma[i]), textcoords="offset points", xytext=(0, 10),
                    ha="center", fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("τ (s)")
    ax.set_ylabel("γ (Ω·cm²)")
    ax.set_title("Distribution of Relaxation Times (DRT)")
    ax.set_xlim(tau_min, tau_max)
    ax.set_ylim(0, gamma.max() * 1.45)  # 顶部留白, 避免峰标注被裁
    ax.grid(alpha=0.3)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    png_path = out_dir / f"{stem}_DRT.png"
    fig.savefig(png_path, dpi=200)
    print(f"已保存: {png_path}")
    print(f"已保存: {drt_csv}")


if __name__ == "__main__":
    main()
