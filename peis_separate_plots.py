# -*- coding: utf-8 -*-
"""EIS 单文件分析: Nyquist / Bode / DRT 三图单独输出 (实测点 + 拟合曲线)

与 RuNiMoHER/scripts/peis_separate_plots.py 对齐, 输入支持 CS Studio
*-EIS.txt 与 galvani 转出 csv (经 echem_process.load_eis 归一化)。

- 等效电路拟合: 先试 Rs+(Rct‖CPE), 再试 Rs+L+(R1‖CPE1)+(R2‖CPE2)
  (双时间常数, 适用于带感抗尾/双弧数据), 双弧模型 GoF 明显更优
  (相对改善 >20%) 时采用之
- DRT 复用 drt_analysis.py 的 Tikhonov + NNLS 实现; λ 默认 1e-2
  (该数据 L-curve 自动值过小会导致过拟合尖峰), 可用 --lambda 或 --lcurve 调整

用法 (在 EchemPro-src/ 下):
  ../RuHeterojunctionHER/.venv/bin/python peis_separate_plots.py \
      ../RuHeterojunctionHER/DATA/EC/20260921/Ru_NiP_annealing/init-EIS.txt
  ../RuHeterojunctionHER/.venv/bin/python peis_separate_plots.py <eis 文件> --lambda 1e-3
  ../RuHeterojunctionHER/.venv/bin/python peis_separate_plots.py <eis 文件> --lcurve

输出 (与输入同目录):
  <stem>_Nyquist.png / <stem>_Bode.png / <stem>_DRT.png / <stem>_DRT.csv
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 脚本直接运行时, 同目录 (scripts/) 已在 sys.path 中, 可直接复用现有工具
sys.path.insert(0, str(Path(__file__).resolve().parent))
from echem_process import fit_eis, circuit_z, load_eis
from drt_analysis import build_matrices, solve_drt, lcurve_lambda, find_peaks

plt.rcParams['font.size'] = 14
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ---------------- 等效电路拟合: Rs + L + (R1‖CPE1) + (R2‖CPE2) ----------------

def circuit_z2(freqs, rs, L, r1, q1, n1, r2, q2, n2):
    """Rs + L + (R1‖CPE1) + (R2‖CPE2) 阻抗; CPE: Z = 1/(Q·(jω)^n)"""
    w = 2 * np.pi * np.asarray(freqs, dtype=float)
    return (rs + 1j * w * L
            + r1 / (1.0 + r1 * q1 * (1j * w) ** n1)
            + r2 / (1.0 + r2 * q2 * (1j * w) ** n2))


def fit_eis2(df: pd.DataFrame):
    """比例加权 (1/|Z|) 最小二乘拟合, 返回 dict(参数, gof, f_fit)"""
    from scipy.optimize import least_squares
    freqs = df['freq/Hz'].to_numpy(dtype=float)
    z_meas = (df['Re(Z)/Ohm'].to_numpy(dtype=float)
              - 1j * df['-Im(Z)/Ohm'].to_numpy(dtype=float))
    span = max(z_meas.real.max() - z_meas.real.min(), 1e-3)

    def residuals(p):
        res = circuit_z2(freqs, *p) - z_meas
        return np.concatenate([res.real, res.imag]) / np.abs(z_meas).repeat(2)

    sol = least_squares(residuals,
                        [z_meas.real[0], 1e-7, 0.3 * span, 0.05, 0.9, 0.7 * span, 0.3, 0.8],
                        bounds=([0, 0, 1e-9, 1e-10, 0.5, 1e-9, 1e-10, 0.5],
                                [np.inf, 1e-3, np.inf, 1e3, 1.0, np.inf, 1e3, 1.0]),
                        max_nfev=20000)
    rel = np.abs(circuit_z2(freqs, *sol.x) - z_meas) / np.abs(z_meas)
    keys = ['rs', 'L', 'r1', 'q1', 'n1', 'r2', 'q2', 'n2']
    return dict(zip(keys, sol.x), gof=100 * np.sqrt(np.mean(rel ** 2)),
                f_fit=np.logspace(np.log10(freqs.min()), np.log10(freqs.max()), 200))


def fit_z(fit, freqs):
    """按所选模型计算拟合阻抗"""
    if fit['model'] == 2:
        return circuit_z2(freqs, fit['rs'], fit['L'], fit['r1'], fit['q1'], fit['n1'],
                          fit['r2'], fit['q2'], fit['n2'])
    return circuit_z(freqs, fit['rs'], fit['rct'], fit['q'], fit['n'])


def fit_label(fit):
    return ('fit: Rs+L+(R1‖CPE1)+(R2‖CPE2)' if fit['model'] == 2
            else 'fit: Rs+(Rct‖CPE)')


def fit_text(fit):
    unit = 'Ω·cm²'
    if fit['model'] == 2:
        return (f"Rs = {fit['rs']:.2f} {unit},  L = {fit['L']:.1e} H\n"
                f"R1 = {fit['r1']:.2f} {unit} (n1 = {fit['n1']:.2f})\n"
                f"R2 = {fit['r2']:.2f} {unit} (n2 = {fit['n2']:.2f})\n"
                f"GoF = {fit['gof']:.1f}%")
    return (f"Rs = {fit['rs']:.2f} {unit}\nRct = {fit['rct']:.2f} {unit}\n"
            f"n = {fit['n']:.2f},  GoF = {fit['gof']:.1f}%")


# ---------------- 绘图 ----------------

def plot_nyquist(f, z_re, z_im, fit, out, stem):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(z_re, -z_im, 'o', color='0.45', markersize=6, label='measured')
    z_fit = fit_z(fit, fit['f_fit'])
    ax.plot(z_fit.real, -z_fit.imag, 'b-', linewidth=2, label=fit_label(fit))
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("Z' (Ω·cm²)")
    ax.set_ylabel("-Z'' (Ω·cm²)")
    ax.set_title(f'Nyquist - {stem}')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    ax.text(0.03, 0.03, fit_text(fit), transform=ax.transAxes, va='bottom', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    fig.tight_layout()
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


def plot_bode(f, z_re, z_im, fit, out, stem):
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax2 = ax1.twinx()
    ax1.semilogx(f, np.abs(z_re + 1j * z_im), 'o', color='0.45', markersize=6,
                 label='|Z| measured')
    ax2.semilogx(f, np.degrees(np.arctan2(z_im, z_re)), 's', color='lightsteelblue',
                 markersize=6, label='Phase measured')
    z_fit = fit_z(fit, fit['f_fit'])
    ax1.semilogx(fit['f_fit'], np.abs(z_fit), 'b-', linewidth=2, label='|Z| fit')
    ax2.semilogx(fit['f_fit'], np.degrees(np.angle(z_fit)), '-', color='navy',
                 linewidth=2, label='Phase fit')
    ax1.set_xlabel('Frequency / Hz')
    ax1.set_ylabel('|Z| / Ω·cm²')
    ax2.set_ylabel('Phase / deg')
    ax1.set_title(f'Bode - {stem}')
    ax1.grid(True, which='both', alpha=0.3)
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=10, loc='best')
    fig.tight_layout()
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


# ---------------- DRT ----------------

def tau_grid(f):
    """τ 网格: 覆盖数据频率对应的时间常数范围, 两端各延伸 1 个数量级 (同 drt_analysis.py)"""
    tau_min = 1 / (2 * np.pi * f.max()) * 0.1
    tau_max = 1 / (2 * np.pi * f.min()) * 10
    n_dec = np.log10(tau_max / tau_min)
    return np.logspace(np.log10(tau_min), np.log10(tau_max), int(n_dec * 10))


def run_drt(f, z_re, z_im, lam, lam_desc, stem, out_dir):
    tau = tau_grid(f)
    tau_min, tau_max = tau[0], tau[-1]

    A, D = build_matrices(f, tau)
    b = np.concatenate([z_re, z_im])
    x = solve_drt(A, D, b, lam)

    r_inf, L, gamma = x[0], x[1], x[2:]
    rel_err = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
    r_pol = gamma.sum() * np.diff(np.log(tau)).mean()
    peaks = find_peaks(tau, gamma)

    print(f'\n===== DRT: λ = {lam:.2e} ({lam_desc}) =====')
    print(f'R∞ = {r_inf:.4f} Ω·cm²,  L = {L:.3e} H,  R_pol(Σγ) ≈ {r_pol:.4f} Ω·cm²')
    print(f'拟合相对残差 = {rel_err:.2%}')
    for i in peaks:
        print(f'DRT 峰: τ = {tau[i]:.3e} s  (f* = {1/(2*np.pi*tau[i]):.3g} Hz),  γ = {gamma[i]:.4f} Ω·cm²')

    drt_csv = out_dir / f'{stem}_DRT.csv'
    pd.DataFrame({'tau/s': tau, 'gamma/Ohm.cm2': gamma}).to_csv(drt_csv, index=False)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(tau, gamma, color='tab:blue', linewidth=1.8)
    ax.fill_between(tau, 0, gamma, color='tab:blue', alpha=0.25)
    for k, i in enumerate(peaks):  # 标注高低交错, 避免相邻峰文字重叠
        ax.annotate(f'τ={tau[i]:.2e} s\n(f*={1/(2*np.pi*tau[i]):.2g} Hz)',
                    (tau[i], gamma[i]), textcoords='offset points',
                    xytext=(0, 10 + 22 * (k % 2)), ha='center', fontsize=9)
    ax.set_xscale('log')
    ax.set_xlabel('τ / s')
    ax.set_ylabel('γ / Ω·cm²')
    ax.set_title(f'DRT - {stem}  (λ = {lam:.1e})')
    ax.set_xlim(tau_min, tau_max)
    ax.set_ylim(0, gamma.max() * 1.55)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    png = out_dir / f'{stem}_DRT.png'
    fig.savefig(png, bbox_inches='tight')
    plt.close(fig)
    print(f'已保存: {png}')
    print(f'已保存: {drt_csv}')


def main():
    ap = argparse.ArgumentParser(description='EIS 分析: Nyquist / Bode / DRT 三图单独输出')
    ap.add_argument('csv', help='EIS 文件路径 (CS Studio *-EIS.txt 或 galvani 转出 csv)')
    ap.add_argument('--lambda', dest='lam', type=float, default=1e-2,
                    help='DRT 正则化参数 (默认 1e-2)')
    ap.add_argument('--lcurve', action='store_true',
                    help='改用 L-curve 自动选取 λ (数据噪声小时适用)')
    args = ap.parse_args()

    csv_path = Path(args.csv)
    out_dir = csv_path.parent
    stem = csv_path.stem

    df = load_eis(csv_path)
    f = df['freq/Hz'].to_numpy(dtype=float)
    z_re = df['Re(Z)/Ohm'].to_numpy(dtype=float)
    z_im = -df['-Im(Z)/Ohm'].to_numpy(dtype=float)
    order = np.argsort(f)
    f, z_re, z_im = f[order], z_re[order], z_im[order]

    # 单弧 vs 双弧模型, 双弧 GoF 相对改善 >20% 才采用
    fit1 = fit_eis(df)
    fit1['model'] = 1
    try:
        fit2 = fit_eis2(df)
        fit2['model'] = 2
    except Exception:
        fit2 = None
        print('WARNING: 双弧模型拟合失败, 使用 Rs+(Rct‖CPE)')
    print('===== 等效电路拟合 =====')
    print(f"Rs+(Rct‖CPE):            Rct = {fit1['rct']:.4f} Ω·cm²,  GoF = {fit1['gof']:.2f}%")
    if fit2 is not None:
        print(f"Rs+L+(R1‖CPE1)+(R2‖CPE2): R1 = {fit2['r1']:.4f} Ω·cm²,  R2 = {fit2['r2']:.4f} Ω·cm²,  "
              f"GoF = {fit2['gof']:.2f}%")
    fit = fit2 if (fit2 is not None and fit2['gof'] < 0.8 * fit1['gof']) else fit1
    print(f'采用模型: {fit_label(fit)[5:]}')
    if fit['model'] == 2:
        print(f"Rs = {fit['rs']:.4f} Ω·cm²,  L = {fit['L']:.3e} H,  "
              f"R1 = {fit['r1']:.4f} Ω·cm² (Q1 = {fit['q1']:.3e}, n1 = {fit['n1']:.3f}),  "
              f"R2 = {fit['r2']:.4f} Ω·cm² (Q2 = {fit['q2']:.3e}, n2 = {fit['n2']:.3f})")
    else:
        print(f"Rs = {fit['rs']:.4f} Ω·cm²,  Rct = {fit['rct']:.4f} Ω·cm²,  "
              f"Q = {fit['q']:.3e} S·s^n,  n = {fit['n']:.3f}")

    plot_nyquist(f, z_re, z_im, fit, out_dir / f'{stem}_Nyquist.png', stem)
    print(f'已保存: {out_dir / f"{stem}_Nyquist.png"}')
    plot_bode(f, z_re, z_im, fit, out_dir / f'{stem}_Bode.png', stem)
    print(f'已保存: {out_dir / f"{stem}_Bode.png"}')

    if args.lcurve:
        tau = tau_grid(f)
        A, D = build_matrices(f, tau)
        lam, _, _, _ = lcurve_lambda(A, D, np.concatenate([z_re, z_im]))
        lam_desc = 'L-curve 自动'
    else:
        lam, lam_desc = args.lam, '手动/默认'
    run_drt(f, z_re, z_im, lam, lam_desc, stem, out_dir)


if __name__ == '__main__':
    main()
