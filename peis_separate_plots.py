# -*- coding: utf-8 -*-
"""EIS 单文件分析: Nyquist / Bode / DRT 三图单独输出 (实测点 + 拟合曲线)

与 RuNiMoHER/scripts/peis_separate_plots.py 对齐, 输入支持 CS Studio
*-EIS.txt 与 galvani 转出 csv (经 echem_process.load_eis 归一化)。

- DRT 复用 drt_analysis.py 的 run_drt (Tikhonov + NNLS)
- λ 默认 1e-2 (该数据 L-curve 自动值过小会导致过拟合尖峰),
  可用 --lambda 手动指定或 --lcurve 改用 L-curve 自动选取

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
from echem_process import load_eis
from drt_analysis import run_drt

plt.rcParams['font.size'] = 14
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def drt_text(drt):
    return (f"R∞ = {drt['r_inf']:.2f} Ω·cm²\n"
            f"R_pol = {drt['r_pol']:.2f} Ω·cm²\n"
            f"λ = {drt['lam']:.1e},  rel. err = {drt['rel_err']:.1%}")


# ---------------- 绘图 ----------------

def plot_nyquist(f, z_re, z_im, drt, out, stem):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(z_re, -z_im, 'o', color='0.45', markersize=6, label='measured')
    ax.plot(drt['z_fit_re'], -drt['z_fit_im'], 'b-', linewidth=2, label='DRT fit')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel("Z' (Ω·cm²)")
    ax.set_ylabel("-Z'' (Ω·cm²)")
    ax.set_title(f'Nyquist - {stem}')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    ax.text(0.03, 0.03, drt_text(drt), transform=ax.transAxes, va='bottom',
            fontsize=11, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    fig.tight_layout()
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


def plot_bode(f, z_re, z_im, drt, out, stem):
    z_fit = drt['z_fit_re'] + 1j * drt['z_fit_im']
    fig, ax1 = plt.subplots(figsize=(8, 6))
    ax2 = ax1.twinx()
    ax1.semilogx(f, np.abs(z_re + 1j * z_im), 'o', color='0.45', markersize=6,
                 label='|Z| measured')
    ax2.semilogx(f, np.degrees(np.arctan2(z_im, z_re)), 's', color='lightsteelblue',
                 markersize=6, label='Phase measured')
    ax1.semilogx(f, np.abs(z_fit), 'b-', linewidth=2, label='|Z| DRT fit')
    ax2.semilogx(f, np.degrees(np.angle(z_fit)), '-', color='navy',
                 linewidth=2, label='Phase DRT fit')
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


def plot_drt(drt, out, stem):
    tau, gamma = drt['tau'], drt['gamma']
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(tau, gamma, color='tab:blue', linewidth=1.8)
    ax.fill_between(tau, 0, gamma, color='tab:blue', alpha=0.25)
    for k, i in enumerate(drt['peaks']):  # 左右交错 + 纵向错位, 避免相邻峰文字重叠
        left = (k % 2 == 0)
        ax.annotate(f'τ={tau[i]:.2e} s\n(f*={1/(2*np.pi*tau[i]):.2g} Hz)',
                    (tau[i], gamma[i]), textcoords='offset points',
                    xytext=(-14 if left else 14, 14 + 12 * (k % 3)),
                    ha='right' if left else 'left', va='bottom', fontsize=9,
                    arrowprops=dict(arrowstyle='-', color='0.5', lw=0.6,
                                    shrinkA=0, shrinkB=3))
    ax.set_xscale('log')
    ax.set_xlabel('τ / s')
    ax.set_ylabel('γ / Ω·cm²')
    ax.set_title(f"DRT - {stem}  (λ = {drt['lam']:.1e})")
    ax.set_xlim(tau[0], tau[-1])
    ax.set_ylim(0, max(gamma.max() * 1.55, 1e-12))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)


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

    drt = run_drt(f, z_re, z_im, lam=None if args.lcurve else args.lam)
    print('===== DRT =====')
    print(f"λ = {drt['lam']:.2e} ({drt['lam_desc']})")
    print(f"R∞ = {drt['r_inf']:.4f} Ω·cm²,  L = {drt['L']:.3e} H,  "
          f"R_pol(Σγ) ≈ {drt['r_pol']:.4f} Ω·cm²")
    print(f"拟合相对残差 = {drt['rel_err']:.2%}")
    for i in drt['peaks']:
        print(f"DRT 峰: τ = {drt['tau'][i]:.3e} s  "
              f"(f* = {1/(2*np.pi*drt['tau'][i]):.3g} Hz),  "
              f"γ = {drt['gamma'][i]:.4f} Ω·cm²")

    drt_csv = out_dir / f'{stem}_DRT.csv'
    pd.DataFrame({'tau/s': drt['tau'], 'gamma/Ohm.cm2': drt['gamma']}).to_csv(
        drt_csv, index=False)

    nyquist_png = out_dir / f'{stem}_Nyquist.png'
    bode_png = out_dir / f'{stem}_Bode.png'
    drt_png = out_dir / f'{stem}_DRT.png'

    plot_nyquist(f, z_re, z_im, drt, nyquist_png, stem)
    plot_bode(f, z_re, z_im, drt, bode_png, stem)
    plot_drt(drt, drt_png, stem)
    for path in (nyquist_png, bode_png, drt_png, drt_csv):
        print(f'已保存: {path}')


if __name__ == '__main__':
    main()
