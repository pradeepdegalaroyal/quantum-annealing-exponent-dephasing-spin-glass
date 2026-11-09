"""Greyscale figures from results/*.csv (dumbbell, small multiples, waterfall, exponent-vs-dephasing)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.color": "0.85", "grid.linewidth": 0.5})
DPI = 200


def dumbbell(res: Path, fig: Path) -> None:
    df = pd.read_csv(res / "coherent_exponents.csv").sort_values(["n", "seed"]).reset_index(drop=True)
    f, ax = plt.subplots(figsize=(6.2, 0.28 * len(df) + 1.4))
    for i, r in df.iterrows():
        ax.plot([r.kappa_sa, r.kappa_qa], [i, i], color="0.6", lw=1.2, zorder=1)
    ax.scatter(df.kappa_sa, range(len(df)), marker="s", s=28, facecolor="white", edgecolor="black", zorder=2,
               label="simulated annealing (SA)")
    ax.scatter(df.kappa_qa, range(len(df)), marker="o", s=28, color="black", zorder=3,
               label="coherent quantum annealing (QA)")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels([f"n={int(r.n)}, seed {int(r.seed)}" for _, r in df.iterrows()], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel(r"fitted residual-energy exponent $\kappa$ (window $4\leq\tau\leq32$)")
    ax.set_title("Residual-energy exponent per instance ensemble:\nsimulated annealing versus coherent quantum annealing")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=8, frameon=False)
    f.tight_layout()
    f.savefig(fig / "fig1_exponent_dumbbell.png", dpi=DPI)
    plt.close(f)


def small_multiples(res: Path, fig: Path) -> None:
    df = pd.read_csv(res / "dephasing_curves.csv")
    n = int(df.n.max())
    d = df[df.n == n]
    gammas = sorted(d.gamma.unique())
    cols = 3
    rows = int(np.ceil(len(gammas) / cols))
    f, axes = plt.subplots(rows, cols, figsize=(8.4, 2.6 * rows), sharex=True, sharey=True, squeeze=False)
    for ax, g in zip(axes.ravel(), gammas):
        dg = d[d.gamma == g].groupby("tau")
        m, s = dg.qa.mean(), dg.qa.std(ddof=1)
        ms, ss = dg.sa.mean(), dg.sa.std(ddof=1)
        ax.errorbar(m.index, m, yerr=s, color="black", marker="o", ms=3.5, lw=1.2, capsize=2, label="QA")
        ax.errorbar(ms.index, ms, yerr=ss, color="0.4", marker="s", ms=3.5, lw=1.2, ls="--", capsize=2,
                    mfc="white", label="SA")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks([2, 4, 8, 16, 32])
        ax.set_xticklabels(["2", "4", "8", "16", "32"])
        ax.minorticks_off()
        ax.set_title(rf"$\gamma$ = {g:g}", fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel(r"anneal time $\tau$ (1/J; SA: sweeps)")
    for ax in axes[:, 0]:
        ax.set_ylabel("residual energy per spin (J)")
    axes[0, 0].legend(fontsize=8)
    f.suptitle(f"Residual energy versus anneal time at n = {n}, mean ± std over seeds, by dephasing rate", y=1.0)
    f.tight_layout()
    f.savefig(fig / "fig2_dephasing_small_multiples.png", dpi=DPI)
    plt.close(f)


def waterfall(res: Path, fig: Path) -> None:
    df = pd.read_csv(res / "ablation_noise_sources.csv")
    summ = json.load(open(res / "summary.json", encoding="utf8"))
    n = summ["ablation"]["n"]
    k = df.groupby("config").kappa.agg(["mean", "std"])
    ideal, ctrl, deph, both = (k.loc[c, "mean"] for c in ("ideal", "control_error", "dephasing", "both"))
    steps = [("ideal\ncoherent QA", ideal, None), ("+ control\nerror", ctrl - ideal, None),
             ("+ dephasing\n(on top)", both - ctrl, None), ("both\nnoise sources", both, None)]
    f, ax = plt.subplots(figsize=(6.0, 3.8))
    level = 0.0
    xs = np.arange(len(steps))
    for i, (lab, v, _) in enumerate(steps):
        if i in (0, 3):
            ax.bar(i, v, color="0.25", edgecolor="black", width=0.55)
            level = v
        else:
            bottom = level if v >= 0 else level + v
            ax.bar(i, abs(v), bottom=bottom, color="white", edgecolor="black", hatch="//", width=0.55)
            level += v
        ax.text(i, (level if i in (0, 3) else max(level, level - v)) + 0.04, f"{v:+.2f}" if i in (1, 2) else f"{v:.2f}",
                ha="center", fontsize=8)
    ksa = pd.read_csv(res / "coherent_exponents.csv")
    ksa = ksa[ksa.n == n].kappa_sa.mean() if (ksa.n == n).any() else ksa.kappa_sa.mean()
    ax.axhline(ksa, color="black", ls=":", lw=1.2)
    ax.text(-0.4, ksa + 0.04, "SA exponent", ha="left", fontsize=8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([s[0] for s in steps], fontsize=8)
    ax.set_ylabel(r"exponent $\kappa$")
    ax.set_title(f"Noise-source ablation of the QA exponent (n = {n})")
    f.tight_layout()
    f.savefig(fig / "fig3_noise_ablation_waterfall.png", dpi=DPI)
    plt.close(f)


def exponent_vs_gamma(res: Path, fig: Path) -> None:
    df = pd.read_csv(res / "dephasing_exponents.csv")
    ns = sorted(df.n.unique())
    f, axes = plt.subplots(1, len(ns), figsize=(2.9 * len(ns) + 0.6, 3.2), sharey=True, squeeze=False)
    for ax, n in zip(axes[0], ns):
        d = df[df.n == n]
        g = d.groupby("gamma").kappa_qa
        m, s = g.mean(), g.std(ddof=1)
        ax.errorbar(m.index, m, yerr=s, color="black", marker="o", ms=4, capsize=2, lw=1.2, label="QA")
        ks = d.groupby("seed").kappa_sa.first()
        ax.axhline(ks.mean(), color="0.35", ls="--", lw=1.2, label="SA (mean)")
        ax.axhspan(ks.mean() - ks.std(ddof=1), ks.mean() + ks.std(ddof=1), color="0.9", zorder=0)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xscale("symlog", linthresh=5e-4)
        ax.set_xticks([0, 5e-4, 1e-3, 2e-3, 4e-3, 8e-3])
        ax.set_xticklabels(["0", "5e-4", "1e-3", "2e-3", "4e-3", "8e-3"], fontsize=7)
        ax.minorticks_off()
        ax.set_title(f"n = {n}", fontsize=9)
        ax.set_xlabel(r"dephasing rate $\gamma$ per qubit (J)")
    axes[0, 0].set_ylabel(r"exponent $\kappa$")
    axes[0, 0].legend(fontsize=8)
    f.suptitle("QA exponent versus dephasing rate; shaded band is SA mean ± std over seeds", y=1.0)
    f.tight_layout()
    f.savefig(fig / "fig4_exponent_vs_dephasing.png", dpi=DPI)
    plt.close(f)


def ring_validation(res: Path, fig: Path) -> None:
    d = pd.read_csv(res / "validation_ring.csv")
    summ = json.load(open(res / "summary.json", encoding="utf8"))
    f, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.loglog(d.tau, d.residual, "ko-", ms=4, lw=1.2, label=f"ring, n = {summ['ring_n']}")
    ref = d.residual.iloc[2] * (d.tau / d.tau.iloc[2]) ** -0.5
    ax.loglog(d.tau, ref, "k:", lw=1.2, label=r"$\tau^{-1/2}$ reference")
    ax.set_xlabel(r"anneal time $\tau$ (1/J)")
    ax.set_ylabel("residual energy per spin (J)")
    ax.set_title(r"Ferromagnetic ring: fitted $\kappa$ = " + f"{summ['ring_kappa']:.2f}")
    ax.legend(fontsize=8)
    f.tight_layout()
    f.savefig(fig / "fig5_ring_validation.png", dpi=DPI)
    plt.close(f)


def make_all(res: Path, fig: Path) -> None:
    dumbbell(res, fig)
    small_multiples(res, fig)
    waterfall(res, fig)
    exponent_vs_gamma(res, fig)
    ring_validation(res, fig)
