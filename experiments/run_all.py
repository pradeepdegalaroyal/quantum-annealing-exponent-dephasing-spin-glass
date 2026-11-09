"""Reproduce every result and figure: ``python experiments/run_all.py`` (add ``--quick`` for a smoke run)."""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from annealexp import (  # noqa: E402
    density_matrix_qa,
    fit_exponent,
    make_ensemble,
    ring_problem,
    simulate_qa,
    simulate_sa,
)

torch.set_num_threads(6)

WINDOW = (4.0, 32.0)  # fit window fixed before any run


def crossover(tau: np.ndarray, qa: np.ndarray, sa: np.ndarray) -> float:
    """Anneal time (log-interpolated) beyond which QA residual energy stays below the SA residual.

    Returns nan when QA never gets below SA on the grid; tau[0] if it is below everywhere.
    """
    below = qa < sa
    if not below.any():
        return float("nan")
    # last index where QA is NOT below SA; crossover is just after it
    bad = np.where(~below)[0]
    if len(bad) == 0:
        return float(tau[0])
    i = bad[-1]
    if i == len(tau) - 1:
        return float("nan")
    r0, r1 = np.log(qa[i] / sa[i]), np.log(qa[i + 1] / sa[i + 1])
    f = r0 / (r0 - r1)
    return float(np.exp(np.log(tau[i]) + f * (np.log(tau[i + 1]) - np.log(tau[i]))))


def gamma_star(gammas: np.ndarray, kappa_qa: np.ndarray, kappa_sa: float) -> float:
    """Dephasing rate at which the QA exponent falls to the SA exponent (linear interpolation).

    nan if the QA exponent never drops below kappa_sa on the grid.
    """
    for i in range(len(gammas) - 1):
        a, b = kappa_qa[i], kappa_qa[i + 1]
        if a >= kappa_sa > b:
            return float(gammas[i] + (a - kappa_sa) / (a - b) * (gammas[i + 1] - gammas[i]))
    return float("nan")


def tune_sa(n: int, taus: np.ndarray, rep: int) -> float:
    """Pick the final inverse temperature of SA on held-out instances (lowest residual at the longest tau)."""
    ens = make_ensemble(n, 6, seed=99)
    best, best_res = 3.0, np.inf
    for bf in (2.0, 3.0, 5.0):
        r = simulate_sa(ens, int(taus[-1]), beta_final=bf, n_rep=rep, seed=1)["residual"].mean()
        if r < best_res:
            best, best_res = bf, r
    return best


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="smoke run (< 1 minute)")
    ap.add_argument("--out", default=str(ROOT))
    args = ap.parse_args()
    out = Path(args.out)
    res_dir, fig_dir = out / "results", out / "figures"
    res_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    q = args.quick

    taus = np.array([2, 4, 8, 16, 32.0]) if q else np.array([2, 3, 4, 6, 8, 12, 16, 24, 32.0])
    seeds = [0, 1] if q else [0, 1, 2, 3, 4]
    n_inst = 3 if q else 6
    sizes_coh = [6, 8] if q else [8, 10, 12]
    sizes_noise = [6] if q else [6, 8, 10]
    gammas = np.array([0.0, 0.002, 0.008]) if q else np.array([0.0, 0.0005, 0.001, 0.002, 0.004, 0.008])
    n_traj = {6: 40 if q else 100, 8: 80, 10: 50}
    rep = 200 if q else 1000
    dt = 0.1 if q else 0.05
    ring_n = 8 if q else 14
    summary: dict = {"window": list(WINDOW), "taus": taus.tolist(), "quick": q}

    # ---- 1. validation -------------------------------------------------
    ring = ring_problem(ring_n)
    ring_res = np.array([simulate_qa([ring], t, dt=dt)["residual"][0] for t in taus])
    pd.DataFrame({"tau": taus, "residual": ring_res}).to_csv(res_dir / "validation_ring.csv", index=False)
    summary["ring_n"] = ring_n
    summary["ring_kappa"] = fit_exponent(taus, ring_res, WINDOW)

    p0 = make_ensemble(6, 1, seed=0)[0]
    conv = {f"dt_{d}": float(simulate_qa([p0], 8.0, dt=d)["residual"][0]) for d in (0.1, 0.05, 0.025)}
    rho_val = density_matrix_qa(p0, 8.0, 0.05, 0.02)
    tr = simulate_qa([p0], 8.0, dt=0.05, gamma=0.02, n_traj=1000 if q else 3000)
    summary["dt_convergence"] = conv
    summary["dm_vs_trajectory"] = {
        "density_matrix": rho_val,
        "trajectory": float(tr["residual"][0]),
        "trajectory_sem": float(tr["residual_sem"][0]),
    }
    print("validation done", round(time.time() - t_start, 1), "s", flush=True)

    # ---- 2. coherent QA vs SA, scaling with n ---------------------------
    rows, exp_rows, tune = [], [], {}
    for n in sizes_coh:
        bf = tune_sa(n, taus, rep)
        tune[n] = bf
        ens = [make_ensemble(n, n_inst, seed=s) for s in seeds]
        flat = [p for e in ens for p in e]
        qa = np.array([simulate_qa(flat, t, dt=dt)["residual"] for t in taus]).T  # (S*B, T)
        sa = np.array([simulate_sa(flat, int(t), beta_final=bf, n_rep=rep, seed=5)["residual"] for t in taus]).T
        for si, s in enumerate(seeds):
            sl = slice(si * n_inst, (si + 1) * n_inst)
            mq, ms = qa[sl].mean(0), sa[sl].mean(0)
            for ti, t in enumerate(taus):
                rows.append({"n": n, "seed": s, "tau": t, "qa": mq[ti], "sa": ms[ti]})
            exp_rows.append(
                {
                    "n": n,
                    "seed": s,
                    "kappa_qa": fit_exponent(taus, mq, WINDOW),
                    "kappa_sa": fit_exponent(taus, ms, WINDOW),
                    "tau_cross": crossover(taus, mq, ms),
                }
            )
        print("coherent n", n, round(time.time() - t_start, 1), "s", flush=True)
    pd.DataFrame(rows).to_csv(res_dir / "coherent_curves.csv", index=False)
    pd.DataFrame(exp_rows).to_csv(res_dir / "coherent_exponents.csv", index=False)
    summary["sa_beta_final"] = {str(k): v for k, v in tune.items()}

    # ---- 3. dephasing sweep (trajectories) ------------------------------
    noisy_g = gammas[gammas > 0]
    nrows, nexp = [], []
    for n in sizes_noise:
        ens = [make_ensemble(n, n_inst, seed=s) for s in seeds]
        flat = [p for e in ens for p in e]
        nb = len(flat)
        big_probs = [p for _ in noisy_g for p in flat]
        big_gam = np.repeat(noisy_g, nb)
        curves = {0.0: np.array([simulate_qa(flat, t, dt=dt)["residual"] for t in taus]).T}
        acc = []
        for t in taus:
            r = simulate_qa(big_probs, t, dt=dt, gamma=big_gam, n_traj=n_traj[n], seed=11)["residual"]
            acc.append(r.reshape(len(noisy_g), nb))
        acc = np.stack(acc, axis=-1)  # (G, S*B, T)
        for gi, g in enumerate(noisy_g):
            curves[float(g)] = acc[gi]
        bf = tune_sa(n, taus, rep)
        sa = np.array([simulate_sa(flat, int(t), beta_final=bf, n_rep=rep, seed=5)["residual"] for t in taus]).T
        for si, s in enumerate(seeds):
            sl = slice(si * n_inst, (si + 1) * n_inst)
            ms = sa[sl].mean(0)
            ks_sa = fit_exponent(taus, ms, WINDOW)
            kq = []
            for g in gammas:
                mq = curves[float(g)][sl].mean(0)
                kq.append(fit_exponent(taus, mq, WINDOW))
                for ti, t in enumerate(taus):
                    nrows.append({"n": n, "seed": s, "gamma": g, "tau": t, "qa": mq[ti], "sa": ms[ti]})
                nexp.append(
                    {
                        "n": n,
                        "seed": s,
                        "gamma": g,
                        "kappa_qa": kq[-1],
                        "kappa_sa": ks_sa,
                        "tau_cross": crossover(taus, mq, ms),
                        "res_at_tau_max": mq[-1],
                        "sa_res_at_tau_max": ms[-1],
                    }
                )
        print("noise n", n, round(time.time() - t_start, 1), "s", flush=True)
    pd.DataFrame(nrows).to_csv(res_dir / "dephasing_curves.csv", index=False)
    nexp_df = pd.DataFrame(nexp)
    nexp_df.to_csv(res_dir / "dephasing_exponents.csv", index=False)

    # gamma* from the seed-mean exponent curve and per seed
    gstar = {}
    for n in sizes_noise:
        d = nexp_df[nexp_df.n == n]
        km = d.groupby("gamma").kappa_qa.mean().reindex(gammas).to_numpy()
        ksa = d.groupby("seed").kappa_sa.first().mean()
        per_seed = []
        for s in seeds:
            ds = d[d.seed == s].sort_values("gamma")
            per_seed.append(gamma_star(gammas, ds.kappa_qa.to_numpy(), float(ds.kappa_sa.iloc[0])))
        gstar[str(n)] = {
            "gamma_star_of_mean": gamma_star(gammas, km, float(ksa)),
            "gamma_star_per_seed": per_seed,
        }
    summary["gamma_star"] = gstar

    # ---- 4. noise-source ablation (n = 8 in quick, 10 otherwise) --------
    n_ab = sizes_noise[-1] if q else 10
    ens = [make_ensemble(n_ab, n_inst, seed=s) for s in seeds]
    flat = [p for e in ens for p in e]
    nb = len(flat)
    g_ref, sig_ref = 0.001, 0.10
    configs = {
        "ideal": dict(gamma=0.0, control_sigma=0.0),
        "control_error": dict(gamma=0.0, control_sigma=sig_ref),
        "dephasing": dict(gamma=g_ref, control_sigma=0.0),
        "both": dict(gamma=g_ref, control_sigma=sig_ref),
    }
    ab_rows = []
    for name, kw in configs.items():
        r = np.array(
            [
                simulate_qa(flat, t, dt=dt, n_traj=n_traj.get(n_ab, 50), seed=21, **kw)["residual"]
                for t in taus
            ]
        ).T
        for si, s in enumerate(seeds):
            mq = r[si * n_inst : (si + 1) * n_inst].mean(0)
            ab_rows.append({"config": name, "seed": s, "kappa": fit_exponent(taus, mq, WINDOW), "res_tau_max": mq[-1]})
    pd.DataFrame(ab_rows).to_csv(res_dir / "ablation_noise_sources.csv", index=False)
    summary["ablation"] = {"n": n_ab, "gamma": g_ref, "control_sigma": sig_ref}

    # ---- summary ---------------------------------------------------------
    cdf = pd.DataFrame(exp_rows)
    summary["coherent"] = {
        str(n): {
            "kappa_qa_mean": float(g.kappa_qa.mean()),
            "kappa_qa_std": float(g.kappa_qa.std(ddof=1)),
            "kappa_sa_mean": float(g.kappa_sa.mean()),
            "kappa_sa_std": float(g.kappa_sa.std(ddof=1)),
            "qa_beats_sa_in_seeds": int((g.kappa_qa > g.kappa_sa).sum()),
            "tau_cross_mean": float(g.tau_cross.mean()),
            "tau_cross_std": float(g.tau_cross.std(ddof=1)),
            "n_seeds": int(len(g)),
        }
        for n, g in cdf.groupby("n")
    }
    summary["runtime_seconds"] = round(time.time() - t_start, 1)
    summary["device"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"
    with open(res_dir / "summary.json", "w", encoding="utf8") as f:
        json.dump(summary, f, indent=2)

    import make_figures

    make_figures.make_all(res_dir, fig_dir)
    print("done in", summary["runtime_seconds"], "s")


if __name__ == "__main__":
    main()
