"""Power-law exponent fits: residual energy epsilon(tau) ~ tau**(-kappa)."""

from __future__ import annotations

import numpy as np


def fit_exponent(tau: np.ndarray, resid: np.ndarray, window: tuple[float, float]) -> float:
    """Least-squares slope kappa of log(resid) against log(tau) inside ``window`` (inclusive).

    Returns ``-slope`` so that a faster decay gives a larger kappa. Residuals are floored at 1e-12.
    """
    tau = np.asarray(tau, dtype=float)
    resid = np.asarray(resid, dtype=float)
    sel = (tau >= window[0]) & (tau <= window[1])
    if sel.sum() < 3:
        raise ValueError("need at least three tau values inside the fit window")
    slope = np.polyfit(np.log(tau[sel]), np.log(np.maximum(resid[sel], 1e-12)), 1)[0]
    return float(-slope)


def bootstrap_exponent(
    tau: np.ndarray,
    resid_per_instance: np.ndarray,
    window: tuple[float, float],
    n_boot: int = 500,
    seed: int = 0,
) -> tuple[float, float]:
    """Exponent of the instance-mean curve and its bootstrap standard deviation over instances.

    Args:
        tau: (T,) anneal times.
        resid_per_instance: (B, T) residual energies per instance.
    """
    rng = np.random.default_rng(seed)
    b = resid_per_instance.shape[0]
    point = fit_exponent(tau, resid_per_instance.mean(0), window)
    boots = [
        fit_exponent(tau, resid_per_instance[rng.integers(0, b, b)].mean(0), window) for _ in range(n_boot)
    ]
    return point, float(np.std(boots, ddof=1))
