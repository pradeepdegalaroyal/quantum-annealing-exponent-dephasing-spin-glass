"""Exact density-matrix reference for dephased annealing (small n), used to validate the trajectories."""

from __future__ import annotations

import numpy as np

from .problems import Problem


def density_matrix_qa(problem: Problem, tau: float, dt: float, gamma: float) -> float:
    """Per-site residual energy from the Lindblad equation with L_k = sqrt(gamma) Z_k.

    Uses the same Strang splitting as the trajectory code; the dephasing superoperator is exactly
    diagonal: rho_ab -> rho_ab exp(-2 gamma dt d_H(a, b)).
    """
    n = problem.n
    dim = 2**n
    idx = np.arange(dim)
    ham = np.bitwise_count(idx[:, None] ^ idx[None, :]) if hasattr(np, "bitwise_count") else None
    if ham is None:  # numpy < 2.0
        x = idx[:, None] ^ idx[None, :]
        ham = np.zeros_like(x)
        for k in range(n):
            ham += (x >> k) & 1
    psi = np.full(dim, 1.0 / np.sqrt(dim), dtype=complex)
    rho = np.outer(psi, psi.conj())
    steps = int(round(tau / dt))
    h = tau / steps
    damp = np.exp(-2.0 * gamma * h * ham)
    c_cache: dict[float, np.ndarray] = {}

    def x_rot(theta: float) -> np.ndarray:
        if theta not in c_cache:
            u1 = np.cos(theta) * np.eye(2) + 1j * np.sin(theta) * np.array([[0, 1], [1, 0]])
            u = np.array([[1.0]])
            for _ in range(n):
                u = np.kron(u1, u)
            c_cache[theta] = u
        return c_cache[theta]

    for k in range(steps):
        s = (k + 0.5) / steps
        d = np.exp(-1j * s * h / 2.0 * problem.energies)
        u = (d[:, None] * x_rot((1.0 - s) * h)) * d[None, :]
        # U = D X D with D diagonal
        rho = u @ rho @ u.conj().T
        rho = rho * damp
    return float((np.real(np.diag(rho)) @ problem.energies - problem.e0) / n)
