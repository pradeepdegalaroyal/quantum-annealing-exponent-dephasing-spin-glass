"""Coherent and dephased quantum annealing by second-order split-step (Strang) evolution.

H(s) = -(1 - s) sum_k X_k + s H_P,  s = t / tau.  Z-dephasing of rate gamma per qubit
(Lindblad operator sqrt(gamma) Z_k) is unravelled exactly as white-noise phase kicks
exp(-i phi Z_k), phi ~ N(0, gamma dt), and averaged over trajectories.
"""

from __future__ import annotations

import numpy as np
import torch

from .problems import Problem, spin_values


def default_device() -> torch.device:
    """CUDA if available, else CPU."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _x_rotation(psi: torch.Tensor, n: int, theta: float) -> torch.Tensor:
    """Apply prod_k exp(i theta X_k) to a (N, 2**n) batch."""
    c, s = float(np.cos(theta)), float(np.sin(theta))
    big = psi.shape[0]
    for k in range(n):
        v = psi.reshape(big, 2 ** (n - 1 - k), 2, 2**k)
        a, b = v[:, :, 0, :], v[:, :, 1, :]
        psi = torch.stack((c * a + 1j * s * b, c * b + 1j * s * a), dim=2).reshape(big, -1)
    return psi


@torch.no_grad()
def simulate_qa(
    problems: list[Problem],
    tau: float,
    dt: float = 0.05,
    gamma: float | np.ndarray = 0.0,
    n_traj: int = 1,
    seed: int = 0,
    control_sigma: float = 0.0,
    device: torch.device | None = None,
) -> dict[str, np.ndarray]:
    """Anneal every problem for total time ``tau`` and return residual energies.

    Args:
        problems: instances with equal n.
        tau: anneal time in units of 1/J.
        dt: split-step size.
        gamma: per-qubit Z-dephasing rate (units of J), a scalar or one value per problem; 0 gives exact
            coherent evolution.
        n_traj: noise trajectories per instance (ignored when gamma == 0).
        seed: noise / control-error seed.
        control_sigma: relative static error J_ij -> J_ij (1 + sigma xi) in the Hamiltonian only
            (the residual energy is still measured against the intended problem).
        device: torch device.

    Returns:
        dict with ``residual`` (per-site residual energy, mean over trajectories, shape (B,)),
        ``p_gs`` (ground-state probability, (B,)) and ``residual_sem`` (trajectory standard error, (B,)).
    """
    dev = device or default_device()
    n = problems[0].n
    nb = len(problems)
    gen = torch.Generator(device=dev).manual_seed(seed)
    rng = np.random.default_rng(seed + 7919)

    e_true = np.stack([p.energies for p in problems])
    e_ham = e_true.copy()
    z = spin_values(n)
    if control_sigma > 0:
        for i, p in enumerate(problems):
            jn = np.triu(p.coupling * (1.0 + control_sigma * rng.standard_normal(p.coupling.shape)), 1)
            jn = jn + jn.T
            e_ham[i] = 0.5 * np.einsum("ik,ij,jk->k", z, jn, z)
    gam = np.broadcast_to(np.asarray(gamma, dtype=float), (nb,)).copy()
    noisy = bool(gam.max() > 0)
    m = n_traj if noisy else 1
    g_vec = torch.tensor(np.repeat(gam, m), dtype=torch.float64, device=dev)[:, None]
    e_h = torch.tensor(np.repeat(e_ham, m, axis=0), dtype=torch.float64, device=dev)
    zmat = torch.tensor(z, dtype=torch.float64, device=dev)  # (n, 2**n)
    big = nb * m
    dim = 2**n
    psi = torch.full((big, dim), 1.0 / np.sqrt(dim), dtype=torch.complex128, device=dev)

    steps = int(round(tau / dt))
    h = tau / steps
    sd = torch.sqrt(g_vec * h / 2.0)  # (big, 1)

    def diag_half(s: float) -> torch.Tensor:
        phase = (s * h / 2.0) * e_h
        if noisy:
            phase = phase + (torch.randn(big, n, generator=gen, dtype=torch.float64, device=dev) * sd) @ zmat
        return torch.exp(-1j * phase)

    for k in range(steps):
        s = (k + 0.5) / steps
        psi = psi * diag_half(s)
        psi = _x_rotation(psi, n, (1.0 - s) * h)
        psi = psi * diag_half(s)

    prob = (psi.conj() * psi).real
    et = torch.tensor(np.repeat(e_true, m, axis=0), dtype=torch.float64, device=dev)
    e0 = torch.tensor(np.repeat([p.e0 for p in problems], m), dtype=torch.float64, device=dev)
    res = (((prob * et).sum(1) - e0) / n).reshape(nb, m)
    pgs = (prob * (et <= e0[:, None] + 1e-9)).sum(1).reshape(nb, m)
    if m > 1:
        sem = res.std(1, unbiased=True) / np.sqrt(m)
    else:
        sem = torch.zeros(nb, dtype=torch.float64, device=dev)
    return {
        "residual": res.mean(1).cpu().numpy(),
        "p_gs": pgs.mean(1).cpu().numpy(),
        "residual_sem": sem.cpu().numpy(),
    }
