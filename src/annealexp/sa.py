"""Classical baseline: Metropolis simulated annealing with a linear inverse-temperature ramp."""

from __future__ import annotations

import numpy as np
import torch

from .problems import Problem
from .qa import default_device


@torch.no_grad()
def simulate_sa(
    problems: list[Problem],
    sweeps: int,
    beta_final: float = 3.0,
    beta_init: float = 0.05,
    n_rep: int = 1000,
    seed: int = 0,
    device: torch.device | None = None,
) -> dict[str, np.ndarray]:
    """Run ``n_rep`` independent annealing replicas per instance for ``sweeps`` sequential sweeps.

    Returns:
        dict with ``residual`` (per-site residual energy, mean over replicas, (B,)) and ``p_gs`` ((B,)).
    """
    dev = device or default_device()
    n = problems[0].n
    nb = len(problems)
    gen = torch.Generator(device=dev).manual_seed(seed)
    jm = torch.tensor(np.stack([p.coupling for p in problems]), dtype=torch.float32, device=dev)  # (B,n,n)
    spins = torch.where(
        torch.rand(nb, n_rep, n, generator=gen, device=dev) < 0.5,
        torch.tensor(1.0, device=dev),
        torch.tensor(-1.0, device=dev),
    )
    for t in range(sweeps):
        beta = beta_init + (beta_final - beta_init) * (t + 0.5) / sweeps
        for i in torch.randperm(n, generator=gen, device=dev).tolist():
            field = (spins * jm[:, None, i, :]).sum(-1)  # (B,R)
            d_e = -2.0 * spins[:, :, i] * field
            acc = torch.rand(nb, n_rep, generator=gen, device=dev) < torch.exp(-beta * d_e).clamp(max=1.0)
            spins[:, :, i] = torch.where(acc, -spins[:, :, i], spins[:, :, i])
    energy = 0.5 * ((spins @ jm) * spins).sum(-1)  # (B,R)
    e0 = torch.tensor([p.e0 for p in problems], dtype=torch.float32, device=dev)[:, None]
    res = ((energy - e0) / n).double()
    pgs = (energy <= e0 + 1e-6).double()
    return {"residual": res.mean(1).cpu().numpy(), "p_gs": pgs.mean(1).cpu().numpy()}
