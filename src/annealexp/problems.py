"""Ising problem instances: frustrated +-J spin glasses on random 3-regular graphs and ferromagnetic rings."""

from __future__ import annotations

from dataclasses import dataclass

import networkx as nx
import numpy as np


@dataclass(frozen=True)
class Problem:
    """Classical Ising cost H_P = sum_{i<j} J_ij z_i z_j on n spins.

    Attributes:
        n: number of spins.
        coupling: dense symmetric (n, n) coupling matrix with zero diagonal.
        energies: (2**n,) energy of every basis state; bit k of the index is spin k (z = 1 - 2 bit).
        e0: ground-state energy.
        label: short description of the instance.
    """

    n: int
    coupling: np.ndarray
    energies: np.ndarray
    e0: float
    label: str


def spin_values(n: int) -> np.ndarray:
    """Return the (n, 2**n) array of z eigenvalues (+1/-1); qubit k is bit k of the basis index."""
    idx = np.arange(2**n)
    return np.stack([1.0 - 2.0 * ((idx >> k) & 1) for k in range(n)]).astype(np.float64)


def energies_from_coupling(coupling: np.ndarray) -> np.ndarray:
    """Energies of all 2**n configurations for a symmetric coupling matrix."""
    n = coupling.shape[0]
    z = spin_values(n)
    return 0.5 * np.einsum("ik,ij,jk->k", z, coupling, z)


def _make(coupling: np.ndarray, label: str) -> Problem:
    e = energies_from_coupling(coupling)
    return Problem(n=coupling.shape[0], coupling=coupling, energies=e, e0=float(e.min()), label=label)


def spin_glass_problem(n: int, seed: int) -> Problem:
    """Random +-J spin glass on a random 3-regular graph (frustrated; n must be even)."""
    rng = np.random.default_rng(seed)
    g = nx.random_regular_graph(3, n, seed=int(rng.integers(2**31 - 1)))
    j = np.zeros((n, n))
    for a, b in g.edges():
        j[a, b] = j[b, a] = rng.choice([-1.0, 1.0])
    return _make(j, f"pmJ-3reg-n{n}-s{seed}")


def ring_problem(n: int) -> Problem:
    """Ferromagnetic ring J=-1 (the transverse-field Ising chain at the end of the anneal)."""
    j = np.zeros((n, n))
    for a in range(n):
        j[a, (a + 1) % n] = j[(a + 1) % n, a] = -1.0
    return _make(j, f"ring-n{n}")


def make_ensemble(n: int, size: int, seed: int) -> list[Problem]:
    """Draw ``size`` frustrated spin-glass instances (rejects instances with no frustration)."""
    out: list[Problem] = []
    k = 0
    while len(out) < size:
        p = spin_glass_problem(n, seed * 1000 + k)
        k += 1
        if p.e0 > -np.abs(np.triu(p.coupling)).sum() + 1e-9:
            out.append(p)
    return out
