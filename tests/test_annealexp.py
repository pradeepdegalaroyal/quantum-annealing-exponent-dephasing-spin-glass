"""Fast tests: physics limits, trajectory-versus-density-matrix agreement, fit and crossover helpers."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))

from annealexp import (  # noqa: E402
    density_matrix_qa,
    fit_exponent,
    make_ensemble,
    ring_problem,
    simulate_qa,
    simulate_sa,
)


def test_ground_energy_and_frustration():
    ens = make_ensemble(8, 3, seed=0)
    for p in ens:
        assert p.energies.shape == (256,)
        assert p.e0 == pytest.approx(p.energies.min())
        assert p.e0 > -np.abs(np.triu(p.coupling)).sum()  # frustrated


def test_qa_long_anneal_reaches_ground_state_and_beats_short():
    p = [ring_problem(6)]
    short = simulate_qa(p, 1.0, dt=0.1)["residual"][0]
    long = simulate_qa(p, 40.0, dt=0.1)["residual"][0]
    assert long < 0.05 < short
    assert simulate_qa(p, 40.0, dt=0.1)["p_gs"][0] > 0.8


def test_trajectories_match_density_matrix():
    p = make_ensemble(4, 1, seed=3)[0]
    exact = density_matrix_qa(p, 4.0, 0.1, 0.05)
    out = simulate_qa([p], 4.0, dt=0.1, gamma=0.05, n_traj=3000, seed=2)
    assert abs(out["residual"][0] - exact) < 5 * out["residual_sem"][0] + 1e-3


def test_dephasing_increases_residual_energy():
    ens = make_ensemble(6, 2, seed=1)
    clean = simulate_qa(ens, 16.0, dt=0.1)["residual"].mean()
    dirty = simulate_qa(ens, 16.0, dt=0.1, gamma=0.02, n_traj=60, seed=0)["residual"].mean()
    assert dirty > clean


def test_sa_improves_with_sweeps():
    ens = make_ensemble(8, 3, seed=2)
    a = simulate_sa(ens, 2, n_rep=300)["residual"].mean()
    b = simulate_sa(ens, 40, n_rep=300)["residual"].mean()
    assert b < a


def test_fit_exponent_recovers_power_law():
    tau = np.array([2, 4, 8, 16, 32.0])
    assert fit_exponent(tau, 3.0 * tau**-0.7, (2, 32)) == pytest.approx(0.7)


def test_crossover_and_gamma_star():
    from run_all import crossover, gamma_star

    tau = np.array([2, 4, 8, 16.0])
    qa = np.array([0.5, 0.3, 0.1, 0.01])
    sa = np.array([0.2, 0.15, 0.12, 0.1])
    x = crossover(tau, qa, sa)
    assert 4 < x < 8
    assert np.isnan(crossover(tau, sa, qa * 0 + 0.001))
    gs = gamma_star(np.array([0.0, 1.0, 2.0]), np.array([3.0, 2.0, 0.0]), 1.0)
    assert gs == pytest.approx(1.5)
