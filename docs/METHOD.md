# Method

## Model
Spins on n qubits with the anneal Hamiltonian H(s) = -(1 - s) sum_k X_k + s H_P, s = t / tau, linear schedule, tau in units of 1/J. The problem H_P = sum_{i<j} J_ij Z_i Z_j is a frustrated +-J spin glass on a random 3-regular graph (6 instances per seed; frustration enforced) or a ferromagnetic ring (validation). The initial state is the uniform superposition.

## Simulation
- **Coherent**: exact state-vector evolution with a second-order split step (half diagonal step, X rotation, half diagonal step), step 0.05 in 1/J; convergence checked at 0.1, 0.05 and 0.025.
- **Dephasing**: Lindblad operator sqrt(gamma) Z_k on every qubit, unravelled exactly as white-noise phase kicks with variance gamma dt and averaged over trajectories (50 to 100 per instance). The trajectory average is validated against an exact density-matrix integration (`dm.py`).
- **Static control error** (ablation only): each coupling is multiplied by 1 + sigma xi with xi standard normal, in the Hamiltonian only; the residual energy is measured against the intended problem.

## Baselines
- **Classical**: Metropolis simulated annealing (SA) with a linear ramp of the inverse temperature from 0.05 to a final value chosen from {2, 3, 5} on held-out instances (lowest residual at the longest time), 1000 replicas per instance. One sweep is identified with one unit of 1/J. This identification sets the horizontal offset of the two curves but not the exponent, so the exponent is the primary comparison.
- **Standard quantum baseline**: coherent (noise-free) quantum annealing (QA) with the same linear schedule, the protocol of the annealing experiment in the literature.

## Metrics
- Residual energy per spin, (<H_P> - E_0) / n, with E_0 from exhaustive enumeration.
- Exponent kappa: least-squares slope of log residual versus log tau over the window 4 <= tau <= 32, fixed before any run, fitted to the instance-mean curve of each seed (9 anneal times).
- tau_cross: anneal time beyond which QA stays below SA (log interpolation).
- gamma*: dephasing rate at which kappa(QA) falls to kappa(SA) (linear interpolation on the grid gamma in {0, 0.0005, 0.001, 0.002, 0.004, 0.008}).

## Experiments
1. Validation: ring exponent, step-size convergence, trajectories versus density matrix.
2. Scaling: coherent QA versus SA at n = 8, 10, 12, 5 seeds.
3. Robustness sweep: dephasing rate at n = 6, 8, 10, 5 seeds.
4. Ablation: noise-source ablation at n = 10 (ideal, 10 % control error, dephasing at gamma = 0.001, both).

Seeds: 0 to 4; each seed draws its own instance ensemble, so the reported spread includes instance-to-instance variation. A fixed random seed controls noise and SA replicas.
