# Research gap

## State of the art (from the supplied literature)

Programmable quantum annealers now reach thousands of qubits. King et al., 2023, doi:10.1038/s41586-023-05867-2 realise quantum-critical dynamics of a three-dimensional spin glass on a superconducting annealer, agree quantitatively with Schrödinger evolution on small instances, and extract critical exponents that separate quantum annealing from the slower stochastic dynamics of Monte Carlo (MC) algorithms. The condensed-matter reading of such experiments is given by Fisher et al., 2023, doi:10.1146/annurev-conmatphys-031720-030658, who frame circuit and annealing dynamics as a tractable setting for universal non-equilibrium physics and digital quantum simulation. On the error-control side, Sivak et al., 2023, doi:10.1038/s41586-023-05782-6 state that whether quantum error correction (QEC) can practically extend coherence is an open question, and Bluvstein et al., 2024, doi:10.1038/s41586-023-06927-3 show that the overhead of logical qubits is the central obstacle to large logical processors. Neutral-atom and erasure-conversion work (Evered et al., 2023, doi:10.1038/s41586-023-06481-y; Ma et al., 2023, doi:10.1038/s41586-023-06438-1) lowers physical error rates, but a residual dephasing rate remains. Recent work stresses that classical simulation is the step that decides feasibility before hardware (Nguyen et al., 2027, doi:10.1007/978-981-92-2587-3_22, on graphics processing units) and that optical Ising machines compete for the same long-range Ising problems (Zhou et al., 2027, doi:10.1016/j.optcom.2026.133796).

## Related work

The role column is my own grouping of the supplied papers: landmark = highly cited, review = survey, recent = dated 2026 to 2027, adjacent = neighbouring method or platform.

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| King et al., doi:10.1038/s41586-023-05867-2 | 2023 | landmark | Quantum-critical spin-glass dynamics on about 5,000 qubits; exponents separate quantum annealing (QA) from MC | A controlled study of how environmental dephasing moves the extracted exponent, or when a tuned classical baseline closes the gap at small n |
| Fisher et al., doi:10.1146/annurev-conmatphys-031720-030658 | 2023 | review | Random circuits as a model of non-equilibrium many-body dynamics | A quantitative noise budget for annealing-type dynamics |
| Bauer et al., doi:10.1103/PRXQuantum.4.027001 | 2023 | review | Quantum simulation for high-energy physics roadmap | Annealing exponents or finite-size effects |
| Sheremet et al., doi:10.1103/RevModPhys.95.015002 | 2023 | review | Collective radiance of emitters in waveguides | Optimisation dynamics on frustrated models |
| Sivak et al., doi:10.1038/s41586-023-05782-6 | 2023 | landmark | Real-time QEC beyond break-even | Whether unprotected annealing tolerates its residual dephasing |
| Bluvstein et al., doi:10.1038/s41586-023-06927-3 | 2024 | landmark | Logical processor on up to 280 physical atoms; scrambling-type circuits | Analogue annealing without encoding |
| Evered et al., doi:10.1038/s41586-023-06481-y | 2023 | landmark | High-fidelity parallel entangling gates on neutral atoms | Dephasing during slow analogue sweeps |
| Ma et al., doi:10.1038/s41586-023-06438-1 | 2023 | landmark | Gates with mid-circuit erasure conversion | Effect of the remaining dephasing on many-body sweeps |
| Shirakawa et al., doi:10.1016/j.future.2026.108731 | 2027 | recent | Closed-loop quantum-classical electronic-structure calculation at full scale | Spin-glass optimisation or exponent scaling |
| Nguyen et al., doi:10.1007/978-981-92-2587-3_22 | 2027 | recent | Graphics processing unit (GPU) acceleration of variational classifier simulation | Many-body dynamics or physical noise models |
| Miao and Jarvis, doi:10.1007/978-3-032-35248-4_1 | 2027 | recent | Tensor-layout optimisation for contractions used in quantum simulation | Any physics benchmark |
| Zhou et al., doi:10.1016/j.optcom.2026.133796 | 2027 | adjacent | Photonic Ising machine encoding long-range couplings | A comparison of exponents with coherent quantum dynamics |

## The gap

The experimental claim of King et al. rests on an exponent: the way the residual energy falls with anneal time. An exponent is only a fair discriminator if (i) the classical competitor is tuned and measured in the same window, and (ii) the quantum exponent survives the dephasing that every platform in the supplied hardware papers still has. Neither the annealing paper nor the error-correction papers quantify the second point: Sivak et al. leave open whether coherence can be extended in practice, and the gate papers report error rates rather than the effect of dephasing on a slow many-body sweep. The finite-size regime also matters, because the exponent of a small system is a crossover quantity.

**Testable gap.** At small, exactly simulable sizes (n = 6 to 12 qubits on frustrated spin glasses), how large a per-qubit dephasing rate can coherent annealing tolerate before its fitted residual-energy exponent falls to that of tuned simulated annealing (SA), and how does this rate depend on n?

**Hypotheses (falsifiable).**
- H1: coherent QA has a larger fitted exponent than tuned SA at n = 8, 10 and 12 in at least 4 of 5 seeds, in a fit window fixed before the runs.
- H2: the break-even dephasing rate gamma* is finite on the swept grid and decreases with n, because dephasing heats every qubit and the heating is extensive.
- H3: at the reference noise levels, Z-dephasing removes more exponent than a 10 % static coupling error.

**Success metrics.** The fitted exponent kappa (mean ± std over 5 seeds), the anneal time at which QA overtakes SA, and gamma* per n. Results are simulation only; no quantum advantage is claimed.

**Novelty in one sentence.** This is a measurement of the dephasing rate at which the annealing exponent of frustrated spin glasses drops to that of a tuned classical annealer, with a trajectory method validated against the exact density matrix and an ablation over two noise sources.
