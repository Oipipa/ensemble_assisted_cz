# Numerical findings

## Model and fidelity

We model the ensemble-assisted CZ gate of [Zhang, Arunseangroj, and Xu (2025)](https://arxiv.org/html/2503.16896v1) with two data atoms and one ancilla. Data $|1\rangle\leftrightarrow|r\rangle$ and ancilla $|A\rangle\leftrightarrow|R\rangle$ are driven in a $\pi$–$2\pi$–$\pi$ sequence. The resonant pair coupling $|rR\rangle\leftrightarrow|r'R'\rangle$ follows Supplementary Eq. (29). At the 10 MHz baseline, the gate lasts 0.20 µs and $V_{Rr}=2\pi\times200$ MHz.

We compute the Haar-average fidelity of the **reduced two-data-qubit operation** to the paper's $-CZ$ target. The ancilla is traced out, and population outside the computational subspace lowers fidelity. Data Rydberg decay at $1/100\,\mu\mathrm{s}^{-1}$ and ancilla decay at $1/200\,\mu\mathrm{s}^{-1}$ are explicit quantum jumps. The calculation in [reproduce.py](reproduce.py) evolves the no-jump amplitudes and integrates all ancilla-jump trajectories. Data jumps enter the non-Hermitian survival evolution and cannot return to the target subspace. This is equivalent to the specified absorbing-loss Lindblad model for this fidelity, rather than a decay subtraction. We assign the same species lifetime to the primed Rydberg levels, whose lifetimes are not separately specified in the paper.

## Baseline reproduction

| Quantity at 10 MHz | Paper | This calculation |
| --- | ---: | ---: |
| Noisy CZ fidelity | 99.80% (Table II) | **99.8008%** |
| Noiseless Förster fidelity | — | 99.9806% |
| Maximum noiseless computational leakage | — | 0.0775% |

Switching on one decay mechanism at a time gives 0.1498 percentage points of data-decay error and 0.0300 percentage points of ancilla-decay error relative to the noiseless Förster result. The remaining noiseless blockade infidelity is 0.0194 percentage points. These round to the paper's Table II estimates of 0.15%, 0.03%, and 0.02%, respectively.

The original diagonal-shift model gave 99.7351% noisy fidelity and a conditional phase error of 0.0784 rad. That interaction is useful as a contrasting approximation, but it is not the resonant Förster coupling used to explain the paper's high fidelity. Combining its decay contributions with the Förster blockade error was not a consistent reproduction. The full noisy Förster result now agrees with the paper's **rounded total fidelity**. The isolated error terms are calculated by switching off other decay channels, so their sum need not equal the combined infidelity exactly.

## Drive-frequency dip and pulse shaping

The old plot smoothed the two data pulses only. During those pulses, the ancilla is predominantly $|A\rangle$, so the resonant pair interaction is inactive. Thus, data smoothing changed the curves negligibly. We infer that smoothing the ancilla $2\pi$ pulse is the relevant test of the paper's mechanism, since $|R\rangle$ and the coupled pair states are populated then. Square and smooth pulses have identical duration and integrated area. The smooth envelope is a truncated Gaussian with zero endpoints, width $0.25$ of the pulse duration, and area normalization. Its peak drive is higher than the square pulse's; the x axis denotes the **area-equivalent** Rabi frequency, not equal peak amplitude.

| Drive | Square noisy fidelity | Smooth ancilla noisy fidelity | Square / smooth maximum noiseless leakage |
| ---: | ---: | ---: | ---: |
| 10 MHz | 99.8008% | 99.8201% | 0.0775% / 0.00009% |
| 30 MHz | 98.9020% | 99.9392% | 1.528% / 0.0014% |
| 37.5 MHz | 98.1468% | 99.9492% | 2.736% / 0.0055% |
| 40 MHz | 99.9052% | 99.9426% | 0.150% / 0.0220% |

![Noisy fidelity and leakage](figures/drive_sweep_noisy.png)

The 30 MHz dip is a coherent finite-coupling effect in this truncated model, not a visible solver artifact. Tightening the solver from relative tolerance $2\times10^{-9}$ to $2\times10^{-11}$ and doubling jump-integration quadrature from 32 to 64 nodes changes the 30 MHz fidelities by less than $2\times10^{-9}$. The square curve oscillates and recovers at 40 MHz; the smooth ancilla curve strongly suppresses those oscillations. This reproduces the qualitative behavior of Supplementary Fig. 4(a).

The paper does not provide enough pulse-shape detail for a unique point-by-point match to that figure. Our model also omits other Rydberg pair states, interaction inhomogeneity, and other experimental errors. The reproduced 99.80% baseline and high-drive trend are therefore results for the stated truncated model and pulse constraint.
