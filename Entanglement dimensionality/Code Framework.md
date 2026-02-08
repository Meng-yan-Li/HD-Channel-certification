# See-Saw Optimization Framework for Certifying Channel's Entanglement Dimensionality

This project implements a numerical see-saw (alternating) optimization framework for studying **high-dimensional quantum channels**. The code is designed for **semi-device-independent (SDI)** analysis, where minimal assumptions are made on internal device structures, while key physical constraint (dimension) are explicitly enforced.

The implementation is modular and MATLAB-based, relying heavily on **semidefinite programming (SDP)** solved via **CVX** with external solvers such as MOSEK, SDPT3, or SCS. **QETLAB** was also used.

---

## Overall Structure

The project consists of:

- A **main script** (`main_ED.m`) that defines numerical tasks and runs simulations.
- Core **see-saw optimization routines** (SN constraint and noisy).
- A collection of **helper functions** for POVM optimization, channel construction, partial traces, and random channel generation.

Each component is briefly described below.

---

## Main Script

### `main_ED.m`

This is the **entry point** of the project. It defines and executes numerical experiments in two main scenarios:

1. **TASK 1 (SN constraint)**  
   Runs the see-saw optimization for a prepare-and-measure task with fixed parameters `(n, d, r)`. The script supports both single runs and batch simulations over different dimensions and Schmidt numbers.

2. **TASK 2 (Noisy scenario)**  
   Studies the effect of noise by introducing a communication channel (depolarizing or dephasing). The objective value is evaluated as a function of the noise parameter `v`, averaged over multiple random initializations.

The script also handles data storage for later analysis and plotting.

---

## Core Optimization Routines

### `see_saw.m`

Implements the **see-saw optimization algorithm** in the noiseless setting. The algorithm alternates between optimizing:

1. **Preparation states** $ \{\rho_x\} $
2. **Measurement POVMs** $ \{M_{b|y}\} $
3. **Bipartite Choi state** $ \Phi $ with a bounded Schmidt number

Each subproblem is convex when the other variables are fixed, allowing efficient optimization via SDPs. The routine iterates until convergence of the objective value.

---

### `see_saw_noisy.m`

A variant of `see_saw.m` that incorporates **noise in the communication channel**. All relevant operators are passed through a noisy channel map before being used in the optimization or objective evaluation.

This function is used to study robustness of certification results against realistic imperfections.

---

## SDP Subroutines

### `optimize_POVM.m`

Solves an SDP to optimize a POVM $ \{M_b\} $​ that maximizes a linear objective of the form
$$
\sum_b \mathrm{Tr}(G_b M_b)
$$
subject to:

- $ M_b \ge 0 $ for all outcomes
- $ \sum_b M_b = I $

This subroutine is called repeatedly inside the see-saw loop.

---

### `optimize_phi.m`

Optimizes a **bipartite quantum state** $ \Phi \in \mathcal{H}_d \otimes \mathcal{H}_d $ by solving an SDP that maximizes
$$
\mathrm{Tr}(A \Phi) 
$$
subject to:

- Positivity and unit trace
- An **entanglement constraint** enforced via the **generalized reduction map**, ensuring that the Schmidt number of $ \Phi $ is bounded by `r`

The file also documents alternative (commented) constraints for generalized DPS/PPT hierarchies.

---

## Channel and Utility Functions

### `commu_channel.m`

Implements simple noisy quantum channels acting on operators:

- **Depolarizing channel**
- **Dephasing channel**

The noise strength is controlled by a parameter `v`. This function is used consistently in both state and measurement updates in the noisy scenario.

---

### `random_SN_channel.m`

Generates a **random CPTP quantum channel** whose Choi state has Schmidt number:

- exactly `r` (if `exactRank = true`), or
- at most `r`

The function explicitly constructs Kraus operators with controlled rank, enforces the CPTP condition, and builds the corresponding Choi state using the Jamiołkowski isomorphism.

---

### `index2code_array.m`

Utility function that converts linear indices into base-`d` strings of length `n`. This is used to systematically enumerate all preparation settings $ x \in \{1,\dots,d\}^n $.

---

## Intended Use

This codebase is intended for:

- Numerical certification of **high-dimensional quantum channels**
- Studying **Schmidt number / entanglement dimensionality** in SDI scenarios
- Exploring **noise robustness** of prepare-and-measure protocols

The modular structure allows easy modification of constraints, noise models, and optimization targets.

---

## Dependencies

- MATLAB
- CVX (with SDP solvers such as MOSEK, SDPT3, or SCS)
- QETLAB

---

## Notes

- All optimization problems are convex at each see-saw step, but the overall algorithm is heuristic and may converge to local optima.
- Multiple random initializations are recommended for reliable results.

---

This document serves as a high-level guide to the codebase. For mathematical details and physical interpretation, see the accompanying manuscript or notes.