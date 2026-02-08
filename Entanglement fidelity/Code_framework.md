# Quantum Channel Entanglement Fidelity 

This folder contains MATLAB implementations for computing **entanglement fidelity bounds of quantum channels** using a hierarchy of semidefinite programming (SDP) relaxations based on localizing matrices.

The framework applies to **semi-device-independent (SDI)** scenarios, where only the system dimension is assumed to be known, while all internal device details remain uncharacterized.

The **main entry file** of this part is:

---

## Overview

The purpose of this code is to compute a **lower bound on the entanglement fidelity** of a quantum channel from an observed **Average Success Probability (ASP)**.

The workflow consists of:
1. Generating symbolic operator lists up to a fixed hierarchy level,
2. Constructing channel-dependent and channel-independent moment matrices,
3. Extracting unique moment variables,
4. Solving a hierarchy of SDPs to certify entanglement preservation.

---

## File Structure and Function Descriptions

### 1. `main_EF.m`

**Main script for entanglement fidelity computation**

This file coordinates the full numerical procedure. It:
- Sets physical parameters such as dimension `d`, number of preparations `n`, and hierarchy level `k`,
- Generates operator lists and corresponding moment matrices,
- Extracts unique moment matrix entries,
- Evaluates the channel hierarchy SDP over a range of ASP values.

The output is a numerical estimate of the **entanglement fidelity bound** as a function of ASP.

---

### 2. `generateOperatorList.m`

**Symbolic operator list generation**

Generates all admissible operator monomials up to hierarchy level `k`, determined by:
- The number of preparations `n`,
- The Hilbert space dimension `d`.

Each operator is encoded symbolically as a pair of index sequences corresponding to preparation and measurement operators.

---

### 3. `generateMomentMatrix.m`

**Moment matrix construction**

Constructs two moment matrices:
- `Γ^Φ` (channel-dependent),
- `Γ^IZ` (channel-independent reference).

For each pair of operator words, the function:
- Applies conjugation,
- Concatenates operator sequences,
- Simplifies them using algebraic constraints,
- Stores both symbolic and string-encoded representations.

The string-encoded matrices are used to identify unique SDP variables.

---

### 4. `simplifyOperator.m`

**Operator word simplification**

Applies algebraic and structural constraints to operator sequences, including:
- Removal of zero operators,
- Consistency checks for operator adjacency,
- Elimination of adjacent duplicates,
- Canonical cyclic reordering (when required).

Invalid operator products are discarded.  
This step ensures a consistent and physically meaningful moment matrix.

---

### 5. `deleteDuplicates.m`

**Duplicate elimination**

Removes duplicate operator entries using string-based hashing while preserving order.
This reduces the number of SDP variables and improves computational efficiency.

---

### 6. `channelHierarchy.m`

**SDP formulation and solution**

Given:
- An observed ASP value,
- Sets of unique moment variables,
- Index mappings of moment matrices,

this function formulates and solves an SDP that lower-bounds the **entanglement fidelity** of the quantum channel.

The SDP enforces:
- Positive semidefiniteness of moment matrices,
- Linear constraints imposed by the observed ASP,
- Normalization and dimension constraints.

---

## Dependencies

- MATLAB
- CVX (with SDP solvers such as MOSEK, SDPT3, or SCS)
- QETLAB

---

## Physical Interpretation

- **Average Success Probability (ASP)** characterizes operational RAC performance.
- **Entanglement fidelity** quantifies how well the channel preserves entanglement with a reference system.
- The hierarchy provides a systematic way to relate ASP to entanglement fidelity in an SDI framework.

---

The resulting fidelity bound is robust against noise and device imperfections, up to the assumed system dimension.

