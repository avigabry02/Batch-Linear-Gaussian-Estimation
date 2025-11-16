# Batch-Linear-Gaussian-Estimation

# Batch Linear Gaussian Estimation (1D Robot Localization)

## Project Overview

This project implements a **Batch Linear-Gaussian Estimator** to determine the one-dimensional position ($\mathbf{x}_k$) of a mobile robot. The primary objective is to demonstrate robust **sensor fusion** and investigate the trade-off between computational efficiency and estimation consistency using a sparse observation framework.

The estimation fuses high-frequency odometry data with sparse range measurements from a laser rangefinder relative to a fixed landmark.

## Key Technical Contributions

* **Model Implementation:** Implemented the linear discrete-time state-space models for motion (odometry input $u_k$) and observation (transformed laser range $y_k$).
* **Weighted Least Squares (WLS):** Derived and solved the system of Normal Equations ($\mathbf{H}\mathbf{x}^* = \mathbf{b}'$) to find the Maximum A Posteriori (MAP) optimal state sequence $\mathbf{x}_{1:K}^*$.
* **Sparse Batch Estimation:** Developed a modified objective function and derived the solution for a sparse system ($\mathbf{H}_{sparse}$) by integrating odometry data over intervals ($\delta$), drastically reducing the matrix dimension from $12709 \times 12709$ to as low as $12 \times 12$.
* **Statistical Validation:** Performed a rigorous **Consistency Check** (Statistical Validation) across varying sparsity intervals ($\delta=1, 10, 100, 1000$). Validated that the theoretical uncertainty ($\mathbf{H}^{-1}$) matches the empirical error distribution.
* **Observability Analysis:** Conducted an analytical derivation confirming the non-observability of the standalone odometry model (due to a singular Information Matrix) and the necessity of sensor fusion.

## Technologies & Tools

* **Language:** MATLAB / Octave
* **Data Handling:** Linear Interpolation, FFT (Implicit in data analysis)
* **Core Concepts:** Linear-Gaussian Systems, Weighted Least Squares (WLS), Information Matrix ($\mathbf{H}$), Sensor Fusion, Observability.

## Key Findings (Consistency vs. Efficiency)

The project concluded that the full solution ($\delta=1$) was **inconsistent** (too optimistic), despite having the lowest error. The sparse solutions ($\delta \ge 10$) achieved **consistency** while maintaining near-optimal accuracy, proving the sparse batch method is the most **computationally practical and robust** approach for long-duration estimation.

---

## Intended Audience

This project is relevant for researchers, master's students, and engineers interested in:

* **State Estimation and Filtering Theory** (Batch vs. Sequential approaches).
* **Simultaneous Localization and Mapping (SLAM)** fundamentals.
* **Computational Efficiency** in large-scale optimization problems.
* **Sensor Fusion** techniques in mobile robotics.
