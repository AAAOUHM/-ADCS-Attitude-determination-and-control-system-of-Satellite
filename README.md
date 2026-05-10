# 🛰️ ADCS — Attitude Determination & Control System of Satellite

> **Quaternion-based Linear MPC for spacecraft attitude control**  
> Developed as part of the AULE Space Inc. GNC Technical Assignment

---

## 📌 Overview

This project implements a full **3-DOF spacecraft Attitude Determination and Control System (ADCS)** using:

- **Quaternion kinematics & dynamics** to avoid gimbal lock
- **Linear Model Predictive Control (LMPC)** as the optimal control law (non-PID)
- **Sinusoidal disturbance torque** in the ECI frame to simulate real space environment
- Numerical linearization via **Jacobian computation** at each timestep

The controller drives the satellite from an initial attitude to a desired final attitude with zero body rates, while rejecting environmental disturbances.

---

## 🎯 Simulation Parameters

| Parameter | Value |
|---|---|
| Initial Attitude (Euler ZYX) | [12°, 31°, −40°] |
| Final Desired Attitude (Euler ZYX) | [40°, −30°, 90°] |
| Initial Body Rates | [10, −20, 12] deg/s |
| Final Desired Body Rates | [0, 0, 0] deg/s |
| Moment of Inertia | [[1.10, 0.05, 0.00], [0.05, 1.90, −0.01], [0.00, −0.01, 1.12]] kg·m² |
| Reaction Wheels | 3 (one per axis) |
| Disturbance Torque | Sinusoidal, 1×10⁻⁶ Nm (ECI frame) |
| Timestep | 0.05 s |
| Prediction Horizon (N) | 5 |

---

## 🧠 Methodology

### 1. Quaternion Representation

The satellite attitude is represented as a unit quaternion **q = [q₀, q₁, q₂, q₃]** where q₀ is the scalar part. This avoids the singularities (gimbal lock) of Euler angle representations.

**Quaternion kinematics:**
```
q̇ = 0.5 × Ω(ω) × q
```

where Ω(ω) is the skew-symmetric matrix built from body angular rates ω.

**Rotational dynamics (Euler's equation, simplified — no gyroscopic term):**
```
ω̇ = I⁻¹ × (τ_control + τ_disturbance)
```

The full state vector is: **x = [q₀, q₁, q₂, q₃, ωx, ωy, ωz]ᵀ** (7 states)

---

### 2. Linear MPC (LMPC)

Since the quaternion dynamics are nonlinear, we **linearize at each timestep** using numerical Jacobians, then discretize using Zero-Order Hold (ZOH):

```
ẋ ≈ A(xₖ)·x + B(xₖ)·u   →   xₖ₊₁ = Ad·xₖ + Bd·uₖ
```

The LMPC solves a **Quadratic Program (QP)** at every step:

```
min   Σ [ (xₖ - x_ref)ᵀ Q (xₖ - x_ref) + uₖᵀ R uₖ ]
U

subject to:  u_min ≤ uₖ ≤ u_max
```

**Prediction matrices** Φ and Γ are built to propagate the state over the horizon N=5, and `quadprog` solves the QP. Only the **first control input** is applied (receding horizon principle).

**Tuning:**
| Parameter | Value | Reasoning |
|---|---|---|
| Q (attitude) | [250, 180, 100] | Aggressive attitude tracking |
| Q (rates) | [10, 15, 5] | Moderate rate damping |
| R | 0.01 | Low control effort penalty |
| τ_max | 10 Nm | Reaction wheel saturation limit |

---

### 3. Disturbance Torque

Environmental disturbances (solar pressure, magnetic, gravity gradient) are modeled as a **sinusoidal torque in the ECI frame**:

```
τ_dist_ECI = A × sin(2π × f × t) × [1, 1, 1]ᵀ
```

Since the dynamics equation lives in the **body frame**, the disturbance is rotated using the current quaternion:

```
τ_dist_body = R(q) × τ_dist_ECI
```

using the quaternion rotation formula:
```
v_body = v + 2·q₀·(qᵥ × v) + 2·(qᵥ × (qᵥ × v))
```

The disturbance enters through the **input channel** of the dynamics, same as control torque — making it possible (in principle) for the controller to reject it.

---

### 4. Can LMPC Reject Disturbance?

**Partially — but not truly.** The LMPC has no internal disturbance model, so:

- It reacts to disturbance effects **after they appear** (one step late)
- Re-linearization every step gives **implicit feedback rejection**
- For tiny disturbances (1e-6 Nm), the effect is negligible
- For larger disturbances, a **steady-state offset** will appear

A proper solution would add **integral augmentation** or an **EKF-based disturbance estimator** — which connects directly to the EKF theory below.

---

## 🔬 EKF for Attitude Estimation (Theory)

In a real AOCS, the true attitude is **never directly measured** — it must be estimated from noisy sensors. The **Extended Kalman Filter (EKF)** is the standard approach.

### Sensor Fusion: IMU + Star Tracker

| Sensor | Provides | Noise |
|---|---|---|
| Gyroscope (IMU) | Angular rates ω (high rate) | Bias drift, white noise |
| Accelerometer (IMU) | Specific force | Low accuracy for attitude |
| Star Tracker | Absolute quaternion attitude (low rate) | High accuracy, slow update |

The EKF fuses these by running two steps every cycle:

**Predict step** (using gyroscope at high rate):
```
x̂ₖ|ₖ₋₁ = f(x̂ₖ₋₁, ωₖ)        % propagate state using dynamics
Pₖ|ₖ₋₁ = Fₖ·Pₖ₋₁·Fₖᵀ + Q    % propagate covariance
```

**Update step** (when Star Tracker measurement arrives):
```
Kₖ = Pₖ|ₖ₋₁·Hᵀ·(H·Pₖ|ₖ₋₁·Hᵀ + R)⁻¹   % Kalman gain
x̂ₖ = x̂ₖ|ₖ₋₁ + Kₖ·(zₖ - h(x̂ₖ|ₖ₋₁))    % state update
Pₖ = (I - Kₖ·H)·Pₖ|ₖ₋₁                  % covariance update
```

After each update, the quaternion estimate **must be renormalized** to stay on the unit sphere.

### Impact of Estimation Errors on Control

If the EKF gives a wrong attitude estimate:
- The **error quaternion** computed by the controller is wrong
- The controller applies torques in the **wrong direction**
- This can cause **limit cycling, overshoot, or divergence**
- Gyro bias drift is especially dangerous — it accumulates over time

### Noise Covariance Tuning Philosophy

| Matrix | Meaning | Tuning Approach |
|---|---|---|
| **Q** (process noise) | How much we trust the dynamics model | Start small; increase if filter diverges |
| **R** (measurement noise) | How much we trust the sensors | Set from sensor datasheet specs |
| **P₀** (initial covariance) | Confidence in initial estimate | Set large if initial attitude is uncertain |

**General rule:** If the filter is **too slow to respond** → increase Q. If it is **too noisy** → increase R. The ratio Q/R matters more than absolute values.

---

## 📁 Repository Structure

```
📦 ADCS-Satellite
 ┣ 📜 LMPC_quaternion.m       # Main simulation script
 ┣ 📜 README.md               # This file
```

---

## ⚙️ Requirements

- MATLAB R2019b or later (R2024b recommended)
- **Toolboxes required:**
  - Optimization Toolbox (`quadprog`)
  - Control System Toolbox (`c2d`, `ss`)
  - Aerospace Toolbox (`eul2quat`, `quat2eul`, `quatmultiply`)

---

## 🚀 How to Run

```matlab
% Simply open and run in MATLAB:
LMPC_quaternion.m
```

The script will produce three figures:
1. **Attitude** — Euler angles converging to target
2. **Angular velocity** — body rates converging to zero
3. **Control torques** — reaction wheel torque history

---

## 📊 Expected Results

- Attitude converges from [12°, 31°, −40°] → [40°, −30°, 90°]
- Body rates converge to [0, 0, 0] rad/s
- Control torques stay within ±10 Nm saturation limits
- Disturbance causes minor steady-state ripple (visible at higher amplitudes)

---

## 👤 Author

**AAAOUHM**  
GNC Technical Assignment — AULE Space Inc.  
*Quaternion-based LMPC for spacecraft attitude control*

---

> *"Control is not about eliminating uncertainty — it's about being robust to it."*
