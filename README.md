# Weyl-Hysteresis-Model

**Non-Local Gravitational Leakage: A Hybrid Braneworld Approach to the Radial Acceleration Relation and Cluster Offsets**

[cite_start]This repository contains the numerical verification and analysis code for a non-local gravitational leakage model derived from a hybrid RSII-DGP braneworld action[cite: 5]. [cite_start]The model accounts for the Radial Acceleration Relation (RAR) and cluster-scale spatial offsets without invoking collisionless dark matter particles[cite: 6].

---

## 🌌 Overview
[cite_start]Traditional $\Lambda$CDM models rely on dark matter particles to explain galactic and cluster dynamics[cite: 15]. [cite_start]This research introduces **Weyl Hysteresis**, a geometric retardation governed by the 5D Bianchi identity, explaining the spatial lag in merging clusters as a gravitational wake effect[cite: 8, 46].

### Key Theoretical Components
* [cite_start]**Effective Propagator**: Derived from the spectral sum of Kaluza-Klein modes, defining an infrared (IR) boost on the 3-brane[cite: 37].
[cite_start]$$G(p) \approx \frac{1}{p^2} \left[ 1 + \left( \frac{a_5}{p^2} \right)^{n/2} \right]^{-1}$$ [cite: 38]
* [cite_start]**Weyl Hysteresis**: A geometric delay ($\Delta t \approx 52$ Myr) between baryonic matter and the gravitational potential[cite: 9, 47].

---

## 📊 Key Results

### 1. Radial Acceleration Relation (RAR)
[cite_start]Analysis of 140 high-quality SPARC galaxies yields a transition index $n = 1.46 \pm 0.05$ and a leakage scale $a_5 \approx 1.2 \times 10^{-10} \text{ m/s}^2$[cite: 7, 64].

### 2. Bullet Cluster Offset (1E 0657-568)
[cite_start]The observed $250$ kpc spatial offset is explained by the finite propagation time of bulk gravitons relative to brane matter[cite: 8, 84].
* [cite_start]**Collision Velocity**: $v_{coll} \approx 4700$ km/s[cite: 48].
* [cite_start]**Retardation Time**: $\Delta t \approx 52$ Myr[cite: 9, 47].
* [cite_start]**Calculated Offset**: $\Delta x \approx 250$ kpc[cite: 49].

---

## 🛠 Usage
The analysis script `analyze.py` verifies the mathematical consistency of the model.

```bash
# Run the verification script
python analyze.py
