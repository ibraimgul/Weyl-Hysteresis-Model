# Weyl-Hysteresis-Model

**Non-Local Gravitational Leakage: A Scalar-Tensor Braneworld Approach to Galactic Dynamics**

This repository contains the numerical verification code and dataset for a non-local gravity model derived from a **5D Einstein-Dilaton action**. The model unifies galactic dynamics (RAR) and cluster-scale spatial offsets (Bullet Cluster) through a single geometric framework involving Lifshitz scaling and metastable bulk gravitons.

---

## 🌌 Abstract & Theory
We propose a UV-complete non-local gravity model where the anomalous scaling of the graviton propagator arises from a **bulk dilaton field** $\phi$. This generates a non-AdS (Lifshitz) background geometry that modifies the Kaluza-Klein spectral density.

### Key Theoretical Mechanisms
1.  **Lifshitz Scaling via Dilaton**: The transition index $n$ is not arbitrary but fixed by the bulk scalar coupling $\gamma$:
    $$ds^2 = \left(\frac{L}{z}\right)^{2\delta} (\eta_{\mu\nu} dx^\mu dx^\nu - dz^2), \quad n = 2\delta - 1$$

2.  **Volcano Potential & Tunneling**: The Bullet Cluster offset is explained not by superluminal travel, but by the **WKB tunneling time** ($\tau$) of metastable gravitons through the effective potential barrier $V_{eff}(z)$ created by the brane.

---

## 📊 Key Results

### 1. Galactic Dynamics (SPARC)
Bayesian MCMC inference on **140 high-quality SPARC galaxies** ($i > 30^\circ, Q \le 2$) yields:
* **Transition Index:** $n = 1.46 \pm 0.05$
* **Implied Lifshitz Exponent:** $\delta \approx 1.23$ (Deviates from standard AdS $\delta=1$)
* **Leakage Scale:** $a_5 \approx 1.2 \times 10^{-10} \text{ m/s}^2$
* **Parsimony:** $\Delta \text{BIC} \approx 18.2$ preference over $\Lambda$CDM+NFW.

![RAR Residuals](Figure_1_RAR_Residuals.png)

### 2. Bullet Cluster Constraint
The observed **250 kpc** offset is reconciled via the decay width of the bulk resonance:
* **Tunneling Time:** $\tau \approx 52 \text{ Myr}$
* **Stability:** The model satisfies the Breitenlohner-Freedman bound ($m^2 \ge -4k^2$) and maintains subluminal signal propagation ($v_g < c$).

![Bullet Offset](Figure_4_Bullet_Offset.png)

---

## 📂 Repository Structure

* `analyze.py`: Python script for verification of Lifshitz scaling parameters, BF stability bounds, and tunneling time consistency.
* `SPARC_Q2.csv`: Processed dataset of 140 galaxies used in the MCMC fit.
* `main.tex`: Source code of the manuscript (RevTeX 4.2).
* `figures/`: High-resolution plots used in the paper.

---

## 🛠 Usage

To reproduce the theoretical consistency checks:

```bash
# Clone the repository
git clone [https://github.com/ibraimgul/Weyl-Hysteresis-Model.git](https://github.com/ibraimgul/Weyl-Hysteresis-Model.git)

# Run the verification script
python analyze.py
