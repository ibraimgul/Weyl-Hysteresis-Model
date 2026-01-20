# Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension

**Author:** Ibrahim Gul  
**Correspondence:** gulpoetika@gmail.com  
**Status:** Under Review (**EPJ C**)

## 📌 Overview
This repository contains the source code, data, and numerical analysis for the paper **"Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension"**. 

The model introduces a massive graviton spectral propagator that facilitates gravitational leakage into a dark dimension at low acceleration scales. This code reproduces the key results, demonstrating that the model resolves galactic and cosmological tensions without particle dark matter.

## 🏆 Key Results
1.  **Galactic Dynamics:** Bayesian inference on the SPARC dataset favors the leakage model ($n > 1$) over standard MOND ($n=1$) with decisive statistical significance.
2.  **Robustness Check:** The analysis consistently processes a "Golden Sample" of **125 high-quality galaxies**. The results ($n \approx 2.15$) provide strong evidence for gravitational leakage, maintaining superior predictive power over larger datasets.
3.  **Cosmology ($S_8$):** The provided `class_mu_patch.diff` implements the model in the **CLASS** Boltzmann code, yielding $S_8 \approx 0.78$ and alleviating the $S_8$ tension.

## ⚠️ Data Selection & Quality Control
To ensure reliable physical constraints, the analysis utilizes a **"Golden Sample"** of 125 galaxies from the SPARC database. Following standard astrophysical practices (Lelli et al. 2016), galaxies are filtered based on:

* **Inclination ($i > 30^\circ$):** To minimize de-projection uncertainties. In face-on galaxies, observational errors are significantly amplified, which can mask the subtle gravitational leakage signal.
* **Quality Flag ($Q=1$):** To ensure regular kinematics and reliable photometry, excluding systems with non-circular motions or tidal disturbances.

## 📂 Repository Structure
* **`analyze.py`**: The main Bayesian analysis script using the Dynesty sampler.
* **`sparc_raw_full_3389_points.csv`**: Raw rotation curve data from the SPARC database.
* **`sparc_subset_125_list.txt`**: The finalized list of 125 high-quality Galaxy IDs used for the analysis.
* **`chain_results.txt` / `analysis_summary.txt`**: Output data from the latest 125-galaxy analysis run.
* **`class_mu_patch.diff`**: Patch for the CLASS Boltzmann code to implement the non-local leakage propagator.
* **`Figure_1` - `Figure_5`**: Official visualizations presented in the manuscript.

## 🚀 Usage
1.  Install dependencies:  
    `pip install -r requirements.txt`
2.  Run analysis:  
    `python analyze.py`

## 📜 Citation
```bibtex
@article{Gul2026Leakage,
  title={Non-Local Gravitational Leakage: A Spectral Propagator Ansatz from the Dark Dimension},
  author={Gul, Ibrahim},
  journal={arXiv preprint / Submitted to EPJ C},
  year={2026}
}
